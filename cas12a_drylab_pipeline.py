#!/usr/bin/env python3
"""
cas12a_drylab_pipeline.py
=========================
Dry-lab "digital twin" pipeline for a CRISPR-Cas12a (+RPA) molecular diagnostic.

Idea: every step that produces a *number* in silico has a matching wet-lab
measurement that can overwrite or recalibrate it. Nothing here is a black box.

  Step 1  Guide discovery      : TTTV-PAM scan on both strands
  Step 2  Target validation    : conservation across variant genomes (position-weighted mismatch model)
  Step 3  Specificity          : cross-reactivity vs background (host / co-pathogens)
  Step 4  RPA primer picking   : heuristic amplicon design around the best guide
  Step 5  Kinetic model (ODE)  : Cas12a trans-cleavage + optional RPA amplification
  Step 6  Wet-lab calibration  : fit kon, kcat/KM to plate-reader curves (+ CI)
  Step 7  Mismatch calibration : learn position weights from a measured mismatch panel
  Step 8  LoD / time-to-result : predicted before any wet-lab run

USAGE
  python cas12a_drylab_pipeline.py --demo                       # synthetic end-to-end run
  python cas12a_drylab_pipeline.py --target ref.fa --variants variants.fa \
        --background host_and_others.fa --wetlab plate.csv       # real data

  plate.csv columns: time_s, conc_nM, rep, fluor_norm   (0-1 normalised to full cleavage)

DEMO NOTE: --demo uses SYNTHETIC sequences and SYNTHETIC "wet-lab" data. It shows the
machinery works (e.g. parameter recovery). It is NOT a claim about any real assay.
Kinetic constants are order-of-magnitude literature-style priors; replace with fits.
"""
import argparse, json, os
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, nnls
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------- utilities
COMP = str.maketrans("ACGT", "TGCA")
ENC = {"A": 0, "C": 1, "G": 2, "T": 3}
L = 20                                    # spacer length
# Prior per-position mismatch penalty (pos 1 = PAM-proximal). Cas12a seed ~ 1-5.
# This is a HYPOTHESIS to be replaced by Step 7 once a mismatch panel is measured.
W_PRIOR = np.array([2.0] * 5 + [1.2] * 5 + [0.5] * 6 + [0.2] * 4)


def revcomp(s):
    return s.translate(COMP)[::-1]


def enc(s):
    return np.array([ENC.get(c, 4) for c in s], dtype=np.int8)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.upper().replace("U", "T"))
    if name:
        seqs[name] = "".join(buf)
    return seqs


# ------------------------------------------------- Step 1: guide discovery
def find_guides(seq):
    out = []
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for i in range(len(s) - 4 - L + 1):
            pam = s[i:i + 4]
            if pam[:3] == "TTT" and pam[3] in "ACG":
                sp = s[i + 4:i + 4 + L]
                start = i + 4 if strand == "+" else len(seq) - (i + 4 + L)
                out.append(dict(strand=strand, pam=pam, spacer=sp, start=start))
    return out


# --------------------------------------- Steps 2-3: mismatch-aware scanning
def scan_penalties(spacer, ref, weights):
    """Position-weighted mismatch penalty of `spacer` against every PAM-valid window in ref."""
    sp, pens = enc(spacer), []
    for s in (ref, revcomp(ref)):
        if len(s) < L + 4:
            continue
        e = enc(s)
        win = sliding_window_view(e, L)[4:]
        j = np.arange(4, len(e) - L + 1)
        pam_ok = (e[j - 4] == 3) & (e[j - 3] == 3) & (e[j - 2] == 3) & (e[j - 1] < 3)
        pens.append(((win != sp)[pam_ok]) @ weights)
    return np.concatenate(pens) if pens else np.array([np.inf])


def best_activity(spacer, ref, weights):
    """Predicted relative trans-cleavage activity (1 = perfect match) at best site."""
    p = scan_penalties(spacer, ref, weights)
    return float(np.exp(-p.min())) if p.size else 0.0


def gc(s):
    return (s.count("G") + s.count("C")) / len(s)


def homopolymer(s, n=5):
    return any(c * n in s for c in "ACGT")


def score_guides(target, variants, background, weights):
    rows = []
    for g in find_guides(target):
        sp = g["spacer"]
        cov = [best_activity(sp, v, weights) for v in variants.values()] or [1.0]
        off = [best_activity(sp, b, weights) for b in background.values()] or [0.0]
        gc_pen = max(0.0, 0.4 - gc(sp), gc(sp) - 0.6)
        seq_term = max(0.0, 1 - 2 * gc_pen) * (0.5 if homopolymer(sp) else 1.0)
        rows.append(dict(**g, gc=round(gc(sp), 2),
                         conservation_mean=round(float(np.mean(cov)), 3),
                         frac_variants_ge50pct=round(float(np.mean(np.array(cov) >= .5)), 3),
                         worst_variant=round(float(np.min(cov)), 3),
                         max_offtarget_activity=round(float(np.max(off)), 3),
                         score=round(float(np.mean(cov) * (1 - np.max(off)) * seq_term), 4)))
    return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)


# ------------------------------------------------ Step 4: RPA primer picking
def primer_ok(p):
    return 0.35 <= gc(p) <= 0.65 and not homopolymer(p, 5) and p[-3:].count("G") + p[-3:].count("C") <= 2


def dimer(a, b):
    return revcomp(a[-5:]) in b or revcomp(b[-5:]) in a or revcomp(a[-5:]) in a or revcomp(b[-5:]) in b


def pick_rpa(target, g, plen=32, amp=(110, 220)):
    lo, hi = g["start"] - 4, g["start"] + L + 4       # PAM + protospacer footprint (fwd coords)
    best = None
    for f in range(max(0, lo - 90), max(0, lo - 10)):
        fp = target[f:f + plen]
        if len(fp) < plen or not primer_ok(fp):
            continue
        for r_end in range(hi + 10, min(len(target), hi + 90)):
            rp = revcomp(target[r_end - plen:r_end])
            size = r_end - f
            if not (amp[0] <= size <= amp[1]) or not primer_ok(rp) or dimer(fp, rp):
                continue
            cost = abs(gc(fp) - .5) + abs(gc(rp) - .5) + abs(size - 150) / 300
            if best is None or cost < best[0]:
                best = (cost, dict(fwd=fp, rev=rp, amplicon_bp=size, fwd_start=f, rev_end=r_end))
    return best[1] if best else None


# ------------------------------------------- Step 5: kinetic model (ODE)
KINETIC_TRUE = dict(kon=1e-4, kcat=5.0, KM=800.0)     # nM^-1 s^-1, s^-1, nM (priors)
R0, S0 = 50.0, 500.0                                   # RNP and reporter, nM
COPIES_TO_NM = 1e6 * 1e9 / 6.022e23                    # copies/uL -> nM


def target_conc(t, T0, amplify, r=0.02, Tmax=100.0, lag=60.0):
    if not amplify:
        return T0
    x = T0 * np.exp(r * np.maximum(t - lag, 0.0))      # logistic RPA growth
    return Tmax * x / (Tmax + x - T0)


def simulate(params, T0, t_eval, amplify=False, mm_factor=1.0):
    kon, kcat, KM = params["kon"] * mm_factor, params["kcat"], params["KM"]

    def rhs(t, y):
        A, P = y
        Tfree = max(target_conc(t, T0, amplify) - A, 0.0)
        return [kon * max(R0 - A, 0.0) * Tfree, kcat * A * max(S0 - P, 0.0) / (KM + max(S0 - P, 0.0))]

    sol = solve_ivp(rhs, (0, t_eval[-1]), [0.0, 0.0], t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-12)
    return sol.y[1] / S0                               # fraction of reporter cleaved


# -------------------------------------------- Step 6: wet-lab calibration
def fit_kinetics(df, KM=KINETIC_TRUE["KM"]):
    """Fit log(kon), log(kcat) on amplification-free curves. KM fixed: with S0 < KM only kcat/KM is identifiable."""
    groups = [(c, g["time_s"].values, g["fluor_norm"].values) for c, g in df.groupby(["conc_nM", "rep"])]

    def resid(theta):
        p = dict(kon=np.exp(theta[0]), kcat=np.exp(theta[1]), KM=KM)
        return np.concatenate([simulate(p, c[0], t) - y for c, t, y in groups])

    res = least_squares(resid, x0=np.log([3e-5, 2.0]), method="lm")
    dof = max(len(res.fun) - 2, 1)
    cov = np.linalg.inv(res.jac.T @ res.jac) * (res.fun @ res.fun / dof)
    sd = np.sqrt(np.diag(cov))
    return dict(kon=float(np.exp(res.x[0])), kcat=float(np.exp(res.x[1])), KM=KM,
                kon_CI95=[float(np.exp(res.x[0] - 1.96 * sd[0])), float(np.exp(res.x[0] + 1.96 * sd[0]))],
                kcat_over_KM=float(np.exp(res.x[1]) / KM), rmse=float(np.sqrt(np.mean(res.fun ** 2))))


def synth_wetlab(rng, concs=(0.05, 0.1, 0.3, 1, 3, 10), reps=3, sd=0.015):
    t = np.arange(0, 3601, 60.0)
    rows = []
    for c in concs:
        clean = simulate(KINETIC_TRUE, c, t)
        for r in range(reps):
            y = clean * (1 + rng.normal(0, .03)) + rng.normal(0, sd, t.size)
            rows += [dict(time_s=ti, conc_nM=c, rep=r, fluor_norm=yi) for ti, yi in zip(t, y)]
    return pd.DataFrame(rows)


# --------------------------------------- Step 7: mismatch-weight calibration
def calibrate_mismatch(M, rel_activity):
    """M: (n_guides, 20) 0/1 mismatch matrix; rel_activity: measured activity / perfect-match activity.
    Model: -ln(activity) = M @ w with w >= 0 (non-negative least squares)."""
    w, _ = nnls(M.astype(float), -np.log(np.clip(rel_activity, 1e-3, 1.0)))
    return w


# ------------------------------------------------- Step 8: LoD prediction
def time_to_positive(params, copies_per_ul, mm_factor=1.0, thresh=0.10, tmax=2700):
    t = np.arange(0, tmax + 1, 15.0)
    f = simulate(params, copies_per_ul * COPIES_TO_NM, t, amplify=True, mm_factor=mm_factor)
    idx = np.argmax(f >= thresh)
    return float(t[idx] / 60) if f[idx] >= thresh else np.nan


# --------------------------------------------------------------- demo data
def demo_inputs(rng):
    rnd = lambda n: "".join(rng.choice(list("ACGT"), n))
    target = rnd(900)
    variants = {}
    for i in range(12):
        v = list(target)
        for pos in rng.choice(len(v), rng.integers(4, 9), replace=False):
            v[pos] = rng.choice([b for b in "ACGT" if b != v[pos]])
        variants[f"variant_{i+1}"] = "".join(v)
    background = {"host_chunk": rnd(30000)}
    g = find_guides(target)
    # plant two decoys so the specificity filter has something to catch
    d1 = list(g[0]["spacer"]); d1[14] = "A" if d1[14] != "A" else "C"; d1[18] = "G" if d1[18] != "G" else "T"
    d2 = list(g[1]["spacer"]); d2[11] = "A" if d2[11] != "A" else "C"
    background["decoy_2mm_distal"] = rnd(200) + "TTTA" + "".join(d1) + rnd(200)
    background["decoy_1mm_mid"] = rnd(200) + "TTTC" + "".join(d2) + rnd(200)
    return target, variants, background


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--target"); ap.add_argument("--variants"); ap.add_argument("--background")
    ap.add_argument("--wetlab"); ap.add_argument("--out", default="drylab_out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(7)

    # ---- inputs
    if a.demo:
        target, variants, background = demo_inputs(rng)
        wet = synth_wetlab(rng)
    else:
        target = next(iter(read_fasta(a.target).values()))
        variants = read_fasta(a.variants) if a.variants else {}
        background = read_fasta(a.background) if a.background else {}
        wet = pd.read_csv(a.wetlab) if a.wetlab else None

    # ---- steps 1-3
    guides = score_guides(target, variants, background, W_PRIOR)
    guides.to_csv(f"{a.out}/guides_ranked.csv", index=False)
    top = guides.iloc[0].to_dict()
    print(f"[1-3] {len(guides)} candidate guides; top: {top['spacer']} ({top['strand']}, start {top['start']}) "
          f"score={top['score']} cons={top['conservation_mean']} offtarget={top['max_offtarget_activity']}")

    # ---- step 4
    primers = pick_rpa(target, top)
    json.dump(primers, open(f"{a.out}/rpa_primers.json", "w"), indent=2)
    print(f"[4]   RPA primers: {primers}")

    # ---- step 6 (calibration) - falls back to priors if no wet-lab data
    if wet is not None:
        fit = fit_kinetics(wet)
        print(f"[6]   fitted kon={fit['kon']:.2e} (95% CI {fit['kon_CI95'][0]:.2e}-{fit['kon_CI95'][1]:.2e}), "
              f"kcat/KM={fit['kcat_over_KM']:.2e}, RMSE={fit['rmse']:.4f}")
        if a.demo:
            print(f"      (synthetic truth: kon={KINETIC_TRUE['kon']:.2e}, kcat/KM={KINETIC_TRUE['kcat']/KINETIC_TRUE['KM']:.2e})")
        params = dict(kon=fit["kon"], kcat=fit["kcat"], KM=fit["KM"])
        json.dump(fit, open(f"{a.out}/calibration.json", "w"), indent=2)
    else:
        params = dict(KINETIC_TRUE)
        print("[6]   no wet-lab data given -> using literature-style priors")

    # ---- hold-out validation (demo): predict an unseen concentration
    t = np.arange(0, 3601, 60.0)
    hold_c = 0.03
    obs = simulate(KINETIC_TRUE, hold_c, t) + rng.normal(0, .015, t.size)
    pred = simulate(params, hold_c, t)
    print(f"[6]   hold-out {hold_c} nM: prediction RMSE = {np.sqrt(np.mean((pred - obs) ** 2)):.4f} (noise sd 0.015)")

    # ---- step 7 (demo: recover hidden mismatch weights from a synthetic panel)
    w_true = np.array([2.6] * 5 + [1.0] * 5 + [0.4] * 6 + [0.1] * 4)
    M = np.zeros((60, L), int)
    for i in range(60):
        M[i, rng.choice(L, rng.integers(1, 3), replace=False)] = 1
    act = np.exp(-(M @ w_true)) * np.exp(rng.normal(0, .12, 60))
    w_fit = calibrate_mismatch(M, np.minimum(act, 1.0))
    print(f"[7]   mismatch-weight recovery: r = {np.corrcoef(w_true, w_fit)[0, 1]:.2f} (synthetic panel of 60 guides)")

    # ---- step 8
    copies = np.array([3, 10, 30, 100, 300, 1000, 1e4, 1e5])
    ttp_wt = [time_to_positive(params, c) for c in copies]
    ttp_mm = [time_to_positive(params, c, mm_factor=float(np.exp(-2.0))) for c in copies]  # one seed mismatch
    pd.DataFrame(dict(copies_per_uL=copies, ttp_min_perfect=ttp_wt, ttp_min_seed_mismatch=ttp_mm)) \
        .to_csv(f"{a.out}/predicted_time_to_positive.csv", index=False)
    print("[8]   predicted time-to-positive (min):", {int(c): round(x, 1) for c, x in zip(copies, ttp_wt)})
    print("      Poisson floor: 95% detection needs ~3 template copies per reaction "
          "(e.g. 2 uL input -> >=1.5 copies/uL), regardless of kinetics.")

    # ---- dashboard figure
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    top10 = guides.head(10)
    ax[0, 0].barh(range(len(top10))[::-1], top10["score"], color="#2b6cb0")
    ax[0, 0].set_yticks(range(len(top10))[::-1]); ax[0, 0].set_yticklabels([s[:10] + "…" for s in top10["spacer"]], fontsize=7)
    ax[0, 0].set_title("A. Guide ranking (conservation × specificity × sequence)"); ax[0, 0].set_xlabel("composite score")
    if wet is not None:
        for c, g in wet.groupby("conc_nM"):
            m = g.groupby("time_s")["fluor_norm"].mean()
            ax[0, 1].plot(m.index / 60, m.values, "o", ms=2.5)
            ax[0, 1].plot(t / 60, simulate(params, c, t), "-", lw=1)
    ax[0, 1].set_title("B. Wet-lab curves (dots) vs calibrated model (lines)")
    ax[0, 1].set_xlabel("min"); ax[0, 1].set_ylabel("fraction reporter cleaved")
    ax[1, 0].semilogx(copies, ttp_wt, "o-", label="perfect match")
    ax[1, 0].semilogx(copies, ttp_mm, "s--", label="1 seed mismatch (kon × e^-2)")
    ax[1, 0].set_title("C. Predicted RPA–Cas12a time-to-positive"); ax[1, 0].set_xlabel("copies / µL"); ax[1, 0].set_ylabel("min to 10% cleavage"); ax[1, 0].legend()
    ax[1, 1].bar(np.arange(1, L + 1) - .2, w_true, .4, label="true (hidden)")
    ax[1, 1].bar(np.arange(1, L + 1) + .2, w_fit, .4, label="recovered from panel")
    ax[1, 1].set_title("D. Mismatch penalty by position (1 = PAM-proximal)"); ax[1, 1].set_xlabel("spacer position"); ax[1, 1].legend()
    plt.tight_layout(); plt.savefig(f"{a.out}/dashboard.png", dpi=130)
    print(f"outputs in ./{a.out}/")


if __name__ == "__main__":
    main()
