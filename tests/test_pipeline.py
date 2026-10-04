"""Basic sanity tests. Run with:  pytest -q"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import cas12a_drylab_pipeline as m


def test_revcomp():
    assert m.revcomp("AACGT") == "ACGTT"


def test_find_guides_requires_tttv_pam():
    seq = "GGGG" + "TTTA" + "ACGTACGTACGTACGTACGT" + "CCCC"
    guides = m.find_guides(seq)
    assert any(g["strand"] == "+" and g["spacer"] == "ACGTACGTACGTACGTACGT" for g in guides)
    assert all(g["pam"][:3] == "TTT" and g["pam"][3] in "ACG" for g in guides)


def test_perfect_match_activity_is_one():
    sp = "ACGTACGTACGTACGTACGT"
    ref = "GGGGTTTA" + sp + "CCCC"
    assert m.best_activity(sp, ref, m.W_PRIOR) == 1.0


def test_seed_mismatch_hurts_more_than_distal():
    sp = "ACGTACGTACGTACGTACGT"
    seed = list(sp); seed[1] = "T"
    distal = list(sp); distal[18] = "A"
    mk = lambda s: "GGGGTTTA" + "".join(s) + "CCCC"
    assert m.best_activity(sp, mk(seed), m.W_PRIOR) < m.best_activity(sp, mk(distal), m.W_PRIOR)


def test_kinetics_monotonic_in_concentration():
    t = np.arange(0, 1801, 60.0)
    lo = m.simulate(m.KINETIC_TRUE, 0.1, t)[-1]
    hi = m.simulate(m.KINETIC_TRUE, 10.0, t)[-1]
    assert hi > lo


def test_calibration_recovers_synthetic_parameters():
    rng = np.random.default_rng(1)
    fit = m.fit_kinetics(m.synth_wetlab(rng))
    assert abs(fit["kon"] / m.KINETIC_TRUE["kon"] - 1) < 0.15


def test_mismatch_calibration_recovers_weights():
    rng = np.random.default_rng(2)
    w_true = np.array([2.6] * 5 + [1.0] * 5 + [0.4] * 6 + [0.1] * 4)
    M = np.zeros((60, m.L), int)
    for i in range(60):
        M[i, rng.choice(m.L, rng.integers(1, 3), replace=False)] = 1
    act = np.minimum(np.exp(-(M @ w_true)) * np.exp(rng.normal(0, .1, 60)), 1.0)
    assert np.corrcoef(w_true, m.calibrate_mismatch(M, act))[0, 1] > 0.9
