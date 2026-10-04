# Cas12a Dry-Lab Pipeline

![tests](https://github.com/aaradhyaaggarwal2701/cas12a-drylab-pipeline/actions/workflows/ci.yml/badge.svg)
![license](https://img.shields.io/badge/license-MIT-blue)

A small, open Python pipeline for designing and modelling **CRISPR-Cas12a (+RPA) molecular diagnostics**, built around one rule: **every number the computer predicts should have a wet-lab measurement that can confirm or overwrite it.**

> **Important:** everything in `examples/` and `sample_output/` is **simulated** (random sequences, simulated plate-reader data). It is a software test, not a result about any real assay. Kinetic constants and mismatch weights are placeholders meant to be replaced by values fitted from real data.

![Demo dashboard](demo_dashboard.png)

*Demo dashboard (simulated data). A: guide ranking. B: plate-reader curves vs calibrated model. C: predicted time-to-positive vs copy number. D: mismatch penalty per guide position, true vs recovered.*

## The idea in one paragraph

A Cas12a test has three steps: RPA copies the target, the guide helps Cas12a find it, and Cas12a cuts a reporter that gives a signal. Each step has a measurable speed. This pipeline describes the steps with simple equations, lets your bench data set the speeds, and then predicts how new guide/primer designs should behave, so fewer designs need to be tried at the bench.

## What it does

| Step | What it does | Wet-lab data that can replace or check it |
|---|---|---|
| 1. Guide discovery | Scans both strands for TTTV-PAM protospacers (20 nt) | - |
| 2. Target validation | Predicts guide activity against every variant genome (position-weighted mismatch model) | Variant sequences, measured activity |
| 3. Specificity | Flags guides with high predicted activity on background sequences (host, other pathogens) | Cross-reactivity tests |
| 4. RPA primers | Picks a primer pair (32 nt, 110-220 bp amplicon) around the best guide | Amplification tests |
| 5. Kinetic model | ODE model of Cas12a trans-cleavage with optional RPA amplification | - |
| 6. Calibration | Fits binding rate (kon) and kcat/KM to plate-reader curves, with 95% CI | Fluorescence time courses |
| 7. Mismatch calibration | Learns per-position mismatch penalties from a mismatch panel (non-negative least squares) | Mismatch panel |
| 8. Prediction | Predicts time-to-positive across target copy numbers | Time-to-positive, LoD |

## Quick start

```bash
git clone https://github.com/aaradhyaaggarwal2701/cas12a-drylab-pipeline.git
cd cas12a-drylab-pipeline
pip install -r requirements.txt

# 1) fully synthetic end-to-end demo
python cas12a_drylab_pipeline.py --demo

# 2) same pipeline through the real-data interface, using the example files
python cas12a_drylab_pipeline.py \
  --target examples/target_reference.fa \
  --variants examples/variants.fa \
  --background examples/background.fa \
  --wetlab examples/plate_synthetic.csv \
  --out my_run
```

Outputs: `guides_ranked.csv`, `rpa_primers.json`, `calibration.json`, `predicted_time_to_positive.csv`, `dashboard.png`. A copy of one run is in `sample_output/`.

## Using real data

| Argument | Format |
|---|---|
| `--target` | FASTA, one reference sequence of the target region |
| `--variants` | FASTA, one entry per variant/lineage genome (aligned region or full) |
| `--background` | FASTA, host/other-pathogen sequences to check cross-reactivity |
| `--wetlab` | CSV with columns `time_s, conc_nM, rep, fluor_norm` (fluorescence normalised 0-1, 1 = all reporter cleaved). See `examples/plate_synthetic.csv` |

Use amplification-free titrations (about 5-6 concentrations, 3 replicates) for calibration first.

## Tests

```bash
pip install -r requirements-dev.txt
pytest -q
```

The tests check PAM detection, that seed mismatches hurt more than distal ones, and that calibration recovers known simulated parameters. Passing them shows the code is internally consistent. It does **not** show the model matches real biology; that needs wet-lab data.

## Assumptions and limits

- Mismatch weights are fixed priors (strong in the PAM-proximal seed, weak far from the PAM) until calibrated (step 7).
- With a single reporter concentration below KM, only kcat/KM is identifiable, so KM is fixed during fitting.
- RPA is modelled as simple logistic growth; primer-dimers and stochastic effects are not modelled.
- No reliable test can detect fewer than a few template copies per reaction (Poisson sampling).
- Guide RNA secondary structure is not yet scored.

## Roadmap

- [ ] Calibrate on real plate-reader data and report prediction error on held-out guides
- [ ] Lateral-flow readout model (band intensity vs cleaved reporter)
- [ ] Secondary-structure scoring (ViennaRNA) and primer checks (Primer3)
- [ ] Cross-checks with CRISPOR / Cas-OFFinder
- [ ] Variant coverage from public genomes (e.g. SARS-CoV-2, M. tuberculosis lineages)

## Status

Prototype, validated on simulated data only. Looking for real wet-lab data to calibrate and blind-test it.

## Author

Aaradhya Aggarwal, M.Sc. Biotechnology (Bioinformatics), TERI SAS.
[LinkedIn](https://www.linkedin.com/in/aaradhya-aggarwal-biotechnology/)
