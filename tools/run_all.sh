#!/bin/sh
# The experiment chain (resumable): cross-check, Shor grid, phase estimation.
cd "$(dirname "$0")/.."
[ -f data/calibration/xcheck.json ] || .venv/bin/python tools/xcheck.py
.venv/bin/python tools/run_shor.py --N 15 21 35 77 143
.venv/bin/python tools/run_qpe.py --mol H2 --factory cultivation 15to1 --shots 600
.venv/bin/python tools/run_qpe.py --mol HeH+ --factory cultivation --shots 600
.venv/bin/python tools/qpe_hist.py
