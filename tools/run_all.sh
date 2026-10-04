#!/bin/sh
# The experiment chain (resumable): cross-check, Shor grid, phase estimation.
cd "$(dirname "$0")/.."
.venv/bin/python tools/xcheck.py
.venv/bin/python tools/run_shor.py
.venv/bin/python tools/run_qpe.py --mol H2 --factory cultivation 15to1 --shots 600
.venv/bin/python tools/run_qpe.py --mol HeH+ --factory cultivation --shots 600
