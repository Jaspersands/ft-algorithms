#!/bin/sh
# Copies the committed fitted models and rotation cache into python/ftalgo/_data, where the next
# wheel built carries them: an installed ftalgo has no data/ directory. The wheels workflow runs it.
set -e
cd "$(dirname "$0")/.."
mkdir -p python/ftalgo/_data
cp data/calibration/model.json python/ftalgo/_data/model.json
cp data/calibration/model_mr.json python/ftalgo/_data/model_mr.json
cp data/synth/cache.json python/ftalgo/_data/synth_cache.json
