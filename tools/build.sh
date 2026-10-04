#!/bin/sh
# Builds the native arm64 module (the host Rust toolchain is x86_64, so cross-target explicitly)
# and installs the wheel into .venv. The "python" profile is release without strip: Rust 1.90's
# strip leaves the arm64 module's symbol table misaligned and macOS refuses to load it.
set -e
cd "$(dirname "$0")/.."
rm -f dist/ftalgo-*.whl
.venv/bin/maturin build --profile python --target aarch64-apple-darwin --out dist -q
.venv/bin/pip uninstall -q -y ftalgo 2>/dev/null || true
.venv/bin/pip install -q --no-deps dist/ftalgo-*.whl
# Also place the module next to the sources, so tests and tools import python/ftalgo directly
# (conftest.py puts python/ first on sys.path) and Python edits need no rebuild.
.venv/bin/python - <<'PY'
import glob, shutil, zipfile
whl = glob.glob("dist/ftalgo-*.whl")[0]
with zipfile.ZipFile(whl) as z:
    name = [n for n in z.namelist() if n.startswith("ftalgo/_ftsim")][0]
    with z.open(name) as src, open("python/" + name, "wb") as dst:
        shutil.copyfileobj(src, dst)
PY
