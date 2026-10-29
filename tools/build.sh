#!/bin/sh
# Builds the native module for this machine and copies it into python/ftalgo. On an arm64 Mac whose
# Rust toolchain is x86_64 (under Rosetta) it cross-targets aarch64-apple-darwin explicitly. The
# "python" profile is release without strip: Rust 1.90's strip leaves the arm64 module's symbol
# table misaligned and macOS refuses to load it.
set -e
cd "$(dirname "$0")/.."
rm -f dist/ftalgo-*.whl
TARGET=""
if [ "$(uname -s)" = Darwin ] && [ "$(uname -m)" = arm64 ] && rustc -vV | grep -q '^host: x86_64'; then
    TARGET="--target aarch64-apple-darwin"
fi
.venv/bin/maturin build --profile python $TARGET --out dist -q
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
# The venv imports python/ftalgo through a redirect package (no installed wheel). Neither a .pth
# file (macOS marks new files here hidden; Python 3.13 skips hidden .pth files) nor
# sitecustomize.py (Anaconda's own shadows it) works here.
SP=$(.venv/bin/python -c "import site; print(site.getsitepackages()[0])")
mkdir -p "$SP/ftalgo"
cat > "$SP/ftalgo/__init__.py" <<PY
import os as _os
__path__ = ["$PWD/python/ftalgo"]
__file__ = _os.path.join(__path__[0], "__init__.py")
with open(__file__) as _f:
    exec(compile(_f.read(), __file__, "exec"))
PY
