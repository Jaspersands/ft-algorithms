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
