# Releasing ftalgo (PyPI, GitHub Release, Zenodo DOI)

A version tag `v*` runs `.github/workflows/wheels.yml`: abi3 wheels for Linux (x86_64 and aarch64,
glibc and musl), macOS (universal2) and Windows (x64), and the sdist, each carrying the fitted model
and the rotation cache (`tools/bundle_data.sh`); a smoke test of the installed wheel on each
operating system (`tools/smoke_installed.py`); then, on the tag only, upload to PyPI by trusted
publishing and a GitHub Release with the files attached. Zenodo archives the release and mints a DOI.

The workflow can also be run by hand (Actions → Wheels → Run workflow): it builds and smoke-tests
without publishing. GitHub creates the `pypi` environment the publishing job names the first time it
runs; no secrets are needed, since trusted publishing uses the workflow's identity.

## Once, before the first release

1. **PyPI trusted publisher.** At [pypi.org/manage/account/publishing](https://pypi.org/manage/account/publishing/),
   under "Add a new pending publisher" → GitHub: PyPI project name `ftalgo`, owner `Jaspersands`,
   repository `ft-algorithms`, workflow name `wheels.yml`, environment name `pypi`.
2. **Zenodo.** At [zenodo.org/account/settings/github](https://zenodo.org/account/settings/github/),
   press **Sync now** and switch on `Jaspersands/ft-algorithms`. `.zenodo.json` supplies the
   record's metadata.

## Each release

1. Set the version in `pyproject.toml` and `Cargo.toml` (they match) and in `CITATION.cff`.
2. Tag and push: `git tag v0.1.0 && git push origin v0.1.0`.
3. When the workflow finishes, the version is on PyPI (`pip install ftalgo`), the GitHub Release
   exists, and within minutes Zenodo shows its DOI and the concept DOI (which always resolves to
   the latest version).

After the first release, put the concept DOI (`10.5281/zenodo.NNNNNNN`) in `CITATION.cff`
(`doi:`) and add its badge to the top of `README.md`:
`[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.NNNNNNN.svg)](https://doi.org/10.5281/zenodo.NNNNNNN)`.
