# Contributing

Contributions are welcome: bug reports, documentation fixes, validation results on new platforms,
and code.

## Development Setup

```bash
git clone https://github.com/Breeding-Insight/MADCATS.git
cd MADCATS

# Editable install: the madcats command runs your working copy
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

No other dependencies are needed to develop or test. The `filter-paralogs` tests that align to a
synthetic genome run only when `mappy` is installed (`pip install -e '.[genome]'`); without it they
are skipped.

## Running Tests

```bash
python3 -m unittest discover -s tests -v
```

The suite covers evidence routes, primer resolution, design checks, read accounting, inference,
indel support, normalization and window freezing, output framing, paired-end merging and tiling,
window-scan, gapped alignment, denoising and partial reads. Tests run on every push and pull
request across Python 3.9–3.13 (`.github/workflows/tests.yml`).

Two tests guard the documentation. `tests/test_docs.py` checks that each page in
`docs/src/commands/` lists exactly the options that command's `--help` shows, with matching
choices and defaults, so update the page whenever you add or change an option.
`tests/test_example.py` runs the example in `examples/madcats-demo` and checks the output against
its simulated truth. If a change alters the example's numbers, update
[Try the Example](../getting-started/example.md) to match. To change the example data itself,
edit and rerun `examples/madcats-demo/make_demo_data.py`.

Changes that touch counting should also be checked against real data. The standard is byte-for-byte:
a refactor regenerates existing tables identically from cache, and a behaviour change is measured
against a comparator before it becomes a default. Record those measurements in
[`docs/DEVELOPMENT.md`](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md).

## Code Style

- Standard library only. A new dependency needs a strong reason; the optional `mappy` extra for
  `filter-paralogs` is imported only inside that command's module, `madcats_paralogs.py`.
- Keep the section layout of `madcats.py` (see [Architecture](architecture.md)).
- Every analysis threshold is a command-line option with a documented default.
- Stop with a clear error on inconsistent input instead of writing a table that may be wrong.

## Building the Documentation

The docs are an [mdBook](https://rust-lang.github.io/mdBook/) in `docs/`, published to Read the
Docs from `.readthedocs.yaml`.

```bash
cargo install mdbook        # once
mdbook serve docs --open    # live preview at http://localhost:3000
mdbook build docs           # static site in docs/book/
```

## Releasing

Releases are driven by a version tag. The version lives in one place, `__version__` in
`madcats.py`; `pyproject.toml` reads it from there.

1. Bump `__version__` in `madcats.py` and add an entry to `docs/src/development/changelog.md`.
2. Commit, then tag and push:

   ```bash
   git tag v0.14.0
   git push origin v0.14.0
   ```

3. `.github/workflows/release.yml` then:
   - checks that the tag matches `__version__`
   - runs the tests
   - builds the sdist and wheel and publishes them to PyPI (trusted publishing)
   - builds standalone binaries with PyInstaller on Linux x86_64, macOS Intel and macOS Apple
     Silicon, smoke-tests each one, and attaches them to a GitHub release

### Building locally

```bash
# Python package (sdist + wheel in dist/)
pip install build twine
python -m build
twine check dist/*

# Standalone binary for the current platform (dist/madcats)
pip install pyinstaller
pyinstaller --onefile --name madcats --paths . packaging/pyinstaller_entry.py
./dist/madcats --version
```

### One-time setup

- PyPI: on [pypi.org](https://pypi.org), add a *pending trusted publisher* for project
  `madcats` with owner `Breeding-Insight`, repository `MADCATS`, workflow `release.yml` and
  environment `pypi`. No API token is stored in GitHub.
- Read the Docs: import the GitHub repository at
  [readthedocs.org](https://readthedocs.org/dashboard/import/) under the project slug
  `madcats`. Builds then run on every push.

## Reporting Issues

Open an issue at [github.com/Breeding-Insight/MADCATS/issues](https://github.com/Breeding-Insight/MADCATS/issues)
with the command you ran, the `madcats --version` output, and the relevant part of the QC JSON.
