# Installation

MADCATS is a pure-Python tool with no third-party dependencies. Install it from PyPI (recommended),
download a standalone binary, or run it straight from a clone of the repository. The one optional
command that aligns to a reference genome, `filter-paralogs`, needs an extra
([below](#optional-genome-aligner)).

## From PyPI (Recommended)

Requires Python 3.9 or newer.

```bash
pip install madcats
```

This puts a `madcats` command on your `PATH`:

```bash
madcats --version
madcats --help
```

To keep MADCATS isolated from other Python packages, install it with
[pipx](https://pipx.pypa.io/) instead:

```bash
pipx install madcats
```

Upgrade later with `pip install --upgrade madcats` (or `pipx upgrade madcats`).

### Optional: Genome Aligner

[`madcats filter-paralogs`](../commands/filter-paralogs.md) aligns haplotypes to a reference genome
with minimap2's Python binding, `mappy`. Install it with the `genome` extra:

```bash
pip install 'madcats[genome]'
```

`mappy` compiles minimap2 when no wheel exists for your platform, which needs a C compiler and zlib
headers. Every other command runs without it, and only `filter-paralogs` imports it. The standalone
binaries do not include it.

## Pre-built Binaries

Standalone executables bundle their own Python interpreter, so they run on machines without Python.
Download the latest release for your platform from
[GitHub Releases](https://github.com/Breeding-Insight/MADCATS/releases):

| Platform | Download |
|----------|----------|
| Linux (x86_64) | `madcats-linux-x86_64.tar.gz` |
| macOS (Intel) | `madcats-macos-x86_64.tar.gz` |
| macOS (Apple Silicon) | `madcats-macos-aarch64.tar.gz` |

### Linux

```bash
# Download the latest release
curl -LO https://github.com/Breeding-Insight/MADCATS/releases/latest/download/madcats-linux-x86_64.tar.gz

# Extract the archive
tar -xzf madcats-linux-x86_64.tar.gz

# Verify the installation
./madcats --help

# (Optional) Move to a directory in your PATH
mkdir -p ~/bin
mv madcats ~/bin/

# Add to PATH if not already present (add to ~/.bashrc or ~/.zshrc)
export PATH="$HOME/bin:$PATH"
```

### macOS

```bash
# For Apple Silicon (M1/M2/M3/M4)
curl -LO https://github.com/Breeding-Insight/MADCATS/releases/latest/download/madcats-macos-aarch64.tar.gz
tar -xzf madcats-macos-aarch64.tar.gz

# For Intel Macs
curl -LO https://github.com/Breeding-Insight/MADCATS/releases/latest/download/madcats-macos-x86_64.tar.gz
tar -xzf madcats-macos-x86_64.tar.gz

# Verify the installation
./madcats --help
```

> **Note for macOS users**: the binaries are not notarized. If macOS reports that the developer
> cannot be verified, remove the quarantine flag:
>
> ```bash
> xattr -d com.apple.quarantine madcats
> ```

The binary takes a couple of seconds to start (it unpacks its bundled interpreter on each run);
processing speed is the same as the pip-installed command.

### Windows

No Windows binary is published yet. Install with `pip install madcats` under Python 3.9+.

## From Source

```bash
git clone https://github.com/Breeding-Insight/MADCATS.git
cd MADCATS

# Install the madcats command from the checkout
pip install .

# ...or run the script directly, with no installation at all
python3 madcats.py --help
```

Every example in this documentation uses `madcats <command>`. When running from a checkout without
installing, substitute `python3 madcats.py <command>`.

## Verifying the Installation

```bash
madcats --version
```

```
madcats 0.14.0
```

To run the test suite from a source checkout:

```bash
python3 -m unittest discover -s tests -v
```

## Next Steps

- [Quick Start](quickstart.md): generate your first MADC
- [Try the Example](example.md): run a small simulated project in seconds
- [Choosing a Workflow](workflows.md): which commands your data needs
