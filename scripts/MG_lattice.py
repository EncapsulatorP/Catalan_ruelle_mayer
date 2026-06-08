#!/usr/bin/env python3
"""Compatibility entrypoint for the historical MG lattice generator."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from catalan_mg.generate import generate_all_artifacts


if __name__ == "__main__":
    generate_all_artifacts(REPO_ROOT, layout="legacy", sections=["mg"])
