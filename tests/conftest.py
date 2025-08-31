#!/usr/bin/env python3
"""Pytest configuration for b3p tests."""

import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Removed import of BuildApp as it's no longer available after refactor
# from b3p.cli.build_app import BuildApp
