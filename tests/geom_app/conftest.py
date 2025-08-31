#!/usr/bin/env python3
"""Pytest configuration for geom_app tests."""

import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))
