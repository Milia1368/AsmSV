#!/usr/bin/env python3
"""Source-checkout launcher for the CPX exporter."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from asmsv_next.cpx_export import main

if __name__ == '__main__':
    raise SystemExit(main())
