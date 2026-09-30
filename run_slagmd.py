#!/usr/bin/env python3
"""Launch the slagmd window.

    python3 run_slagmd.py                       open with no project
    python3 run_slagmd.py /path/to/project      open a project folder straight away
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from slagmd.gui import main                      # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
