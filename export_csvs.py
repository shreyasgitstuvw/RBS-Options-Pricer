"""
export_csvs.py — Generate static CSV snapshots from rel_bs.exe for all 6 modes.

Usage:
    python export_csvs.py

Writes to data/ directory:
    data/iv_surface.csv
    data/density.csv
    data/greeks.csv
    data/american.csv
    data/barrier.csv
    data/calibration.csv

Use cases:
  1. Seed Lovable's design session with real data schemas
  2. Generate static assets for GitHub Pages deployment
"""

import subprocess, os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(SCRIPT_DIR, "rel_bs.exe")
if not os.path.exists(EXE):
    EXE = os.path.join(SCRIPT_DIR, "rel_bs")

MODES = ["iv_surface", "density", "greeks", "american", "barrier", "calibration"]
OUT_DIR = os.path.join(SCRIPT_DIR, "data")

os.makedirs(OUT_DIR, exist_ok=True)

for mode in MODES:
    print(f"  exporting {mode} ...", end=" ", flush=True)
    result = subprocess.run([EXE, "--csv", mode], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"FAILED\n{result.stderr}")
        sys.exit(1)
    out_path = os.path.join(OUT_DIR, f"{mode}.csv")
    with open(out_path, "w", newline="") as f:
        f.write(result.stdout)
    lines = result.stdout.strip().count("\n")  # rows excluding header
    print(f"{lines} rows -> {out_path}")

print("\nAll CSVs written to data/")
