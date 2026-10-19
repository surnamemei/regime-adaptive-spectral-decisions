"""Check release integrity and stored simulation counts without new draws."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
for line in (ROOT / "SHA256SUMS.txt").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
    assert actual == expected, relative

study = ROOT / "noise_stochastic_transfer"
original = json.loads((study / "adaptive_selection/results/raw.json").read_text())
assert len(original["configs"]) == 196
assert sum(len(c["records"]) for c in original["configs"]) == 196000
meta = json.loads((study / "targeted_validation/results/fresh_meta.json").read_text())
assert len(meta["cells"]) == 72
with np.load(study / "targeted_validation/results/fresh_records.npz") as records:
    assert len(records.files) == 72
    assert sum(records[k].shape[0] for k in records.files) == 36000
    assert all(records[k].shape == (500, 9) for k in records.files)
    assert all(np.isfinite(records[k]).all() for k in records.files)
with (study / "targeted_validation/results/fresh_cells.csv").open() as stream:
    rows = [r for r in csv.DictReader(stream) if r["rule"] == "S"]
assert len(rows) == 72
sums = {k: sum(int(r[k]) for r in rows) for k in ["n", "correct", "wrong", "abstain", "failures"]}
assert sums == {"n": 36000, "correct": 14149, "wrong": 741, "abstain": 21110, "failures": 0}, sums
print(json.dumps({"manifest": "PASS", "original_cells": 196, "original_draws": 196000,
                  "fresh_cells": 72, "fresh_totals": sums}, indent=2))
