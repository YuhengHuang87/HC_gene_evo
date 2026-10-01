#!/usr/bin/env python3
"""Collect lnL / omega from cerv_X branch-model runs and do the LRTs.
Run from cerv_X_Mel_group/:  python3 branch_models/grab_branch_models.py > branch_models/cerv_X_branch_models_LRT.txt
"""
import math
import os
import re

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "branch_models")


def chi2_sf(x, df):
    if x <= 0:
        return 1.0
    if df % 2 == 0:
        k = df // 2
        return math.exp(-x / 2) * sum((x / 2) ** i / math.factorial(i) for i in range(k))
    # odd df: erfc term plus series
    s = math.erfc(math.sqrt(x / 2))
    term = math.sqrt(2 * x / math.pi) * math.exp(-x / 2)
    for i in range(1, (df - 1) // 2 + 1):
        s += term
        term *= x / (2 * i + 1)
    return s


def read(run):
    path = os.path.join(OUT, run, f"cerv_X_{run}_output")
    if not os.path.exists(path):
        return None
    text = open(path).read()
    m = re.findall(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s+(-?\d+\.\d+)", text)
    if not m:
        return None
    np_, lnl = m[-1]
    res = {"np": int(np_), "lnL": float(lnl), "text": text}
    w = re.search(r"w \(dN/dS\) for branches:\s+([\d.\s]+)", text)
    if w:
        res["w_branches"] = w.group(1).split()
    return res


runs = ["M0", "two_ratio", "three_ratio", "free_ratio", "BSA_null", "BSA_alt", "M2a_rel", "CmC"]
R = {r: read(r) for r in runs}

print("run\tnp\tlnL\tomega")
for r in runs:
    if R[r] is None:
        print(f"{r}\tNA\tNA\tnot finished")
        continue
    print(f"{r}\t{R[r]['np']}\t{R[r]['lnL']}\t{' '.join(R[r].get('w_branches', []))}")

for r in ("BSA_alt", "CmC"):
    if R[r]:
        block = re.search(r"MLEs of dN/dS \(w\) for site classes.*?(?=\n\n\n|\nnaive|\nBayes|\Z)", R[r]["text"], re.S)
        if block:
            print(f"\n## {r}\n{block.group(0).strip()}")

tests = [
    ("M0 vs two_ratio (omega differs on mel complex)", "M0", "two_ratio", False),
    ("two_ratio vs three_ratio (single-copy tips differ)", "two_ratio", "three_ratio", False),
    ("M0 vs free_ratio (descriptive)", "M0", "free_ratio", False),
    ("BSA_null vs BSA_alt (positive selection on mel complex; 50:50 mixture)", "BSA_null", "BSA_alt", True),
    ("M2a_rel vs CmC (divergent class differs on mel complex)", "M2a_rel", "CmC", False),
]
print("\ntest\t2dlnL\tdf\tp")
for label, null, alt, mixture in tests:
    if not (R[null] and R[alt]):
        print(f"{label}\tNA\tNA\tNA")
        continue
    lrt = 2 * (R[alt]["lnL"] - R[null]["lnL"])
    df = R[alt]["np"] - R[null]["np"]
    p = chi2_sf(lrt, df) / (2 if mixture else 1)
    flag = "  <- alt lnL below null: rerun alt from other starting omega" if lrt < 0 else ""
    print(f"{label}\t{lrt:.4f}\t{df}\t{p:.4g}{flag}")
