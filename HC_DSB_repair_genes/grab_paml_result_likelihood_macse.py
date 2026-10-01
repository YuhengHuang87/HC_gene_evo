#!/usr/bin/env python3
"""Collect M7 / M8 / M8a likelihoods for each gene in a gene list (MACSE PAML runs).

Python version of grab_paml_result_likelihood_update.pl for the folders made by
macse_paml_input_each_gene.py:
    <input_dir>/M7_8/sa<N>/<gene>_M8_output_macse          (NSsites = 7 8)
    <input_dir>/M8a/sa<N>/<gene>_fix_omega_output_macse    (NSsites = 8, fix_omega = 1, omega = 1)

The gene in each sa<N> folder is read from the outfile line of its codeml.ctl, so gene names
containing '_' (e.g. Su_var_2-10) are handled.

Output: one tab-separated row per gene, in gene-list order. Genes with no finished run are marked
"check" in the note column. Columns:
    gene, n_species, n_codons,
    lnL_M7, np_M7, lnL_M8, np_M8, M8_p0, M8_p, M8_q, M8_p1, M8_w, M8_BEB_sites_P95, M8_BEB_sites_P99,
    lnL_M8a, np_M8a,
    LRT_M7_M8 (2*dlnL), p_M7_M8 (chi2, df = 2),
    LRT_M8a_M8 (2*dlnL), p_M8a_M8 (chi2, df = 1), p_M8a_M8_mix (50:50 mixture of 0 and chi2 df = 1),
    note

Usage:
    python3 grab_paml_result_likelihood_macse.py
    python3 grab_paml_result_likelihood_macse.py --input_dir <input_files_ms> --gene_list <list> --out <file>
"""
import argparse
import glob
import math
import os
import re
import sys

BASE = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species"

RUNS = {"M7_8": "M8_output", "M8a": "fix_omega_output"}

RE_MODEL = re.compile(r"^Model (\d+):")
RE_LNL = re.compile(r"^lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s+(-?\d+\.\d+)")
RE_NSLS = re.compile(r"^ns\s*=\s*(\d+)\s+ls\s*=\s*(\d+)")
RE_M8_P0 = re.compile(r"p0\s*=\s*([\d.]+)\s+p\s*=\s*([\d.]+)\s+q\s*=\s*([\d.]+)")
RE_M8_P1 = re.compile(r"\(p1\s*=\s*([\d.]+)\)\s+w\s*=\s*([\d.]+)")
RE_BEB_SITE = re.compile(r"^\s*\d+\s+\S\s+[\d.]+(\*{1,2})")


def parse_codeml(path, default_model=None):
    """Return {'ns', 'ls', 'models': {model_no: {...}}} from a codeml main output file.

    codeml writes a "Model N:" header only when NSsites lists several models (e.g. "7 8").
    With a single NSsites value there is no header, so results go to default_model.
    """
    res = {"ns": None, "ls": None, "models": {}}
    model = default_model
    if model is not None:
        res["models"][model] = {"beb95": 0, "beb99": 0}
    in_beb = False
    with open(path) as f:
        for line in f:
            m = RE_NSLS.match(line)
            if m and res["ns"] is None:
                res["ns"], res["ls"] = int(m.group(1)), int(m.group(2))
                continue
            m = RE_MODEL.match(line)
            if m:
                model = int(m.group(1))
                res["models"][model] = {"beb95": 0, "beb99": 0}
                in_beb = False
                continue
            if model is None:
                continue
            cur = res["models"][model]
            m = RE_LNL.match(line)
            if m:
                cur["np"], cur["lnL"] = int(m.group(1)), float(m.group(2))
                continue
            m = RE_M8_P0.search(line)
            if m:
                cur["p0"], cur["p"], cur["q"] = m.groups()
                continue
            m = RE_M8_P1.search(line)
            if m:
                cur["p1"], cur["w"] = m.groups()
                continue
            if line.startswith("Bayes Empirical Bayes"):
                in_beb = True
                continue
            if in_beb:
                if line.startswith("The grid"):
                    in_beb = False
                    continue
                m = RE_BEB_SITE.match(line)
                if m:
                    cur["beb95"] += 1
                    cur["beb99"] += m.group(1) == "**"
    return res


def read_ctl(ctl, tag):
    """Return (gene, outfile, single NSsites model or None) from codeml.ctl."""
    gene = out = single_model = None
    with open(ctl) as f:
        for line in f:
            line = line.split("*", 1)[0]
            m = re.match(r"\s*outfile\s*=\s*(\S+)", line)
            if m:
                out = m.group(1)
                gene = re.sub(r"_" + tag + r".*$", "", out)
            m = re.match(r"\s*NSsites\s*=\s*([\d\s]+)", line)
            if m:
                models = m.group(1).split()
                single_model = int(models[0]) if len(models) == 1 else None
    return gene, out, single_model


def collect(run_dir, tag):
    """Map gene -> (sa folder, parsed output or None)."""
    found = {}
    for sa in sorted(glob.glob(os.path.join(run_dir, "sa*"))):
        ctl = os.path.join(sa, "codeml.ctl")
        if not os.path.isfile(ctl):
            continue
        gene, out, single_model = read_ctl(ctl, tag)
        if gene is None or out is None:
            continue
        out_path = os.path.join(sa, out)
        parsed = parse_codeml(out_path, single_model) if os.path.isfile(out_path) else None
        found[gene] = (os.path.basename(sa), parsed)
    return found


def chi2_sf(x, df):
    """Upper tail of chi-square for df 1 or 2."""
    if x <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(x / 2))
    if df == 2:
        return math.exp(-x / 2)
    raise ValueError(df)


def fmt(v, nd=6):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.{nd}g}" if abs(v) < 1e-3 and v != 0 else f"{v:.{nd}f}".rstrip("0").rstrip(".")
    return str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input_dir", default=os.path.join(BASE, "MACSE_PRANK_23species", "input_files_ms"),
                    help="folder holding M7_8/ and M8a/")
    ap.add_argument("--gene_list", default=os.path.join(BASE, "Aniek_gene_list_2026.txt"),
                    help="gene names in the first (tab-separated) column")
    ap.add_argument("--out", default="paml_likelihood_result_MACSE_Aniek_gene_list_2026.tsv")
    ap.add_argument("--m8_rerun_dirs", nargs="+", default=[],
                    help="extra folders of M8-only runs (e.g. <input_dir>_rerun/M8_w1.5 from "
                         "rerun_paml_M8_initial_omega.py); the best M8 lnL over all runs is reported")
    args = ap.parse_args()

    runs = {}
    for run, tag in RUNS.items():
        run_dir = os.path.join(args.input_dir, run)
        if not os.path.isdir(run_dir):
            print(f"Warning: {run_dir} not found; {run} columns will be NA", file=sys.stderr)
            runs[run] = {}
        else:
            runs[run] = collect(run_dir, tag)

    reruns = {}
    for d in args.m8_rerun_dirs:
        d = d if os.path.isabs(d) else os.path.join(args.input_dir, d)
        if os.path.isdir(d):
            reruns[os.path.basename(d.rstrip("/"))] = collect(d, "M8_output")
        else:
            print(f"Warning: rerun folder {d} not found", file=sys.stderr)

    with open(args.gene_list) as f:
        genes = [line.rstrip("\n").split("\t")[0].strip() for line in f]
    genes = [g for g in genes if g]

    header = ["gene", "n_species", "n_codons",
              "lnL_M7", "np_M7", "lnL_M8", "np_M8", "M8_p0", "M8_p", "M8_q", "M8_p1", "M8_w",
              "M8_BEB_sites_P95", "M8_BEB_sites_P99", "lnL_M8a", "np_M8a",
              "LRT_M7_M8", "p_M7_M8", "LRT_M8a_M8", "p_M8a_M8", "p_M8a_M8_mix", "M8_run", "note"]
    n_done = 0
    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for gene in genes:
            notes = []
            m78 = runs["M7_8"].get(gene)
            m8a = runs["M8a"].get(gene)
            r78 = m78[1] if m78 else None
            r8a = m8a[1] if m8a else None
            if m78 is None:
                notes.append("no M7_8 folder")
            elif r78 is None:
                notes.append(f"M7_8 {m78[0]}: no output file")
            if m8a is None:
                notes.append("no M8a folder")
            elif r8a is None:
                notes.append(f"M8a {m8a[0]}: no output file")

            M7 = (r78 or {}).get("models", {}).get(7, {})
            M8 = (r78 or {}).get("models", {}).get(8, {})
            M8a = (r8a or {}).get("models", {}).get(8, {})

            # best M8 over the original run and any re-runs from other initial omega values
            m8_run = "M7_8" if "lnL" in M8 else None
            for name, found in reruns.items():
                hit = found.get(gene)
                parsed = hit[1] if hit else None
                if parsed is None:
                    continue
                cand = parsed["models"].get(8, {})
                if "lnL" not in cand:
                    notes.append(f"{name} {hit[0]}: unfinished")
                elif "lnL" not in M8 or cand["lnL"] > M8["lnL"]:
                    M8, m8_run = cand, name
            for name, mod, parsed in (("M7", M7, r78), ("M8", M8, r78), ("M8a", M8a, r8a)):
                if parsed is not None and "lnL" not in mod:
                    notes.append(f"{name} unfinished")

            ns = (r78 or r8a or {}).get("ns")
            ls = (r78 or r8a or {}).get("ls")

            lrt78 = p78 = lrt8a = p8a = p8a_mix = None
            if "lnL" in M7 and "lnL" in M8:
                lrt78 = 2 * (M8["lnL"] - M7["lnL"])
                p78 = chi2_sf(lrt78, 2)
                if lrt78 < -1e-3:
                    notes.append("lnL M8 < M7: rerun with other initial omega")
            if "lnL" in M8a and "lnL" in M8:
                lrt8a = 2 * (M8["lnL"] - M8a["lnL"])
                p8a = chi2_sf(lrt8a, 1)
                p8a_mix = 0.5 * p8a if lrt8a > 0 else 1.0
                if lrt8a < -1e-3:
                    notes.append("lnL M8 < M8a: rerun with other initial omega")

            row = [gene, ns, ls,
                   M7.get("lnL"), M7.get("np"), M8.get("lnL"), M8.get("np"),
                   M8.get("p0"), M8.get("p"), M8.get("q"), M8.get("p1"), M8.get("w"),
                   M8.get("beb95") if "lnL" in M8 else None, M8.get("beb99") if "lnL" in M8 else None,
                   M8a.get("lnL"), M8a.get("np"),
                   lrt78, p78, lrt8a, p8a, p8a_mix, m8_run]
            if lrt78 is not None and lrt8a is not None and not notes:
                n_done += 1
            note = "check: " + "; ".join(notes) if notes else ""
            out.write("\t".join(fmt(v) for v in row) + "\t" + note + "\n")

    print(f"{n_done}/{len(genes)} genes with complete M7/M8/M8a results -> {args.out}")


if __name__ == "__main__":
    main()
