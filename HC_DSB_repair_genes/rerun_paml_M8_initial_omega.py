#!/usr/bin/env python3
"""Re-run codeml M8 from other initial omega values for genes where M8 did not converge.

M8 must fit at least as well as M7 (nested, same beta) and as M8a (same model with w fixed at 1).
When grab_paml_result_likelihood_macse.py reports "lnL M8 < M7" or "lnL M8 < M8a", M8 stopped at a
worse optimum - typically stuck with p1 near zero and a large w - and the fit has to be repeated from
other starting values.

This makes one folder per initial omega, reusing the alignment and tree already in the M7_8 sa folders:
    <out_dir>/M8_w<omega>/sa<N>/   same sa<N> numbering and gene as <input_dir>/M7_8/sa<N>
with codeml.ctl copied from the original M7_8 run but NSsites = 8 (M8 only, M7 is not repeated),
fix_omega = 0 and omega = <omega>, writing <gene>_M8_output_<tag>_w<omega>.

Then:
    sbatch MACSE_PAML_submit.sub "M8_w0.5 M8_w1.5 M8_w3" input_files_unmask_rerun
    sbatch MACSE_PAML_grab.sub input_files_unmask        # picks the best M8 across all runs

Usage:
    python3 rerun_paml_M8_initial_omega.py --input_dir <input_files_unmask> \\
        --result paml_likelihood_result_MACSE_input_files_unmask_Aniek_gene_list_2026.tsv
    python3 rerun_paml_M8_initial_omega.py --input_dir <...> --genes His2Av mei-41 --omegas 0.2 2
"""
import argparse
import csv
import os
import re
import shutil
import sys

FLAG = "rerun with other initial omega"


def read_ctl(ctl):
    values = {}
    with open(ctl) as f:
        text = f.read()
    for key in ("seqfile", "outfile", "treefile"):
        m = re.search(rf"^\s*{key}\s*=\s*(\S+)", text, re.M)
        if m:
            values[key] = m.group(1)
    return values, text


def make_ctl(text, outfile, omega):
    """Original M7_8 control file -> M8 only, started from the given omega."""
    text = re.sub(r"^\s*outfile\s*=.*$", f"outfile = {outfile}", text, flags=re.M)
    text = re.sub(r"^(\s*NSsites\s*=\s*)[\d\s]+", r"\g<1>8 ", text, flags=re.M)
    text = re.sub(r"^(\s*fix_omega\s*=\s*)\d+", r"\g<1>0", text, flags=re.M)
    text = re.sub(r"^(\s*omega\s*=\s*)\S+", rf"\g<1>{omega}", text, flags=re.M)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input_dir", required=True, help="folder with the M7_8/ runs to repeat")
    ap.add_argument("--result", help="result table from grab_paml_result_likelihood_macse.py; "
                                     "genes whose note asks for a re-run are selected")
    ap.add_argument("--genes", nargs="+", help="gene names to re-run (instead of --result)")
    ap.add_argument("--omegas", nargs="+", default=["0.5", "1.5", "3"], help="initial omega values")
    ap.add_argument("--out_dir", default=None, help="default: <input_dir>_rerun")
    ap.add_argument("--tag", default=None, help="outfile suffix (default: taken from the original run)")
    args = ap.parse_args()

    src_dir = os.path.join(args.input_dir, "M7_8")
    if not os.path.isdir(src_dir):
        sys.exit(f"Not found: {src_dir}")
    out_dir = args.out_dir or args.input_dir.rstrip("/") + "_rerun"

    genes = set(args.genes or [])
    if args.result:
        with open(args.result) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if FLAG in row.get("note", ""):
                    genes.add(row["gene"])
    if not genes:
        sys.exit("No genes selected (use --result and/or --genes)")

    # gene -> sa folder, from the map written at setup, else from each codeml.ctl
    sa_of = {}
    gene_map = os.path.join(src_dir, "sa_gene_map.tsv")
    if os.path.isfile(gene_map):
        with open(gene_map) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                sa_of[row["gene"]] = row["sa_dir"]
    else:
        for sa in os.listdir(src_dir):
            ctl = os.path.join(src_dir, sa, "codeml.ctl")
            if os.path.isfile(ctl):
                values, _ = read_ctl(ctl)
                gene = re.sub(r"_M8_output.*$", "", values.get("outfile", ""))
                if gene:
                    sa_of[gene] = sa

    missing = sorted(g for g in genes if g not in sa_of)
    todo = sorted(g for g in genes if g in sa_of)

    for omega in args.omegas:
        model_dir = os.path.join(out_dir, f"M8_w{omega}")
        os.makedirs(model_dir, exist_ok=True)
        with open(os.path.join(model_dir, "sa_gene_map.tsv"), "w") as gm:
            gm.write("sa_dir\tgene\tinitial_omega\toutfile\n")
            for gene in todo:
                sa = sa_of[gene]
                src = os.path.join(src_dir, sa)
                dst = os.path.join(model_dir, sa)
                os.makedirs(dst, exist_ok=True)
                values, text = read_ctl(os.path.join(src, "codeml.ctl"))
                tag = args.tag or re.sub(r"^.*_M8_output_?", "", values["outfile"]) or "macse"
                outfile = f"{gene}_M8_output_{tag}_w{omega}"
                for key in ("seqfile", "treefile"):
                    shutil.copyfile(os.path.join(src, values[key]), os.path.join(dst, values[key]))
                with open(os.path.join(dst, "codeml.ctl"), "w") as ctl:
                    ctl.write(make_ctl(text, outfile, omega))
                gm.write(f"{sa}\t{gene}\t{omega}\t{outfile}\n")
        print(f"[M8_w{omega}] {len(todo)} genes -> {model_dir}")

    if missing:
        print(f"\nNot found in {src_dir}: {', '.join(missing)}", file=sys.stderr)
    print("\nNext: sbatch MACSE_PAML_submit.sub \"" +
          " ".join(f"M8_w{w}" for w in args.omegas) + f"\" {os.path.basename(out_dir)}")


if __name__ == "__main__":
    main()
