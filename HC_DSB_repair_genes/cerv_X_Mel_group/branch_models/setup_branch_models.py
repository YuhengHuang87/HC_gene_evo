#!/usr/bin/env python3
"""Write codeml inputs for cerv_X branch / branch-site / clade model C runs.

Foreground (#1) = D. melanogaster complex (mel, sim, sec, mau) tips, internal branches, and the stem
branch on which CG42299_X and qjt_3L arose, i.e. the lineage where cerv_X has 3-4 paralogs.
Background (#0) = lineages with 2 copies (cerv_X + one autosomal copy).
Three-ratio tree adds #2 = single-copy tips. Only D. ficusphila: tblastn finds a single cerv hit there,
while D. pseudotakahashii has a second, unannotated cerv36B_2L-like copy on scaffold84
(JAJJHV010000913.1:1099897-1100708), so it is scored as 2 copies and stays in the background.

Trees are unrooted (the rooted tree in species_num/ makes codeml fit 32 instead of 31 branches).
Run from cerv_X_Mel_group/:  python3 branch_models/setup_branch_models.py
"""
import os
import re
import shutil

WD = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # cerv_X_Mel_group/
OUT = os.path.join(WD, "branch_models")
CODON = os.path.join(WD, "input_files_unmask/M7_8/sa1/cerv_X.macse.codon")
TEMPLATE = os.path.join(WD, "codeml_update_M7_8.txt")
NSPECIES = 17

MEL_SUBGROUP = "(D_erecta,(D_teissieri_273_3,(D_santomea,D_yakuba)))"
TAKA_SUZ = "((D_biarmipes,(D_subpulchrella,D_suzukii)),(D_pseudotakahashii,D_takahashii))"
ELEG = "(D_rhopaloa,(D_gunungcola,D_elegans))"

TREES = {
    "unlabeled": f"(({TAKA_SUZ},({MEL_SUBGROUP},(((D_mauritiana,D_simulans),D_sechellia),D_melanogaster))),{ELEG},D_ficusphila);",
    "fg_melcomplex": f"(({TAKA_SUZ},({MEL_SUBGROUP},(((D_mauritiana #1,D_simulans #1) #1,D_sechellia #1) #1,D_melanogaster #1) #1)),{ELEG},D_ficusphila);",
    "three_ratio": f"(({TAKA_SUZ},({MEL_SUBGROUP},(((D_mauritiana #1,D_simulans #1) #1,D_sechellia #1) #1,D_melanogaster #1) #1)),{ELEG},D_ficusphila #2);",
}

# name: (tree, model, NSsites, fix_omega, omega)
RUNS = {
    "M0":          ("unlabeled",     0, 0,  0, 0.4),   # one omega for all branches
    "two_ratio":   ("fg_melcomplex", 2, 0,  0, 0.4),   # omega_bg, omega_fg(3-4 copies)
    "three_ratio": ("three_ratio",   2, 0,  0, 0.4),   # omega_2copies, omega_3-4copies, omega_1copy
    "free_ratio":  ("unlabeled",     1, 0,  0, 0.4),   # descriptive only: one omega per branch
    "BSA_null":    ("fg_melcomplex", 2, 2,  1, 1.0),   # branch-site model A, fg omega2 fixed = 1
    "BSA_alt":     ("fg_melcomplex", 2, 2,  0, 1.5),   # branch-site model A, fg omega2 estimated
    "M2a_rel":     ("unlabeled",     0, 22, 0, 0.4),   # null for clade model C
    "CmC":         ("fg_melcomplex", 3, 2,  0, 0.4),   # clade model C, divergent class differs fg vs bg
}


def check_tree(newick):
    tips = re.findall(r"D_[A-Za-z0-9_]+", newick)
    assert len(tips) == NSPECIES == len(set(tips)), tips
    assert newick.count("(") == newick.count(")")
    # unrooted: trifurcation at the top level
    depth, commas = 0, 0
    for ch in newick:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 1:
            commas += 1
    assert commas == 2, "tree is not unrooted"


def set_option(text, key, value):
    # replace the whole value up to the '*' comment (template has e.g. "NSsites = 7 8 * ...")
    new, n = re.subn(rf"(?m)^(\s*{key} = )[^*\n]*", rf"\g<1>{value} ", text)
    assert n == 1, key
    return new


template = open(TEMPLATE).read()
codon_species = {l.strip() for l in open(CODON) if l.startswith("D_")}
for name, newick in TREES.items():
    check_tree(newick)
    assert set(re.findall(r"D_[A-Za-z0-9_]+", newick)) == codon_species, name

for run, (tree, model, nssites, fix_omega, omega) in RUNS.items():
    d = os.path.join(OUT, run)
    os.makedirs(d, exist_ok=True)
    shutil.copy(CODON, os.path.join(d, "cerv_X.macse.codon"))
    with open(os.path.join(d, "cerv_X_tree.txt"), "w") as f:
        f.write(f"{NSPECIES} 1\n{TREES[tree]}\n")
    ctl = template
    for key, value in (("model", model), ("NSsites", nssites), ("fix_omega", fix_omega), ("omega", omega)):
        ctl = set_option(ctl, key, value)
    header = f"seqfile = cerv_X.macse.codon\noutfile = cerv_X_{run}_output\ntreefile = cerv_X_tree.txt\n"
    with open(os.path.join(d, "codeml.ctl"), "w") as f:
        f.write(header + ctl)
    print(f"{run:12s} tree={tree:14s} model={model} NSsites={nssites} fix_omega={fix_omega} omega={omega}")
