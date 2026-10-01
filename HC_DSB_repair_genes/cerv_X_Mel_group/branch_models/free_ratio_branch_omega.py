#!/usr/bin/env python3
"""Per-branch dN/dS from a codeml free-ratio run (model = 1) as a table for ggplot2 / ggtree.

Writes
  <prefix>_branch_omega.tsv   one row per branch (the branch above node `label`)
  <prefix>_codeml_nodes.nwk   the tree with species tip labels, internal nodes labelled n<codeml id>
                              and branch lengths t, so ggtree(...) %<+% table joins on `label`

codeml numbers tips 1..s in the order of the sequence file, and internal nodes s+1.. with the root
first. Every parsed omega is checked against the "w ratios as labels for TreeView" tree in the same
output, so a wrong tip mapping stops the script instead of producing a mislabelled table.

Copy number per species comes from the headers of cerv-rename-header.fa (species|copy), which
includes the unannotated D_pseudotakahashii cerv36B copy found by tblastn on scaffold84.

Usage:
  python3 free_ratio_branch_omega.py <codeml_output> <codon_alignment> <copies.fa> <prefix> [fg_tips]
"""
import re
import sys

OUT_FILE, CODON, COPIES_FA, PREFIX = sys.argv[1:5]
# optional 5th argument: comma-separated foreground tips (default: the melanogaster complex)
FOREGROUND = set(sys.argv[5].split(",")) if len(sys.argv) > 5 else None

MEL_COMPLEX = FOREGROUND or {"D_melanogaster", "D_simulans", "D_sechellia", "D_mauritiana"}  # foreground (#1)
RENAME = {"D_teissieri": "D_teissieri_273_3"}                                   # fasta -> tree name
LOW_INFO_SDS = 2.0   # branches with fewer inferred synonymous substitutions than this: omega unstable


def tip_order(codon_path):
    """species in sequence-file order = codeml tip numbers 1..s"""
    lines = [l.strip() for l in open(codon_path) if l.strip()]
    n = int(lines[0].split()[0])
    names = [l for l in lines[1:] if l.startswith("D_")]
    assert len(names) == n, (len(names), n)
    return names


def copies_per_species(fa):
    copies = {}
    for line in open(fa):
        if line.startswith(">"):
            sp, copy = (x.strip() for x in line[1:].strip().split("|", 1))
            sp = sp if sp.startswith("D_") else "D_" + sp          # headers may be "species | copy"
            copies.setdefault(RENAME.get(sp, sp), []).append(copy)
    return copies


def parse_branch_table(text):
    block = text.split(" branch          t       N       S   dN/dS")[-1]
    rows = []
    for line in block.splitlines()[1:]:
        m = re.match(r"\s+(\d+)\.\.(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
        if m:
            p, c = int(m.group(1)), int(m.group(2))
            t, N, S, w, dN, dS, NdN, SdS = (float(x) for x in m.groups()[2:])
            rows.append(dict(parent=p, child=c, t=t, N=N, S=S, omega=w, dN=dN, dS=dS, NdN=NdN, SdS=SdS))
        elif rows and not line.strip():
            break
    assert rows, "no free-ratio branch table found (is this a model = 1 run?)"
    return rows


def treeview_omegas(text):
    """{frozenset(descendant tips): omega} from the 'w ratios as labels for TreeView' tree"""
    nwk = text.split("w ratios as labels for TreeView:")[1].strip().splitlines()[0]
    out, stack, i = {}, [], 0
    toks = re.findall(r"\(|\)|,|;|#[0-9.eE+-]+|D_[A-Za-z0-9_]+", nwk)
    last = None
    for tok in toks:
        if tok == "(":
            stack.append(set())
        elif tok == ")":
            last = frozenset(stack.pop())
            if stack:
                stack[-1] |= last
        elif tok.startswith("D_"):
            last = frozenset([tok])
            stack[-1].add(tok)
        elif tok.startswith("#"):
            out[last] = float(tok[1:])
    return out


text = open(OUT_FILE).read()
tips = tip_order(CODON)
ntip = len(tips)
branches = parse_branch_table(text)
children = {}
for b in branches:
    children.setdefault(b["parent"], []).append(b["child"])
root = ({b["parent"] for b in branches} - {b["child"] for b in branches}).pop()


def desc(node):
    if node <= ntip:
        return {tips[node - 1]}
    return set().union(*(desc(c) for c in children[node]))


def label(node):
    return tips[node - 1] if node <= ntip else f"n{node}"


copies = copies_per_species(COPIES_FA)
missing = [sp for sp in tips if sp not in copies]
assert not missing, f"no copy information for {missing}"


def copy_class(n):
    return "1" if n == 1 else "2" if n == 2 else "3-4"


tv = treeview_omegas(text)
rows = []
for b in branches:
    d = desc(b["child"])
    # the TreeView tree rounds differently from the table: compare at the table's 4 decimals
    assert abs(round(tv[frozenset(d)], 4) - b["omega"]) < 1.5e-4, (b, tv[frozenset(d)])
    ncopies = sorted({len(copies[sp]) for sp in d})
    classes = sorted({copy_class(len(copies[sp])) for sp in d})
    rows.append(dict(
        label=label(b["child"]),
        branch=f"{b['parent']}..{b['child']}",
        node_type="tip" if b["child"] <= ntip else "internal",
        descendants=";".join(sorted(d)),
        n_descendants=len(d),
        copy_number=str(ncopies[0]) if len(ncopies) == 1 else "-".join(map(str, ncopies)),
        copy_class=classes[0] if len(classes) == 1 else "mixed",
        copies=";".join(copies[next(iter(d))]) if len(d) == 1 else "",
        foreground=("two_copies" if d <= MEL_COMPLEX else "background"),
        t=b["t"], N=b["N"], S=b["S"], omega=b["omega"], dN=b["dN"], dS=b["dS"],
        NdN=b["NdN"], SdS=b["SdS"],
        low_info=b["SdS"] < LOW_INFO_SDS,
    ))

cols = ["label", "branch", "node_type", "descendants", "n_descendants", "copy_number", "copy_class",
        "copies", "foreground", "t", "N", "S", "omega", "dN", "dS", "NdN", "SdS", "low_info"]
with open(PREFIX + "_branch_omega.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join("TRUE" if r[c] is True else "FALSE" if r[c] is False else str(r[c])
                          for c in cols) + "\n")

blen = {b["child"]: b["t"] for b in branches}


def newick(node):
    if node <= ntip:
        return f"{label(node)}:{blen[node]}"
    inner = ",".join(newick(c) for c in children[node])
    return f"({inner}){label(node)}" + (f":{blen[node]}" if node in blen else "")


with open(PREFIX + "_codeml_nodes.nwk", "w") as f:
    f.write(newick(root) + ";\n")

print(f"{PREFIX}: {len(rows)} branches ({sum(r['node_type'] == 'tip' for r in rows)} tips), root n{root}; "
      f"{sum(r['low_info'] for r in rows)} low-info branches (S*dS < {LOW_INFO_SDS}); all omegas match the TreeView tree")
