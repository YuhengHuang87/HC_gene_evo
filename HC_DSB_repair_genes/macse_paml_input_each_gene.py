#!/usr/bin/env python3
"""Set up codeml input folders for each gene from the OMM_MACSE (PRANK) alignments.

Python version of aa_paml_input_each_gene_species_update.pl, for MACSE output.

For each gene in MACSE_PRANK_23species/<gene>/<gene>_final_align_NT.aln this writes
    <out_dir>/<model>/sa<N>/<gene>.macse.codon   (PHYLIP, same layout as the clustalo .codon files)
    <out_dir>/<model>/sa<N>/<gene>_tree.txt      (copied from species_num/)
    <out_dir>/<model>/sa<N>/codeml.ctl           (seqfile/outfile/treefile + template)
plus <out_dir>/<model>/sa_gene_map.tsv linking sa<N> to gene names.

Models:
    M7_8 -> codeml_update_M7_8.txt        outfile <gene>_M8_output_<tag>
    M8a  -> codeml_update_fix_omerga.txt  outfile <gene>_fix_omega_output_<tag>

Before writing, each alignment is checked against its tree (same species, same count) and for
codeml problems (length not a multiple of 3, stop codons, MACSE '!'/'*', invalid characters).
Genes that fail are skipped and reported.

--export_frameshifts converts a MACSE alignment that was not exported (e.g. *_final_unmask_align_NT.aln)
the same way OMM_MACSE makes *_final_align_NT.aln: every codon containing a frameshift '!' and every
stop codon (internal or terminal) becomes NNN. Applied to *_final_mask_align_NT.aln this reproduces
*_final_align_NT.aln exactly. The number of converted codons per gene goes to sa_gene_map.tsv.

Usage (submit through MACSE_PAML_setup.sub on HPC):
    python3 macse_paml_input_each_gene.py --template_dir /path/to/templates
    # unmasked alignments, kept apart from the masked runs:
    python3 macse_paml_input_each_gene.py --template_dir /path/to/templates \\
        --aln_suffix _final_unmask_align_NT.aln --export_frameshifts \\
        --out_dir <macse_dir>/input_files_unmask --tag macse_unmask
"""
import argparse
import os
import re
import sys

BASE = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/MACSE_PRANK_23species"

MODELS = {
    "M7_8": ("codeml_update_M7_8.txt", "M8_output"),
    "M8a": ("codeml_update_fix_omerga.txt", "fix_omega_output"),
}
STOPS = {"TAA", "TAG", "TGA"}
VALID = set("ACGTUNRYSWKMBDHV-?")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in seqs:
                    raise ValueError(f"duplicate sequence name {name}")
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def export_frameshifts(seqs):
    """MACSE-style export: codons with '!' and stop codons -> NNN. Returns (seqs, n_fs, n_stop)."""
    exported, n_fs, n_stop = {}, 0, 0
    for name, s in seqs.items():
        codons = []
        for i in range(0, len(s), 3):
            c = s[i:i + 3]
            if "!" in c:
                n_fs += 1
                c = "NNN"
            elif c in STOPS:
                n_stop += 1
                c = "NNN"
            codons.append(c)
        exported[name] = "".join(codons)
    return exported, n_fs, n_stop


def read_tree(path):
    with open(path) as f:
        header = f.readline().split()
        newick = f.read()
    tips = re.findall(r"[(,]\s*([^(),:;\s]+)", newick)
    return int(header[0]), tips


def check_gene(seqs, n_tree, tips):
    problems = []
    names = set(seqs)
    if len(tips) != len(set(tips)):
        problems.append("duplicate tip names in tree")
    if n_tree != len(tips):
        problems.append(f"tree header says {n_tree} species but tree has {len(tips)} tips")
    if names != set(tips):
        if names - set(tips):
            problems.append(f"in alignment but not tree: {sorted(names - set(tips))}")
        if set(tips) - names:
            problems.append(f"in tree but not alignment: {sorted(set(tips) - names)}")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        problems.append(f"unequal sequence lengths {sorted(lengths)}")
    elif lengths.pop() % 3:
        problems.append("alignment length is not a multiple of 3")
    for name, s in seqs.items():
        bad = set(s) - VALID
        if bad:
            problems.append(f"{name}: invalid characters {''.join(sorted(bad))} (MACSE '!'/'*' not exported?)")
        stops = [i // 3 + 1 for i in range(0, len(s) - 2, 3) if s[i:i + 3] in STOPS]
        if stops:
            problems.append(f"{name}: stop codon(s) at codon {stops}")
    return problems


def write_phylip(seqs, order, path, width=60):
    length = len(next(iter(seqs.values())))
    with open(path, "w") as out:
        out.write(f"  {len(seqs)}   {length}\n")
        for name in order:
            out.write(f"{name}\n")
            s = seqs[name]
            for i in range(0, len(s), width):
                out.write(s[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--macse_dir", default=BASE, help="folder with one MACSE output folder per gene")
    ap.add_argument("--tree_dir", default=None, help="folder with <gene>_tree.txt (default: <macse_dir>/species_num)")
    ap.add_argument("--out_dir", default=None, help="output root (default: <macse_dir>/input_files_ms)")
    ap.add_argument("--template_dir", default=".", help="folder with the codeml_update_*.txt templates")
    ap.add_argument("--models", nargs="+", choices=sorted(MODELS), default=sorted(MODELS))
    ap.add_argument("--aln_suffix", default="_final_align_NT.aln")
    ap.add_argument("--export_frameshifts", action="store_true",
                    help="convert codons with '!' and stop codons to NNN (for non-exported MACSE alignments)")
    ap.add_argument("--tag", default="macse", help="suffix of codeml outfile names")
    args = ap.parse_args()

    tree_dir = args.tree_dir or os.path.join(args.macse_dir, "species_num")
    out_dir = args.out_dir or os.path.join(args.macse_dir, "input_files_ms")

    templates = {}
    for model in args.models:
        tpath = os.path.join(args.template_dir, MODELS[model][0])
        if not os.path.isfile(tpath):
            sys.exit(f"Template not found: {tpath} (set --template_dir)")
        with open(tpath) as f:
            templates[model] = f.read()

    genes = sorted(
        g for g in os.listdir(args.macse_dir)
        if os.path.isfile(os.path.join(args.macse_dir, g, g + args.aln_suffix))
    )
    if not genes:
        sys.exit(f"No *{args.aln_suffix} alignments found under {args.macse_dir}")

    ready, skipped = [], []
    for gene in genes:
        tree_file = os.path.join(tree_dir, f"{gene}_tree.txt")
        if not os.path.isfile(tree_file):
            skipped.append((gene, ["no tree file " + tree_file]))
            continue
        seqs = read_fasta(os.path.join(args.macse_dir, gene, gene + args.aln_suffix))
        n_fs = n_stop = 0
        if args.export_frameshifts:
            seqs, n_fs, n_stop = export_frameshifts(seqs)
        n_tree, tips = read_tree(tree_file)
        problems = check_gene(seqs, n_tree, tips)
        if problems:
            skipped.append((gene, problems))
        else:
            ready.append((gene, seqs, tips, tree_file, n_fs, n_stop))

    for model in args.models:
        model_dir = os.path.join(out_dir, model)
        os.makedirs(model_dir, exist_ok=True)
        outfile_tag = f"{MODELS[model][1]}_{args.tag}"
        with open(os.path.join(model_dir, "sa_gene_map.tsv"), "w") as gene_map:
            gene_map.write("sa_dir\tgene\tn_species\tn_codons\tn_frameshift_codons_to_NNN\tn_stop_codons_to_NNN\n")
            for num, (gene, seqs, tips, tree_file, n_fs, n_stop) in enumerate(ready, start=1):
                sa = os.path.join(model_dir, f"sa{num}")
                os.makedirs(sa, exist_ok=True)
                seqfile = f"{gene}.macse.codon"
                treefile = f"{gene}_tree.txt"

                write_phylip(seqs, tips, os.path.join(sa, seqfile))
                with open(tree_file) as src, open(os.path.join(sa, treefile), "w") as dst:
                    dst.write(src.read())
                with open(os.path.join(sa, "codeml.ctl"), "w") as ctl:
                    ctl.write(f"seqfile = {seqfile}\n")
                    ctl.write(f"outfile = {gene}_{outfile_tag}\n")
                    ctl.write(f"treefile = {treefile}\n")
                    ctl.write(templates[model])

                n_codons = len(next(iter(seqs.values()))) // 3
                gene_map.write(f"sa{num}\t{gene}\t{len(seqs)}\t{n_codons}\t{n_fs}\t{n_stop}\n")
        print(f"[{model}] wrote {len(ready)} genes to {model_dir}/sa1..sa{len(ready)}")

    if skipped:
        print(f"\nSkipped {len(skipped)} gene(s):", file=sys.stderr)
        for gene, problems in skipped:
            print(f"  {gene}", file=sys.stderr)
            for p in problems:
                print(f"    - {p}", file=sys.stderr)


if __name__ == "__main__":
    main()
