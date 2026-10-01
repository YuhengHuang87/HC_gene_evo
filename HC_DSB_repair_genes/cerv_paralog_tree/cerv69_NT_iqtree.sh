#!/bin/bash
# IQ-TREE 3.1.3 trees of the 69 cerv copies from the nucleotide alignment made by cerv69_MACSE_alignment.sh.
# Run from the folder holding cerv69_NT_clean.fa, codonpos.nex and constraint_cervL4.tre:
#   bash cerv69_NT_iqtree.sh
# Regenerated on 2026-10-01 from the .iqtree reports and au_con.log of the original runs (IQ-TREE 3.1.3, bioconda).
# The seeds are the ones those runs recorded. Re-running gives the same trees; log-likelihoods and
# SH-aLRT/UFBoot values can differ slightly.
# Node labels in the .treefile are SH-aLRT/UFBoot support (read by circular_tree.R).
set -euo pipefail

IQTREE=iqtree3

# 1. one model for all sites (ModelFinder: TIM2+F+I+R3)
$IQTREE -s cerv69_NT_clean.fa -B 1000 -alrt 1000 -seed 669497 --prefix cerv69_NT

# 2. partitioned by codon position, edge-linked proportional branch lengths (-p);
#    ModelFinder: GTR+F+I+G4, TPM3u+F+I+G4, TPM3u+I+G4 for positions 1-3 (cerv69_NTpart.best_scheme.nex).
#    This is the tree in Fig. S6 (circular_tree.R)
$IQTREE -s cerv69_NT_clean.fa -p codonpos.nex -B 1000 -alrt 1000 -seed 674912 --prefix cerv69_NTpart

# 3. same partition models, best tree with cerv-L4 (eugracilis__cerv83B_3R) forced into the clade of the
#    melanogaster-subgroup copies (constraint_cervL4.tre); in the tree from step 2 it groups with ficusphila__cerv_X
$IQTREE -s cerv69_NT_clean.fa -p cerv69_NTpart.best_scheme.nex -g constraint_cervL4.tre -B 1000 -alrt 1000 -T 2 \
    -seed 583650 --prefix cerv69_NTcon

# 4. AU test (10000 RELL replicates): unconstrained tree (step 2) vs constrained tree (step 3)
cat cerv69_NTpart.treefile cerv69_NTcon.treefile > au_con.tre
$IQTREE -s cerv69_NT_clean.fa -p cerv69_NTpart.best_scheme.nex -z au_con.tre -n 0 -zb 10000 -au -T 2 \
    -seed 710330 --prefix au_con -redo
