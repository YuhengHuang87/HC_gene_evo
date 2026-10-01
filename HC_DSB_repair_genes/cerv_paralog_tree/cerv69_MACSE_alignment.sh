#!/bin/bash
# Codon-aware alignment of the 69 cerv copies (cerv69.fa, headers species__copy) with MACSE v2.03,
# input for the IQ-TREE runs in cerv69_NT_iqtree.sh.
# Run from the folder holding cerv69.fa:  bash cerv69_MACSE_alignment.sh
# Regenerated on 2026-10-01 from the output files: both steps reproduce cerv69_NT.fa, cerv69_AA.fa,
# cerv69_NT_clean.fa and cerv69_AA_clean.fa byte for byte.
set -euo pipefail

MACSE=/Applications/macse_v2.03.jar

# 1. align with default settings ('!' = frameshift, '*' = stop codon in the output)
java -jar $MACSE -prog alignSequences -seq cerv69.fa -out_NT cerv69_NT.fa -out_AA cerv69_AA.fa

# 2. export for IQ-TREE: final stop codons and codons with internal frameshifts become '---';
#    MACSE then drops the last codon column, which is gap-only after this
java -jar $MACSE -prog exportAlignment -align cerv69_NT.fa \
    -codonForFinalStop --- -codonForInternalFS --- \
    -out_NT cerv69_NT_clean.fa -out_AA cerv69_AA_clean.fa
