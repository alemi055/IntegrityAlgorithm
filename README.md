# IntegrityAlgorithm
Each provirus undergoes analysis to identify the presence of (i) inversions, (ii) hypermutations, (iii) large internal deletions, (iv) stop codons, (v) $\Psi$ defects, and (vi) small internal deletions in *gag*, *pol*, *vif*, *vpr*, *tat*, *rev*, *vpu*, *nef*, *env*, and RRE.

1. <span style="color:#7EB6D6;">**Inversions** (HEX #7EB6D6)</span>: manually identified in the alignment using Geneious
2. <span style="color:#989895;">**Hypermutations** (HEX #989895)</span>: identified with QCTool
3. <span style="color:#CA5F2E;">**Large internal deletions** (HEX #CA5F2E)</span>: size $<$ 8800 bp without primers, can be identified using ProSeq-IT
4. <span style="color:#E6AB48;">**Stop codons** (HEX #E6AB48)</span>: identified with GeneCutter (and QCTool); stop codons in Tat2 and/or Nef are not considered defects
5. <span style="color:#B0B2CB;">**$\Psi$ defects** (HEX #B0B2CB)</span>: manually identified in the alignment using Geneious
    a) Psi deletion
    b) SL2 deletion
    c) MSD point mutation
6. <span style="color:#941100;">**Small internal deletions** (HEX #941100)</span>: identified using ProSeq-IT with the values provided by Sarah Palmer (see [Appendix A, Criteria for defining defects][A. Criteria for defining defects])
    a) $>$  50 nt deleted in *gag*
    b) $>$ 50 nt deleted in *pol*
    c) $>$ 30 nt deleted in *vif*
    d) $>$ 15 nt deleted in *vpr*
    e) $>$ 11 nt deleted in *tat* exon 1 only, or $>$ 15 nt deleted in complete *tat*
    f) $>$ 18 nt deleted in *rev*
    g) $>$ 12 nt deleted in *vpu*
    h) $>$ 31 nt deleted in *nef* (a defective *nef* is not considered a defect for the integrity)
    i) $>$ 100 nt deleted in *env*
7. <span style="color:#76974C;">**Intact** (HEX #76974C)</span>
