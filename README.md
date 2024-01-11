# IntegrityAlgorithm
Each provirus undergoes analysis to identify the presence of (i) inversions, (ii) hypermutations, (iii) large internal deletions, (iv) stop codons, (v) $\Psi$ defects, and (vi) small internal deletions in *gag*, *pol*, *vif*, *vpr*, *tat*, *rev*, *vpu*, *nef*, *env*, and RRE.

1. **Inversions**: manually identified in the alignment using Geneious
2. **Hypermutations**: identified with QCTool
3. **Large internal deletions**: size $<$ 8800 bp without primers, can be identified using ProSeq-IT
4. **Stop codons**: identified with GeneCutter (and QCTool); stop codons in Tat2 and/or Nef are not considered defects for the intactness
5. **$\Psi$ defects**: identified using ProSeq-IT and visually confirmed in the alignment using Geneious
    a) Psi deletion
    b) SL2 deletion
    c) MSD point mutation
6. **Small internal deletions**: identified using ProSeq-IT;
           - $>$  50 nt deleted in *gag*
           - $>$ 50 nt deleted in *pol*
           - $>$ 30 nt deleted in *vif*
           - $>$ 15 nt deleted in *vpr*
           - $>$ 11 nt deleted in *tat* exon 1 only, or $>$ 15 nt deleted in complete *tat*
           - $>$ 18 nt deleted in *rev*
           - $>$ 12 nt deleted in *vpu*
           - $>$ 31 nt deleted in *nef* (a defective *nef* is not considered a defect for the integrity)
           - $>$ 100 nt deleted in *env*
8. **Intact**
