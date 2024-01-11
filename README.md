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
6. **Small internal deletions**: identified using ProSeq-IT;<br>
   to be considered intact, they must have a start codon, no premature stop codons, no frameshifts, and no deletion >5% of ORF<br>
     - \> 50 nt deleted in *gag*
     - \> 50 nt deleted in *pol*
     - \> 30 nt deleted in *vif*
     - \> 15 nt deleted in *vpr*
     - \> 11 nt deleted in *tat* exon 1 only, or \> 15 nt deleted in complete *tat*
     - \> 18 nt deleted in *rev*
     - \> 12 nt deleted in *vpu*
     - \> 31 nt deleted in *nef* (a defective *nef* is not considered a defect for the integrity)
     - \> 100 nt deleted in *env*
8. **Intact**


# Running the algorithm

The algorithm requires **3 files** to run:
1. **Excel template that contains 3 tabs**
   - _ProseqIT_criteria_: used for the analysis of ProSeq-IT
   - _Manual_assessment_: contains the results of the manual assessment with Geneious (to identify inversions)
   - _Hyperlinks_: contains the hyperlinks of HIV Database QCTool and Gene Cutter's results
2. Text file containing the **summary of QCTool's results** (_downloaded directly from QCTool's output_)
3. Excel file containing **ProSeq-IT's results** (_downloaded directly from ProSeq-IT's output_)

The Excel template can be found in this repo.

# Installation (from R or RStudio)

Install `IntegrityAlgorithm` from Github

```
install.packages("devtools")
install_github("alemi055/IntegrityAlgorithm")
```

# Usage (within R)

### HIV_IntegrityAnalysis()

To analyze the intactness of HIV proviruses, run the function `HIV_IntegrityAnalysis()` <br>
1. The algorithm will first analyse the results from QCTool to identify ORFs containing stop codons and hypermutated sequences.
2. It will then analyse the results from Gene Cutter to obtain the start and stop codons, if applicable, of and within each ORF.
3. Next, it will analyse the results from ProSeq-IT for large and small internal deletions.
4. Finally, all results from QCTool, Gene Cutter, ProSeq-IT, and the manual assessment with Geneious (through the Excel template) will be combined to assess the intactness of the HIV proviruses.

```
HIV_IntegrityAnalysis(template_filename, QCTool_summary, ProseqIT_rx, ProseqIT_RefSeq = FALSE, RefSeq = TRUE, analyses = 4)
```

| Argument | Required/Optional | Note |
| --- | --- | --- |
| `template_filename` | Required | Excel template containing the inversions manual assessment and the hyperlinks for QCTool and Gene Cutter |
|`QCTool_summary` | Required | Text file containing the results from QCTool |
| `ProseqIT_rx` | Required | Excel file containing the results from ProSeq-IT |
| `ProseqIT_RefSeq` | Optional; default if `FALSE` | Logical. If TRUE, the reference sequence (HXB2) is included in ProSeq-IT's results |
| `RefSeq` | Optional argument; default is `TRUE` | Logical. If TRUE, the reference sequence (HXB2) is included in QCTool and GeneCutter's results |
| `analyses` | optional argument; default is `4` | Specifies the analyses to be done. <br> 1: QCTool only, 2: Gene Cutter and ProSeq-IT only, 3: IntegrateInfo only, 4: All |

### Clonality_Analysis()

To analyze the clonality of the HIV proviruses, run the function `Split_files()`, then `Clonality_Analysis()` <br>
1. Files in the FASTA file input will be splitted according to the list of donors supplied. Donors with single sequences will not be analyzed.
2. The clonality will be assessed on each of the splitted files. A list of clones (0 different nt between two sequences) and potential clones, established by a threshold, will be output.

```
Split_files(FASTA_file, donors)
Clonality_Analysis(threshold)
```
| Argument | Required/Optional | Note |
| --- | --- | --- |
| `FASTA_file` | Required | FASTA file containing all sequences, including the reference sequence (HXB2) |
| `donors` | Required | Str vector containing the list of donors. If more than one, use the `c()` function |
| `threshold` | Optional; default is `5` | Threshold number of different nucleotides to consider two sequences as "potential clones" |

## Output

**tmp folder**
| Output file | Note |
| --- | --- |
| `Analyzed_GeneCutter.csv` | Analyzed results from Gene Cutter. <br> Details, for each sequence, the list of ORFs that have start and premature stop codons |
| `Analyzed_ProseqIT.csv` | Analyzed results from ProSeq-IT. <br> Details, for each sequence, the sequence length, the presence (binary: 0 [absence] and 1 [presence]) of large internal deletions, $\Psi$ mutations, and small internal deletions in each of the ORFs. For the ORFs, the list of detailed defects is also included |
| `Analyzed_QCTool.csv` | Analyzed results from QCTool. <br> Details, for each sequence, the number of stop codons, the list of stop codons, the number of incomplete codons, and the presence or absence of hypermutations |
| `intactness_detailedsummary.csv` | Detailed summary of defects for all sequences. <br> Details, for each sequence, the inferred intactness, the number of "main defects" (see list 1-7 above), the "main defect", and all defects in the sequence, including in each of the ORFs. Note that this summary **does not** hierarchize the defects |

**FINAL_OUTPUT folder**
| Output file | Note |
| --- | --- |
| `intactness_summary.csv` | Summary of **hierarchized** defects for all sequences |
| `[donor]_ClonalityAnalysis.csv` | Can be found in the **Clonality** subfolder. List of clones and potential clones for each of the donor's sequences |

# Citation

If you use the IntegrityAlgorithm, please cite this repository.

# Contact information
Audrée Lemieux <br>
lemieux.audree@gmail.com

