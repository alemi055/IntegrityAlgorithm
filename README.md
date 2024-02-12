# IntegrityAlgorithm: a R package to analyze the intactness and clonality of near full-length HIV proviruses

**The algorithm automatizes the near full-length (NFL) pipeline developed by Dr. Caroline Dufour:**
- Sannier, G., Dubé, M., Dufour, C., ..., Fromentin, R., Chomont, N., & Kaufmann, D.E. [Combined single-cell transcriptional, translational, and genomic profiling reveals HIV-1 reservoir diversity](https://doi.org/10.1016/j.celrep.2021.109643). *Cell Rep* **36**, 109643 (2021).
- Dufour, C., ..., Fromentin, R., & Chomont, N. [Phenotypic characterization of single CD4+ T cells harboring genetically intact and inducible HIV genomes](https://doi.org/10.1038/s41467-023-36772-x). *Nat Commun* **14**, 1115 (2023).
- Dubé, M., Tastet, O., Dufour, C., ..., Fromentin, R., Chomont, N., & Kaufmann, D.E. [Spontaneous HIV expression during suppressive ART is associated with the magnitude and function of HIV-specific CD4+ and CD8+ T cells](https://doi.org/10.1016/j.chom.2023.08.006). *Cell Host Microbe* **13**, 1507-1522 (2023).
- Dufour, C., ..., & Chomont, N. [Near full-length HIV sequencing in multiple tissues collected postmortem reveals shared clonal expansions across distinct reservoirs during ART](https://doi.org/10.1016/j.celrep.2023.113053). *Cell Rep* **42**, 113053 (2023).

### Defects

Each provirus undergoes analysis to identify the presence of (i) inversions, (ii) hypermutations, (iii) large internal deletions, (iv) stop codons, (v) $\Psi$ defects, and (vi) small internal deletions in *gag*, *pol*, *vif*, *vpr*, *tat*, *rev*, *vpu*, *nef*, *env*, and RRE.

1. **Inversions**: manually identified in the alignment using Geneious
2. **Hypermutations**: identified with QCTool
3. **Large internal deletions**: size $<$ 8800 bp without primers, can be identified using ProSeq-IT
4. **Stop codons**: identified with GeneCutter; stop codons in Tat2 and/or Nef are not considered defects for the intactness
5. **$\Psi$ defects**: identified using ProSeq-IT and visually confirmed in the alignment using Geneious <br>
    - Psi deletion <br>
    - SL2 deletion <br>
    - MSD point mutation <br>
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


# Before running the algorithm

To analyze the HIV proviral sequences, you should have a FASTA file containing sequences:
- where inversions have been corrected (to do so, replace the inversions with their reverse-complements)
- where primers have been removed

**1. Align the inversions-free, primers-free, sequences with the reference sequence (HXB2) using MAFFT (E-INS-i method for more accuracy)** <br>

**Online tool:** [https://mafft.cbrc.jp/alignment/server/](https://mafft.cbrc.jp/alignment/server/)
- [x] Direction of nucleotide sequences: ``Adjust direction according to the first sequence``
- [x] Output order: ``Same as input``
- [x] Strategy: ``E-INS-i``
- [x] Scoring matrix for nucleotide sequences: ``1PAM/k=2``

**Command line:**
```
mafft --thread 8 --threadtb 5 --threadit 0 --inputorder --adjustdirection --kimura 1 --maxiterate 2 --retree 1 --genafpair INPUT_FILE.fasta > OUTPUT_FILE.fasta
```

**2. Rename your sequences** <br>
- Sequence names cannot contain spaces. Rename and simplify your sequence names using a text editor.
- Make sure that the reference sequence (e.g., either "KO3455" or "HXB2" for subtype B) is not in any of the sequences names, except the reference sequence.

**3. Submit your FASTA file to [HIV Database QCTool](https://www.hiv.lanl.gov/content/sequence/QC/index.html)** <br>
- Submit the FASTA file of aligned, primers-free, proviral sequences (including the reference sequence, usually HXB2). The results will be returned by email.
- Download the results as a text file: ``Download Summary``

**4. Submit your FASTA file to [HIV Database Gene Cutter](https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html)** <br>
- [x] ``Input sequences are pre-aligned``.
- Submit the FASTA file of aligned, primers-free, proviral sequences (including the reference sequence, usually HXB2). The results will be returned by email.

**5. Prepare your FASTA file for submission to ProSeq-IT** <br>
As per [ProSeq-IT's instructions](https://psd.cancer.gov/tools/tool_index.php), if the reference sequence is included in the FASTA file, it has to be renamed. <br>
- Duplicate your FASTA file.
- In this duplicated file, rename the reference sequence (usually HXB2) to **"Reference_sequence"**.

**6. Submit your FASTA file (containing the renamed reference sequence) to [ProSeq-IT](https://psd.cancer.gov/tools/pvs_annot.php)** <br>
- [x] Subtype: ``B``
- Make sure that the number of sequences in ProSeq-IT's output corresponds to the number of sequences in your FASTA file.
- Download the results as an Excel file: ``Download Result``.

**7. Complete the [Excel template](https://github.com/alemi055/IntegrityAlgorithm/blob/main/Template_IntegrityAlgorithm.xlsx)** <br>
The Excel templates contains 3 tabs:
- ProseqIT_criteria [locked]: used for the analysis of ProSeq-IT
- Manual_assessment: contains the results of the manual assessment with Geneious (to identify inversions) <br>
    - Fill the "Name column" with the names of the sequences.
    - Fill the "Inversions columns" with `Y` for inversions, or `N` otherwise.
- Hyperlinks: contains the hyperlinks of HIV Database QCTool and Gene Cutter's results
    - Fill the "Hyperink" column with the URL links of the results from QCTool and Gene Cutter (that were sent to your email).
 
### Files in directory before running the algorithm
- [x] The FASTA file containing the aligned sequences + reference sequence (HXB2)
- [x] The Excel template
- [x] The text file containing QCTool's results
- [x] The Excel file containing ProSeq-IT's results

# Installation

Install the R package `IntegrityAlgorithm` from Github

```
install.packages("devtools")
install_github("alemi055/IntegrityAlgorithm")
```

# Usage

### HIV_IntegrityAnalysis()

To analyze the intactness of HIV proviruses, run the function `HIV_IntegrityAnalysis()`

![Figure 1. IntegrityAnalysis](/Images/240212_IntegrityImage_GitHub.jpg)

1. The algorithm will first analyse the results from QCTool to identify ORFs containing stop codons and hypermutated sequences.
2. It will then analyse the results from Gene Cutter to obtain the start and stop codons, if applicable, of and within each ORF.
3. Next, it will analyse the results from ProSeq-IT for large and small internal deletions.
4. Finally, all results from QCTool, Gene Cutter, ProSeq-IT, and the manual assessment with Geneious (through the Excel template) will be combined to assess the intactness of the HIV proviruses.

```
HIV_IntegrityAnalysis(template_filename, QCTool_summary, ProseqIT_rx, ProseqIT_RefSeq = TRUE, RefSeq = TRUE, analyzes = 4)
```

| Argument | Required/Optional | Note |
| --- | --- | --- |
| `template_filename` | Required | Excel template containing the inversions manual assessment and the hyperlinks for QCTool and Gene Cutter |
|`QCTool_summary` | Required | Text file containing the results from QCTool |
| `ProseqIT_rx` | Required | Excel file containing the results from ProSeq-IT |
| `ProseqIT_RefSeq` | Optional; default if `TRUE` | Logical. If TRUE, the reference sequence (HXB2) is included in ProSeq-IT's results. It should be renamed to "Reference_sequence" (see [instructions](#before-running-the-algorithm))  |
| `RefSeq` | Optional argument; default is `TRUE` | Logical. If TRUE, the reference sequence (HXB2) is included in QCTool and GeneCutter's results |
| `analyzes` | optional argument; default is `4` | Specifies the analyzes to be done: <br> ``1`` QCTool only <br>``2`` Gene Cutter and ProSeq-IT only <br> ``3`` IntegrateInfo only <br> ``4`` All |

### Clonality_Analysis()

To analyze the clonality of the HIV proviruses, run the function `Clonality_Analysis()`

![Figure 2. ClonalityAnalysis](/Images/240212_ClonalityImage_GitHub.jpg)

1. Files in the FASTA file input will be splitted according to the list of donors supplied. Donors with single sequences will not be analyzed.
2. The clonality will be assessed on each of the splitted files. A list of clones (0 different nt between two sequences) and potential clones, established by a threshold, will be output.

```
Clonality_Analysis(FASTA_file, donors, threshold)
```
| Argument | Required/Optional | Note |
| --- | --- | --- |
| `FASTA_file` | Required | FASTA file containing all sequences, including the reference sequence (HXB2) |
| `donors` | Required | Str vector containing the list of donors. If more than one, use the `c()` function |
| `threshold` | Optional; default is `1` | Threshold number of different nucleotides to consider two sequences as "potential clones" |

## To confirm clones and potential clones in Geneious
1. Import the "donor1_forClonality.fasta" file into Geneious. In the imported alignment, select all sequences.
3. ``Extract regions`` :arrow_right: ``Extract region as a list of sequences``.
4. Select the newly created list.
5. ``Workflows`` :arrow_right: ``Run Workflow`` :arrow_right: ``Extract Sequences By Name``
6. Write the name of the unique sequence in the ``Names to Extract`` box. Then, write ``, `` (comma followed by a space) and the names of the clones/potential clones [copied-pasted from the **clones** and **potential clones** columns, in the "donor1_ClonalityAnalysis.csv" output files.]
7. From that, it is going to create a new sublist of sequences. Reselect all the sequences from this sublist and ``Extract regions``.
8. Select all the extracted sequences from the sublist.
9. ``Multiple Align`` :arrow_right: ``MAFFT Alignment`` <br>
If it is a clone, a message will appear: "All sequences are identical. Options will not affect the result of the alignment". Click ``Ok``.
10. Make sure that the **distances** between sequences are all 0 nucleotides for clones, or below the threshold value for potential clones.
11. Repeat for all clones/potential clones.

# Output

**tmp folder**
| Output file | Note |
| --- | --- |
| `Analyzed_GeneCutter.csv` | Analyzed results from Gene Cutter. <br> Details, for each sequence, the list of ORFs that have start and premature stop codons |
| `Analyzed_ProseqIT.csv` | Analyzed results from ProSeq-IT. <br> Details, for each sequence, the sequence length, the presence (binary: ``0`` [absence] and ``1`` [presence]) of large internal deletions, $\Psi$ mutations, and small internal deletions in each of the ORFs. For the ORFs, the list of detailed defects is also included |
| `Analyzed_QCTool.csv` | Analyzed results from QCTool. <br> Details, for each sequence, the number of stop codons, the list of stop codons, the number of incomplete codons, and the presence or absence of hypermutations |
| `intactness_detailedsummary.csv` | Detailed summary of defects for all sequences. <br> Details, for each sequence, the inferred intactness, the number of "main defects" (see [list 1-7 above](#defects)), the "main defect", and all defects in the sequence, including in each of the ORFs. <br> Note that this summary **does not** hierarchize the defects |

**FINAL_OUTPUT folder**
| Output file | Note |
| --- | --- |
| `intactness_summary.csv` | Summary of **hierarchized** defects for all sequences |
| `[donor]_ClonalityAnalysis.csv` | Can be found in the **Clonality** subfolder. List of clones and potential clones for each of the donor's sequences |

# Datasets
| Paper | GenBank accession numbers | PopSet |
| --- | --- | --- |
| Sannier et al., Cell Rep 2021 | MZ662560 - MZ662755 (n = 196) | [2083716276](https://www.ncbi.nlm.nih.gov/popset/2083716276)
| Dufour et al., Nat Commun 2023 | ON816029 - ON816663 (n = 635) | [2306699925](https://www.ncbi.nlm.nih.gov/popset/2306699925) |
| Dubé et al., Cell Host Microbe 2023 | OR105517 - OR105556 (n = 40) | [2569148931](https://www.ncbi.nlm.nih.gov/popset/2569148931) |
| Dufour et al., Cell Rep 2023 | ON816664 - ON817104 (n = 441) | [2306704384](https://www.ncbi.nlm.nih.gov/popset/2306704384)|

# Citation
If you use the IntegrityAlgorithm, please cite this repository.

# Contact information
If there are any questions, please contact: <br>

**Audrée Lemieux** <br>
lemieux.audree@gmail.com

