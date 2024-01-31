\name{HIV_IntegrityAnalysis}
\alias{HIV_IntegrityAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Infer the intactness of HIV proviruses
}
\description{
This function integrates the analyzes from all tools (Geneious, QCTool, Gene Cutter, and ProSeq-IT) to infer the intactness of the provirus ("intact" or "defective"). For each sequence, if the sequence is defective, it will detail all the defects in the "intactness_detailedsummary.csv" file (in the "tmp" folder), or only the main defect (defects are hierarchized) in the "intactness_summary.csv" (in the "FINAL_OUTPUT" folder). It also outputs, in the "tmp" folder, one .csv file per tool (QCTool: "Analyzed_QCTool.csv"; Gene Cutter: "Analyzed_GeneCutter.csv"; ProSeq-IT: "Analyzed_ProseqIT.csv").
}
\usage{
HIV_IntegrityAnalysis(template_filename, QCTool_summary, ProseqIT_rx, ProseqIT_RefSeq = TRUE, RefSeq = TRUE, analyzes = 4)
}
\arguments{
  \item{template_filename}{
the name of the Excel template file
}
  \item{QCTool_summary}{
the name of the summary text (.txt) file downloaded from QCTool's output
}
  \item{ProseqIT_rx}{
the name of the summary Excel (.xls) file downloaded from ProSeq-IT
}
  \item{ProseqIT_RefSeq}{
logical. If TRUE, the reference sequence is included in ProSeq-IT's results. Default is TRUE
}
  \item{RefSeq}{
logical. If TRUE, the reference sequence is included in QCTool's and Gene Cutter's results. Default is TRUE
}
  \item{analyzes}{
integer. The functions to run. 1: QCTool only; 2: GeneCutter and ProSeq-IT; 3: IntegrateInfo only; 4: All. Default is 4
}
}
\details{
To analyze the data, the function assesses the results from all tools: manual assessment with Geneious (through the Excel template), QCTool, Gene Cutter, and Proseq-IT. Each provirus undergoes analysis to identify the presence of, in order, (i) inversions, (ii) hypermutations, (iii) large internal deletions, (iv) stop codons, (v) \eqn{\Psi} defects, and (vi) small internal deletions in gag, pol, vif, vpr, tat, rev, vpu, nef, env, and RRE. The defects are either hierarchized (i.e., only one defect per sequence, see "intactness_summary.csv" file) or not (exhaustive list of all defects, see "intactness_detailedsummary.csv").
}
\value{
None
}
\references{
QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
ProSeq-IT website: \url{https://psd.cancer.gov/tools/pvs_annot.php}
}
\author{
Audrée Lemieux
}
\note{
Troubleshooting:
"Warning in file(file, ifelse(append, "a", "w")): cannot open file 'tmp/Analyzed_QCTool.csv': No such file or directory
Error in file(file, ifelse(append, "a", "w")): cannot open the connection

This error occurs when the user does not have administrative rights. To fix this error, manually create the folders: "tmp" and "FINAL_OUTPUT".

~~~

Error in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 11, 10

This type of error occurs when the number of sequences between all files and tools (Excel template, QCTool, Gene Cutter, and ProSeq-IT) is not the same and/or in the same order. To fix this error, check the temporary .csv files (in the "tmp" folder) and verify at which tool the number of sequences varies. It can sometimes happen when the reference sequence is present (or absent) in one of the tools, and the \'RefSeq\' and \'ProseqIT_RefSeq\' options have not been modified.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
%
% ## The function is currently defined as
% function (x)
% {
%   }
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory (show via RShowDoc("KEYWORDS")):
% % \keyword{ ~kwd1 }
% % \keyword{ ~kwd2 }
% % Use only one keyword per line.
% % For non-standard keywords, use \concept instead of \keyword:
% % \concept{ ~cpt1 }
% % \concept{ ~cpt2 }
% % Use only one concept per line.
