\name{IntegrateInfo}
\alias{IntegrateInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrate information of all tools
}
\description{
This function integrates the analyzes from all tools (Geneious, QCTool, Gene Cutter, and ProSeq-IT) to infer the intactness of the provirus ("intact" or "defective"). For each sequence, if the sequence is defective, it will detail all the defects in the "intactness_detailedsummary.csv" file (in the "tmp" folder), or only the main defect in the "intactness_summary.csv" (in the "FINAL_OUTPUT" folder).
}
\usage{
IntegrateInfo(filename)
}
\arguments{
  \item{filename}{
the name of the Excel template file
}
}
\details{
To analyze the data, the function assesses the results from all tools: manual assessment with Geneious (through the Excel template), QCTool, Gene Cutter, and Proseq-IT. This function is used in the main function HIV_IntegrityAnalysis().
}
\value{
None
}
% \references{
% ProSeq-IT website: \url{https://psd.cancer.gov/tools/pvs_annot.php}
% }
\author{
Audrée Lemieux
}
\note{
Troubleshooting:\cr
# "Warning in file(file, ifelse(append, "a", "w")): cannot open file 'tmp/Analyzed_QCTool.csv': No such file or directory\cr
# Error in file(file, ifelse(append, "a", "w")): cannot open the connection

This error occurs when the user does not have administrative rights. To fix this error, manually create the folders: "tmp" and "FINAL_OUTPUT".

~~~

# Error in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 11, 10

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
