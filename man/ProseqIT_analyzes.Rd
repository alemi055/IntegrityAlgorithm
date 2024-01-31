\name{ProseqIT_analyzes}
\alias{ProseqIT_analyzes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyze the results from ProSeq-IT
}
\description{
This function analyzes the results from ProSeq-IT. Verifies the sequence length and the intactness of \eqn{\Psi}. For each ORF, reads each annotation to identify defects (insertions, deletions, frameshifts). It then infers the intactness of the ORF (0: intact, 1: defect) and details, if the ORF is defective, the type of defect. This function creates a temporary file named "Analyzed_ProseqIT.csv" in the "tmp" folder.
}
\usage{
ProseqIT_analyzes(filename, ProseqIT_filename, ProseqIT_RefSeq = TRUE)
}
\arguments{
  \item{filename}{
the name of the Excel template file
}
  \item{ProseqIT_filename}{
the name of the summary Excel (.xls) file downloaded from ProSeq-IT
}
  \item{ProseqIT_RefSeq}{
logical. If TRUE, the reference sequence is included in ProSeq-IT's results. Default is TRUE
}
}
\details{
To analyze the data, the function first assesses the results from GeneCutter to obtain the start and stop codons. It then assesses the Excel (.xls) file downloaded from ProSeq-IT, using the criteria in the "ProseqIT_criteria" tab of the Excel template file. This function is used in the main function HIV_IntegrityAnalysis().
}
\value{
None
}
\references{
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
