\name{GeneCutter_analyzes}
\alias{GeneCutter_analyzes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyze the results from Gene Cutter
}
\description{
This function analyzes the results from Gene Cutter. For each ORF, it opens the amino acid result (aa), identifies the presence of a start codon (presence of start codon: will be written in the temporary file), and then identifies the sequences with stop codons (presence of stop codon: will be written in the temporary file). This function creates a temporary file named "Analyzed_GeneCutter.csv" in the "tmp" folder.
}
\usage{
GeneCutter_analyzes(filename, RefSeq = TRUE)
}
\arguments{
  \item{filename}{
the name of the Excel template file
}
  \item{RefSeq}{
logical. If TRUE, the reference sequence is included in QCTool's and Gene Cutter's results. Default is TRUE
}
}
\details{
To analyze the data, the function assesses the URL link of the results (to access the start and stop codons). This function is used in the main function HIV_IntegrityAnalysis().
}
\value{
None
}
\references{
Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
}
\author{
Audrée Lemieux
}
\note{
Troubleshooting:\cr
# "Warning in file(file, ifelse(append, "a", "w")): cannot open file 'tmp/Analyzed_QCTool.csv': No such file or directory\cr
# Error in file(file, ifelse(append, "a", "w")): cannot open the connection

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
