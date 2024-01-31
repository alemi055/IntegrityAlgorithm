\name{IntegrateInfo}
\alias{IntegrateInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Integrate information of all tools
}
\description{
This function integrates the analyzes from all tools (\strong{Geneious}, \strong{QCTool}, \strong{Gene Cutter}, and \strong{Proseq-IT}) to infer the intactness of the provirus (\strong{intact} or \strong{defective}).

For each sequence, if the sequence is defective, it will detail all the defects in the \file{intactness_detailedsummary.csv} file (in the \strong{tmp} folder), or only the main defect (defects are hierarchized) in the \file{intactness_summary.csv} (in the \strong{FINAL_OUTPUT} folder).
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
\strong{Troubleshooting}

~~~~~~~~~~~~~~~~~~~~

\code{Warning in file(file, ifelse(append, "a", "w")): cannot open file 'tmp/Analyzed_QCTool.csv': No such file or directory}\cr
\code{Error in file(file, ifelse(append, "a", "w")): cannot open the connection}

This error occurs when the user does not have administrative rights. To fix this error, manually create the folders: \strong{tmp} and \strong{FINAL_OUTPUT}.

~~~~~~~~~~~~~~~~~~~~

\code{Error in data.frame(..., check.names = FALSE): arguments imply differing number of rows: 11, 10}

This type of error occurs when the number of sequences between all files and tools (Excel template, QCTool, Gene Cutter, and ProSeq-IT) is not the same and/or in the same order.

To fix this error, check the temporary .csv files (in the \strong{tmp} folder) and verify at which tool the number of sequences varies. It could be in the downloaded files (QCTool's summary .txt file and/or ProSeq-IT's .xls file) as well.

It can sometimes happen when the reference sequence is present (or absent) in one of the tools, and the \samp{RefSeq} and \samp{ProseqIT_RefSeq} options have not been modified.
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
