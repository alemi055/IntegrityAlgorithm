\name{psi_defects}
\alias{psi_defects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify \eqn{\Psi} defects
}
\description{
This function identifies \eqn{\Psi} defects in sequences based on ProSeq-IT's results. It returns, for each sequence, a data frame summarizing the presence or absence (0 or 1) of \eqn{\Psi} defects. Further details are added in the 'Psi_defects_comments' if there is a defect.
}
\usage{
psi_defects(ProseqIT_rx, ProseqIT_criteria)
}
\arguments{
  \item{ProseqIT_rx}{
Data frame containing the results from ProSeq-IT (the Excel downloaded file)
}
  \item{ProseqIT_criteria}{
Data frame containing, for each ProSeq-IT annotation, the criteria used to infer proviral intactness
}
}
\details{
For this function to work, ProSeq-IT's output Excel file must be in your current directory. This function is used in the main function ProseqIT_analyzes().
}
\value{
A data frame summarizing the {\eqn{\Psi} defects.
}
\references{
ProSeq-IT website: \url{https://psd.cancer.gov/tools/pvs_annot.php}
}
\author{
Audrée Lemieux
}
\note{
%%  ~~further notes~~
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
