\name{check_ProseqIT}
\alias{check_ProseqIT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the ProSeq-IT's Excel file provided is usable
}
\description{
This function verifies if the ProSeq-IT's Excel file provided is indeed an Excel file, if it can be found in your current directory, and if it contains the appropriate header.
}
\usage{
check_ProseqIT(ProseqIT_rx)
}
\arguments{
  \item{ProseqIT_rx}{
the name of the summary Excel (.xls) file downloaded from ProSeq-IT
}
}
\details{
This function is used at the beginning of the main function \preformatted{HIV_IntegrityAnalysis()}.
}
\value{
\preformatted{TRUE} if the ProSeq-IT's Excel file is usable, or \preformatted{FALSE} otherwise.
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
