\name{check_logical}
\alias{check_logical}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Verify if the \preformatted{RefSeq} and \preformatted{ProseqIT_RefSeq} values provided are logical
}
\description{
This function verifies if \preformatted{RefSeq} and \preformatted{ProseqIT_RefSeq} values provided are logical (either \preformatted{TRUE} or \preformatted{FALSE})
}
\usage{
check_logical(RefSeq, ProseqIT_RefSeq)
}
\arguments{
  \item{RefSeq}{
logical. If \preformatted{TRUE}, the reference sequence is included in QCTool's and Gene Cutter's results
}
  \item{ProseqIT_RefSeq}{
logical. If \preformatted{TRUE}, the reference sequence is included in ProSeq-IT's results
}
}
\details{
This function is used at the beginning of the main function \preformatted{HIV_IntegrityAnalysis()}.
}
\value{
None
}
% \references{
% QCTool website: \url{https://www.hiv.lanl.gov/content/sequence/QC/index.html}
% Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
% }
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
