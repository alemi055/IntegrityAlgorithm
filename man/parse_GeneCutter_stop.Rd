\name{parse_GeneCutter_stop}
\alias{parse_GeneCutter_stop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse Gene Cutter's stop codons results
}
\description{
This function reads and parses Gene Cutter's HIV ORF's amino acid 'stop codons' section. It reads the lines detailing the stop codons and returns a data frame with the sequence names and the position of the stop codon.
}
\usage{
parse_GeneCutter_stop(lines)
}
\arguments{
  \item{lines}{
lines of Gene Cutter's HIV ORF's amino acid 'stop codons' section
}
}
\details{
For this function to work, the URL link to Gene Cutter's results should be input in the Excel Template file. This function is used in the main function GeneCutter_analyzes().
}
\value{
A data frame containing, for a specific HIV ORF, the sequence names and the position of the stop codon.
}
\references{
Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
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
