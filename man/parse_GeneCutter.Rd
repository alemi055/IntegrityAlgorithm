\name{parse_GeneCutter}
\alias{parse_GeneCutter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parse Gene Cutter's results
}
\description{
This function reads the HIV protein's amino acid results page, obtained from Gene Cutter. It reads the lines of the web page and returns a data frame with the sequence names, the base positions, and the corresponding amino acid sequences.
}
\usage{
parse_GeneCutter(lines)
}
\arguments{
  \item{lines}{
lines of the HIV protein's amino acid results page, obtained from Gene Cutter
}
}
\details{
The URL link to Gene Cutter's results should be input in the Excel Template file for this function to work. This function is used in the main function GeneCutter_analyzes.
}
\value{
A data frame containing, for a specific HIV protein, the sequence names, the base positions, and the corresponding amino acid sequences.
}
\references{
Gene Cutter website: https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html
}
\author{
Audrée Lemieux
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
