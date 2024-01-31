\name{vif_vpr_tat_rev_vpu_nef_small_delet}
\alias{vif_vpr_tat_rev_vpu_nef_small_delet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify defects in Vif, Vpr, Tat, Rev, Vpu, and Nef
}
\description{
This function identifies HIV Vif, Vpr, Tat, Rev, Vpu, and Nef defects in sequences based on ProSeq-IT's results, except for the stop codons (Gene Cutter).

It returns, for each sequence, a data frame summarizing the presence or absence (\samp{0} or \samp{1}) of Vif, Vpr, Tat, Rev, Vpu, and Nef defects.

Further details are added in the \emph{[protein_name]_defects_comments} if there is a defect.
}
\usage{
vif_vpr_tat_rev_vpu_nef_small_delet(ProseqIT_rx, ProseqIT_criteria, Analyzed_GeneCutter)
}
\arguments{
  \item{ProseqIT_rx}{
Data frame containing the results from ProSeq-IT (the Excel downloaded file)
}
  \item{ProseqIT_criteria}{
Data frame containing, for each ProSeq-IT annotation, the criteria used to infer proviral intactness
}
  \item{Analyzed_GeneCutter}{
Data frame containing the analyzed results from Gene Cutter
}
}
\details{
Based on ProSeq-IT's values \emph{[protein_name]_insertion}, \emph{[protein_name]_deletion}, \emph{[protein_name]_frameshift} (to see values, see the \strong{ProseqIT_criteria} tab of the Excel template).

To be considered \strong{intact}, the ORFs must have a start codon, no frameshift, and insertions/deletions representing <5\% of their length.

For this function to work, Gene Cutter's results must have been previously analyzed and ProSeq-IT's output Excel file must be in the user's current directory.

This function is used in the main function \samp{ProseqIT_analyzes()}.
}
\value{
A data frame summarizing the Pol defects.
}
\references{
Gene Cutter website: \url{https://www.hiv.lanl.gov/content/sequence/GENE_CUTTER/cutter.html}
\cr ProSeq-IT website: \url{https://psd.cancer.gov/tools/pvs_annot.php}
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
