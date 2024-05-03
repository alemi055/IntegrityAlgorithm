
########################################################################################
########################################################################################
########################################################################################

                                    ################## 
                                    # 1 Installation #
                                    ##################

# Only do this once

# # Install the "Integrity Algorithm package" (available on Github)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
#
# BiocManager::install("GenomicAlignments")
# 
# install.packages("devtools")
#
# # If you need to update the "IntegrityAlgorithm" package, rerun these 3 lines only:
# library(devtools)
# install_github("alemi055/IntegrityAlgorithm")

# Load package
library(IntegrityAlgorithm)

########################################################################################

                                  #######################
                                  # 2 Run the Algorithm #
                                  #######################

# Replace the "TO_FILL" with the names of your files/participants

# Assess the intactness of your sequences
HIV_IntegrityAnalysis(
  template_filename = "TO_FILL.xlsx", # Name of the Excel template
  QCTool_summary = "TO_FILL.txt", # Name of the txt file containing QCTool's results
  ProseqIT_rx = "TO_FILL.annot.xls", # Name of the Excel file containing ProSeq-IT's results
)


# Assess the Clonality of your sequences
Clonality_Analysis(
  FASTA_file = "TO_FILL.fasta", # Name of the FASTA file containing all seqs
  participants = c("TO_FILL", "TO_FILL", "TO_FILL") # Name of the participants
  # Default threshold is 1
)

########################################################################################
########################################################################################
########################################################################################
