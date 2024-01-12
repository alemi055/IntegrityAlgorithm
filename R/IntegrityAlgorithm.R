
########################################################################################
########################################################################################
########################################################################################

# Audrée Lemieux
# Updated on January 11, 2023
# Command version

########################################################################################

library(readxl)
library(RCurl)
library(XML)
library(stringr)
library(ape)
library(insect)
library(alakazam)

########################################################################################

                                    ##########################
                                    # Declare MAIN functions #
                                    ##########################

# MAIN FUNCTION 1 - HIV INTEGRITY ANALYSIS
HIV_IntegrityAnalysis <- function(template_filename, QCTool_summary, ProseqIT_rx, ProseqIT_RefSeq = TRUE, RefSeq = TRUE, analyzes = 4){
  # (str, str, str, log, log) - > None
  #
  # Input:
  #   - template_filename: the name of the Template file
  #   - QCTool_summary: the name of the summary .txt file downloaded from QCTool's output
  #   - ProseqIT_rx: the name of the summary .xls file downloaded from ProSeq-IT
  #   - ProseqIT_RefSeq: logical. If TRUE, the reference sequence is included in ProSeq-IT's results.
  #   - RefSeq: logical. If TRUE, the reference sequence is included in QCTool's and Gene Cutter's results.
  #   - analyzes: the functions to run. 1: QCTool only; 2: GeneCutter and ProSeq-IT; 3: IntegrateInfo only; 4: All
  #
  # Analyzes the results from QCTool, GeneCutter, ProseqIT, as well as our manual
  # assessment with Geneious, and combine them. Exports a CSV file summarizing
  # the integrity of all queried sequences, including their defects (if applicable).
  
  # Double check the parameters input
  check_template(template_filename)
  check_QCTool(QCTool_summary)
  check_ProseqIT(ProseqIT_rx)
  check_both_links(template_filename)
  check_logical(RefSeq, ProseqIT_RefSeq)
  check_integer(analyzes)
  
  if (analyzes == 1){
    QCTool_analyzes(template_filename, QCTool_summary, RefSeq) # QCTool
  }else if (analyzes == 2){
    GeneCutter_analyzes(template_filename, RefSeq) # GeneCutter
    ProseqIT_analyzes(template_filename, ProseqIT_rx, ProseqIT_RefSeq) # ProseqIT
  }else if (analyzes == 3){
    IntegrateInfo(template_filename) # Integrate info
  }else if (analyzes == 4){
    QCTool_analyzes(template_filename, QCTool_summary, RefSeq) # QCTool
    GeneCutter_analyzes(template_filename, RefSeq) # GeneCutter
    ProseqIT_analyzes(template_filename, ProseqIT_rx, ProseqIT_RefSeq) # ProseqIT
    IntegrateInfo(template_filename) # Integrate info
  }
}


# MAIN FUNCTION 2 - QCTool
QCTool_analyzes <- function(filename, summary, RefSeq = TRUE){
  # (str, str, log) -> None
  #
  # Input:
  #   - filename: the name of the Template file
  #   - summary: the name of the summary .txt file downloaded from QCTool's output
  #   - RefSeq: logical. If TRUE, the reference sequence is included in QCTool's and Gene Cutter's results.
  #
  # Analyzes the results from QCTool
  
  hyperlinks <- as.data.frame(read_xlsx(filename, sheet = "Hyperlinks"))
  QC_hyperlink <- hyperlinks$Hyperlink[which(hyperlinks$Tool == "QCTool")]
  
  ##################################################
  
  # Extract the stop codons #
  
  # Clean QCTool summary file
  lines <- readLines(summary)
  lines <- gsub("Cannot determine|Cannot Determine", "Cannotdetermine", lines)
  QCTool_excel <- read.table(textConnection(lines), header = T, sep = " ", row.names = NULL)
  QCTool_excel[,ncol(QCTool_excel)-1] <- paste0(QCTool_excel[,ncol(QCTool_excel)-1], " ", QCTool_excel[,ncol(QCTool_excel)]) # Resolve the last column (% of ACGT)
  colnames(QCTool_excel) <- c(colnames(QCTool_excel[2:ncol(QCTool_excel)]), "X")
  QCTool_excel <- QCTool_excel[,-(ncol(QCTool_excel))]
  rm(lines) # Clean space
  
  # Find the stop codons hyperlinks
  webpage <- getURL(QC_hyperlink)
  lines <- strsplit(webpage, "\n")[[1]]
  stop_hyperlinks <- lines[grep(".*/tmp/download/QC/.*/GENE_CUTTER/.*qcstop.html", lines)]
  stop_hyperlinks <- as.character(sapply(stop_hyperlinks, function(x){strsplit(x, " href=")[[1]][2]}))
  stop_hyperlinks <- gsub("</A></td>", "", stop_hyperlinks)
  
  # Create a df with the sequence name, the number of stop codons, and the hyperlink on the right
  QCstop_df <- NULL
  for (i in 1:length(stop_hyperlinks)){
    tmp <- strsplit(stop_hyperlinks, ">")[[i]]
    seqn <- as.numeric(strsplit(strsplit(tmp[1], "GENE_CUTTER/")[[1]][2], "/qcstop.html")[[1]][1])
    seqname <- QCTool_excel$SeqName[seqn + 1]
    QCstop_df <- rbind(QCstop_df, cbind(SeqName = seqname, nstop = tmp[2], URL = paste0("https://www.hiv.lanl.gov/", tmp[1])))
  }
  QCstop_df <- as.data.frame(QCstop_df)
  class(QCstop_df$nstop) <- "numeric"
  rm(i, lines, seqn, seqname, tmp, webpage, stop_hyperlinks) # Clean space
  
  # Add column: stop codons comments
  QCTool_excel <- cbind(QCTool_excel[,1:which(colnames(QCTool_excel) == "StopCodons")], stop_comments = NA, QCTool_excel[,(((which(colnames(QCTool_excel) == "StopCodons"))+1):ncol(QCTool_excel))])
  
  # Read webpages
  pb <- txtProgressBar(min = 1, max = nrow(QCstop_df), style = 3, width = 50, char = "=") # Add progress bar
  cat("\nNow analyzing the results from QCTool...\n")
  
  for (i in 1:nrow(QCstop_df)){
    # Update progress bar
    setTxtProgressBar(pb, i)
    
    seqn <- QCstop_df$SeqName[i]
    URL <- QCstop_df$URL[i]
    webpage <- getURL(URL)
    lines <- strsplit(webpage, "\n")[[1]]
    
    # Find the regions with stop codons in the middle
    pos <- grep("Regions with stop codons in the middle", lines) # Position of this little sentence
    n_stopcodons <- QCstop_df$nstop[i] # How many stop codons are there for this sequence?
    stopurl <- lines[(pos+1):(pos+n_stopcodons)] # URL links to the stop codons
    
    stop_codons <- NULL # Find the names of the genes where there are stop codons
    for (j in stopurl){
      tmp <- strsplit(j, "</a>")[[1]]
      tmp <- strsplit(tmp, ".aa>")[[1]][2]
      stop_codons <- c(stop_codons, tmp)
    }
    pos_QCTool_excel <- which(QCTool_excel$SeqName == seqn)
    QCTool_excel$stop_comments[pos_QCTool_excel] <- paste0(stop_codons, collapse = ", ")
  }
  rm(pb, i, j, lines, n_stopcodons, pos, seqn, stop_codons, stopurl, tmp, URL, webpage, pos_QCTool_excel) # Clean space
  
  # Remove weird characters in the column "Blast"
  QCTool_excel$Blast <- gsub("</a>", "", QCTool_excel$Blast)
  
  ##################################################
  
  # Export as CSV
  cols_to_keep <- c("SeqName", "StopCodons", "stop_comments", "IncompleteCodons", "Hypermutation")
  QCTool_excel <- QCTool_excel[,(match(cols_to_keep, colnames(QCTool_excel)))]
  if (RefSeq){ # If RefSeq is present - by default = TRUE
    QCTool_excel <- QCTool_excel[-1,] # Remove RefSeq (first line)
  }
  rm(cols_to_keep) # Clean space
  
  if (!file.exists("tmp")){system("mkdir tmp")}
  # if (!file.exists("tmp/RData")){system("mkdir tmp/RData")}
  write.table(QCTool_excel, "tmp/Analyzed_QCTool.csv", row.names = F, quote = F, sep = "\t", na = "")
  # save.image("tmp/RData/Analyzed_QCTools.RData")
  
  cat("\nDone.\n")
}


# MAIN FUNCTION 3 - GENE CUTTER
GeneCutter_analyzes <- function(filename, RefSeq = TRUE){
  # (str) -> None
  #
  # Input:
  #   - filename: the name of the Template file
  #
  # Analyzes the results from GeneCutter
  
  hyperlinks <- as.data.frame(read_xlsx(filename, sheet = "Hyperlinks"))
  GC_hyperlink <- hyperlinks$Hyperlink[which(hyperlinks$Tool == "GeneCutter")]
  
  # Extract the specific genes hyperlinks
  webpage <- getURL(GC_hyperlink)
  lines <- strsplit(webpage, "\n")[[1]]
  line <- lines[grep("aa.", lines)]
  line <- strsplit(line, "href=")[[1]]
  line <- line[grep("\\.aa\\.", line)]
  gene_hyperlinks <- as.character(sapply(line, function(x){strsplit(x, ">aa<")[[1]][1]}))
  gene_hyperlinks <- paste0("https://www.hiv.lanl.gov", gene_hyperlinks)
  
  # # Add columns: "start_codon"
  # # If there is a start codon: will write the name in start_codon
  # AnalyzedQCTool <- read.table("tmp/Analyzed_QCTool.csv", sep = "\t", header = T) # Read Analyzed_QCTool file
  # # AnalyzedQCTool <- AnalyzedQCTool[order(AnalyzedQCTool$SeqName),]
  # GC_excel <- as.data.frame(cbind(Name = AnalyzedQCTool$SeqName)) # Should be in the same order than the QCTool file
  
  # Read webpages
  genes <- c("Gag", "Pol", "Vif", "Vpr", "Tat", "Rev", "Vpu", "Env", "Nef")
  
  # data_list <- vector(mode = 'list', length = nrow(AnalyzedQCTool))
  # data_list_stop <- vector(mode = 'list', length = nrow(AnalyzedQCTool))
  # # AnalyzedQCTool <- AnalyzedQCTool[order(AnalyzedQCTool$SeqName),]
  # names(data_list) <- AnalyzedQCTool$SeqName
  # names(data_list_stop) <- AnalyzedQCTool$SeqName
  
  pb <- txtProgressBar(min = 1, max = length(gene_hyperlinks), style = 3, width = 50, char = "=") # Add progress bar
  cat("\nNow analyzing the results from Gene Cutter...\n")
  
  for (i in 1:length(gene_hyperlinks)){
    # Update progress bar
    setTxtProgressBar(pb, i)
    
    tmp <- strsplit(gene_hyperlinks[i], ".aa.html")[[1]]
    tmp <- strsplit(tmp, "/")[[1]]
    tmp <- tmp[length(tmp)]
    
    if (length(grep(tmp, genes)) > 0){
      
      # Find the start codons
      webpage_gene <- getURL(gene_hyperlinks[i])
      lines_gene <- strsplit(webpage_gene, "\n")[[1]]
      lines_gene <- strsplit(lines_gene, "<br>")[[1]]
      
      tmp1 <- grep("---------- List of .*", lines_gene) # Only keep everything above "List of Stop Codons"
      if (length(tmp1) != 0){ # If this sentence can't be found in the webpage, keep everything
        lines_start <- lines_gene[-((tmp1[1]-1):length(lines_gene))]
      }else{
        lines_start <- lines_gene
      }
      lines_start <- gsub("<a.*.html><font color=red>.*</font></a>", "", lines_start) # Remove when there is a hyperlink in red
      
      if (tmp != "Pol"){ # Only exception for start codon
        # Does the gene have a start codon?
        parsed_df <- parse_GeneCutter(lines_start)
        
        if (i == 1){ # Create the data frame and data list to store the analyzed data
          if (RefSeq){
            names <- parsed_df$seqn[2:(which(is.na(parsed_df$seqn))[1]-1)] # First one is the RefSeq (if in the Results)
          }else{
            names <- parsed_df$seqn[1:(which(is.na(parsed_df$seqn))[1]-1)]
          }
          GC_excel <- as.data.frame(cbind(Name = names)) # Should be in the same order than the QCTool file
          data_list <- vector(mode = 'list', length = length(names))
          data_list_stop <- vector(mode = 'list', length = length(names))
          names(data_list) <- names
          names(data_list_stop) <- names
        }
        
        parsed_df <- parsed_df[1:((which(is.na(parsed_df[1]) & is.na(parsed_df[2]))[1])-1),]  # Only keep first positions
        if (RefSeq){ # If RefSeq is present - by default = TRUE
          parsed_df <- parsed_df[-1,] # First one is the reference sequence
        }
        # parsed_df <- parsed_df[order(parsed_df$seqn),]
        # if (!identical(names(data_list), parsed_df$seqn)){
        #   last_c <- as.character(sapply(parsed_df$seqn, function(x){str_sub(x, -1)}))
        #   pos <- which(last_c == "_")
        #   if (length(pos) >= 1){
        #     for (i in pos){
        #       parsed_df$seqn[i] <- str_sub(parsed_df$seqn[i], 1, -2)
        #     }
        #   }
        # }
        # parsed_df <- parsed_df[match(names(data_list), parsed_df$seqn),]
        
        for (j in 1:nrow(parsed_df)){
          tmp2 <- strsplit(parsed_df$seq[j], "")[[1]]
          if (tmp2[1] == "M"){
            data_list[[j]] <- c(data_list[[j]], tmp)
          }
        }
      }

      # Find the stop codons
      lines_gene <- gsub("</pre>", "", lines_gene)
      tmp3 <- grep("---------- List of Stop Codons Within Sequences.*", lines_gene) # Find stop codons
      if (length(tmp3) != 0){
        tmp4 <- grep("^$", lines_gene[tmp3:length(lines_gene)])[1]-1 # -1 to get the line above
        lines_stop <- lines_gene[(tmp3+1):(tmp3+tmp4-1)]
        parsed_df_stop <- parse_GeneCutter_stop(lines_stop)
        
        for (k in 1:nrow(parsed_df_stop)){
          if (!grepl("^K03455", parsed_df_stop$seqn[k])){
            data_list_stop[[parsed_df_stop$seqn[k]]] <- c(data_list_stop[[parsed_df_stop$seqn[k]]], tmp)
          }
        }
      }
    }
  }
  rm(parsed_df, pb, i, j, k, line, lines, lines_gene, lines_start, lines_stop, tmp, tmp1, tmp2, tmp3, tmp4, webpage, webpage_gene, names) # Clean space
  
  # Add to Excel
  newcol <- NULL # Start codons
  newcol <- as.data.frame(sapply(data_list, function(x){
    tmp <- paste0(x, collapse = ", ")
    ncol <- rbind(newcol, tmp)
  }))
  colnames(newcol) <- "start_codon"
  
  newcol2 <- NULL # Start codons
  newcol2 <- as.data.frame(sapply(data_list_stop, function(x){
    tmp <- paste0(unique(x), collapse = ", ")
    ncol <- rbind(newcol2, tmp)
  }))
  colnames(newcol2) <- "stop_codon"
  GC_excel <- cbind(GC_excel, newcol, newcol2)
  rm(newcol, newcol2) # Clean space
  
  ##################################################
  
  # Export as CSV
  write.table(GC_excel, "tmp/Analyzed_GeneCutter.csv", row.names = F, quote = F, sep = "\t", na = "")
  # save.image("tmp/RData/Analyzed_GeneCutter.RData")
  cat("\nDone.\n")
}


# MAIN FUNCTION 4 - PROSEQIT
ProseqIT_analyzes <- function(filename, ProseqIT_filename, ProseqIT_RefSeq = TRUE){
  # (str) -> None
  #
  # Input:
  #   - filename: the name of the Template file
  #   - ProseqIT_filename: the name of the summary .xls file downloaded from ProSeq-IT
  #   - ProseqIT_RefSeq: logical. If TRUE, the reference sequence is included in ProSeq-IT's results.
  #
  # Analyzes the results from ProseqIT
  
  cat("\nNow analyzing the results from ProSeq-IT...\n")
  
  directory <- getwd()
  
  ProseqIT_summary <- NULL
  criteria <- c("Large_inter_delet", "psi_defects", "gag_small_delet", "pol_small_delet", "vif_small_delet", "vpr_small_delet", "tat_small_delet", "rev_small_delet", "vpu_small_delet", "nef_small_delet")
  
  # Also read the "Analyzed GeneCutter" Excel file
  AnalyzedGeneCutter <- read.table("tmp/Analyzed_GeneCutter.csv", header = TRUE, sep = "\t")
  
  ##################################################
  
  # Analyze ProseqIT results #
  
  # Clean ProseqIT file
  ProseqIT_excel <- read.table(ProseqIT_filename, header = T, sep = "\t", row.names = NULL, fill = TRUE)
  pos <- which(ProseqIT_excel[,1] ==  "ID") # Find the header of the table
  colnames(ProseqIT_excel) <- ProseqIT_excel[pos,]
  ProseqIT_excel <- ProseqIT_excel[-1,]
  # ProseqIT_excel <- ProseqIT_excel[-which(is.na(ProseqIT_excel[,2])),] # Removes empty rows
  ProseqIT_excel <- as.data.frame(ProseqIT_excel)
  # coln <- c(2:11, 13:16, 19, 20, 22:24, 26:28, 30:32, 34:36, 38:40, 42:44, 46:48, 50:52, 54, 56)
  # ProseqIT_excel[,coln] <- apply(ProseqIT_excel[,coln], 2, function(x) as.numeric(as.character(x))) # Change variables with numbers as "numerical" variables

  
  # # Is the RefSeq present in the file? In theory, should not be.
  if (ProseqIT_RefSeq){
    pos2 <- which(ProseqIT_excel$ID == "Reference_sequence")
    if (length(pos2) > 0){
      ProseqIT_excel <- ProseqIT_excel[-pos2,]
    }
  }
  
  # Load the criteria for ProseqIT
  ProseqIT_criteria <- as.data.frame(read_excel(filename, sheet = "ProseqIT_criteria"))
  
  # Create summary object
  ProseqIT_summary <- as.data.frame(ProseqIT_excel$seq_length)
  row.names(ProseqIT_summary) <- ProseqIT_excel$ID
  ProseqIT_summary <- cbind(ProseqIT_summary, ProseqIT_excel$seq_length)
  colnames(ProseqIT_summary) <- c("seq_length", criteria[1])
  
  # Large internal deletions
  pos <- which(as.numeric(ProseqIT_excel$seq_length) <= as.numeric(ProseqIT_criteria[ProseqIT_criteria$Column_name == "seq_length",3]))
  ProseqIT_summary$Large_inter_delet[pos] <- 1 # If large internal deletion
  ProseqIT_summary$Large_inter_delet[setdiff(1:nrow(ProseqIT_excel), pos)] <- 0
  
  # Psi defects
  ProseqIT_summary <- cbind(ProseqIT_summary, psi_defects(ProseqIT_excel, ProseqIT_criteria)) # Psi
  
  # Small internal deletions
  ProseqIT_summary <- cbind(ProseqIT_summary, gag_small_delet(ProseqIT_excel, ProseqIT_criteria, AnalyzedGeneCutter)) # Gag
  ProseqIT_summary <- cbind(ProseqIT_summary, pol_small_delet(ProseqIT_excel, ProseqIT_criteria, ProseqIT_summary, AnalyzedGeneCutter)) # Pol
  ProseqIT_summary <- cbind(ProseqIT_summary, vif_vpr_tat_rev_vpu_nef_small_delet(ProseqIT_excel, ProseqIT_criteria, AnalyzedGeneCutter)) # Vif, Vpr, Tat, Rev, Vpu, Nef
  ProseqIT_summary <- cbind(ProseqIT_summary, env_small_delet(ProseqIT_excel, ProseqIT_criteria, AnalyzedGeneCutter)) # Env
  
  # Invert Env and Nef
  ProseqIT_summary <- ProseqIT_summary[,c(1:18, 21, 22, 19, 20)]
  
  # RRE status
  pos <- which(ProseqIT_excel$rre_status != "good")
  ProseqIT_summary$rre_status <- 0
  if (length(pos) > 0){
    ProseqIT_summary$rre_status[pos] <- 1 # If large internal deletion
  }
  
  ##################################################
  
  # Export as csv
  ProseqIT_summary <- cbind(Name = row.names(ProseqIT_summary), ProseqIT_summary)
  write.table(ProseqIT_summary, file = "tmp/Analyzed_ProseqIT.csv", quote = F, sep = "\t", na = "", row.names = F)
  
  # save.image("tmp/RData/Analyzed_ProseqIT.RData")
  cat("Done.\n")
}


# MAIN FUNCTION 5 - Integrate all info
IntegrateInfo <- function(filename){
  # (str) -> None
  #
  # Input:
  #   - filename: the name of the Template file
  #
  # Analyzes the results from QCTool
  
  # cat("\n\n############################################\n# Now integrating information of all tools #\n############################################\n\n")
  
  directory <- getwd()
  
  # Load analyzed QCTool, GeneCutter, and ProseqIT files
  QC_filen <- "Analyzed_QCTool.csv"
  GC_filen <- "Analyzed_GeneCutter.csv"
  ProseqIT_filen <- "Analyzed_ProseqIT.csv"
  QCTool_excel <- as.data.frame(read.table(paste0(directory, "/tmp/", QC_filen), sep = "\t", header = TRUE))
  GC_excel <- as.data.frame(read.table(paste0(directory, "/tmp/", GC_filen), sep = "\t", header = TRUE))
  ProseqIT_excel <- as.data.frame(read.table(paste0(directory, "/tmp/", ProseqIT_filen), sep = "\t", header = TRUE))
  
  # Load template for manual assessment
  manual_excel <- read_xlsx(paste0(directory, "/", filename), sheet = "Manual_assessment")
  
  # Put all sequences in alphabetical order
  manual_excel <- manual_excel[order(manual_excel$Name),]
  ProseqIT_excel <- ProseqIT_excel[order(ProseqIT_excel$Name),]
  QCTool_excel <- QCTool_excel[order(QCTool_excel$SeqName),]
  GC_excel <- GC_excel[order(GC_excel$Name),]
  
  # Check if the seq names of all analyzed tools are the same. If not, compute error
  if (!(all(sapply(list(manual_excel$Name, ProseqIT_excel$Name, QCTool_excel$SeqName, GC_excel$Name), FUN = identical, manual_excel$Name)))){
    stop("\nThe number of sequences or the sequence names are not the same (or are not in the same order) in all files.")
  }
  
  # Intactness summary + sequence length + number of main defects + empty column
  intact_summary <- as.data.frame(cbind(Name = ProseqIT_excel[,1], seq_length = ProseqIT_excel[,2]))
  
  ##################################################
  
  # 1 Extract the inversions #
  
  # Inversion = defective
  intact_summary <- as.data.frame(cbind(intact_summary, inversions = manual_excel$Inversions))
  intact_summary$inversions <- gsub("Y|y", 1, intact_summary$inversions)
  intact_summary$inversions <- gsub("N|n", 0, intact_summary$inversions)
  
  ##################################################
  
  # 2 Hypermutations #
  
  # Hypermutation = defective
  intact_summary <- cbind(intact_summary, hypermutations = QCTool_excel$Hypermutation)
  intact_summary$hypermutations <- gsub("Possible", 1, intact_summary$hypermutations)
  intact_summary$hypermutations <- gsub("Not Detected|NotDetected", 0, intact_summary$hypermutations)
  
  ##################################################
  
  # 3 Large internal deletions #
  
  # Large internal deletion (<8800 bp without primers) = defective
  intact_summary <- cbind(intact_summary, large_intern_delet = ProseqIT_excel$Large_inter_delet)
  
  ########################################################################################
  
  # 4 Stop codons #
  
  # Stop codons in all proteins, except for Nef and Tat2 = defective
  # Look at GeneCutter first. If there is a Tat stop codon, look at QCTool results to see
  # if it's only in Tat2
  stop_codon <- NULL
  
  for (i in 1:nrow(GC_excel)){
    stops <- strsplit(GC_excel$stop_codon[i], ", ")[[1]]
    stops <- gsub("Nef", "", stops) # Nef
    if (length(grep("Tat", stops)) > 0){
      if ((length(grep("Tat2", QCTool_excel$stop_comments[i])) > 0) & (length(grep("Tat1", QCTool_excel$stop_comments[i])) == 0)){ # Is Tat2 in QCTool results?
        stops <- gsub("Tat", "", stops) # Tat is not a def
      }
    }
    stops <- str_subset(stops, ".+") # Keep only those with at least one character
    
    if (length(stops) > 0){
      stop_codon <- c(stop_codon, 1)
    }else{
      stop_codon <- c(stop_codon, 0)
    }
  }
  intact_summary <- cbind(intact_summary, stop_codon = stop_codon)
  rm(i, stop_codon, stops) # Clean space
  
  # # Add which stop codon
  # intact_summary <- cbind(intact_summary, stopcodon_comments = GC_excel$stop_codon)
  
  ##################################################
  
  # 5 Psi mutations #
  
  intact_summary <- cbind(intact_summary, psi_defects = ProseqIT_excel$Psi_defects)
  
  # # Psi deletion, SL2 deletion, or MDS point mutation = defective
  # intact_summary <- cbind(intact_summary, psi_defects = manual_excel$`Psi defects`)
  # intact_summary$psi_defects <- gsub("Y|y", 1, intact_summary$psi_defects)
  # intact_summary$psi_defects <- gsub("N|n", 0, intact_summary$psi_defects)
  
  ##################################################
  
  # 6 Small Internal Deletions #
  
  # Deletion in Gag, Pol, Vif, Vpr, Tat, Rev, Vpu, or Env = defective
  # Deletion in Nef only = not defective
  small_intern_delet <- NULL
  colnames <- c("Gag_defects", "Pol_defects", "Vif_defects", "Vpr_defects", "Tat_defects", "Rev_defects", "Vpu_defects", "Env_defects", "Nef_defects", "rre_status")

  tmp <- ProseqIT_excel[,as.numeric(sapply(colnames, function(x){which(colnames(ProseqIT_excel) == x)}))]
  tmp <- tmp[,-9] # Defect in Nef doesn't count as a small internal deletion
  sum_sid <- as.numeric(apply(tmp, 1, sum))
  for (i in 1:nrow(intact_summary)){
    if (sum_sid[i] > 0){
      small_intern_delet <- c(small_intern_delet, 1)
    }else{
      small_intern_delet <- c(small_intern_delet, 0)
    }
  }
  intact_summary <- cbind(intact_summary, small_intern_delet)
  
  # Add back the defects in Nef. Calculate the number of small internal deletions
  tmp <- ProseqIT_excel[,as.numeric(sapply(colnames, function(x){which(colnames(ProseqIT_excel) == x)}))]
  sum_sid <- as.numeric(apply(tmp, 1, sum))
  intact_summary <- cbind(intact_summary, n_small_inter_delet = sum_sid, tmp)
  rm(tmp, i, small_intern_delet, sum_sid, colnames) # Clean space
  
  # # Invert Env and Nef for small internal deletions
  # intact_summary <- intact_summary[,c(1:17, 19, 18, 20)]

  ##################################################
  
  # 7 Summarize the information #
  
  colnames <- c("inversions", "hypermutations", "large_intern_delet", "stop_codon", "psi_defects", "small_intern_delet") # A defect in small intern delet = defect in at least one of the individual (Gag, Pol, etc.) gene
  intactness <- NULL
  defects_comments <- vector(mode = 'list', length = nrow(ProseqIT_excel))
  names(defects_comments) <- ProseqIT_excel$Name
  ndefects <- NULL
  tmp <- NULL
  
  for (i in 1:nrow(intact_summary)){
    for (j in colnames){
      pos <- which(colnames(intact_summary) == j)
      if (intact_summary[i,pos] == 1){
        defects_comments[[i]] <- c(defects_comments[[i]], j)
      }
    }
    tmp <- rbind(tmp, paste0(defects_comments[[i]], collapse = ", "))
    
    # Intact?
    if (length(defects_comments[[i]]) == 0){
      intactness <- c(intactness, "intact")
    }else{
      intactness <- c(intactness, "defective")
    }
    
    # Number of "main" defects (1+ small internal deletions count for one)
    ndefects <- c(ndefects, length(defects_comments[[i]]))
  }
  # intact_summary$number_of_defects <- ndefects
  # intact_summary <- cbind(name = intact_summary$Name, nseq = intact_summary$nseq, intactness = intactness, defects_comments = tmp, empty_column = NA, intact_summary[,3:ncol(intact_summary)])
  intact_summary <- cbind(Name = intact_summary$Name, intactness = intactness, n_main_defects = ndefects, defects_comments = tmp, empty_column = NA, intact_summary[,2:ncol(intact_summary)])
  rm(defects_comments, i, j, intactness, ndefects, pos) # Clean space
  
  
  # Add main defect (if any) for each sequence
  main_defect <- NULL
  for (i in 1:nrow(intact_summary)){
    tmp <- NULL
    if (intact_summary$intactness[i] == "defective"){
      tmp <- strsplit(intact_summary$defects_comments, ", ")[[i]][1]
      main_defect <- c(main_defect, tmp)
    }else{
      tmp <- NA
      main_defect <- c(main_defect, tmp)
    }
  }
  intact_summary <- cbind(intact_summary[,c(1:3)], main_defect, intact_summary[,4:ncol(intact_summary)])
  rm(i, tmp) # Clean space
  
  # Remove title of "empty column"
  colnames(intact_summary) <- gsub("empty_column", "", colnames(intact_summary))
  
  # Export as CSV
  write.table(intact_summary, "tmp/intactness_detailedsummary.csv", quote = F, row.names = F, sep = "\t", na = "")
  
  ##################################################
  ##################################################

  # Hierarchize the defects
  # Everytime a sequence has a defect, it cannot be considered for "lower" defects
  intact_main <- cbind(intact_summary[1:nrow(QCTool_excel),c(1:2)])
  defects_ord <- c("inversions", "hypermutations", "large_intern_delet", "stop_codon", "psi_defects", "small_intern_delet") 
  
  intact_main <- cbind(intact_main, as.data.frame(matrix(data = 0, nrow = nrow(intact_main), ncol = length(defects_ord))))
  colnames(intact_main)[(ncol(intact_main)-length(defects_ord)+1):(ncol(intact_main))] <- defects_ord
  for (i in 1:nrow(intact_main)){
    pos <- which(colnames(intact_main) == main_defect[i])
    intact_main[i,pos] <- 1
  }
  
  # Add columns for all the small internal deletions
  intact_main <- cbind(intact_main, empty_column = NA, intact_summary[1:nrow(QCTool_excel),(which(colnames(intact_summary) == "Gag_defects")):ncol(intact_summary)])
  
  # Add column for sequence length
  intact_main <- cbind(intact_main[,c(1:2)], seq_length = ProseqIT_excel$seq_length, intact_main[,c(3:ncol(intact_main))])
  
  # Add column for number of main defects
  intact_main <- cbind(intact_main[,c(1:3)], n_main_defects = intact_summary$n_main_defects, empty_column = NA, intact_main[,c(4:ncol(intact_main))])
  
  # Remove title of "empty column"
  colnames(intact_main) <- gsub("empty_column", "", colnames(intact_main))
  
  ########################################################################################
  
  # Export as CSV
  # system("mkdir FINAL_OUTPUT")
  if (!file.exists("FINAL_OUTPUT")){system("mkdir FINAL_OUTPUT")}
  write.table(intact_main, "FINAL_OUTPUT/intactness_summary.csv", quote = F, row.names = F, sep = "\t", na = "")
  
  # save.image("tmp/RData/integrateinfo.RData")
  
  cat("\n############################################\n")
  cat("\nJOB COMPLETED.\n\nThe summary of intactness can be found as a CSV file in the folder \'FINAL OUTPUT\'\nAll other temporary files (i.e., outputs of individual analyzes) can be found in the \'tmp\' folder.\n\nNote: Manually confirm the Psi defects by looking at the alignment in Geneious.")
}


# MAIN FUNCTION 6
Split_files <- function(FASTA_file, donors){
  # (str, str) -> None
  #
  # Input:
  #   - FASTA_file: name of the FASTA files containing all aligned sequences.
  #   - donors: name of the donors. Use c() if more than one donor.
  #
  # Split the sequences in individual files (per donor)
  
  all_seqs <- read.FASTA(FASTA_file)
  
  # pb <- txtProgressBar(min = 1, max = length(donors), style = 3, width = 50, char = "=") # Add progress bar
  cat("\nNow splitting sequences into individual files...\n")
  
  for (i in 1:length(donors)){
    # setTxtProgressBar(pb, i)
    pos <- grep(donors[i], names(all_seqs))
    if (length(pos) > 0){
      # cat(paste0("Now doing donor \'", donors[i], "\'\n"))
      tmp <- all_seqs[pos]
      if (length(tmp) > 1){
        write.FASTA(tmp, paste0(donors[i], "_forClonality.fasta"))
      }else{
        cat(paste0("   - The donor with the name \'", donors[i], "\' only had 1 sequence. No file created.\n"))
      }
      
    }else{
      cat(paste0("   - The donor with the name \'", donors[i], "\' could not be found.\n"))
    }
  }
  cat("Done.\n")
}


# MAIN FUNCTION 7
Clonality_Analysis <- function(threshold = 5){
  # (int) -> None
  #
  # Input:
  #   - threshold: the threshold of different nucleotides between two sequences to consider them as "potential clones"
  #
  # Analyzes the clonality of the sequences
  
  directory <- getwd()
  fasta <- list.files(paste0(directory, "/"), pattern = "_forClonality.fasta")
  
  
  cat("\nNow analyzing the clonality...\n")
  
  suppressWarnings(
    for (i in fasta){
      # print(paste0("--- Now doing ", i, " ---"))
      seqs <- read.FASTA(i)
      # if (RefSeq){ # Removes the RefSeq
      #   seqs <- seqs[-1]
      # }
      
      # Convert DNAbin to a character vector
      seqs_vector <- dna2char(seqs)
      
      # Look at the differences between sequences
      # Gaps are treated as universally non-matching characters
      # Look at seqs with less than 5 (or threshold) differences - if 0 = clone, but if between 1 and threshold = potential clones
      pairwise_matrix <- pairwiseDist(seqs_vector, dist_mat = getDNAMatrix(gap = 1))
      
      unique_seqs <- list()
      seqs_names <- row.names(pairwise_matrix)
      
      for (j in 1:length(seqs_names)){
        name <- seqs_names[j]
        pos <- as.numeric(which(pairwise_matrix[j,] <= threshold))
        
        if (length(pos) > 1){ # If length(pos) = 1, it is only the pairwise comparison with itself. It needs to be at least 2 for a duplicate.
          pos <- pos[-which(pos == j)] # Remove self-comparison
          tmp <- NULL
          tmp2 <- NULL
          for (k in pos){
            if (pairwise_matrix[j,k] == 0){ # identical seq
              tmp <- c(tmp, seqs_names[k])
            }else if (pairwise_matrix[j,k] >= 1 & pairwise_matrix[j,k] <= threshold){ # potential seq
              tmp2 <- c(tmp2, seqs_names[k])
            }
          }
          # tmp <- paste0(tmp, collapse = ", ")
          # tmp2 <- paste0(tmp2, collapse = ", ")
          tmp <- as.data.frame(cbind(tmp))
          tmp2 <- as.data.frame(cbind(tmp2))
          
          if (nrow(tmp2) == 0){ # Only clones
            unique_seqs[[name]] <- list(tmp, data.frame())
          }else if (nrow(tmp) == 0){ # Only potential clones
            unique_seqs[[name]] <- list(data.frame(), tmp2)
          }else{
            unique_seqs[[name]] <- list(tmp, tmp2)
          }
          
        }else{
          unique_seqs <- c(unique_seqs, list(NULL))
          names(unique_seqs)[length(unique_seqs)] <- name
        }
      }
      rm(j, k, name, pos, tmp, tmp2) # Clean space
      
      
      # Create a df from the list
      flag <- TRUE
      j <- 1
      n <- length(unique_seqs)
      clonality_df <- NULL
      
      while (flag){
        if (j <= n){
          tmp <- unique_seqs[[j]]
          
          if (length(tmp) > 0){ # If there are clones/potential clones
            if (nrow(tmp[[1]]) > nrow(tmp[[2]])){
              nrow <- nrow(tmp[[1]])
            }else{
              nrow <- nrow(tmp[[2]])
            }
            
            tmp_df <- data.frame(name = rep(names(unique_seqs)[j], nrow), clones = NA, potential_clones = NA) # Create "empty" df
            if (nrow(tmp[[1]]) > 0){
              tmp_df$clones[1:nrow(tmp[[1]])] <- as.character(unlist(tmp[[1]]))
              
              for (k in as.character(unlist(tmp[[1]]))){ # Remove seqs (downstream) to avoid duplicates
                pos <- which(names(unique_seqs) == k)
                if (length(pos) > 0){
                  unique_seqs <- unique_seqs[-pos]
                }
              }
            }
            if (nrow(tmp[[2]]) > 0){
              tmp_df$potential_clones[1:nrow(tmp[[2]])] <- as.character(unlist(tmp[[2]]))
            }
            
            clonality_df <- rbind(clonality_df, tmp_df)
            j <- j + 1
            n <- length(unique_seqs)
            
          }else{ # If there are no clones
            clonality_df <- rbind(clonality_df, data.frame(name = names(unique_seqs)[j], clones = NA, potential_clones = NA))
            j <- j + 1
            n <- length(unique_seqs)
          }
          
        }else{
          flag <- FALSE
        }
      }
      rm(flag, j, n, k, nrow, tmp, tmp_df) # Clean space
      
      
      # Number of clones and potential clones
      clonality_df <- cbind(name = clonality_df$name, n_clones = NA, clones = clonality_df$clones, n_potential_clones = NA, potential_clones = clonality_df$potential_clones)
      clonality_df <- as.data.frame(clonality_df)
      unique_seq_names <- as.character(unique(na.omit(clonality_df$name)))
      
      for (j in unique_seq_names){
        pos <- which(clonality_df$name == j)
        clonality_df$n_clones[pos[1]] <- length(na.omit(clonality_df$clones[pos]))
        clonality_df$n_potential_clones[pos[1]] <- length(na.omit(clonality_df$potential_clones[pos]))
      }
      
      
      # Write the name of the unique sequence (with clones) only once
      pos <- table(clonality_df$name)
      pos <- names(pos)[which(pos > 1)]
      
      if (length(pos) > 0){
        for (j in pos){
          pos2 <- which(clonality_df$name == j)
          clonality_df$name[pos2[2:length(pos2)]] <- NA
        }
      }
      
      
      # Export
      if (!file.exists("FINAL_OUTPUT/Clonality")){system("mkdir FINAL_OUTPUT/Clonality")}
      write.table(clonality_df, paste0("FINAL_OUTPUT/Clonality/", gsub("forClonality.fasta", "ClonalityAnalysis.csv", i)), na = "", row.names = F, sep = "\t", quote = F)
    }
  )

  cat("Done.\n")
}


########################################################################################

                                ###########################
                                # Declare other functions #
                                ###########################

# Function 1
parse_GeneCutter <- function(lines){
  # (str) -> df
  #
  # Input:
  #   - lines: lines of sequence in GeneCutter
  #
  # Returns a df containing the parsed lines
  
  df <- NULL
  
  for (i in 2:length(lines)){
    newline <- NULL
    tmp <- str_squish(lines[i])
    tmp <- strsplit(tmp, " ")[[1]]
    
    newline <- cbind(seqn = tmp[1], pos1 = as.numeric(tmp[2]), seq = tmp[3], pos2 = as.numeric(tmp[4]))
    df <- rbind(df, newline)
  }
  df <- as.data.frame(df)
  return(df)
}


# Function 2
parse_GeneCutter_stop <- function(lines){
  # (str) -> df
  #
  # Input:
  #   - lines: lines of sequence in GeneCutter
  #
  # Returns a df containing the parsed lines
  
  df <- NULL
  
  for (i in 1:length(lines)){
    newline <- NULL
    tmp <- str_squish(lines[i])
    tmp <- strsplit(tmp, " ")[[1]]
    
    newline <- cbind(seqn = tmp[1], pos = as.numeric(tmp[3]))
    df <- rbind(df, newline)
  }
  df <- as.data.frame(df)
  return(df)
}


# Function 3
psi_defects <- function(ProseqIT_rx, ProseqIT_criteria){
  # (df, df) -> df
  #
  # Input:
  #   - ProseqIT_rx: df containing the results from ProseqIT
  #   - ProseqIT_criteria: df containing the criteria for each ProseqIT variable
  #
  # Returns a df containing (1) the presence/absence of Psi defects and (2) the defects themselves
  
  cols <- c("msd_status", "package_deletion") # Cols to look at
  pos <- match(cols, colnames(ProseqIT_rx)) # Indexes of columns
  
  tmp <- as.data.frame(cbind(Psi_defects = rep(0, nrow(ProseqIT_rx)), Psi_defects_comments = rep("", nrow(ProseqIT_rx)))) # Summarized object with comments
  row.names(tmp) <- ProseqIT_rx$ID
  if (length(cols) != length(pos)){
    cat("Psi defects: You have a problem with the column names of ProseqIT results.\n\n")
  }
  
  # Note each defect in the "comments" section
  comments <- vector(mode = 'list', length = nrow(ProseqIT_rx))
  
  # Assess the intactness based on each criteria
  for (j in 1:length(cols)){
    name <- cols[j]
    row_crit <- which(ProseqIT_criteria[,1] == name)
    
    if (name == "msd_status"){
      tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "correct")]
    }else if(name == "package_deletion"){
      tmp_seqn <- ProseqIT_rx$ID[which(as.numeric(ProseqIT_rx[,pos[j]]) >= as.numeric(ProseqIT_criteria[row_crit,3]))]
    }
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
  }
  
  # Add comments to tmp object
  for (i in 1:length(comments)){
    if (length(comments[[i]]) != 0){
      tmp[i,2] <- paste0(comments[[i]], collapse = ", ")
    }
  }
  return(tmp)
}


# Function 4
gag_small_delet <- function(ProseqIT_rx, ProseqIT_criteria, Analyzed_GeneCutter){
  # (df, df) -> df
  #
  # Input:
  #   - ProseqIT_rx: df containing the results from ProseqIT
  #   - ProseqIT_criteria: df containing the criteria for each ProseqIT variable
  #   - Analyzed_GeneCutter: df containing the analyzed results from GeneCutter
  #
  # Returns a df containing (1) the presence/absence of Gag defects and (2) the defects themselves
  
  # cols <- c("U5_gag_pair_R2_deletion", "2base_before_gag_status", "gag_start_codon", "gag_insertion", "gag_deletion", "gag_frameshift", "gag_stop_codon") # Cols to look at
  cols <- c("U5_gag_pair_R2_deletion", "2base_before_gag_status", "gag_insertion", "gag_deletion", "gag_frameshift") # Cols to look at
  pos <- match(cols, colnames(ProseqIT_rx)) # Indexes of columns
  
  tmp <- as.data.frame(cbind(Gag_defects = rep(0, nrow(ProseqIT_rx)), Gag_defects_comments = rep("", nrow(ProseqIT_rx)))) # Summarized object with comments
  row.names(tmp) <- ProseqIT_rx$ID
  if (length(cols) != length(pos)){
    cat("Gag small internal deletion: You have a problem with the column names of ProseqIT results.\n\n")
  }
  
  # Note each defect in the "comments" section
  comments <- vector(mode = 'list', length = nrow(ProseqIT_rx))
  
  
  # Assess the intactness based on each criteria
  # First, look if there is a start codon
  name <- "StartCodon"
  tmp_seqn <- Analyzed_GeneCutter$Name[setdiff(1:nrow(Analyzed_GeneCutter), grep("Gag", Analyzed_GeneCutter$start_codon))] # Those who don't have a Start codon
  if (length(tmp_seqn)){
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
  }
  
  # Assess the intactness based on each criteria
  for (j in 1:length(cols)){
    name <- cols[j]
    row_crit <- which(ProseqIT_criteria[,1] == name)
    
    if (name == "U5_gag_pair_R2_deletion"){
      tmp_seqn <- ProseqIT_rx$ID[which(as.numeric(ProseqIT_rx[,pos[j]]) >= as.numeric(ProseqIT_criteria[row_crit,3]))]
    # }else if(name == "2base_before_gag_status" | name == "gag_start_codon"){
    }else if(name == "2base_before_gag_status"){
      tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "correct")]
    }else if (name == "gag_frameshift"){
      tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "no")]
    # }else if(name == "gag_insertion" | name == "gag_deletion" | name == "gag_stop_codon"){
    }else if(name == "gag_insertion" | name == "gag_deletion"){
      tmp_seqn <- ProseqIT_rx$ID[which(as.numeric(ProseqIT_rx[,pos[j]]) > as.numeric(ProseqIT_criteria[row_crit,3]))]
    }
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
  }
  
  # Look at stop codons
  pos2 <- match(Analyzed_GeneCutter$Name[grep("Gag", Analyzed_GeneCutter$stop_codon)], row.names(tmp))
  tmp[pos2,1] <- 1
  comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], "gag_stop_codon")}else{comments[[x]] <- comments[[x]]}})
  
  # Add comments to tmp object
  for (i in 1:length(comments)){
    if (length(comments[[i]]) != 0){
      tmp[i,2] <- paste0(comments[[i]], collapse = ", ")
    }
  }
  return(tmp)
}


# Function 5
pol_small_delet <- function(ProseqIT_rx, ProseqIT_criteria, ProseqIT_summary, Analyzed_GeneCutter){
  # (df, df) -> df
  #
  # Input:
  #   - ProseqIT_rx: df containing the results from ProseqIT
  #   - ProseqIT_criteria: df containing the criteria for each ProseqIT variable
  #   - ProseqIT_summary: df that will be the summary of the results from ProseqIT
  #   - Analyzed_GeneCutter: df containing the analyzed results from GeneCutter
  #
  # Returns a df containing (1) the presence/absence of Pol defects and (2) the defects themselves
  
  # cols <- c("gag_start_codon", "pol_insertion", "pol_deletion", "pol_frameshift", "pol_stop_codon") # Cols to look at
  cols <- c("pol_insertion", "pol_deletion", "pol_frameshift") # Cols to look at
  pos <- match(cols, colnames(ProseqIT_rx)) # Indexes of columns
  
  tmp <- as.data.frame(cbind(Pol_defects = rep(0, nrow(ProseqIT_rx)), Pol_defects_comments = rep("", nrow(ProseqIT_rx)))) # Summarized object with comments
  row.names(tmp) <- ProseqIT_rx$ID
  if (length(cols) != length(pos)){
    cat("Pol small internal deletion: You have a problem with the column names of ProseqIT results.\n\n")
  }
  
  # Note each defect in the "comments" section
  comments <- vector(mode = 'list', length = nrow(ProseqIT_rx))
  
  
  # Assess the intactness based on each criteria
  # First, look if there is a start codon for Gag
  name <- "StartCodon"
  tmp_seqn <- row.names(ProseqIT_summary)[grep("StartCodon", ProseqIT_summary$Gag_defects_comments)] # Which seqs don't have a StartCodon?
  if (length(tmp_seqn) > 0){
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], "gag_start_codon")}else{comments[[x]] <- comments[[x]]}})
  }
  
  for (j in 1:length(cols)){
    name <- cols[j]
    row_crit <- which(ProseqIT_criteria[,1] == name)
    
    # if(name == "gag_start_codon"){
    #   tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "correct")]
    if (name == "pol_frameshift"){
      tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "no")]
    # }else if(name == "pol_insertion" | name == "pol_deletion" | name == "pol_stop_codon"){
    }else if(name == "pol_insertion" | name == "pol_deletion"){
      tmp_seqn <- ProseqIT_rx$ID[which(as.numeric(ProseqIT_rx[,pos[j]]) > as.numeric(ProseqIT_criteria[row_crit,3]))]
    }
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
  }
  
  # Look at stop codons
  pos2 <- match(Analyzed_GeneCutter$Name[grep("Pol", Analyzed_GeneCutter$stop_codon)], row.names(tmp))
  tmp[pos2,1] <- 1
  comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], "pol_stop_codon")}else{comments[[x]] <- comments[[x]]}})
  
  # Add comments to tmp object
  for (i in 1:length(comments)){
    if (length(comments[[i]]) != 0){
      tmp[i,2] <- paste0(comments[[i]], collapse = ", ")
    }
  }
  return(tmp)
}


# Function 6
vif_vpr_tat_rev_vpu_nef_small_delet <- function(ProseqIT_rx, ProseqIT_criteria, Analyzed_GeneCutter){
  # (df, df, df) -> df
  #
  # Input:
  #   - ProseqIT_rx: df containing the results from ProseqIT
  #   - ProseqIT_criteria: df containing the criteria for each ProseqIT variable
  #   - Analyzed_GeneCutter: df containing the analyzed results from GeneCutter
  #
  # Returns a df containing (1) the presence of Vif/Vpr/Tat/Rev/Vpu/Nef defects and (2) the defects themselves
  
  proteins <- c("Vif", "Vpr", "Tat", "Rev", "Vpu", "Nef")
  final_tmp <- as.data.frame(cbind(Name = ProseqIT_rx$ID))
  
  for (i in proteins){
    lower_protein <- tolower(i)
    # cols <- paste0("vif", c("_deletion", "_frameshift", "stop_codon")) # Cols to look at in ProseqIT
    cols <- paste0(lower_protein, c("_deletion", "_frameshift")) # Cols to look at in ProseqIT
    pos <- match(cols, colnames(ProseqIT_rx)) # Indexes of columns
    
    tmp <- as.data.frame(cbind(rep(0, nrow(ProseqIT_rx)), rep("", nrow(ProseqIT_rx)))) # Summarized object with comments
    colnames(tmp) <- paste0(i, c("_defects", "_defects_comments"))
    row.names(tmp) <- ProseqIT_rx$ID
    if (length(cols) != length(pos)){
      cat(i, "small internal deletion: You have a problem with the column names of ProseqIT results.\n\n")
    }
    
    # Note each defect in the "comments" section
    comments <- vector(mode = 'list', length = nrow(ProseqIT_rx))
    
    
    # Assess the intactness based on each criteria
    # First, look if there is a start codon
    name <- "StartCodon"
    tmp_seqn <- Analyzed_GeneCutter$Name[setdiff(1:nrow(Analyzed_GeneCutter), grep(i, Analyzed_GeneCutter$start_codon))] # Those who don't have a Start codon
    if (length(tmp_seqn)){
      pos2 <- match(tmp_seqn, row.names(tmp))
      tmp[pos2,1] <- 1
      comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
    }
    
    for (j in 1:length(cols)){
      name <- cols[j]
      row_crit <- which(ProseqIT_criteria[,1] == name)
      
      if (name == paste0(lower_protein, "_frameshift")){
        tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "no")]
        # }else if(name == paste0(lower_protein, "_deletion") | name == paste0(lower_protein, "_stop_codon")){
      }else if(name == paste0(lower_protein, "_deletion")){
        tmp_seqn <- ProseqIT_rx$ID[which(as.numeric(ProseqIT_rx[,pos[j]]) > as.numeric(ProseqIT_criteria[row_crit,3]))]
      }
      pos2 <- match(tmp_seqn, row.names(tmp))
      tmp[pos2,1] <- 1
      comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
    }
    
    # Look at stop codons
    pos2 <- match(Analyzed_GeneCutter$Name[grep(i, Analyzed_GeneCutter$stop_codon)], row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], paste0(lower_protein, "_stop_codon"))}else{comments[[x]] <- comments[[x]]}})
    
    # Add comments to tmp object
    for (k in 1:length(comments)){
      if (length(comments[[k]]) != 0){
        tmp[k,2] <- paste0(comments[[k]], collapse = ", ")
      }
    }
    final_tmp <- cbind(final_tmp, tmp)
  }
  final_tmp <- final_tmp[,-1]
  return(final_tmp)
}

# Function 7
env_small_delet <- function(ProseqIT_rx, ProseqIT_criteria, Analyzed_GeneCutter){
  # (df, df, df) -> df
  #
  # Input:
  #   - ProseqIT_rx: df containing the results from ProseqIT
  #   - ProseqIT_criteria: df containing the criteria for each ProseqIT variable
  #   - Analyzed_GeneCutter: df containing the analyzed results from GeneCutter
  #
  # Returns a df containing (1) the presence/absence of Env defects and (2) the defects themselves
  
  # cols <- c("env_insertion", "env_deletion", "env_frameshift", "env_stop_codon") # Cols to look at
  cols <- c("env_insertion", "env_deletion", "env_frameshift") # Cols to look at
  pos <- match(cols, colnames(ProseqIT_rx)) # Indexes of columns
  
  tmp <- as.data.frame(cbind(Env_defects = rep(0, nrow(ProseqIT_rx)), Env_defects_comments = rep("", nrow(ProseqIT_rx)))) # Summarized object with comments
  row.names(tmp) <- ProseqIT_rx$ID
  if (length(cols) != length(pos)){
    cat("Env small internal deletion: You have a problem with the column names of ProseqIT results.\n\n")
  }
  
  # Note each defect in the "comments" section
  comments <- vector(mode = 'list', length = nrow(ProseqIT_rx))
  
  
  # Assess the intactness based on each criteria
  # First, look if there is a start codon
  name <- "StartCodon"
  tmp_seqn <- Analyzed_GeneCutter$Name[setdiff(1:nrow(Analyzed_GeneCutter), grep("Env", Analyzed_GeneCutter$start_codon))]
  if (length(tmp_seqn)){
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
  }
  
  # Assess the intactness based on each criteria
  for (j in 1:length(cols)){
    name <- cols[j]
    row_crit <- which(ProseqIT_criteria[,1] == name)
    
    if (name == "env_frameshift"){
      tmp_seqn <- ProseqIT_rx$ID[which(ProseqIT_rx[,pos[j]] != "no")]
    # }else if(name == "env_insertion" | name == "env_deletion" | name == "env_stop_codon"){
    }else if(name == "env_insertion" | name == "env_deletion"){
      tmp_seqn <- ProseqIT_rx$ID[which(as.numeric(ProseqIT_rx[,pos[j]]) > as.numeric(ProseqIT_criteria[row_crit,3]))]
    }
    pos2 <- match(tmp_seqn, row.names(tmp))
    tmp[pos2,1] <- 1
    comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], name)}else{comments[[x]] <- comments[[x]]}})
  }
  
  # Look at stop codons
  pos2 <- match(Analyzed_GeneCutter$Name[grep("Env", Analyzed_GeneCutter$stop_codon)], row.names(tmp))
  tmp[pos2,1] <- 1
  comments <- lapply(1:length(comments), function(x){if (x %in% pos2){comments[[x]] <- c(comments[[x]], "env_stop_codon")}else{comments[[x]] <- comments[[x]]}})
  
  # Add comments to tmp object
  for (i in 1:length(comments)){
    if (length(comments[[i]]) != 0){
      tmp[i,2] <- paste0(comments[[i]], collapse = ", ")
    }
  }
  return(tmp)
}


# Function 8
check_template <- function(template_filename){
  # (str) -> bool
  #
  # Input:
  #   - template_filename: the name of the Template file
  #
  # Returns TRUE if the template file provided is good, or FALSE otherwise
  
  files_directory <- list.files()
  
  if (!grepl(".xlsx$", template_filename)){
    stop("\n  The template file provided is not an Excel file (.xlsx).")
  }else if (length(grep(template_filename, files_directory)) == 0){
    stop("\n  The template file provided could not be found in your current directory.")
  }else if (!all(excel_sheets(template_filename) == c("ProseqIT_criteria", "Manual_assessment", "Hyperlinks"))){
    stop("\n  One or more sheets are missing from the template file provided.")
  }
}


# Function 9
check_QCTool <- function(QCTool_summary){
  # (str) -> bool
  #
  # Input:
  #   - QCTool_summary: the name of the summary .txt file downloaded from QCTool's output
  #
  # Returns TRUE if the file provided is good, or FALSE otherwise
  
  files_directory <- list.files()
  
  if (!grepl(".txt$", QCTool_summary)){
    stop("\n  The QCTool summary file provided is not a text file (.txt).")
  }else if (length(which(files_directory == QCTool_summary)) == 0){
    stop("\n  The QCTool summary file provided could not be found in your current directory.")
  }else{
    tmp <- read.table(QCTool_summary, header = T, sep = " ")
    if (length(grep("SeqName", colnames(tmp))) == 0){
      stop("\n  The column \'SeqName\' is at least missing from your file.")
    }
  }
}


# Function 10
check_link_QC <- function(hyperlink){
  # (str) -> bool
  #
  # Input:
  #   - hyperlink: URL link of LANL's HIV Sequence Database's QCTool's results
  #
  # Returns TRUE if the URL link provided is good, or FALSE otherwise
  
  lines <- getURL(hyperlink)
  
  if (!grepl("https://", hyperlink)){
    stop("\n  The URL link provided does not start with \'https://'.")
  }else if (!grepl("www.hiv.lanl.gov", hyperlink)){
    stop("\n  The URL link provided is not from LANL's HIV Sequence Database.")
  }else if (!grepl("/download/QC/.*/summary.html", hyperlink)){
    stop("\n  The URL link provided does not contain the results from QCTool.")
  }else if (grepl("The document has moved", lines)){
    stop("\n  The URL link provided for QCTool has expired.")
  }
}


# Function 11
check_link_GC <- function(hyperlink){
  # (str) -> bool
  #
  # Input:
  #   - hyperlink: URL link of LANL's HIV Sequence Database's Gene Cutter's results
  #
  # Returns TRUE if the URL link provided is good, or FALSE otherwise
  
  lines <- getURL(hyperlink)
  
  if (!grepl("https://", hyperlink)){
    stop("\n  The URL link provided does not start with \'https://'.")
  }else if (!grepl("www.hiv.lanl.gov", hyperlink)){
    stop("\n  The URL link provided is not from LANL's HIV Sequence Database.")
  }else if (!grepl("/tmp/GENE_CUTTER/.*/out.html", hyperlink)){
    stop("\n  The URL link provided does not contain the results from Gene Cutter.")
  }else if (grepl("The document has moved", lines)){
    stop("\n  The URL link provided for Gene Cutter has expired.")
  }
}


# Function 12
check_both_links <- function(template_filename){
  # (str) -> bool
  #
  # Input:
  #   - template_filename: the name of the Template file
  #
  # Returns TRUE if the URL link provided are good, or FALSE otherwise
  
  hyperlinks <- as.data.frame(read_xlsx(template_filename, sheet = "Hyperlinks"))
  
  # QCTool
  tmp1 <- hyperlinks$Hyperlink[which(hyperlinks$Tool == "QCTool")]
  if (is.na(tmp1) | tmp1 == ""){
    stop("\n  No URL link was provided for QCTool.")
  }else{
    check_link_QC(tmp1)
  }
  
  # GeneCutter
  tmp2 <- hyperlinks$Hyperlink[which(hyperlinks$Tool == "GeneCutter")]
  if (is.na(tmp2) | tmp1 == ""){
    stop("\n  No URL link was provided for GeneCutter.")
  }else{
    check_link_GC(tmp2)
  }
}


# Function 13
check_logical <- function(RefSeq, ProseqIT_RefSeq){
  # (bool, bool) -> None
  #
  # Input:
  #   - RefSeq: logical. If TRUE, the reference sequence is included in QCTool's and Gene Cutter's results.
  #   - ProseqIT_RefSeq: logical. If TRUE, the reference sequence is included in ProSeq-IT's results.
  #
  # Returns TRUE if the values are logical, or FALSE otherwise
  
  if (!is.logical(RefSeq)){
    stop("\n  The value provided for the argument \'RefSeq\' is not logical.")
  }else if (!is.logical(ProseqIT_RefSeq)){
    stop("\n  The value provided for the argument \'ProseqIT_RefSeq\' is not logical.")
  }
}


# Function 14
check_integer <- function(analyzes){
  # (int) -> None
  #
  # Input:
  #   - RefSeq: logical. If TRUE, the reference sequence is included in QCTool's and Gene Cutter's results.
  #   - ProseqIT_RefSeq: logical. If TRUE, the reference sequence is included in ProSeq-IT's results.
  #   - analyzes: the functions to run. 1: QCTool only; 2: GeneCutter and ProSeq-IT; 3: IntegrateInfo only; 4: All
  #
  # Returns TRUE if the values are logical, or FALSE otherwise
  
  if (analyzes != 1 & analyzes != 2 & analyzes != 3 & analyzes != 4){
    stop("\n  The value provided for the argument \'analyzes\' is not an integer between 1 and 4.")
  }
}


# Function 15
check_ProseqIT <- function(ProseqIT_rx){
  # (str) -> bool
  #
  # Input:
  #   - ProseqIT_rx: the name of the summary .xls file downloaded from ProSeq-IT
  #
  # Returns TRUE if the file provided is good, or FALSE otherwise
  
  files_directory <- list.files()
  
  if (!grepl(".xls$", ProseqIT_rx)){
    stop("\n  The ProSeq-IT results file provided is not an Excel file (.xls and not .xlsx).")
  }else if (length(which(files_directory == ProseqIT_rx)) == 0){
    stop("\n  The ProSeq-IT results file provided could not be found in your current directory.")
  }else{
    tmp <- read.table(ProseqIT_rx, header = T, sep = "\t", row.names = NULL, fill = TRUE)
    if (length(grep("ID", tmp[1])) == 0){
      stop("\n  The column \'ID\' is at least missing from your file.")
    }
  }
}

########################################################################################
########################################################################################
########################################################################################

