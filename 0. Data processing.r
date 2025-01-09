create_directory_structure <- function(base_dir = ".") {
  dirs <- c(
    "./1_basic_analysis/",
    
  
  for(dir in dirs) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Created directory: %s", dir))
  }
}


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("mixOmics", "Amelia", "pcaMethods"))

install.packages(c("svglite"))


install.packages(c("extrafont"))

library(readxl)
library(writexl)
library(stringr)
library(devtools)
library(Biobase)
library(preprocessCore)
library(ggplot2)
library(gplots)
library(NMF)
library(impute)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggthemes)
library(ggbreak)
library(tidyverse) 
library(ggvenn)
library(extrafont)
library(dplyr)
library(rstatix)
library(ggpubr)
library(vsn)
library(stats)
library(modeest)
library(mixtools)
library(genefu)
library(ggalluvial)
library(seqinr)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rentrez)
library(GenomicRanges)
library(pbapply)
library(progress)
library(gtable)
library(svglite)


required_packages <- c("ggplot2", "viridis", "gridExtra", "scales", "dplyr", "tidyr")
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}





create_directory_structure()

sample.info <- read.table("./sample.info.txt", sep="\t", header=TRUE, quote="")

xls_list <- list.files(pattern=".xlsx")
pep_list <- xls_list[grepl("pep", xls_list)]

message(sprintf("\nProcessing %d TMT batch files...", length(pep_list)))

pb <- progress_bar$new(
  format = "Processing TMT batches [:bar] :percent eta: :eta",
  total = length(pep_list),
  clear = FALSE,
  width = 60
)

###Total Global data
total_list <- list()

for(i in 1:length(pep_list)){
  pb$tick()
  message(sprintf("\nProcessing batch %d/%d: %s", i, length(pep_list), pep_list[i]))
  
  # Data load
  tmp_data <- read_excel(pep_list[i], sheet=1)
  
  # Confidence Filtering
  confidence_pos <- grepl("Confidence ", colnames(tmp_data))
  tmp_data <- tmp_data[as.data.frame(tmp_data)[, confidence_pos] == "High",]
  
  # XCorr Filtering
  xcorr_col <- grep("XCorr .* Sequest HT", colnames(tmp_data))
  xcorr_values <- as.numeric(as.character(unlist(tmp_data[,xcorr_col])))
  tmp_data <- tmp_data[xcorr_values >= 2.0,]
  
  col_pos <- grepl("^Abundances [[:punct:]]Grouped[[:punct:]][[:punct:]]|Annotated Sequence|Master Protein Accessions|XCorr (by Search Engine): Sequest HT|PSM", 
                   colnames(tmp_data))
  
  colnames(tmp_data)[grepl("^Abundances [[:punct:]]Grouped[[:punct:]][[:punct:]]", 
                           colnames(tmp_data))] <- 
    gsub("^Abundances [[:punct:]]Grouped[[:punct:]][[:punct:]] ", "", 
         colnames(tmp_data)[grepl("^Abundances [[:punct:]]Grouped[[:punct:]][[:punct:]]", 
                                  colnames(tmp_data))])
  
  tmp_col_data <- tmp_data[,col_pos]
  
  batch_name <- as.numeric(gsub("B", "", unlist(strsplit(pep_list[i], "_"))[4]))
  report_ion <- sample.info$Reporter_ion[sample.info$Batch==batch_name]
  sample_id <- sample.info$Sample_ID[sample.info$Batch==batch_name]
  
  # Abundance Column Matching
  tmp_abundance <- colnames(tmp_col_data)[grepl("^[0-9][0-9][0-9].*", colnames(tmp_col_data))]
  col_match_pos <- match(tmp_abundance, report_ion, nomatch=0)
  colnames(tmp_col_data)[grepl("^[0-9][0-9][0-9].*", colnames(tmp_col_data))] <- sample_id[col_match_pos]
  
  colnames(tmp_col_data)[grepl("Annotated Sequence|PSM|Master Protein Accessions", 
                               colnames(tmp_col_data))] <- c("Sequence", "PSM", "Accessions")
  
  # Sequence chech
  tmp_col_data$Sequence <- gsub("^[[:punct:]](.)+[[:punct:]][.]|[.][[:punct:]](.)+$", "", 
                                tmp_col_data$Sequence)
  
  total_list[[i]] <- tmp_col_data
}

names(total_list) <- unlist(lapply(pep_list, function(x){
  unlist(strsplit(x, "_"))[4]
}))

message("\nCreating total peptide list...")

# total_peptide
total_peptide <- as.vector(na.omit(unique(unlist(lapply(pep_list, function(x){
  tmp_data <- read_excel(x, sheet=1)
  
  # Confidence Filtering
  confidence_pos <- grepl("Confidence ", colnames(tmp_data))
  tmp_data <- tmp_data[as.data.frame(tmp_data)[, confidence_pos] == "High",]
  
  # XCorr Filtering 
  xcorr_col <- grep("XCorr .* Sequest HT", colnames(tmp_data))
  xcorr_values <- as.numeric(as.character(unlist(tmp_data[,xcorr_col])))
  tmp_data <- tmp_data[xcorr_values >= 2.0,]
  
  unique(tmp_data$'Annotated Sequence')
})))))

total_peptide <- gsub("^[[:punct:]](.)+[[:punct:]][.]|[.][[:punct:]](.)+$", "", total_peptide)

message("\nIntegrating multiple batches data...")
pb <- progress_bar$new(
  format = "Processing batch [:bar] :percent eta: :eta",
  total = length(total_list),
  clear = FALSE,
  width = 60
)

for(x in 1:length(total_list)){
  pb$tick()
  
  tmp_df <- data.frame(Sequence = total_peptide, test = rep(NA, length(total_peptide)))
  tmp_data <- total_list[[x]]
  tmp_pos2 <- match(tmp_data$Sequence, total_peptide, nomatch=0)
  quant_pos <- !grepl("Sequence|PSM|Accessions|Ref", colnames(tmp_data))
  
  if(sum(tmp_data$Sequence == total_peptide[tmp_pos2]) == dim(tmp_data)[1]){
  }else{
    message(sprintf("\nWarning: Position error in batch %s", names(total_list)[x]))
    next
  }
  
  Ref_pos <- grepl("^Ref", colnames(tmp_data))
  
  # CR ratio
  tmp_data_new <- tmp_data
  CR_ratio <- t(matrix(unlist(apply(tmp_data_new, 1, function(x){
    as.numeric(x[quant_pos])/mean(as.numeric(x[Ref_pos]), na.rm=T)
  })), ncol = dim(tmp_data_new)[1]))
  
  tmp_df[tmp_pos2, c(2:(sum(quant_pos)+1))] <- CR_ratio
  colnames(tmp_df)[c(2:(sum(quant_pos)+1))] <- colnames(tmp_data)[quant_pos]
  quant_df <- tmp_df[,-1]
  
  tmp_seq_ac <- tmp_data[, c("Sequence", "Accessions")]
  
  if(x==1){
    total_df <- tmp_df
    total_seq_ac <- tmp_seq_ac
  }else{
    total_df <- cbind(total_df, quant_df)
    total_seq_ac <- rbind(total_seq_ac, tmp_seq_ac)
  }
  
  message(sprintf("Batch %s processed", names(total_list)[x]))
}

# Accessions Matching
total_seq_ac <- unique(total_seq_ac)
total_df$Accessions <- total_seq_ac$Accessions[match(total_df$Sequence, 
                                                     total_seq_ac$Sequence, nomatch=0)]

save(total_list, file="./1_basic_analysis/raw_data/total_list.RData")
write.table(total_df, "./1_basic_analysis/raw_data/total_df.txt", 
            sep="\t", row.names=FALSE, quote=FALSE)
writexl::write_xlsx(total_df, "./1_basic_analysis/raw_data/total_df.xlsx")

message("\nInitial TMT data processing completed.")



################# 4. RefSeq/UniProt Filtering #################

# RefSeq Database Loading
message("\nInitiating database filtering process...")

read_refseq_database <- function() {
  message("Loading RefSeq database...")
  
  # .faa file
  fasta_files <- list.files(pattern = "human.*\\.protein\\.faa$")
  message(sprintf("Found %d FASTA files", length(fasta_files)))
  
  # Progress bar 
  pb <- progress_bar$new(
    format = "Processing FASTA files [:bar] :percent eta: :eta",
    total = length(fasta_files),
    clear = FALSE,
    width = 60
  )
  
  # Sequence list
  all_sequences <- lapply(fasta_files, function(file) {
    pb$tick()
    message(sprintf("Reading file: %s", file))
    read.fasta(file, seqtype="AA", as.string=TRUE, strip.desc=TRUE)
  })
  
  message("Processing sequences for tryptic peptides...")
  
  pb <- progress_bar$new(
    format = "Processing sequences [:bar] :percent",
    total = length(unlist(lapply(all_sequences, unlist))),
    clear = FALSE,
    width = 60
  )
  
  #RefSeq sequence to tryptic peptide
  refseq_tryptic <- unlist(lapply(unlist(lapply(all_sequences, unlist)), function(x) {
    pb$tick()
    tmp_k <- gsub("K", "K@", x)
    tmp_r <- gsub("R", "R@", tmp_k)
    tmp_r <- gsub("K@P", "KP", tmp_r)
    tmp_r <- gsub("R@P", "RP", tmp_r)
    peptides <- strsplit(tmp_r, "@")[[1]]
    peptides[nchar(peptides) >= 7 & nchar(peptides) <= 100]
  }))
  
  message(sprintf("Total unique RefSeq tryptic peptides: %d", length(unique(refseq_tryptic))))
  return(unique(refseq_tryptic))
}

# RefSeq Filtering
message("\nStarting RefSeq filtering...")
refseq_tryptic <- read_refseq_database()

# Matching & Filtering RefSeq
match_pos_refseq <- match(total_df$Sequence, refseq_tryptic, nomatch=0)
filtered_df <- total_df[match_pos_refseq == 0,]

message(sprintf("Peptides after RefSeq filtering: %d", nrow(filtered_df)))

write.table(filtered_df, file="./1_basic_analysis/filtered_data/identification_non_cano.txt", 
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
writexl::write_xlsx(filtered_df, "./1_basic_analysis/filtered_data/identification_non_cano.xlsx")

# UniProt Filtering
message("\nStarting UniProt filtering...")
message("Loading UniProt database...")
Uniprot <- read_excel("./Uni_2411_20428.xlsx", sheet=1)

# UniProt sequence to tryptic peptide
message("Processing UniProt sequences...")
pb <- progress_bar$new(
  format = "Processing UniProt sequences [:bar] :percent",
  total = length(Uniprot$Sequence),
  clear = FALSE,
  width = 60
)

Uniprot_seq <- unlist(lapply(Uniprot$Sequence, function(x){
  pb$tick()
  tmp_k <- gsub("K", "K@", x)
  tmp_r <- gsub("R", "R@", tmp_k)
  return(tmp_r)
}))

message("Processing proline sites...")
Uniprot_seq_proline <- lapply(Uniprot_seq, function(x){
  tmp_k <- gsub("K@P", "KP", x)
  tmp_r <- gsub("R@P", "RP", tmp_k)
})

message("Generating tryptic peptides...")
Uniprot_tryptic <- unlist(lapply(Uniprot_seq_proline, function(x){
  strsplit(x, "@")
}))

Uniprot_tryptic <- unique(Uniprot_tryptic)

# RefSeq data + UniProt Filtering
message("\nApplying UniProt filtering...")
match_pos_uniprot <- match(filtered_df$Sequence, Uniprot_tryptic, nomatch=0)
dual_filtered_df <- filtered_df[match_pos_uniprot == 0,]

message(sprintf("Final peptides after dual filtering: %d", nrow(dual_filtered_df)))

message("\nSaving filtered results...")
write.table(dual_filtered_df, file="./1_basic_analysis/filtered_data/identification_non_cano_dual.txt", 
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
writexl::write_xlsx(dual_filtered_df, 
                    "./1_basic_analysis/filtered_data/identification_non_cano_dual.xlsx")

# Filtering Summary
filtering_summary <- data.frame(
  Stage = c("Initial", "After RefSeq", "After UniProt"),


  
  Peptides = c(nrow(total_df), nrow(filtered_df), nrow(dual_filtered_df)),
  Filtered_Out = c(0, nrow(total_df) - nrow(filtered_df), 
                   nrow(filtered_df) - nrow(dual_filtered_df))
)

write.csv(filtering_summary, 
          "./1_basic_analysis/filtered_data/filtering_summary.csv", 
          row.names=FALSE)

message("\nDatabase filtering completed.")
