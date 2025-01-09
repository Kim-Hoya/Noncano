# =============================================
# 1. Setup and Package Installation
# =============================================
# Function to install and load required packages
install_and_load_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  invisible(lapply(packages, library, character.only = TRUE))
}

# Required packages
required_packages <- c(
  "tidyverse", "readxl", "writexl", "tximport", "biomaRt", 
  "stringr", "BiocManager"
)

# Install BioConductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bioc_packages <- c("tximport", "biomaRt")
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc)) BiocManager::install(new_bioc)

# Load all packages
install_and_load_packages(required_packages)

# =============================================
# 2. Create Directory Structure
# =============================================
dirs <- c(
  "results/2_mapping_analysis/ensembl",
  "results/2_mapping_analysis/ncbi"
)

for(dir in dirs) {
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# =============================================
# 3. Data Import Functions
# =============================================
import_initial_data <- function(peptide_file, peptide_info_file) {
  # Read data files
  peptides <- read_excel(peptide_file)
  peptide_info <- read_excel(peptide_info_file)
  
  # Initial data validation
  if(nrow(peptides) != 2448) {
    warning("Warning: Unexpected number of peptides")
  }
  
  list(
    peptides = peptides,
    peptide_info = peptide_info
  )
}

# =============================================
# 4. Peptide Data Processing
# =============================================
process_peptides <- function(peptides, peptide_info) {
  processed <- peptides %>%
    mutate(
      accession_clean = gsub("\\..*$", "", accession),
      refseq_type = case_when(
        grepl("^NM_", accession) ~ "NM",
        grepl("^NR_", accession) ~ "NR",
        grepl("^XM_", accession) ~ "XM",
        grepl("^XR_", accession) ~ "XR",
        TRUE ~ "other"
      )
    ) %>%
    left_join(
      peptide_info %>% 
        select(Sequence, position, locations, peptide_type, peptide_subtype),
      by = "Sequence"
    )
  
  # Validation
  refseq_counts <- table(processed$refseq_type)
  if(sum(refseq_counts) != 2448) {
    warning("Warning: Lost peptides during processing")
  }
  
  processed
}

# =============================================
# 5. Ensembl Mapping Function
# =============================================
get_ensembl_mapping <- function(refseq_ids) {
  # Connect to Ensembl with multiple mirror attempts
  mirrors <- c("asia", "useast", "uswest")
  mart <- NULL
  
  for(mirror in mirrors) {
    tryCatch({
      mart <- useEnsembl(biomart = "ensembl", 
                         dataset = "hsapiens_gene_ensembl",
                         mirror = mirror)
      break
    }, error = function(e) {
      warning(paste("Mirror", mirror, "failed:", e$message))
    })
  }
  
  if(is.null(mart)) {
    stop("Failed to connect to any Ensembl mirror")
  }
  
  # Get different types of RefSeq mappings
  get_mapping <- function(ids, type) {
    attributes <- c(
      "ensembl_transcript_id",
      paste0("refseq_", type),
      "ensembl_gene_id",
      "external_gene_name"
    )
    
    getBM(
      attributes = attributes,
      filters = paste0("refseq_", type),
      values = ids,
      mart = mart
    ) %>%
      rename(refseq_mrna = paste0("refseq_", type))
  }
  
  # Get mappings for each RefSeq type
  mappings <- list(
    nm = get_mapping(refseq_ids[grep("^NM_", refseq_ids)], "mrna"),
    nr = get_mapping(refseq_ids[grep("^NR_", refseq_ids)], "ncrna"),
    xm = get_mapping(refseq_ids[grep("^XM_", refseq_ids)], "mrna_predicted"),
    xr = get_mapping(refseq_ids[grep("^XR_", refseq_ids)], "ncrna_predicted")
  )
  
  # Combine all mappings
  all_mapping <- bind_rows(
    mappings$nm %>% mutate(refseq_type = "NM"),
    mappings$nr %>% mutate(refseq_type = "NR"),
    mappings$xm %>% mutate(refseq_type = "XM"),
    mappings$xr %>% mutate(refseq_type = "XR")
  )
  
  all_mapping
}

# =============================================
# 6. RNA-seq Data Processing
# =============================================
process_rna_seq <- function(salmon_dir, tx2gene) {
  # Get Salmon file paths
  salmon_files <- list.files(salmon_dir, 
                             pattern = "quant.sf$", 
                             recursive = TRUE, 
                             full.names = TRUE)
  
  if(length(salmon_files) == 0) {
    stop("No Salmon quantification files found")
  }
  
  # Name samples based on directory names
  names(salmon_files) <- basename(dirname(salmon_files))
  
  # Import RNA-seq data
  txi <- tximport(salmon_files, 
                  type = "salmon",
                  tx2gene = tx2gene %>% select(-refseq_type),
                  ignoreTxVersion = TRUE)
  
  # Convert to data frame
  as.data.frame(txi$abundance) %>%
    rownames_to_column("gene_id")
}

# =============================================
# 7. Fix Mapping Issues
# =============================================
fix_peptide_mapping <- function(peptide_gene_mapping) {
  # Check for duplicates
  duplicate_check <- peptide_gene_mapping %>%
    group_by(Sequence) %>%
    summarise(
      mapping_count = n(),
      gene_count = n_distinct(ensembl_gene_id[!is.na(ensembl_gene_id)])
    ) %>%
    filter(mapping_count > 1)
  
  # Fix mapping by keeping one mapping per peptide per refseq type
  fixed_mapping <- peptide_gene_mapping %>%
    group_by(Sequence, refseq_type) %>%
    slice(1) %>%
    ungroup()
  
  # Generate mapping summary
  mapping_summary <- fixed_mapping %>%
    group_by(refseq_type) %>%
    summarise(
      total_peptides = n(),
      mapped_to_ensembl = sum(!is.na(ensembl_gene_id)),
      mapping_rate = sprintf("%.1f%%", mean(!is.na(ensembl_gene_id)) * 100)
    )
  
  list(
    fixed_mapping = fixed_mapping,
    mapping_summary = mapping_summary,
    duplicate_info = duplicate_check
  )
}

# =============================================
# 8. Main Pipeline Function
# =============================================
run_mapping_pipeline <- function(peptide_file, peptide_info_file, salmon_dir) {
  # 1. Import data
  data <- import_initial_data(peptide_file, peptide_info_file)
  
  # 2. Process peptides
  processed_peptides <- process_peptides(data$peptides, data$peptide_info)
  
  # 3. Get Ensembl mapping
  ensembl_mapping <- get_ensembl_mapping(processed_peptides$accession_clean)
  
  # 4. Create tx2gene mapping
  tx2gene <- ensembl_mapping %>%
    select(ensembl_transcript_id, ensembl_gene_id, refseq_type) %>%
    distinct() %>%
    group_by(ensembl_transcript_id) %>%
    slice(1)
  
  # 5. Process RNA-seq data
  rna_expression <- process_rna_seq(salmon_dir, tx2gene)
  
  # 6. Create and fix peptide-gene mapping
  peptide_gene_mapping <- processed_peptides %>%
    select(Sequence, accession_clean, refseq_type) %>%
    left_join(
      ensembl_mapping %>% 
        select(refseq_mrna, ensembl_gene_id, refseq_type),
      by = c("accession_clean" = "refseq_mrna")
    )
  
  mapping_results <- fix_peptide_mapping(peptide_gene_mapping)
  
  # Return results
  return(list(
    processed_peptides = processed_peptides,
    ensembl_mapping = ensembl_mapping,
    rna_expression = rna_expression,
    mapping_results = mapping_results,
    tx2gene = tx2gene
  ))
}


# =============================================
# Export Mapped Results
# =============================================
# 1. Export mapped peptide information
mapped_peptide_info <- fixed_peptide_mapping %>%
  left_join(peptide_info, by = "Sequence") %>%
  select(
    Sequence, accession_clean, refseq_type, ensembl_gene_id,
    position, locations, peptide_type, peptide_subtype
  )

# 2. Export RNA-seq expression for mapped transcripts
mapped_rna_expression <- rna_expression %>%
  filter(gene_id %in% mapped_peptide_info$ensembl_gene_id) %>%
  arrange(gene_id)

# 3. Export full results with mapping information
mapped_results <- list(
  "Mapping_Summary" = updated_mapping_summary,
  "Mapped_Peptides" = mapped_peptide_info,
  "RNA_Expression" = mapped_rna_expression
)

# Create directory if it doesn't exist
dir.create("results/2_mapping_analysis", recursive = TRUE, showWarnings = FALSE)

# Save results
write_xlsx(mapped_results, "results/2_mapping_analysis/mapped_results.xlsx")

print("Mapping results exported successfully.")
print(paste("Number of mapped peptides:", nrow(mapped_peptide_info)))
print(paste("Number of mapped genes with RNA-seq data:", nrow(mapped_rna_expression)))











# =============================================
# Library Loading and Setup
# =============================================
library(tidyverse)
library(readxl)
library(writexl)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(scales)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(survival)
library(survminer)

# Set working directory
setwd("C:/Users/KDH/Desktop/Noncano/241202_rna")

# Create directories
dirs <- c(
  "results/mapping_analysis",
  "results/expression_analysis",
  "results/distribution_analysis",
  "results/correlation_analysis",
  "results/differential_analysis",
  "results/clinical_analysis",
  "results/tables"
)

for(dir in dirs) {
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# =============================================
# Data Import
# =============================================
# Load main data
mapped_results <- read_excel("results/2_mapping_analysis/mapped_results.xlsx", 
                             sheet = "RNA_Expression")
mapped_peptides <- read_excel("results/2_mapping_analysis/mapped_results.xlsx", 
                              sheet = "Mapped_Peptides")
mapping_summary <- read_excel("results/2_mapping_analysis/mapped_results.xlsx", 
                              sheet = "Mapping_Summary") %>%
  mutate(mapping_rate = as.numeric(sub("%", "", mapping_rate)))

# Load clinical data
clinical_info <- read_excel("clinical_info_with_clusters.xlsx")

# =============================================
# Color Settings
# =============================================
# Location colors
location_colors <- c(
  "3' UTR" = "#663399",         
  "N-terminal extension" = "#FFD700",
  "5' UTR" = "#90EE90",        
  "CDS-internal" = "#000080",  
  "lncRNA" = "#20B2AA"        
)

# RefSeq colors
refseq_colors <- c(
  "NM" = "#66C2A5",
  "NR" = "#FC8D62",
  "XM" = "#8DA0CB",
  "XR" = "#E78AC3"
)

# Expression level colors
expression_colors <- c(
  "Not expressed" = "#D3D3D3",
  "Low" = "#FEE0D2",
  "Medium" = "#FC9272",
  "High" = "#DE2D26"
)

# Clinical annotation colors
clinical_colors <- list(
  DX = c("AD" = "#c49c94", "SC" = "#f7b6d2", "MA" = "#c7c7c7", 
         "NC" = "#dbdb8d", "Others" = "#9edae5"),
  Smoking = c("Heavy" = "#1f77b4", "Light" = "#ff7f0e", "None" = "#2ca02c", 
              "Unknown" = "#cccccc"),
  OS_type = c("Short" = "#9467bd", "Alive" = "#8c564b", "Long" = "#e377c2"),
  Celltype.based.subtype = c("Cold_Immunogram" = "#98df8a", 
                             "Hot_Immunogram" = "#ff9896", 
                             "Unknown" = "#cccccc"),
  MO.subtype = c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3", 
                 "4" = "#e7298a", "5" = "#66a61e")
)

# =============================================
# Data Processing
# =============================================
# Process RNA expression data
rna_expression_long <- mapped_results %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  mutate(
    sample_type = ifelse(grepl("\\.T$", sample), "Tumor", "Normal"),
    patient = sub("\\.(T|N)$", "", sample),
    log2_expression = log2(expression + 1)
  )

# Calculate location mapping statistics
location_mapping_stats <- mapped_peptides %>%
  group_by(locations) %>%
  summarise(
    total_count = n(),
    mapped_count = sum(!is.na(ensembl_gene_id)),
    mapping_rate = (mapped_count/total_count) * 100
  )

# Calculate T/N fold changes
tn_fold_changes <- mapped_results %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  mutate(
    sample_type = ifelse(grepl("\\.T$", sample), "Tumor", "Normal"),
    patient = sub("\\.(T|N)$", "", sample)
  ) %>%
  group_by(patient) %>%
  filter(n_distinct(sample_type) == 2) %>%  # Only keep paired samples
  ungroup() %>%
  pivot_wider(
    id_cols = c(gene_id, patient),
    names_from = sample_type,
    values_from = expression
  ) %>%
  mutate(log2FC = log2((Tumor + 0.1)/(Normal + 0.1))) %>%
  # Now join with peptide information
  left_join(
    mapped_peptides %>% select(gene_id = ensembl_gene_id, locations),
    by = "gene_id"
  )

# Process expression data for analysis
expression_summary <- rna_expression_long %>%
  left_join(mapped_peptides %>% 
              select(gene_id = ensembl_gene_id, locations),
            by = "gene_id") %>%
  mutate(
    expr_level = case_when(
      expression == 0 ~ "Not expressed",
      expression < quantile(expression[expression > 0], 0.25) ~ "Low",
      expression < quantile(expression[expression > 0], 0.75) ~ "Medium",
      TRUE ~ "High"
    )
  )

# =============================================
# Visualization Functions
# =============================================
# 1. Mapping Analysis Functions
create_mapping_rate_plot <- function(data, colors, title) {
  ggplot(data, aes(x = reorder(refseq_type, -mapping_rate), 
                   y = mapping_rate, 
                   fill = refseq_type)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = title,
      x = "RefSeq Type",
      y = "Mapping Rate (%)"
    ) +
    scale_y_continuous(limits = c(0, 100))
}

create_location_mapping_plot <- function(data, colors) {
  ggplot(data, aes(x = reorder(locations, -mapping_rate), 
                   y = mapping_rate, 
                   fill = locations)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = "RNA-Seq Mapping Rate by Location",
      x = "Location",
      y = "Mapping Rate (%)"
    ) +
    scale_y_continuous(limits = c(0, 100))
}

# 2. Expression Analysis Functions
create_expression_boxplot <- function(data) {
  ggplot(data, aes(x = locations, y = log2_expression, fill = locations)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    scale_fill_manual(values = location_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = "RNA-Seq Expression Distribution by Location",
      x = "Location",
      y = "log2(Expression + 1)"
    )
}

create_expression_density_plot <- function(data) {
  ggplot(data, aes(x = log2_expression, fill = locations)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = location_colors) +
    theme_minimal() +
    labs(
      title = "Expression Density by Location",
      x = "log2(Expression + 1)",
      y = "Density"
    )
}

create_expression_category_plot <- function(data) {
  ggplot(data, aes(x = locations, fill = expr_level)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = expression_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = "Expression Level Distribution by Location",
      x = "Location",
      y = "Proportion",
      fill = "Expression Level"
    )
}

# 3. T/N Analysis Functions
create_tn_boxplot <- function(data) {
  ggplot(data, aes(x = locations, y = log2FC, fill = locations)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    scale_fill_manual(values = location_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"
    ) +
    labs(
      title = "Tumor/Normal Expression Changes by Location",
      x = "Location",
      y = "log2(Fold Change)"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")
}

# Create distribution plot for expressed genes
create_expression_distribution_plots <- function(data, output_dir) {
  # Calculate percentage of expressed genes per sample
  sample_stats <- data %>%
    group_by(sample) %>%
    summarise(
      expressed = sum(expression > 0),
      total = n(),
      percent = (expressed/total) * 100,
      sample_type = unique(sample_type)
    )
  
  # 1. All samples with color
  p1 <- ggplot(sample_stats, aes(x = sample, y = percent, fill = sample_type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Normal" = "#4575B4", "Tumor" = "#D73027")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - All Samples",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # 2. All samples without color
  p2 <- ggplot(sample_stats, aes(x = sample, y = percent)) +
    geom_bar(stat = "identity", fill = "grey50") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - All Samples (No Color)",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # 3. Tumor samples only
  p3 <- ggplot(filter(sample_stats, sample_type == "Tumor"), 
               aes(x = sample, y = percent)) +
    geom_bar(stat = "identity", fill = "#D73027") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - Tumor Samples",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # 4. Normal samples only
  p4 <- ggplot(filter(sample_stats, sample_type == "Normal"), 
               aes(x = sample, y = percent)) +
    geom_bar(stat = "identity", fill = "#4575B4") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - Normal Samples",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # Save all plots
  plots <- list(
    all_color = p1,
    all_plain = p2,
    tumor = p3,
    normal = p4
  )
  
  for(name in names(plots)) {
    ggsave(paste0(output_dir, "/distribution_", name, ".pdf"),
           plots[[name]], width = 12, height = 6)
    ggsave(paste0(output_dir, "/distribution_", name, ".png"),
           plots[[name]], width = 12, height = 6, dpi = 300)
    ggsave(paste0(output_dir, "/distribution_", name, ".svg"),
           plots[[name]], width = 12, height = 6)
  }
  
  return(plots)
}

# =============================================
# Create and Save Visualizations
# =============================================

# 1. Mapping Analysis
# ------------------
# Create mapping rate plots
refseq_mapping_plot <- create_mapping_rate_plot(
  mapping_summary,
  refseq_colors,
  "RNA-Seq Mapping Rate by RefSeq Type"
)

location_mapping_plot <- create_location_mapping_plot(
  location_mapping_stats,
  location_colors
)

# Save mapping plots with tables
mapping_plots <- list(
  refseq_mapping = list(plot = refseq_mapping_plot, data = mapping_summary),
  location_mapping = list(plot = location_mapping_plot, data = location_mapping_stats)
)

for(name in names(mapping_plots)) {
  p <- mapping_plots[[name]]$plot
  table_grob <- tableGrob(mapping_plots[[name]]$data %>%
                            select(-total_peptides, -mapped_to_ensembl), 
                          rows = NULL)
  combined_plot <- plot_grid(p, table_grob, ncol = 1, rel_heights = c(0.7, 0.3))
  
  ggsave(paste0("results/mapping_analysis/", name, ".pdf"),
         combined_plot, width = 8, height = 10, device = cairo_pdf)
  ggsave(paste0("results/mapping_analysis/", name, ".png"),
         combined_plot, width = 8, height = 10, dpi = 300)
  ggsave(paste0("results/mapping_analysis/", name, ".svg"),
         combined_plot, width = 8, height = 10)
}

# 2. Expression Analysis
# ---------------------
# Process expression data by sample type
expression_by_type <- expression_summary %>%
  filter(!is.na(locations)) %>%
  group_by(locations, sample_type) %>%
  summarise(
    mean_expression = mean(log2_expression, na.rm = TRUE),
    median_expression = median(log2_expression, na.rm = TRUE),
    sd_expression = sd(log2_expression, na.rm = TRUE),
    .groups = 'drop'
  )

# Create expression plots
expression_plots <- list(
  boxplot = create_expression_boxplot(expression_summary),
  density = create_expression_density_plot(expression_summary),
  category = create_expression_category_plot(expression_summary)
)

# Save expression plots
for(name in names(expression_plots)) {
  ggsave(paste0("results/expression_analysis/expression_", name, ".pdf"),
         expression_plots[[name]], width = 10, height = 8, device = cairo_pdf)
  ggsave(paste0("results/expression_analysis/expression_", name, ".png"),
         expression_plots[[name]], width = 10, height = 8, dpi = 300)
  ggsave(paste0("results/expression_analysis/expression_", name, ".svg"),
         expression_plots[[name]], width = 10, height = 8)
}

# 3. Distribution Analysis
# -----------------------
# Create and save distribution plots
distribution_plots <- create_expression_distribution_plots(
  expression_summary,
  "results/distribution_analysis"
)

# 4. T/N Analysis
# --------------
# Create and save T/N fold change plot
tn_plot <- create_tn_boxplot(tn_fold_changes)

ggsave("results/differential_analysis/tn_foldchange.pdf",
       tn_plot, width = 10, height = 8, device = cairo_pdf)
ggsave("results/differential_analysis/tn_foldchange.png",
       tn_plot, width = 10, height = 8, dpi = 300)
ggsave("results/differential_analysis/tn_foldchange.svg",
       tn_plot, width = 10, height = 8)

# Calculate summary statistics for T/N analysis
tn_stats <- tn_fold_changes %>%
  group_by(locations) %>%
  summarise(
    mean_fc = mean(log2FC, na.rm = TRUE),
    median_fc = median(log2FC, na.rm = TRUE),
    sd_fc = sd(log2FC, na.rm = TRUE),
    up_regulated = sum(log2FC > 1, na.rm = TRUE),
    down_regulated = sum(log2FC < -1, na.rm = TRUE),
    total = n()
  ) %>%
  mutate(
    percent_up = (up_regulated/total) * 100,
    percent_down = (down_regulated/total) * 100
  )

# 5. Save Summary Tables
# --------------------
write_xlsx(list(
  Mapping_Summary = mapping_summary,
  Location_Mapping = location_mapping_stats,
  Expression_Summary = expression_by_type,
  TN_Analysis = tn_stats
), "results/tables/analysis_summary.xlsx")

print("Analysis completed. Check the results directory for all outputs.")

# Print summary of outputs
print("\nSummary of created files:")
for(dir in dirs) {
  if(dir.exists(dir)) {
    files <- list.files(dir, pattern = "\\.pdf|\\.png|\\.svg|\\.xlsx$")
    print(paste("\nFiles in", dir, ":"))
    print(files)
  }
}

