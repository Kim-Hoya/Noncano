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
setwd("C:/Users/starh/OneDrive/Desktop/Noncano/코드정리/Figure 2")

# Create directories
dirs <- c(
  "results/mapping_analysis",
  "results/expression_analysis",
  "results/distribution_analysis",
  "results/differential_analysis",
  "results/tables"
)

for(dir in dirs) {
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Create a consistent theme for all plots
create_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica", size = base_size),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size + 2),
      axis.text.y = element_text(size = base_size + 2),
      axis.title = element_text(size = base_size + 4, face = "bold"),
      plot.title = element_text(size = base_size + 6, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 2, face = "bold"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
}

# Save function for high-quality output
save_plot <- function(plot, filename, width = 8, height = 10) {
  # TIFF version
  ggsave(paste0(filename, ".tiff"),
         plot,
         width = width,
         height = height,
         dpi = 300,
         device = "tiff",
         compression = "lzw")
  
  # PDF version
  ggsave(paste0(filename, ".pdf"),
         plot,
         width = width,
         height = height,
         device = cairo_pdf)
  
  # PNG version
  ggsave(paste0(filename, ".png"),
         plot,
         width = width,
         height = height,
         dpi = 300)
}

# =============================================
# Data Import
# =============================================
# Load main data
mapped_results <- read_excel("mapped_results.xlsx", 
                             sheet = "RNA_Expression")
mapped_peptides <- read_excel("mapped_results.xlsx", 
                              sheet = "Mapped_Peptides")
mapping_summary <- read_excel("mapped_results.xlsx", 
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
  group_by(peptide_type) %>%
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
  filter(n_distinct(sample_type) == 2) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(gene_id, patient),
    names_from = sample_type,
    values_from = expression
  ) %>%
  mutate(log2FC = log2((Tumor + 0.1)/(Normal + 0.1))) %>%
  left_join(
    mapped_peptides %>% 
      select(gene_id = ensembl_gene_id, peptide_type) %>%
      distinct(),  # 중복 제거
    by = "gene_id",
    relationship = "many-to-many"  # 관계 명시
  )

# Process expression data for analysis
expression_summary <- rna_expression_long %>%
  left_join(mapped_peptides %>% 
              select(gene_id = ensembl_gene_id, peptide_type),
            by = "gene_id") %>%
  mutate(
    expr_level = case_when(
      expression == 0 ~ "Not expressed",
      expression < quantile(expression[expression > 0], 0.25) ~ "Low",
      expression < quantile(expression[expression > 0], 0.75) ~ "Medium",
      TRUE ~ "High",
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
    create_theme(base_size = 14) +
    theme(legend.position = "none") +
    labs(
      title = title,
      x = "RefSeq Type",
      y = "Mapping Rate (%)"
    ) +
    scale_y_continuous(limits = c(0, 100))
}

create_location_mapping_plot <- function(data, colors) {
  ggplot(data, aes(x = reorder(peptide_type, -mapping_rate), 
                   y = mapping_rate, 
                   fill = peptide_type)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = colors) +
    create_theme(base_size = 14) +
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
  ggplot(data, aes(x = peptide_type, y = log2_expression, fill = peptide_type)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    scale_fill_manual(values = location_colors) +
    create_theme(base_size = 14) +
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
  ggplot(data, aes(x = log2_expression, fill = peptide_type)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = location_colors) +
    create_theme(base_size = 14) +
    labs(
      title = "Expression Density by Location",
      x = "log2(Expression + 1)",
      y = "Density"
    )
}

create_expression_category_plot <- function(data) {
  ggplot(data, aes(x = peptide_type, fill = expr_level)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = expression_colors) +
    create_theme(base_size = 14) +
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
  ggplot(data, aes(x = peptide_type, y = log2FC, fill = peptide_type)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
    scale_fill_manual(values = location_colors) +
    create_theme(base_size = 14) +
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
    create_theme(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - All Samples",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # 2. All samples without color
  p2 <- ggplot(sample_stats, aes(x = sample, y = percent)) +
    geom_bar(stat = "identity", fill = "grey50") +
    create_theme(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - All Samples (No Color)",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # 3. Tumor samples only
  p3 <- ggplot(filter(sample_stats, sample_type == "Tumor"), 
               aes(x = sample, y = percent)) +
    geom_bar(stat = "identity", fill = "#D73027") +
    create_theme(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Expressed Genes Distribution - Tumor Samples",
         x = "Sample", y = "Percentage of Expressed Genes (%)")
  
  # 4. Normal samples only
  p4 <- ggplot(filter(sample_stats, sample_type == "Normal"), 
               aes(x = sample, y = percent)) +
    geom_bar(stat = "identity", fill = "#4575B4") +
    create_theme(base_size = 14) +
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


# location_mapping_stats 데이터 확인을 위한 코드 추가
print("Location mapping stats:")
print(str(location_mapping_stats))
print(head(location_mapping_stats))

# Save plots
# 수정된 저장 코드
for(name in names(mapping_plots)) {
  p <- mapping_plots[[name]]$plot
  # 데이터프레임 구조에 따라 select 부분 수정
  table_data <- mapping_plots[[name]]$data %>%
    select(everything())  # 모든 컬럼 선택
  table_grob <- tableGrob(table_data, rows = NULL)
  combined_plot <- plot_grid(p, table_grob, ncol = 1, rel_heights = c(0.7, 0.3))
  
  save_plot(combined_plot, 
            paste0("results/mapping_analysis/", name),
            width = 8, 
            height = 10)
}


# 2. Expression Analysis
# ---------------------
# Process expression data by sample type
expression_by_type <- expression_summary %>%
  filter(!is.na(peptide_type)) %>%
  group_by(peptide_type, sample_type) %>%
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


# 단순화된 저장 함수
save_plot <- function(plot, filename, width = 10, height = 8) {
  # 각 형식별로 순차적 저장
  ggsave(paste0(filename, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(filename, ".png"), plot, width = width, height = height, dpi = 300)
  ggsave(paste0(filename, ".tiff"), plot, width = width, height = height, 
         device = "tiff", compression = "lzw", dpi = 300)
}

# 사용
for(name in names(expression_plots)) {
  save_plot(
    expression_plots[[name]], 
    paste0("results/expression_analysis/expression_", name)
  )
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

# 사용
for(name in names(tn_plot)) {
  save_plot(
    expression_plots[[name]], 
    paste0("results/differential_analysis/tn_foldchange_", name)
  )
}

# Calculate summary statistics for T/N analysis
tn_stats <- tn_fold_changes %>%
  group_by(peptide_type) %>%
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


# Save Summary Tables
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
    files <- list.files(dir, pattern = "\\.pdf|\\.png|\\.svg|\\.tiff|\\.xlsx$")
    print(paste("\nFiles in", dir, ":"))
    print(files)
  }
}



