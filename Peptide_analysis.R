# Part 1: Data Analysis and Export
library(tidyverse)
library(dplyr)
library(writexl)
library(impute)
library(mixtools)
library(stats)

setwd("C:/Users/KDH/Desktop/Noncano/Code뭉치/Figure 3")

# 데이터 로드 함수
load_data <- function() {
  translation_results <- read_excel("Enhanced_Peptide_analysis_TMT(2448).xlsx")
  
  tmt_data <- read.table("identification_non_cano_dual.txt", 
                         sep="\t", header=TRUE, stringsAsFactors = FALSE)
  
  sample_info <- read.table("sample.info.txt", 
                            sep="\t", header=TRUE, stringsAsFactors = FALSE)
  
  return(list(
    translation = translation_results,
    tmt = tmt_data,
    samples = sample_info
  ))
}

# NA 필터링 함수
filter_na_threshold <- function(data, threshold=70) {
  quant_cols <- c(
    grep("^RE\\.P\\d+\\.T$", colnames(data), value = TRUE),
    grep("^RE\\.P\\d+\\.N$", colnames(data), value = TRUE)
  )
  
  if(length(quant_cols) == 0) {
    stop("No quantitative data columns found")
  }
  
  na_ratio <- rowSums(is.na(data[, quant_cols])) / length(quant_cols) * 100
  filtered_data <- data[na_ratio < threshold, ]
  
  message(sprintf("Found %.0f T/N sample pairs", length(quant_cols)/2))
  message(sprintf("Retained %d out of %d rows after %d%% NA filtering", 
                  nrow(filtered_data), nrow(data), threshold))
  
  return(filtered_data)
}

# 정규화 함수
normalize_tmt_data <- function(data) {
  quant_cols <- c(
    grep("^RE\\.P\\d+\\.T$", colnames(data), value = TRUE),
    grep("^RE\\.P\\d+\\.N$", colnames(data), value = TRUE)
  )
  
  # Log2 변환된 데이터 생성
  log_data <- data
  log_matrix <- log2(as.matrix(data[, quant_cols]))
  log_matrix[is.infinite(log_matrix)] <- NA
  
  # 각 컬럼별 정규화
  norm_matrix <- matrix(NA, nrow=nrow(log_matrix), ncol=ncol(log_matrix))
  colnames(norm_matrix) <- colnames(log_matrix)
  
  for(k in 1:ncol(log_matrix)) {
    a <- na.omit(log_matrix[,k])
    b <- bw.SJ(a)
    d <- stats::density(a, bw=b, kernel="gaussian")
    mode <- d$x[which.max(d$y)]
    
    tryCatch({
      set.seed(1234)
      model <- normalmixEM(a, mu=c(mode,mode), sigma=NULL, k=2, maxit = 10000)
      j = which.min(model$sigma)
      sigma = model$sigma[j]
      mu = mode
    }, error = function(e) {
      sigma = sd(a, na.rm=TRUE)
      mu = mean(a, na.rm=TRUE)
    })
    
    norm_matrix[,k] <- (log_matrix[,k] - mu) / sigma
  }
  
  # 결과 데이터프레임 생성
  norm_data <- data
  norm_data[,quant_cols] <- norm_matrix
  
  return(list(
    log_data = log_data,
    norm_data = norm_data
  ))
}

# KNN imputation 함수
perform_knn_imputation <- function(data) {
  quant_cols <- c(
    grep("^RE\\.P\\d+\\.T$", colnames(data), value = TRUE),
    grep("^RE\\.P\\d+\\.N$", colnames(data), value = TRUE)
  )
  
  # 정량 데이터 행렬 추출
  expr_matrix <- as.matrix(data[, quant_cols])
  
  # KNN imputation 수행 (k=10)
  message("Performing KNN imputation...")
  expr_matrix_imputed <- impute.knn(expr_matrix, k=10)$data
  
  # 결과 데이터프레임 생성
  imputed_data <- data
  imputed_data[, quant_cols] <- expr_matrix_imputed
  
  return(imputed_data)
}

# 통계 분석 함수
perform_statistical_analysis <- function(data) {
  tumor_cols <- grep("^RE\\.P\\d+\\.T$", colnames(data), value = TRUE)
  normal_cols <- grep("^RE\\.P\\d+\\.N$", colnames(data), value = TRUE)
  
  stats_results <- data.frame(
    Sequence = data$Sequence,
    peptide_type = data$peptide_type,
    locations = data$locations,
    p_value = NA,
    fold_change = NA,
    stringsAsFactors = FALSE
  )
  
  for(i in 1:nrow(data)) {
    tumor_vals <- as.numeric(unlist(data[i, tumor_cols]))
    normal_vals <- as.numeric(unlist(data[i, normal_cols]))
    
    if(sum(!is.na(normal_vals)) > 1 && sum(!is.na(tumor_vals)) > 1) {
      t_test <- t.test(tumor_vals, normal_vals, paired = FALSE)
      stats_results$p_value[i] <- t_test$p.value
      
      stats_results$fold_change[i] <- 
        (2^mean(tumor_vals, na.rm=TRUE)) / (2^mean(normal_vals, na.rm=TRUE))
    }
  }
  
  return(stats_results)
}

# 메인 분석 실행
main_analysis <- function() {
  # 데이터 로드
  data <- load_data()
  
  # Noncanonical peptide 필터링
  noncanonical_data <- data$tmt %>%
    inner_join(
      data$translation %>%
        filter(!peptide_type %in% c("Canonical", "Not found", "Error")) %>%
        select(Sequence, peptide_type, locations) %>%
        distinct(Sequence, .keep_all = TRUE),
      by = "Sequence"
    )
  
  # NA 필터링 (70%)
  filtered_data <- filter_na_threshold(noncanonical_data)
  
  # 정규화
  norm_results <- normalize_tmt_data(filtered_data)
  
  # 정규화된 데이터에 대해 imputation 수행
  imputed_norm_data <- perform_knn_imputation(norm_results$norm_data)
  
  # 통계 분석 (정규화 데이터 사용, imputation 전)
  stats_results <- perform_statistical_analysis(norm_results$norm_data)
  
  # 결과 저장
  output_dir <- "3_figure_data"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # DEP 데이터 저장
  dep_dir <- file.path(output_dir, "DEP")
  dir.create(file.path(dep_dir, "Up"), recursive = TRUE)
  dir.create(file.path(dep_dir, "Down"), recursive = TRUE)
  
  # Up/Down regulated peptides 분류 및 저장
  up_peptides <- subset(stats_results, 
                        p_value < 0.05 & fold_change > 1.5) %>%
    arrange(p_value)
  
  down_peptides <- subset(stats_results, 
                          p_value < 0.05 & fold_change < 1/1.5) %>%
    arrange(p_value)
  
  write_xlsx(up_peptides, path = file.path(dep_dir, "Up", "Up.xlsx"))
  write_xlsx(down_peptides, path = file.path(dep_dir, "Down", "Down.xlsx"))
  
  # 전체 결과 저장
  write_xlsx(list(
    Raw_Log2_Data = norm_results$log_data,
    Normalized_Data = norm_results$norm_data,
    Normalized_Imputed_Data = imputed_norm_data,
    Statistics = stats_results
  ), path = file.path(output_dir, "filtered_peptide_data.xlsx"))
  
  return(list(
    log_data = norm_results$log_data,
    norm_data = norm_results$norm_data,
    imputed_data = imputed_norm_data,
    stats = stats_results
  ))
}

# 분석 실행
results <- main_analysis()


# Part 2: Figure Generation
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)
library(dplyr)
library(readxl)

# Base theme 설정
base_theme <- theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    legend.key.size = unit(0.5, "lines")
  )

# PCA plot 생성
create_pca_plot <- function(imputed_data) {
  quant_cols <- c(
    grep("^RE\\.P\\d+\\.T$", colnames(imputed_data), value = TRUE),
    grep("^RE\\.P\\d+\\.N$", colnames(imputed_data), value = TRUE)
  )
  
  pca_result <- prcomp(t(imputed_data[, quant_cols]), scale. = TRUE)
  
  pca_data <- data.frame(
    PC1 = pca_result$x[,1],
    PC2 = pca_result$x[,2],
    Type = factor(ifelse(grepl("\\.T$", rownames(pca_result$x)), "Tumor", "NAT"))
  )
  
  var_explained <- round(pca_result$sdev^2/sum(pca_result$sdev^2)*100, 1)
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Type)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95) +
    scale_color_manual(values = c("NAT" = "#00BA38", "Tumor" = "#F8766D")) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = "PCA Plot (Tumor vs NAT)") +
    base_theme
  
  return(p)
}

# Heatmap 생성
create_heatmap <- function(imputed_data, stats_results) {
  significant_peptides <- stats_results$Sequence[
    stats_results$p_value < 0.05 & 
      (stats_results$fold_change > 1.5 | stats_results$fold_change < 1/1.5)
  ]
  
  quant_cols <- c(
    grep("^RE\\.P\\d+\\.T$", colnames(imputed_data), value = TRUE),
    grep("^RE\\.P\\d+\\.N$", colnames(imputed_data), value = TRUE)
  )
  
  mat <- as.matrix(imputed_data[imputed_data$Sequence %in% significant_peptides, quant_cols])
  rownames(mat) <- imputed_data$Sequence[imputed_data$Sequence %in% significant_peptides]
  
  mat_scaled <- t(scale(t(mat)))
  
  column_anno <- HeatmapAnnotation(
    Type = ifelse(grepl("\\.T$", colnames(mat)), "Tumor", "NAT"),
    col = list(Type = c("Tumor" = "#F8766D", "NAT" = "#00BA38")),
    show_legend = TRUE
  )
  
  ht <- Heatmap(mat_scaled,
                name = "Z-score",
                top_annotation = column_anno,
                show_row_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_column_names = FALSE,
                column_names_gp = gpar(fontsize = 8),
                col = colorRamp2(c(-2, 0, 2), 
                                 c("#4575B4", "white", "#D73027")))
  
  return(ht)
}

# Volcano plot 생성
create_volcano_plot <- function(stats_results) {
  stats_results$log2FC <- log2(stats_results$fold_change)
  
  up_count <- sum(stats_results$p_value < 0.05 & stats_results$fold_change > 1.5, 
                  na.rm = TRUE)
  down_count <- sum(stats_results$p_value < 0.05 & stats_results$fold_change < 1/1.5, 
                    na.rm = TRUE)
  
  p <- ggplot(stats_results, aes(x = log2FC, y = -log10(p_value))) +
    geom_point(color = "grey80", alpha = 0.5, size = 0.8) +
    geom_point(
      data = subset(stats_results, 
                    p_value < 0.05 & fold_change < 1/1.5),
      color = "#4575B4", 
      alpha = 0.8,
      size = 0.8
    ) +
    geom_point(
      data = subset(stats_results, 
                    p_value < 0.05 & fold_change > 1.5),
      color = "#D73027", 
      alpha = 0.8,
      size = 0.8
    ) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)),
               linetype = "dashed", 
               color = "grey50", 
               alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", 
               color = "grey50", 
               alpha = 0.5) +
    annotate("text",
             x = max(stats_results$log2FC, na.rm=TRUE) - 0.5,
             y = max(-log10(stats_results$p_value), na.rm=TRUE),
             label = paste0("Up: ", up_count),
             color = "#D73027",
             size = 2) +
    annotate("text",
             x = min(stats_results$log2FC, na.rm=TRUE) + 0.5,
             y = max(-log10(stats_results$p_value), na.rm=TRUE),
             label = paste0("Down: ", down_count),
             color = "#4575B4",
             size = 2) +
    labs(x = "log2(Fold Change)",
         y = "-log10(p-value)") +
    base_theme
  
  return(p)
}

# Distribution plot 생성
create_distribution_plot <- function(data, group) {
  tryCatch({
    # NA 제거
    data_clean <- data %>%
      filter(!is.na(Peptide_Type))
    
    p <- ggplot(data_clean, aes(x = Peptide_Type, fill = Peptide_Type)) +
      geom_bar() +
      scale_fill_manual(values = c(
        "3' UTR" = "#663399",
        "N-terminal extension" = "#FFD700",
        "5' UTR" = "#90EE90",
        "CDS-internal" = "#000080",
        "lncRNA" = "#20B2AA"
      )) +
      labs(title = paste(group, "- Translation Location Distribution"),
           x = "Peptide Type", 
           y = "Count") +
      base_theme
    return(p)
  }, error = function(e) {
    message(sprintf("Error in distribution plot for %s: %s", group, e$message))
    return(NULL)
  })
}

# Survival plot 생성
create_survival_plot <- function(data, group) {
  tryCatch({
    p <- ggplot(filter(data, !is.na(Levels_P)), 
                aes(x = Peptide_Type, y = -log10(Levels_P), fill = Peptide_Type)) +
      geom_boxplot() +
      scale_fill_manual(values = c(
        "3' UTR" = "#663399",
        "N-terminal extension" = "#FFD700",
        "5' UTR" = "#90EE90",
        "CDS-internal" = "#000080",
        "lncRNA" = "#20B2AA"
      )) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      labs(title = paste(group, "- Survival Analysis P-values"),
           x = "Peptide Type",
           y = "-log10(P-value)") +
      base_theme
    return(p)
  }, error = function(e) {
    message(sprintf("Error in survival plot for %s: %s", group, e$message))
    return(NULL)
  })
}

# Figure 3 생성
create_figure3 <- function(log_data, norm_data, imputed_data, stats_results) {
  dir.create("figures", showWarnings = FALSE)
  
  # 개별 플롯 생성
  pca_plot <- create_pca_plot(imputed_data)
  heatmap_plot <- create_heatmap(imputed_data, stats_results)
  volcano_plot <- create_volcano_plot(stats_results)
  
  # 개별 그림 저장
  ggsave("figures/figure3a_pca.pdf", pca_plot,
         width = 3.5, height = 3.5, dpi = 300)
  
  pdf("figures/figure3b_heatmap.pdf", width = 5, height = 6)
  draw(heatmap_plot)
  dev.off()
  
  ggsave("figures/figure3c_volcano.pdf", volcano_plot,
         width = 3.5, height = 3.5, dpi = 300)
  
  # Combined figure
  pdf("figures/figure3_combined.pdf", width = 8.5, height = 11)
  lay <- rbind(c(1,2,2),
               c(3,2,2))
  
  grid.arrange(
    pca_plot,
    grid.grabExpr(draw(heatmap_plot, annotation_legend_side = "right")),
    volcano_plot,
    layout_matrix = lay,
    heights = c(1, 1),
    widths = c(1, 1, 0.1)
  )
  dev.off()
}

# Figure 4 생성
create_figure4 <- function(up_analysis, down_analysis) {
  dir.create("figures", showWarnings = FALSE)
  
  # Distribution과 Survival plot만 생성
  plots <- list(
    up_dist = create_distribution_plot(up_analysis, "Up"),
    up_survival = create_survival_plot(up_analysis, "Up"),
    down_dist = create_distribution_plot(down_analysis, "Down"),
    down_survival = create_survival_plot(down_analysis, "Down")
  )
  
  # Save individual plots
  plot_names <- c("dist", "survival")
  groups <- c("up", "down")
  
  for(group in groups) {
    for(plot_type in plot_names) {
      plot_obj <- plots[[paste0(group, "_", plot_type)]]
      if(!is.null(plot_obj)) {
        filename <- sprintf("figures/figure4_%s_%s.pdf", group, plot_type)
        ggsave(filename, plot_obj,
               width = 5, height = 4, dpi = 300)
      }
    }
  }
  
  # Combined figure
  valid_plots <- plots[!sapply(plots, is.null)]
  if(length(valid_plots) > 0) {
    pdf("figures/figure4_combined.pdf", width = 8.5, height = 11)
    grid.arrange(grobs = valid_plots, 
                 ncol = 2,
                 nrow = 2)
    dev.off()
  }
}

# Save significant survival results
save_survival_results <- function(data, group) {
  tryCatch({
    sig_survival <- data %>%
      filter(Levels_P < 0.05) %>%
      select(Sequence, Levels_P, Cox_HR, Cox_HR_Lower, Cox_HR_Upper, 
             Peptide_Type, kozak_score) %>%
      arrange(Levels_P)
    
    write.csv(sig_survival, 
              sprintf("figures/significant_survival_%s.csv", tolower(group)), 
              row.names = FALSE)
  }, error = function(e) {
    message(sprintf("Error saving survival results for %s: %s", group, e$message))
  })
}

# Main visualization 함수
run_visualization <- function() {
  # Figure 3 생성
  create_figure3(results$log_data, results$norm_data, 
                 results$imputed_data, results$stats)
  
  # Figure 4를 위한 데이터 로드
  peptide_info <- read_excel("Enhanced_Peptide_analysis_TMT(2448).xlsx")
  up_peptides <- read_excel("DEP/Up/Up.xlsx")
  down_peptides <- read_excel("DEP/Down/Down.xlsx")
  survival_data <- read_excel("survival_analysis_results.xlsx")
  
  # Up/Down analysis 데이터 준비
  up_analysis <- peptide_info %>%
    inner_join(up_peptides, by="Sequence") %>%
    left_join(survival_data, by=c("Sequence"="Peptide"))
  
  down_analysis <- peptide_info %>%
    inner_join(down_peptides, by="Sequence") %>%
    left_join(survival_data, by=c("Sequence"="Peptide"))
  
  # Figure 4 생성
  create_figure4(up_analysis, down_analysis)
}

# Significant survival results 저장
# Save significant survival results
save_survival_results <- function(data, group) {
  tryCatch({
    sig_survival <- data %>%
      filter(Levels_P < 0.05) %>%
      select(Sequence, Levels_P, Cox_HR, Cox_HR_Lower, Cox_HR_Upper, 
             Peptide_Type, kozak_score) %>%
      arrange(Levels_P)
    
    write.csv(sig_survival, 
              sprintf("figures/significant_survival_%s.csv", tolower(group)), 
              row.names = FALSE)
  }, error = function(e) {
    message(sprintf("Error saving survival results for %s: %s", group, e$message))
  })
}

# Main visualization 함수
run_visualization <- function() {
  # Figure 3 생성
  create_figure3(results$log_data, results$norm_data, 
                 results$imputed_data, results$stats)
  
  # Figure 4를 위한 데이터 로드
  peptide_info <- read_excel("Enhanced_Peptide_analysis_TMT(2448).xlsx")
  up_peptides <- read_excel("3_figure_data/DEP/Up/Up.xlsx")
  down_peptides <- read_excel("3_figure_data/DEP/Down/Down.xlsx")
  survival_data <- read_excel("survival_analysis_results.xlsx")
  
  # Up/Down analysis 데이터 준비
  up_analysis <- peptide_info %>%
    inner_join(up_peptides, by="Sequence") %>%
    left_join(survival_data, by=c("Sequence"="Peptide"))
  
  down_analysis <- peptide_info %>%
    inner_join(down_peptides, by="Sequence") %>%
    left_join(survival_data, by=c("Sequence"="Peptide"))
  
  # Figure 4 생성
  create_figure4(up_analysis, down_analysis)
  
  # Significant survival results 저장
  save_survival_results(up_analysis, "Up")
  save_survival_results(down_analysis, "Down")
}

# 시각화 실행
run_visualization()
