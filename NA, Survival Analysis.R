# Part 1: Data Analysis and Table Generation
# Part 2: High-quality Visualization

#################### Global Setup and Functions ####################

# Required packages installation and loading
required_packages <- c(
  "tidyverse", "openxlsx", "survival", "survminer", 
  "viridis", "scales", "gridExtra", "Cairo", "ComplexHeatmap","circlize", "readxl"
)

for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Global theme for visualization
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

# High-quality plot saving function
save_plot <- function(plot, filename, width = 8, height = 10) {
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  base_path <- "results/figures/"
  
  # TIFF version
  ggsave(paste0(base_path, filename, ".tiff"),
         plot,
         width = width,
         height = height,
         dpi = 300,
         device = "tiff",
         compression = "lzw")
  
  # PDF version
  ggsave(paste0(base_path, filename, ".pdf"),
         plot,
         width = width,
         height = height,
         device = cairo_pdf)
  
  # PNG version
  ggsave(paste0(base_path, filename, ".png"),
         plot,
         width = width,
         height = height,
         dpi = 300)
}

#################### Part 1: Data Analysis Functions ####################

# Data loading function
load_data <- function() {
  expression_data <- read.xlsx("identification_non_cano_dual_2448.xlsx")
  peptide_info <- read.xlsx("Enhanced_Peptide_analysis_TMT(2448).xlsx")
  clinical_data <- read.xlsx("clinical_info_with_clusters.xlsx")
  
  return(list(
    expression_data = expression_data,
    peptide_info = peptide_info,
    clinical_data = clinical_data
  ))
}

#################### Part 1: NA Analysis ####################

run_na_analysis <- function() {
  # 데이터 로드
  data <- load_data()
  expression_data <- data$expression_data
  clinical_data <- data$clinical_data
  peptide_info <- data$peptide_info
  
  # Expression 열 식별
  expression_cols <- grep("RE.P", colnames(expression_data), value = TRUE)
  
  # NA 패턴 분석
  na_stats <- expression_data %>%
    select(Sequence, all_of(expression_cols)) %>%
    pivot_longer(cols = all_of(expression_cols),
                names_to = "Sample",
                values_to = "Expression") %>%
    group_by(Sequence) %>%
    summarise(
      na_count = sum(is.na(Expression)),
      total_samples = n(),
      na_ratio = na_count / total_samples
    ) %>%
    left_join(peptide_info[, c("Sequence", "peptide_type")], by = "Sequence")
  
  # 샘플 타입별 NA 패턴
  sample_patterns <- expression_data %>%
    select(Sequence, all_of(expression_cols)) %>%
    pivot_longer(cols = all_of(expression_cols),
                names_to = "Sample",
                values_to = "Expression") %>%
    mutate(Type = ifelse(grepl("\\.T$", Sample), "Tumor", "NAT")) %>%
    group_by(Sequence, Type) %>%
    summarise(
      na_ratio = mean(is.na(Expression)),
      .groups = "drop"
    ) %>%
    left_join(peptide_info[, c("Sequence", "peptide_type")], by = "Sequence")
  
  # 결과 저장
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  write.xlsx(list(
    NA_Statistics = na_stats,
    Sample_Patterns = sample_patterns
  ), "results/tables/na_pattern_analysis.xlsx")
  
  return(list(
    na_stats = na_stats,
    sample_patterns = sample_patterns
  ))
}

#################### Part 1: Survival Analysis ####################

run_survival_analysis <- function() {
  # 데이터 로드
  data <- load_data()
  expression_data <- data$expression_data
  clinical_data <- data$clinical_data
  peptide_info <- data$peptide_info
  
  # Tumor 샘플만 선택
  expression_cols <- grep("RE.P.*\\.T$", colnames(expression_data), value = TRUE)
  tumor_clinical <- clinical_data %>%
    filter(!is.na(OS_time), !is.na(OS_event))
  
  survival_results <- data.frame()
  
  # 각 펩타이드별 분석
  for(peptide in expression_data$Sequence) {
    peptide_expr <- expression_data[expression_data$Sequence == peptide, expression_cols]
    
    if(ncol(peptide_expr) >= 10) {
      expr_data <- data.frame(
        Sample = names(peptide_expr),
        Expression = as.numeric(peptide_expr[1,]),
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          has_expression = !is.na(Expression),
          expr_level = case_when(
            is.na(Expression) ~ "Not Detected",
            Expression > median(Expression, na.rm = TRUE) ~ "High",
            TRUE ~ "Low"
          )
        )
      
      analysis_data <- expr_data %>%
        left_join(tumor_clinical, by = "Sample")
      
      # Cox 분석
      if(sum(!is.na(analysis_data$Expression)) >= 10) {
        cox_model <- tryCatch({
          coxph(Surv(OS_time, OS_event == "X") ~ Expression, data = analysis_data)
        }, error = function(e) NULL)
        
        if(!is.null(cox_model)) {
          cox_summary <- summary(cox_model)
          
          result_row <- data.frame(
            Peptide = peptide,
            Peptide_Type = peptide_info$peptide_type[peptide_info$Sequence == peptide],
            HR = cox_summary$conf.int[1],
            HR_lower = cox_summary$conf.int[3],
            HR_upper = cox_summary$conf.int[4],
            P_value = cox_summary$coefficients[5],
            N_samples = nrow(analysis_data),
            N_expressed = sum(!is.na(analysis_data$Expression))
          )
          
          survival_results <- rbind(survival_results, result_row)
        }
      }
    }
  }
  
  # 결과 저장
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  write.xlsx(survival_results, "results/tables/survival_analysis.xlsx")
  
  return(survival_results)
}

#################### Part 1: NA Analysis ####################

run_na_analysis <- function() {
  # 데이터 로드
  data <- load_data()
  expression_data <- data$expression_data
  clinical_data <- data$clinical_data
  peptide_info <- data$peptide_info
  
  # Expression 열 식별
  expression_cols <- grep("RE.P", colnames(expression_data), value = TRUE)
  
  # NA 패턴 분석
  na_stats <- expression_data %>%
    select(Sequence, all_of(expression_cols)) %>%
    pivot_longer(cols = all_of(expression_cols),
                 names_to = "Sample",
                 values_to = "Expression") %>%
    group_by(Sequence) %>%
    summarise(
      na_count = sum(is.na(Expression)),
      total_samples = n(),
      na_ratio = na_count / total_samples
    ) %>%
    left_join(peptide_info[, c("Sequence", "peptide_type")], by = "Sequence")
  
  # 샘플 타입별 NA 패턴
  sample_patterns <- expression_data %>%
    select(Sequence, all_of(expression_cols)) %>%
    pivot_longer(cols = all_of(expression_cols),
                 names_to = "Sample",
                 values_to = "Expression") %>%
    mutate(Type = ifelse(grepl("\\.T$", Sample), "Tumor", "NAT")) %>%
    group_by(Sequence, Type) %>%
    summarise(
      na_ratio = mean(is.na(Expression)),
      .groups = "drop"
    ) %>%
    left_join(peptide_info[, c("Sequence", "peptide_type")], by = "Sequence")
  
  # 결과 저장
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  write.xlsx(list(
    NA_Statistics = na_stats,
    Sample_Patterns = sample_patterns
  ), "results/tables/na_pattern_analysis.xlsx")
  
  return(list(
    na_stats = na_stats,
    sample_patterns = sample_patterns
  ))
}

#################### Part 1: Survival Analysis ####################

run_survival_analysis <- function() {
  # 데이터 로드
  data <- load_data()
  expression_data <- data$expression_data
  clinical_data <- data$clinical_data
  peptide_info <- data$peptide_info
  
  # Tumor 샘플만 선택
  expression_cols <- grep("RE.P.*\\.T$", colnames(expression_data), value = TRUE)
  tumor_clinical <- clinical_data %>%
    filter(!is.na(OS_time), !is.na(OS_event))
  
  survival_results <- data.frame()
  
  for(peptide in expression_data$Sequence) {
    peptide_expr <- expression_data[expression_data$Sequence == peptide, expression_cols]
    
    if(ncol(peptide_expr) >= 10) {
      expr_data <- data.frame(
        Sample = names(peptide_expr),
        Expression = as.numeric(peptide_expr[1,]),
        stringsAsFactors = FALSE
      ) %>%
        mutate(
          has_expression = !is.na(Expression),
          expr_level = case_when(
            is.na(Expression) ~ "Not Detected",
            Expression > median(Expression, na.rm = TRUE) ~ "High",
            TRUE ~ "Low"
          )
        )
      
      analysis_data <- expr_data %>%
        left_join(tumor_clinical, by = "Sample")
      
      # Presence/Absence 분석
      if(sum(!is.na(analysis_data$Expression)) >= 10) {
        # Log-rank test for presence/absence
        lr_presence <- survdiff(Surv(OS_time, OS_event == "X") ~ has_expression, 
                                data = analysis_data)
        presence_p <- 1 - pchisq(lr_presence$chisq, df = 1)
        
        # Log-rank test for expression levels
        expr_data_detected <- analysis_data[!is.na(analysis_data$Expression),]
        if(nrow(expr_data_detected) >= 10) {
          lr_levels <- survdiff(Surv(OS_time, OS_event == "X") ~ expr_level,
                                data = expr_data_detected)
          levels_p <- 1 - pchisq(lr_levels$chisq, df = 1)
        } else {
          levels_p <- NA
        }
        
        # Cox analysis
        cox_model <- coxph(Surv(OS_time, OS_event == "X") ~ Expression, 
                           data = analysis_data)
        cox_summary <- summary(cox_model)
        
        result_row <- data.frame(
          Peptide = peptide,
          Peptide_Type = peptide_info$peptide_type[peptide_info$Sequence == peptide],
          Presence_P = presence_p,
          Levels_P = levels_p,
          Cox_P = cox_summary$coefficients[5],
          HR = cox_summary$conf.int[1],
          HR_lower = cox_summary$conf.int[3],
          HR_upper = cox_summary$conf.int[4],
          N_samples = nrow(analysis_data),
          N_expressed = sum(!is.na(analysis_data$Expression))
        )
        
        survival_results <- rbind(survival_results, result_row)
      }
    }
  }
  
  write.xlsx(survival_results, "results/tables/survival_analysis.xlsx")
  return(survival_results)
}
#################### Part 2: Visualization ####################
visualize_na_patterns <- function() {
  library(readxl)
  
  na_results_file <- "na_pattern_analysis.xlsx"
  
  # readxl로만 시트 읽기
  na_stats <- read_excel(na_results_file, sheet = "NA_Statistics")
  sample_patterns <- read_excel(na_results_file, sheet = "Sample_Patterns")
  
  print("Data loaded successfully")
  print(paste("NA_stats dimensions:", paste(dim(na_stats), collapse = "x")))
  print(paste("sample_patterns dimensions:", paste(dim(sample_patterns), collapse = "x")))
  
  # Type 열 처리
  sample_patterns$Type[sample_patterns$Type == ""] <- "NAT"
  sample_patterns$Type <- factor(sample_patterns$Type,
                                 levels = c("Tumor", "NAT"))
  
  # NA 비율 분포 플롯
  p1 <- ggplot(na_stats, aes(x = na_ratio * 100)) +
    geom_histogram(bins = 50, fill = "#4A90E2", alpha = 0.7) +
    labs(
      title = "Distribution of Missing Values",
      x = "Percentage of Missing Values (%)",
      y = "Number of Peptides"
    ) +
    create_theme()
  
  # 펩타이드 타입별 NA 비율 박스플롯
  p2 <- ggplot(sample_patterns, 
               aes(x = peptide_type, y = na_ratio * 100, fill = Type)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("Tumor" = "#E74C3C", "NAT" = "#2ECC71")) +
    labs(
      title = "Missing Values by Peptide Type and Sample Type",
      x = "Peptide Type",
      y = "Percentage of Missing Values (%)",
      fill = "Sample Type"
    ) +
    create_theme()
  
  # 결과 저장
  save_plot(p1, "na_distribution", width = 10, height = 8)
  save_plot(p2, "na_by_peptide_type", width = 12, height = 8)
  
  return(list(p1 = p1, p2 = p2))
}

visualize_survival_results <- function() {
  # 생존 분석 결과 불러오기
  survival_results <- read.xlsx("survival_analysis.xlsx")
  
  # P-value 분포 플롯
  p_value_long <- survival_results %>%
    select(Peptide, Presence_P, Levels_P, Cox_P) %>%
    pivot_longer(
      cols = c(Cox_P, Levels_P, Presence_P),
      names_to = "Analysis_Type",
      values_to = "P_value"
    )
  
  p1 <- ggplot(p_value_long, aes(x = P_value, fill = Analysis_Type)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.7) +
    facet_wrap(~Analysis_Type) +
    scale_fill_manual(values = c(
      "Cox_P" = "#98FB98",      
      "Levels_P" = "#FFB6C1",   
      "Presence_P" = "#B0C4DE"  
    )) +
    labs(
      title = "Distribution of P-values Across Different Analyses",
      x = "P-value",
      y = "Count"
    ) +
    create_theme()
  
  # Peptide type별 유의미한 결과 분포
  significant_summary <- survival_results %>%
    group_by(Peptide_Type) %>%
    summarise(
      Total = n(),
      Significant_Presence = sum(Presence_P < 0.05, na.rm = TRUE),
      Significant_Levels = sum(Levels_P < 0.05, na.rm = TRUE),
      Significant_Cox = sum(Cox_P < 0.05, na.rm = TRUE)
    ) %>%
    pivot_longer(
      cols = starts_with("Significant"),
      names_to = "Analysis_Type",
      values_to = "Significant_Count"
    ) %>%
    mutate(
      Percent_Significant = (Significant_Count/Total) * 100,
      Analysis_Type = gsub("Significant_", "", Analysis_Type)
    )
  
  p2 <- ggplot(significant_summary, 
               aes(x = Peptide_Type, y = Percent_Significant, fill = Analysis_Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    scale_fill_manual(values = c(
      "Presence" = "#B0C4DE",
      "Levels" = "#FFB6C1",
      "Cox" = "#98FB98"
    )) +
    geom_text(aes(label = sprintf("%.1f%%", Percent_Significant)),
              position = position_dodge(width = 0.9),
              vjust = -0.5) +
    labs(
      title = "Percentage of Significant Results by Peptide Type",
      x = "Peptide Type",
      y = "Percentage of Significant Results (%)",
      fill = "Analysis Type"
    ) +
    create_theme()
  
  # 결과 저장
  save_plot(p1, "survival_pvalue_distribution", width = 12, height = 6)
  save_plot(p2, "survival_significant_by_type", width = 12, height = 8)
  
  return(list(p1 = p1, p2 = p2))
}
#################### Main Execution Functions ####################

# Part 1 실행 함수
run_part1_analysis <- function() {
  cat("Running NA analysis...\n")
  na_results <- run_na_analysis()
  
  cat("Running survival analysis...\n")
  survival_results <- run_survival_analysis()
  
  cat("Part 1 analysis completed. Results saved in 'results/tables' directory.\n")
}

# Part 2 실행 함수
run_part2_visualization <- function() {
  cat("Generating NA analysis visualizations...\n")
  na_plots <- visualize_na_patterns()
  
  cat("Generating survival analysis visualizations...\n")
  survival_plots <- visualize_survival_results()
  
  cat("Part 2 visualization completed. Plots saved in 'results/figures' directory.\n")
}

setwd("C:/Users/samsung/OneDrive/바탕 화면/프로그래밍/KDH/Figure 6")

library(openxlsx)

######실행 예시###
# Part 1만 실행
run_part1_analysis()

# Part 2만 실행 (Part 1의 결과가 있어야 함)
run_part2_visualization()

# 또는 개별 분석만 실행
na_results <- run_na_analysis()
survival_results <- run_survival_analysis()

# 개별 시각화만 실행
na_plots <- visualize_na_patterns()
survival_plots <- visualize_survival_results()

