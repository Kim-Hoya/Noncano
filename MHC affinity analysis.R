#################### Global Setup ####################

# Required packages
required_packages <- c(
  "dplyr", "ggplot2", "tidyr", "viridis", "gridExtra", 
  "scales", "writexl", "Biostrings", "purrr", "stringr",
  "reticulate", "Cairo"
)

for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}


# Environment setup
setwd("C:/Users/starh/OneDrive/Desktop/Noncano/MHC")
Sys.setenv(RETICULATE_PYTHON = "C:/Users/starh/anaconda3/envs/mhcflurry/python.exe")
library(reticulate)
use_condaenv("mhcflurry", required = TRUE)
mhcflurry <- import("mhcflurry")
predictor <- mhcflurry$Class1AffinityPredictor$load()

# Create consistent theme for all plots
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
      axis.line = element_line(color = "black", linewidth = 0.5),
      strip.text = element_text(size = base_size + 2, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = NA)
    )
}

# Save high-quality plots
save_plot <- function(plot, filename, width = 8, height = 10) {
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  base_path <- "results/figures/"
  
  # TIFF version
  ggsave(paste0(base_path, filename, ".tiff"),
         plot, width = width, height = height, dpi = 300,
         device = "tiff", compression = "lzw")
  
  # PDF version
  ggsave(paste0(base_path, filename, ".pdf"),
         plot, width = width, height = height,
         device = cairo_pdf)
  
  # PNG version
  ggsave(paste0(base_path, filename, ".png"),
         plot, width = width, height = height, dpi = 300)
}

# Progress printing function
print_progress <- function(message, appendNewLine = TRUE) {
  cat(paste0("\n[", format(Sys.time(), "%H:%M:%S"), "] ", message))
  if (appendNewLine) cat("\n")
  flush.console()
}

#################### Part 1: Analysis Functions ####################

# Peptide splitting function
split_peptide <- function(peptide, window_size = 9, overlap = 8) {
  if (nchar(peptide) <= window_size) {
    return(peptide)
  }
  
  fragments <- character()
  step_size <- window_size - overlap
  for (i in seq(1, nchar(peptide) - window_size + 1, by = step_size)) {
    fragments <- c(fragments, substr(peptide, i, i + window_size - 1))
  }
  return(unique(fragments))
}

# FASTA processing function
process_fasta_with_fragments <- function(file_path, min_length = 5, max_length = 15, 
                                         window_size = 9, overlap = 8) {
  print_progress("Reading FASTA file...")
  sequences <- readAAStringSet(file_path)
  
  result_df <- data.frame(
    original_name = character(),
    original_sequence = character(),
    fragment = character(),
    fragment_position = integer(),
    fragment_length = integer(),
    original_length = integer(),
    is_fragment = logical(),
    stringsAsFactors = FALSE
  )
  
  total_seqs <- length(sequences)
  print_progress(sprintf("Processing %d sequences...", total_seqs))
  pb <- txtProgressBar(min = 0, max = total_seqs, style = 3)
  
  for (i in 1:total_seqs) {
    seq <- as.character(sequences[i])
    name <- names(sequences[i])
    
    if (nchar(seq) > max_length) {
      fragments <- split_peptide(seq, window_size, overlap)
      
      for (fragment in fragments) {
        pos <- regexpr(fragment, seq)[1]
        result_df <- rbind(result_df, data.frame(
          original_name = name,
          original_sequence = seq,
          fragment = fragment,
          fragment_position = pos,
          fragment_length = nchar(fragment),
          original_length = nchar(seq),
          is_fragment = TRUE,
          stringsAsFactors = FALSE
        ))
      }
    } else if (nchar(seq) >= min_length && nchar(seq) <= max_length) {
      result_df <- rbind(result_df, data.frame(
        original_name = name,
        original_sequence = seq,
        fragment = seq,
        fragment_position = 1,
        fragment_length = nchar(seq),
        original_length = nchar(seq),
        is_fragment = FALSE,
        stringsAsFactors = FALSE
      ))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(result_df)
}

# Main analysis function
run_mhc_analysis <- function(fasta_path, alleles, output_dir = "results") {
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  
  # Process peptides
  peptides_df <- process_fasta_with_fragments(fasta_path)
  
  # Run predictions for each allele
  all_predictions <- data.frame()
  pb <- txtProgressBar(min = 0, max = length(alleles), style = 3)
  
  for (i in seq_along(alleles)) {
    current_allele <- alleles[i]
    print_progress(sprintf("Processing allele %s...", current_allele), FALSE)
    
    predictions <- predictor$predict_to_dataframe(
      peptides = reticulate::r_to_py(peptides_df$fragment),
      alleles = reticulate::r_to_py(rep(current_allele, nrow(peptides_df)))
    )
    
    temp_df <- cbind(
      peptides_df,
      prediction = predictions$prediction,
      allele = current_allele
    )
    
    all_predictions <- rbind(all_predictions, temp_df)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Add HLA groups and binding levels
  all_predictions <- all_predictions %>%
    mutate(
      hla_group = ifelse(grepl("HLA-A", allele), "HLA-A", "HLA-B"),
      binding_level = factor(case_when(
        prediction <= 50 ~ "Strong",
        prediction <= 500 ~ "Moderate",
        TRUE ~ "Weak"
      ), levels = c("Strong", "Moderate", "Weak"))
    )
  
  # Generate fragment summary
  fragment_summary <- all_predictions %>%
    group_by(original_name, original_sequence, allele) %>%
    summarise(
      total_fragments = n(),
      strong_fragments = sum(binding_level == "Strong"),
      moderate_fragments = sum(binding_level == "Moderate"),
      best_fragment = fragment[which.min(prediction)],
      best_position = fragment_position[which.min(prediction)],
      best_affinity = min(prediction),
      .groups = 'drop'
    )
  
  # Generate best binders by allele
  best_by_allele <- all_predictions %>%
    group_by(allele) %>%
    slice_min(prediction, n = 20) %>%
    arrange(allele, prediction)
  
  # Save results
  write_xlsx(list(
    "All_Predictions" = all_predictions,
    "Fragment_Summary" = fragment_summary,
    "Best_by_Allele" = best_by_allele,
    "Strong_Binders" = filter(all_predictions, binding_level == "Strong") %>% 
      arrange(prediction),
    "Moderate_Binders" = filter(all_predictions, binding_level == "Moderate") %>% 
      arrange(prediction)
  ), path = file.path(output_dir, "tables", "mhc_analysis_results.xlsx"))
  
  return(list(
    predictions = all_predictions,
    fragment_summary = fragment_summary,
    best_by_allele = best_by_allele
  ))
}

#################### Part 3: Visualization Functions ####################

# MHC 분석의 Binding Affinity Distribution plot 추가
create_mhc_plots <- function(results_file) {
  results <- read_xlsx(results_file, sheet = "All_Predictions")
  
  if(!"hla_group" %in% colnames(results)) {
    results$hla_group <- ifelse(grepl("HLA-A", results$allele), "HLA-A", "HLA-B")
  }
  if(!"binding_level" %in% colnames(results)) {
    # binding_level의 순서 수정
    results$binding_level <- factor(case_when(
      results$prediction <= 50 ~ "Strong",
      results$prediction <= 500 ~ "Moderate",
      TRUE ~ "Weak"
    ), levels = c("Strong", "Moderate", "Weak"))  # 레벨 순서 명시적 지정
  }
  
  # 1. Overall Distribution Plot by HLA Group
  p1 <- ggplot(results, aes(x = allele, y = prediction)) +
    geom_violin(aes(fill = hla_group), alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.5) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = c("HLA-A" = "#FFB6C1", "HLA-B" = "#87CEEB")) +
    facet_wrap(~hla_group, scales = "free_x") +
    create_theme() +
    labs(title = "Distribution of HLA Binding Affinities",
         x = "HLA Allele", y = "Predicted Affinity (nM)",
         fill = "HLA Group")
  
  # 2. Binding Level Distribution
  p2 <- ggplot(results, aes(x = prediction)) +
    geom_histogram(aes(fill = hla_group), bins = 100) +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(values = c("HLA-A" = "#FFB6C1", "HLA-B" = "#87CEEB")) +
    facet_wrap(~binding_level) +
    create_theme() +
    labs(title = "Distribution of Binding Affinities by Level",
         x = "Predicted Affinity (nM)", y = "Count",
         fill = "HLA Group")
  
  # 3. Density Plot
  p3 <- ggplot(results, aes(x = prediction, fill = hla_group)) +
    geom_density(alpha = 0.7) +
    scale_x_log10() +
    facet_wrap(~binding_level, nrow = 1) +  # nrow=1 추가해서 가로로 한 줄로 배치
    scale_fill_manual(values = c("HLA-A" = "#FFB6C1", "HLA-B" = "#87CEEB")) +
    create_theme() +
    labs(title = "Binding Affinity Distribution by HLA Group",
         x = "Predicted Affinity (nM)",
         y = "Density",
         fill = "HLA Group")
  # 4. Heatmap of Binding Levels
  heatmap_data <- results %>%
    group_by(allele, binding_level) %>%
    summarise(count = n(), .groups = 'drop')
  
  p4 <- ggplot(heatmap_data, aes(x = allele, y = binding_level, fill = count)) +
    geom_tile() +
    scale_fill_viridis_c() +
    create_theme() +
    labs(title = "Binding Level Distribution Across HLA Alleles",
         x = "HLA Allele", y = "Binding Level",
         fill = "Count")
  
  # 5. Fragment Length Distribution
  p5 <- ggplot(results, aes(x = fragment_length, fill = hla_group)) +
    geom_bar(position = "dodge") +
    scale_fill_manual(values = c("HLA-A" = "#FFB6C1", "HLA-B" = "#87CEEB")) +
    create_theme() +
    labs(title = "Peptide Length Distribution by HLA Group",
         x = "Peptide Length", y = "Count",
         fill = "HLA Group")
  
  # 6. Top Binders Distribution
  top_binders <- results %>%
    group_by(allele) %>%
    slice_min(order_by = prediction, n = 10)
  
  p6 <- ggplot(top_binders, aes(x = reorder(allele, prediction), y = prediction)) +
    geom_boxplot(aes(fill = hla_group), alpha = 0.7) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = c("HLA-A" = "#FFB6C1", "HLA-B" = "#87CEEB")) +
    create_theme() +
    labs(title = "Top 10 Binders Distribution by HLA Allele",
         x = "HLA Allele", y = "Predicted Affinity (nM)",
         fill = "HLA Group")
  
  # 7. Fragment Position vs Affinity
  p7 <- ggplot(filter(results, is_fragment), 
               aes(x = fragment_position, y = prediction, color = hla_group)) +
    geom_point(alpha = 0.5) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c("HLA-A" = "#FFB6C1", "HLA-B" = "#87CEEB")) +
    facet_wrap(~allele) +
    create_theme() +
    labs(title = "Fragment Position vs Binding Affinity",
         x = "Fragment Position", y = "Predicted Affinity (nM)",
         color = "HLA Group")
  
  # Save all plots
  plot_list <- list(
    distribution = p1,
    binding_level = p2,
    density = p3,
    heatmap = p4,
    length_dist = p5,
    top_binders = p6,
    position_affinity = p7
  )
  
  for(name in names(plot_list)) {
    save_plot(plot_list[[name]], paste0("mhc_", name), 
              width = ifelse(name == "position_affinity", 15, 10),
              height = ifelse(name == "position_affinity", 12, 8))
  }
  
  # Summary panels
  summary_panel1 <- grid.arrange(p1, p2, p3, p4, ncol = 2)
  summary_panel2 <- grid.arrange(p5, p6, p7, ncol = 2)
  
  save_plot(summary_panel1, "mhc_summary_panel1", width = 20, height = 16)
  save_plot(summary_panel2, "mhc_summary_panel2", width = 20, height = 16)
  
  return(plot_list)
}

#################### Main Execution Functions ####################

# Part 1 실행 함수
run_part1_analysis <- function(fasta_path, alleles) {
  print_progress("Starting Part 1: Analysis and Table Generation")
  results <- run_mhc_analysis(fasta_path, alleles)
  print_progress("Part 1 complete. Results saved in results/tables directory.")
  return(results)
}

# Part 2 실행 함수
run_part2_visualization <- function(results_file = "mhc_analysis_results.xlsx") {
  print_progress("Starting Part 2: Visualization")
  plots <- create_mhc_plots(results_file)
  print_progress("Part 2 complete. Plots saved in results/figures directory.")
  return(plots)
}


#################### Execution Code ####################

# HLA 알렐 목록
alleles <- c(
  # HLA-A alleles
  "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06",
  "HLA-A*03:01", "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02",
  "HLA-A*26:01", "HLA-A*30:01", "HLA-A*30:02", "HLA-A*31:01",
  "HLA-A*32:01", "HLA-A*33:01", "HLA-A*68:01", "HLA-A*68:02",
  
  # HLA-B alleles
  "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01",
  "HLA-B*40:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01",
  "HLA-B*53:01", "HLA-B*57:01", "HLA-B*58:01"
)

# FASTA 파일 경로
fasta_path <- "C:/Users/starh/OneDrive/Desktop/Noncano/MHC/MHC_check.fasta"

# 실행 예시:
# 전체 분석 실행
results <- run_part1_analysis(fasta_path, alleles)
plots <- run_part2_visualization()

# 또는 개별적으로 실행
# Part 1만 실행
# results <- run_part1_analysis(fasta_path, alleles)

# Part 2만 실행 (Part 1의 결과가 있어야 함)
# plots <- run_part2_visualization()