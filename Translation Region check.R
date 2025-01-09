################################# Part 1: 데이터 분석 #################################
# 필요한 패키지 설치 및 로드
required_packages <- c(
  "BiocManager",
  "readxl", "writexl", "dplyr", "tidyr",
  "Biostrings", "GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38", 
  "rentrez", "GenomicRanges", "IRanges"
)

# 패키지 설치 함수
install_required_packages <- function() {
  for(pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if(pkg == "BiocManager") {
        install.packages(pkg)
      } else if(pkg %in% c("Biostrings", "GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38", "rentrez", "GenomicRanges")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
  }
}

# 패키지 로드 함수
load_packages <- function() {
  lapply(required_packages, library, character.only = TRUE)
}

# 초기 설정
install_required_packages()
load_packages()

# NCBI 요청 함수 개선
fetch_from_ncbi <- function(accession, max_attempts = 3, delay = 1) {
  for (attempt in 1:max_attempts) {
    tryCatch({
      Sys.sleep(delay)  # API 요청 간 딜레이
      cat(sprintf("\nFetching data for %s (Attempt %d/%d)\n", accession, attempt, max_attempts))
      
      search_result <- entrez_search(db="nucleotide", term=accession)
      if (search_result$count == 0) {
        stop(paste("No sequence found for accession:", accession))
      }
      seq_record <- entrez_fetch(db="nucleotide", id=search_result$ids[1], rettype="gb", retmode="text")
      
      cat("Successfully retrieved data\n")
      return(seq_record)
      
    }, error = function(e) {
      cat(sprintf("Attempt %d failed: %s\n", attempt, as.character(e)))
      if (attempt == max_attempts) stop("Maximum attempts reached")
      delay <- delay * 2  # 지수 백오프
      cat(sprintf("Waiting %d seconds before retry...\n", delay))
      Sys.sleep(delay)
    })
  }
}

# Kozak sequence 강도 평가 함수
evaluate_kozak_strength <- function(sequence) {
  # Kozak consensus: (gcc/gccrccATGg)
  if (nchar(sequence) != 10) return(0)
  
  score <- 0
  # 핵심 위치 점수
  if (substr(sequence, 4, 6) == "ATG") score <- score + 3
  if (substr(sequence, 1, 1) == "G") score <- score + 2
  if (substr(sequence, 7, 7) == "G") score <- score + 2
  
  # 부가 위치 점수
  if (substr(sequence, 2, 3) %in% c("CC", "GC")) score <- score + 1
  
  return(score)
}

# Translation Initiation Site 분석 함수
find_potential_TIS <- function(mRNA_seq, peptide_pos, search_range = 300) {
  search_region_start <- max(1, peptide_pos - search_range)
  search_region <- subseq(mRNA_seq, start=search_region_start, end=peptide_pos)
  
  # ATG 위치 찾기 - 여기를 수정
  potential_starts <- matchPattern("ATG", search_region)
  
  # 각 ATG에 대한 분석
  tis_info <- lapply(start(potential_starts), function(pos) {
    if (pos >= 4 && pos <= length(search_region)-6) {
      context <- toString(subseq(search_region, pos-3, pos+6))
      kozak_score <- evaluate_kozak_strength(context)
      
      list(
        position = pos + search_region_start - 1,
        context = context,
        kozak_score = kozak_score,
        distance_to_peptide = peptide_pos - (pos + search_region_start - 1)
      )
    }
  })
  
  return(tis_info[!sapply(tis_info, is.null)])
}

# 서열 특성 분석 함수
analyze_sequence_features <- function(mRNA_seq, start_pos, end_pos, window_size = 50) {
  region_start <- max(1, start_pos - window_size)
  region_end <- min(length(mRNA_seq), end_pos + window_size)
  
  region_seq <- subseq(mRNA_seq, region_start, region_end)
  
  # GC content
  gc_content <- letterFrequency(region_seq, "GC") / length(region_seq)
  
  # Sequence complexity (단순 추정)
  dinucleotides <- dinucleotideFrequency(region_seq)
  complexity <- sd(as.vector(dinucleotides)) / mean(as.vector(dinucleotides))
  
  return(list(
    sequence = toString(region_seq),
    gc_content = gc_content,
    complexity = complexity,
    region = paste(region_start, "-", region_end)
  ))
}

# N-terminal extension length 계산 추가
classify_noncanonical_peptide <- function(mRNA_data, peptide_pos, features, frame_info, peptide_seq) {
  classification <- list()
  
  # Frame 정보 파싱
  frame_parts <- strsplit(frame_info, "-")[[1]]
  frame_num <- as.numeric(substr(frame_parts[1], 6, 6))
  stop_num <- as.numeric(frame_parts[2])
  
  # CDS 특성 확인
  cds_features <- Filter(function(x) x$type == "CDS", features)
  
  if (length(cds_features) > 0) {
    cds <- cds_features[[1]]
    canonical_frame <- (cds$start) %% 3 + 1
    actual_frame <- (peptide_pos) %% 3 + 1
    
    if (peptide_pos < cds$start) {
      if (actual_frame == canonical_frame) {
        classification$type <- "N-terminal extension"
        classification$subtype <- "Upstream start"
        classification$extension_length <- cds$start - peptide_pos  # 추가
      } else {
        classification$type <- "5' UTR"
        classification$subtype <- "upstream ORF (uORF)"
      }
      classification$distance_to_cds <- cds$start - peptide_pos
    } 
    else if (peptide_pos > cds$end) {
      classification$type <- "3' UTR"
      classification$subtype <- "downstream ORF (dORF)"
      classification$distance_from_cds <- peptide_pos - cds$end
    }
    else {
      if (actual_frame == canonical_frame) {
        classification$type <- "Canonical"
        classification$subtype <- "Main ORF peptide"
      } else {
        classification$type <- "CDS-internal"
        classification$subtype <- "Frame shift"
      }
    }
    
    classification$frame_offset <- (peptide_pos - cds$start) %% 3
    classification$canonical_frame <- canonical_frame
    classification$actual_frame <- actual_frame
  } else {
    if (any(grepl("NR_|XR_", mRNA_data$accession))) {
      classification$type <- "lncRNA"
      classification$subtype <- "Novel coding region in lncRNA"
    } else {
      classification$type <- "Other"
      classification$subtype <- "Unknown context"
    }
    classification$frame_offset <- (peptide_pos) %% 3
  }
  
  classification$stop_codon_number <- stop_num
  return(classification)
}

# mRNA 구조 정보 가져오기 함수 수정
get_mRNA_structure <- function(accession) {
  cat(sprintf("Fetching data for accession: %s\n", accession))
  
  # 버전 번호 제거한 기본 accession
  base_accession <- sub("\\.[0-9]+$", "", accession)
  
  # 먼저 주어진 정확한 버전으로 시도
  tryCatch({
    search_result <- entrez_search(db="nucleotide", term=accession)
    if (search_result$count > 0) {
      seq_record <- entrez_fetch(db="nucleotide", id=search_result$ids[1], rettype="gb", retmode="text")
      cat(sprintf("Found exact version: %s\n", accession))
    } else {
      # 정확한 버전이 없으면 모든 버전 검색
      search_term <- paste0(base_accession, "[Accession]")
      search_result <- entrez_search(db="nucleotide", term=search_term)
      
      if (search_result$count == 0) {
        stop(paste("No sequence found for accession:", accession))
      }
      
      # 가장 최신 버전 사용
      seq_record <- entrez_fetch(db="nucleotide", id=search_result$ids[1], rettype="gb", retmode="text")
      cat(sprintf("Using latest version for %s\n", base_accession))
    }
    
    # 서열 추출 및 정제
    seq_lines <- strsplit(seq_record, "\n")[[1]]
    seq_string <- paste(seq_lines[grep("^\\s+\\d+\\s+[a-z]+", seq_lines, ignore.case=TRUE)], collapse="")
    seq_string <- gsub("\\s+\\d+\\s+", "", seq_string)
    seq_string <- gsub("[^ACGT]", "", toupper(seq_string))
    
    # 구조 정보 추출
    features <- list()
    current_feature <- NULL
    for (line in seq_lines) {
      if (grepl("^     [a-zA-Z]", line)) {
        if (!is.null(current_feature)) {
          features[[length(features) + 1]] <- current_feature
        }
        current_feature <- list(type = sub("^     ([a-zA-Z]+).*", "\\1", line))
        if (grepl("\\d+\\.+\\d+", line)) {
          positions <- strsplit(sub(".*?(\\d+\\.+\\d+).*", "\\1", line), "\\.+")[[1]]
          current_feature$start <- as.numeric(positions[1])
          current_feature$end <- as.numeric(positions[2])
        }
      } else if (!is.null(current_feature) && grepl("^                     /", line)) {
        key_value <- strsplit(sub("^                     /", "", line), "=")[[1]]
        if (length(key_value) > 1) {
          current_feature[[key_value[1]]] <- gsub('"', '', key_value[2])
        }
      }
    }
    if (!is.null(current_feature)) {
      features[[length(features) + 1]] <- current_feature
    }
    
    return(list(seq=DNAString(seq_string), features=features))
    
  }, error = function(e) {
    cat(sprintf("Error fetching %s: %s\n", accession, as.character(e)))
    stop(paste("Failed to fetch sequence:", as.character(e)))
  })
}

# find_peptide_in_frame 함수 수정
find_peptide_in_frame <- function(mRNA_seq, peptide_seq, frame_info) {
  frame_parts <- strsplit(frame_info, "-")[[1]]
  frame_num <- as.numeric(substr(frame_parts[1], 6, 6))
  
  # 모든 가능한 프레임에서 검색
  all_frames <- 1:3
  for (current_frame in all_frames) {
    seq_length <- length(mRNA_seq)
    start_pos <- current_frame
    translate_length <- floor((seq_length - start_pos + 1) / 3) * 3
    
    translated_seq <- as.character(translate(subseq(mRNA_seq, start=start_pos, width=translate_length)))
    peptide_match <- gregexpr(peptide_seq, translated_seq, fixed = TRUE)[[1]]
    
    if (peptide_match[1] != -1) {
      actual_pos <- peptide_match[1] * 3 + start_pos - 1
      cat(sprintf("Found peptide in frame %d (requested frame %d)\n", 
                  current_frame, frame_num))
      return(actual_pos)
    }
  }
  
  return(NA)
}

# 구조 정보를 GRanges 객체로 변환하는 함수
convert_to_granges <- function(features) {
  gr_list <- lapply(features, function(f) {
    if ("start" %in% names(f) && "end" %in% names(f)) {
      GRanges(seqnames = "1",
              ranges = IRanges(start = f$start, end = f$end),
              type = f$type)
    }
  })
  gr_list <- gr_list[!sapply(gr_list, is.null)]
  do.call(c, gr_list)
}

# find_peptide_location 함수 끝부분의 result 반환 부분 수정
find_peptide_location <- function(peptide_seq, accession, frame_info) {
  cat(sprintf("\nProcessing peptide: %s, Accession: %s, Frame: %s\n", 
              peptide_seq, accession, frame_info))
  
  tryCatch({
    mRNA_data <- get_mRNA_structure(accession)
    mRNA_data$accession <- accession
    peptide_pos <- find_peptide_in_frame(mRNA_data$seq, peptide_seq, frame_info)
    
    if (is.na(peptide_pos)) {
      return(list(
        peptide = peptide_seq,
        accession = accession,
        frame = frame_info,
        position = NA,
        locations = "Not found",
        peptide_type = "Not found",
        peptide_subtype = "Not found",
        frame_offset = NA,
        stop_codon_number = NA,
        distance_to_feature = NA,
        gc_content = NA,
        sequence_complexity = NA,
        surrounding_sequence = NA,
        nearest_TIS = NA,
        kozak_context = NA,
        kozak_score = NA
      ))
    }
    
    # 상세 분석 수행
    tis_info <- find_potential_TIS(mRNA_data$seq, peptide_pos)
    seq_features <- analyze_sequence_features(mRNA_data$seq, 
                                              peptide_pos, 
                                              peptide_pos + nchar(peptide_seq) * 3 - 1)
    classification <- classify_noncanonical_peptide(mRNA_data, 
                                                    peptide_pos, 
                                                    mRNA_data$features, 
                                                    frame_info,
                                                    peptide_seq)
    
    # 위치 정보 생성
    peptide_range <- GRanges(seqnames = "1",
                             ranges = IRanges(start = peptide_pos, 
                                              end = peptide_pos + nchar(peptide_seq) * 3 - 1))
    
    gr_features <- convert_to_granges(mRNA_data$features)
    overlaps <- findOverlaps(peptide_range, gr_features)
    
    if (length(overlaps) == 0) {
      locations <- "Intergenic"
    } else {
      overlapping_regions <- data.frame(
        start = start(ranges(peptide_range))[queryHits(overlaps)],
        end = end(ranges(peptide_range))[queryHits(overlaps)],
        feature_type = gr_features$type[subjectHits(overlaps)]
      )
      
      total_length <- nchar(peptide_seq) * 3
      location_summary <- overlapping_regions %>%
        mutate(length = end - start + 1) %>%
        group_by(feature_type) %>%
        summarise(
          total_bases = sum(length),
          percentage = round(sum(length) / total_length * 100, 1)
        )
      
      locations <- paste(
        sprintf("%s(%g%%)", 
                location_summary$feature_type, 
                location_summary$percentage),
        collapse = ", "
      )
    }
    
    result <- list(
      peptide = peptide_seq,
      accession = accession,
      frame = frame_info,
      position = peptide_pos,
      locations = locations,
      peptide_type = classification$type,
      peptide_subtype = classification$subtype,
      frame_offset = classification$frame_offset,
      stop_codon_number = classification$stop_codon_number,
      distance_to_feature = if(!is.null(classification$distance_to_cds)) 
        classification$distance_to_cds else 
          if(!is.null(classification$distance_from_cds)) 
            classification$distance_from_cds else NA,
      gc_content = seq_features$gc_content,
      sequence_complexity = seq_features$complexity,
      surrounding_sequence = seq_features$sequence,
      nearest_TIS = if(length(tis_info) > 0) 
        min(sapply(tis_info, function(x) x$distance_to_peptide)) else NA,
      kozak_context = if(length(tis_info) > 0) 
        tis_info[[which.min(sapply(tis_info, function(x) x$distance_to_peptide))]]$context else NA,
      kozak_score = if(length(tis_info) > 0) 
        tis_info[[which.min(sapply(tis_info, function(x) x$distance_to_peptide))]]$kozak_score else NA,
      canonical_frame = classification$canonical_frame,
      actual_frame = classification$actual_frame
    )
    
    return(result)
  }, error = function(e) {
    cat(sprintf("Error: %s\n", as.character(e)))
    return(list(
      peptide = peptide_seq,
      accession = accession,
      frame = frame_info,
      position = NA,
      locations = "Error",
      peptide_type = "Error",
      peptide_subtype = "Error",
      frame_offset = NA,
      stop_codon_number = NA,
      distance_to_feature = NA,
      gc_content = NA,
      sequence_complexity = NA,
      surrounding_sequence = NA,
      nearest_TIS = NA,
      kozak_context = NA,
      kozak_score = NA
    ))
  })
}

#결과 처리를 위한 함수
summarize_peptide_analysis <- function(results) {
  total_peptides <- nrow(results)
  not_found <- sum(results$peptide_type == "Not found", na.rm = TRUE)
  canonical <- sum(results$peptide_type == "Canonical", na.rm = TRUE)
  
  cat("분석 결과 요약:\n")
  cat("총 펩타이드 수:", total_peptides, "\n")
  cat("찾지 못한 펩타이드 수:", not_found, "\n")
  cat("Canonical 펩타이드 수:", canonical, "\n")
  cat("비율:\n")
  cat("- Not found:", round(not_found/total_peptides * 100, 2), "%\n")
  cat("- Canonical:", round(canonical/total_peptides * 100, 2), "%\n")
  
  # Canonical이나 Not found가 아닌 펩타이드만 필터링
  noncanonical_results <- results %>%
    filter(peptide_type != "Canonical" & peptide_type != "Not found")
  
  return(noncanonical_results)
}

# 데이터 처리 개선
process_excel_file <- function(input_file, output_file, save_interim = TRUE) {
  cat("\nStarting data processing...\n")
  
  # 진행 상황 표시 함수
  progress <- function(current, total) {
    cat(sprintf("\rProcessing: %d/%d (%.1f%%)", current, total, current/total*100))
    if(current == total) cat("\n")
  }
  
  # 데이터 읽기
  cat("Reading input file...\n")
  data <- read_excel(input_file)
  total_rows <- nrow(data)
  
  # 중간 결과 저장 함수
  save_interim_results <- function(results, step) {
    if(save_interim) {
      interim_file <- sprintf("interim_results_%s.xlsx", step)
      write_xlsx(results, interim_file)
      cat(sprintf("Saved interim results to %s\n", interim_file))
    }
  }
  
  # 결과 저장을 위한 새로운 열 추가
  new_cols <- c(
    "position", "location", "peptide_type", "peptide_subtype",
    "frame_offset", "stop_codon_number", "distance_to_feature",
    "gc_content", "sequence_complexity", "nearest_TIS",
    "kozak_context", "kozak_score", "surrounding_sequence",
    "error_message", "canonical_frame", "actual_frame"
  )
  
  for(col in new_cols) {
    data[[col]] <- NA
  }
  
  # 데이터 처리
  for (i in 1:nrow(data)) {
    tryCatch({
      result <- find_peptide_location(
        data$peptide_seq[i],
        data$accession_num[i],
        data$frame_info[i]
      )
      
      # 결과 저장
      data[i, names(result)] <- result
      
    }, error = function(e) {
      data$error_message[i] <- as.character(e)
      cat(sprintf("\nError in row %d: %s\n", i, as.character(e)))
    })
    
    if(i %% 10 == 0) {
      progress(i, total_rows)
      if(save_interim) save_interim_results(data[1:i,], sprintf("batch_%d", i))
    }
  }
  
  # 결과 요약 및 필터링
  cat("\nSummarizing results...\n")
  filtered_results <- summarize_peptide_analysis(data)
  
  # 최종 결과 저장
  write_xlsx(filtered_results, output_file)
  cat(sprintf("\nAnalysis complete. Results saved to %s\n", output_file))
  
  return(filtered_results)
}

## Part 1: 데이터 분석
# input_file <- "identification_non_cano_dual.xlsx"
# output_file <- "Enhanced_Peptide_analysis_TMT(2448).xlsx"
# results_tmt <- process_excel_file(input_file, output_file)


# 데이터 확인
data <- read_excel("Enhanced_Peptide_analysis_TMT(2448).xlsx")
head(data)  # 데이터의 처음 몇 행을 확인
str(data)   # 데이터 구조 확인

# 필터링된 데이터 확인
filtered_data <- data %>%
  filter(!is.na(peptide_type),
         !peptide_type %in% c("Canonical", "Not found", "Error"))
unique(filtered_data$peptide_type)  # 어떤 peptide_type이 있는지 확인

# Distance data 확인
distance_data <- filtered_data %>%
  filter(peptide_type %in% c("5' UTR", "3' UTR", "N-terminal extension"),
         !is.na(distance_to_feature))
nrow(distance_data)  # 행 수 확인

# Kozak data 확인
kozak_data <- filtered_data %>%
  filter(!is.na(kozak_score))
nrow(kozak_data)  # 행 수 확인


################################# Part 2: 시각화 #################################
library(ggplot2)
library(gridExtra)
library(cowplot)
library(magick)
library(tidyr)      # filter 함수와 %>% 연산자를 위해 필요
library(readxl)     # read_excel 함수를 위해 필요
library(dplyr)      # 데이터 처리를 위해 필요
library(rsvg)      # SVG 파일 처리를 위해 필요


setwd("C:/Users/starh/OneDrive/Desktop/Noncano/코드정리/Figure1")

create_figure1 <- function(excel_file = "Enhanced_Peptide_analysis_TMT(2448).xlsx") {
  # 데이터 로드 및 전처리
  data <- read_excel(excel_file) %>%
    filter(!is.na(peptide_type),
           !peptide_type %in% c("Canonical", "Not found", "Error")) %>%
    mutate(
      distance_to_feature = as.numeric(distance_to_feature),
      kozak_score = as.numeric(kozak_score)
    )
  
  # Color schemes (color blind friendly)
  location_colors <- c(
    "3' UTR" = "#663399",         
    "N-terminal extension" = "#FFD700",
    "5' UTR" = "#90EE90",        
    "CDS-internal" = "#000080",  
    "lncRNA" = "#20B2AA"        
  )
  
  subtype_colors <- c(
    "downstream ORF (dORF)" = "#663399",  
    "Frame shift" = "#000080",           
    "Novel coding region in lncRNA" = "#20B2AA", 
    "upstream ORF (uORF)" = "#90EE90",    
    "Upstream start" = "#FFD700"          
  )
  
  # Create base theme
  base_theme <- theme_minimal() +
    theme(
      text = element_text(family = "Helvetica", size = 12),  # 기본 텍스트 크기 증가
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # x축 텍스트 크기 증가
      axis.text.y = element_text(size = 14),  # y축 텍스트 크기 증가
      axis.title = element_text(size = 16, face = "bold"),  # 축 제목 크기 증가 및 볼드체
      plot.title = element_text(size = 18, face = "bold"),  # 플롯 제목 크기 증가
      legend.text = element_text(size = 12),  # 범례 텍스트 크기 증가
      legend.title = element_text(size = 14, face = "bold"),  # 범례 제목 크기 증가
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
  
  # A: Distribution plot
  p1 <- ggplot(data, aes(x = peptide_type, fill = peptide_subtype)) +
    geom_bar(position = "stack") +
    scale_fill_manual(values = subtype_colors) +
    base_theme +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    ) +
    labs(
      title = "Distribution of Noncanonical Peptide Types",
      x = "Peptide Type",
      y = "Count",
      fill = "Peptide Subtype"
    )
  
  # C: Distance plot
  if("distance_to_feature" %in% colnames(data)) {
    distance_data <- data %>%
      filter(peptide_type %in% c("5' UTR", "3' UTR", "N-terminal extension"),
             !is.na(distance_to_feature))
    
    if(nrow(distance_data) > 0) {
      p3 <- ggplot(distance_data, aes(x = peptide_type, y = distance_to_feature)) +
        geom_boxplot(aes(fill = peptide_type)) +
        scale_fill_manual(values = location_colors) +
        base_theme +
        theme(legend.position = "none") +
        labs(
          title = "Distance from Known Features",
          x = "Peptide Type",
          y = "Distance (bp)"
        )
    }
  }
  
  # D: Kozak score distribution
  if("kozak_score" %in% colnames(data)) {
    kozak_data <- data %>%
      filter(!is.na(kozak_score))
    
    if(nrow(kozak_data) > 0) {
      p4 <- ggplot(kozak_data, aes(x = peptide_type, y = kozak_score)) +
        geom_violin(aes(fill = peptide_type)) +
        scale_fill_manual(values = location_colors) +
        base_theme +
        theme(legend.position = "none") +
        labs(
          title = "Kozak Score Distribution",
          x = "Peptide Type",
          y = "Kozak Score"
        )
    }
  }
  
  # Save individual panels with consistent settings
  save_plot <- function(plot, filename) {
    ggsave(filename,
           plot,
           width = 7,      # Increased size
           height = 8,
           dpi = 300,
           device = "tiff",
           compression = "lzw")
    cat(sprintf("Saved %s\n", filename))
  }
  
  # Save individual plots
  save_plot(p1, "figure1A.tiff")
  if(exists("p3")) save_plot(p3, "figure1C.tiff")
  if(exists("p4")) save_plot(p4, "figure1D.tiff")
}

# 실행
create_figure1("Enhanced_Peptide_analysis_TMT(2448).xlsx")


