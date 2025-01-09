#################### Global Setup & Color Scheme ####################
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(grid)
library(Cairo)

colors <- list(
  Cluster = c(
    "1" = "#E41A1C",
    "2" = "#4DAF4A",
    "3" = "#984EA3",
    "4" = "#FF7F00"
  ),
  Smoking = c(
    "O" = "#D81B60",
    "X" = "yellow",
    "Unknown" = "gray"
  ),
  Death_Status = c(
    "N" = "#00C853",
    "Y" = "#FF5252"
  ),
  Celltype.based.subtype = c(
    "Cold_Immunogram" = "skyblue",
    "Hot_Immunogram" = "pink",
    "Unknown" = "gray"
  ),
  DX = c(
    "AD" = "#0288D1",
    "SC" = "#FFA700",
    "MA" = "green",
    "NC" = "#E53935",
    "Others" = "#5E35B1"
  )
)

create_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica", size = base_size),
      axis.text = element_text(size = base_size + 2),
      axis.title = element_text(size = base_size + 4, face = "bold"),
      plot.title = element_text(size = base_size + 6, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 2, face = "bold"),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
}

#################### Clinical Annotation Function ####################

create_clinical_plots <- function(clinical_data) {
  clinical_data <- clinical_data[order(clinical_data$Cluster),]
  
  ha <- HeatmapAnnotation(
    "Cluster" = clinical_data$Cluster,
    "Smoking" = clinical_data$Smoking,
    "Death Status" = clinical_data$Death_Status,
    "Celltype-based subtype" = clinical_data$Celltype.based.subtype,
    "DX" = clinical_data$DX,
    col = list(
      "Cluster" = colors$Cluster,
      "Smoking" = c("O" = "#D81B60", "X" = "yellow", "Unknown" = "gray"),
      "Death Status" = c("N" = "#00C853", "Y" = "#FF5252"),
      "Celltype-based subtype" = c("Cold_Immunogram" = "skyblue", "Hot_Immunogram" = "pink", "Unknown" = "gray"),
      "DX" = colors$DX
    ),
    annotation_name_side = "left",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 14, fontface = "bold"),
    gap = unit(1, "mm"),
    height = unit(2, "cm")
  )
  
  mat <- matrix(nrow=0, ncol=length(clinical_data$Cluster))
  
  for(format in c("pdf", "png", "tiff")) {
    if(format == "pdf") {
      pdf("annotation_heatmap_with_legend.pdf", width = 14, height = 6)
    } else if(format == "png") {
      png("annotation_heatmap_with_legend.png", width = 14*300, height = 6*300, res = 300)
    } else {
      tiff("annotation_heatmap_with_legend.tiff", width = 14*300, height = 6*300, res = 300, compression = "lzw")
    }
    
    ht <- Heatmap(mat, 
                  name = "Sample groups",
                  top_annotation = ha,
                  show_heatmap_legend = FALSE,
                  column_title = NULL,
                  cluster_columns = FALSE,
                  show_row_names = FALSE,
                  height = unit(0.5, "cm"))
    
    draw(ht, show_annotation_legend = TRUE,
         annotation_legend_side = "right",
         padding = unit(c(2, 2, 2, 20), "mm"))
    dev.off()
  }
  
  for(format in c("pdf", "png", "tiff")) {
    if(format == "pdf") {
      pdf("annotation_legends.pdf", width = 12, height = 2)
    } else if(format == "png") {
      png("annotation_legends.png", width = 12*300, height = 2*300, res = 300)
    } else {
      tiff("annotation_legends.tiff", width = 12*300, height = 2*300, res = 300, compression = "lzw")
    }
    
    layout(matrix(1:5, nrow = 1, ncol = 5))
    par(mar = c(1, 1, 1, 1))
    
    legend_info <- list(
      list(title="Cluster", labels=c("1", "2", "3", "4")),
      list(title="Smoking", labels=c("O", "X", "Unknown")),
      list(title="Death Status", labels=c("N", "Y")),
      list(title="Celltype-based subtype", 
           labels=c("Cold_Immunogram", "Hot_Immunogram", "Unknown"),
           colors=c("skyblue", "pink", "gray")),
      list(title="DX", labels=c("AD", "SC", "MA", "NC", "Others"))
    )
    
    for(info in legend_info) {
      plot.new()
      title(info$title, line = -1, cex.main = 1.2)
      legend("center", 
             legend = info$labels,
             fill = if(!is.null(info$colors)) {
               info$colors
             } else {
               unname(colors[[gsub(" ", "_", info$title)]][info$labels])
             },
             bty = "n",
             cex = 1)
    }
    dev.off()
  }
  
  return(list(heatmap = ha))
}

#################### Survival Analysis Function ####################

create_survival_plot <- function(clinical_data) {
  fit <- survfit(Surv(OS_time, Death_Status == "Y") ~ Cluster, 
                 data = clinical_data)
  
  surv_plot <- ggsurvplot(
    fit,
    data = clinical_data,
    pval = TRUE,
    conf.int = TRUE,
    palette = unname(colors$Cluster),
    ggtheme = create_theme(),
    risk.table = TRUE,
    risk.table.height = 0.3,
    xlab = "Time (months)",
    ylab = "Survival probability",
    legend.title = "Cluster",
    legend.labs = paste("Cluster", names(colors$Cluster)),
    font.legend = 12,
    font.x = 12,
    font.y = 12,
    font.title = 14
  )
  
  for(format in c("pdf", "png", "tiff")) {
    if(format == "pdf") {
      ggsave(paste0("survival_plot.", format),
             surv_plot$plot, 
             width = 10, height = 8,
             device = cairo_pdf)
    } else if(format == "png") {
      ggsave(paste0("survival_plot.", format),
             surv_plot$plot,
             width = 10, height = 8, 
             dpi = 300)
    } else {
      ggsave(paste0("survival_plot.", format),
             surv_plot$plot,
             width = 10, height = 8, 
             dpi = 300,
             device = "tiff",
             compression = "lzw")
    }
  }
  
  return(surv_plot)
}

#################### Main Execution ####################

run_analysis <- function(clinical_file) {
  clinical_data <- read.csv(clinical_file)
  
  print("Creating clinical annotation plots...")
  annotation_plots <- create_clinical_plots(clinical_data)
  
  print("Creating survival plot...")
  survival_plot <- create_survival_plot(clinical_data)
  
  print("Analysis complete. All plots have been saved.")
  
  return(list(
    annotation = annotation_plots,
    survival = survival_plot
  ))
}

# Example execution
setwd("C:/Users/samsung/OneDrive/바탕 화면/프로그래밍/KDH/Figure 5")
results <- run_analysis("clinical_info_with_clusters.csv")

# 또는 개별적으로 실행
clinical_data <- read.csv("clinical_info_with_clusters.csv")
annotation_plots <- create_clinical_plots(clinical_data)
survival_plot <- create_survival_plot(clinical_data)
