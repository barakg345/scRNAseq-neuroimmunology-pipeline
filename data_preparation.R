**Overview:** This script documents the granular, exploratory development and step-by-step testing of the pairwise differential expression (PDE) pipeline. 
#Subsequent analyses for other cell types utilize the fully automated, streamlined pipeline derived from this foundational workflow.

# 1. Environment Setup & Custom Signatures
## 1.1 Package Loading & Cell-Type Signature Definitions
#Libraries and Signatures

{
  
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)   
    library(data.table)
    library(cowplot)
    library(ggplot2)
    library(scran)
    library(celldex)
    library(SingleR)
    library(tidyr)
    library(scales)
    library(corrplot)
    library(magrittr)
    library(dittoSeq)
    library(viridis)
    library(grid)
    library(magrittr)
    library(ggrepel)
    library(glmGamPoi)
    library(qs)
    library(patchwork)
    library(RColorBrewer)
    library(Matrix)
    library(parallel)
    library(future)
    library(GEOquery)
    library(harmony)
    library(ggprism)
    library(scCustomize)
    library(UCell)
    library(colorspace)
    library(SignatuR)
    library(SCPA)
    library(msigdbr)
    library(ComplexHeatmap)
    library(tidyverse)
    library(ggpubr)
    library(purrr)
    library(knitr)
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(reshape2)
    library(knitr)
    
  })
  
  
  #Prepare the sub-cell types signatures for the upcoming analyses
  
  signatures <- list(
    classical_monocytes = c("CD14", "CCR2", "FCGR1A", "CD36", "SELL"),
    intermediate_monocytes = c("CD14", "FCGR3A", "FCGR3B", "HLA-DRA", "HLA-DRB3", "HLA-DRB5", "HLA-DRB1", "CCR5", "CD86", "ITGAX", "CD74"),
    non_classical_monocytes = c("FCGR3A", "FCGR3B", "CX3CR1","CXCR4", "HLA-DRA", "HLA-DRB3", "HLA-DRB5", "HLA-DRB1", "SECISBP2L"),
    
    CD56_bright_nk = c("NCAM1", "CCR7", "CSF2", "CXCR3", "IFNG","IL2RB", "IL7R", "KIT", "KLRC1", "KLRD1", "NCR1", "SELL"),
    CD56_dim_nk = c("FCGR3A", "FCGR3B", "NCAM1", "CX3CR1", "CXCR1", "ITGB2", "KIR2DL1", "KIR3DL1", "KIR2DL3", "KIR2DL2",  "KLRC2", "KLRG1", "PRF1"),
    Adaptive_NK = c("B3GAT1", "KLRC2", "LILRB1", "Klra8"),
    
    
    naive_b_cells = c("CD19", "IGHD", "CD38", "CD24", "CR2", "IGHM", "FCER2"),
    non_switched_memory_b_cells = c("IGHD", "CD27", "CD38", "CD24", "CR2", "IGHM", "CD1c"),
    exhausted_b_cells = c("CD19", "CD38", "CD24", "ITGAX", "TBX21", "FcRL4", "FcRL5"),
    switched_memory_b_cells = c("CD27", "CD38", "CD24", "CR2", "IGHG1", "IGHG2","IGHG3", "IGHG4", "IGHA1", "IGHA2", "FAS"),
    Plasmablasts=c("CD19", "IGHD", "CD27", "CD38", "MS4A1", "HLA-DRA", "HLA-DRB3", "HLA-DRB5", "HLA-DRB1"),
    
    follicular_helper_t_cells = c("CXCR5", "PDCD1", "BCL6", "ICOS", "SH2D1A", "CD40LG", "IL21", "IL10", "B3GAT1"),
    regulatory_tcell = c("IL2RA", "FOXP3", "IL7R", "ICOS", "TNFRSF18", "CTLA4", "TIGIT", "ENTPD1", "TNFRSF9", "TNFRSF1B"),
    th1 = c("TBX21", "STAT4", "IFNG", "IL12RB1", "IL12RB2", "CXCR3", "CCR5", "CD44", "Fas", "	IFNGR1", "HAVCR2"),
    th17 = c("IL17A", "IL17F", "IL21", "IL22", "KLRB1", "CCR4", "CCR6", "APOD", "C1QL1", "CTSL"),
    th2 = c("GATA3", "STAT6", "IL4", "IL5", "IL13", "CCR4", "CCR8"),
    naive_cd4_t_cells=c("LCK", "ZAP70", "KLF2", "CD3G", "CD3D", "CD3GE", "CD28", "SELL", "PTPRC", "CD5", "CCR7", "CD27", "IL7R"),
    terminal_effector_cd4_t_cells= c("PTPRC", "CD44", "SELL", "IL2RB", "IL2RA", "IL2RG"),
    
    
    naive_cd8_t_cells = c("CD27", "CD28", "PTPRC", "SELL", "CCR7", "IL7R", "CD5", "CD8A", "CD8B", "CXCR3"),
    central_memory_cd8_t_cells = c("SELL", "CCR7", "IL7R", "CD27", "EOMES", "HNF1A", "BCL6", "ID3", "STAT3"),
    effector_memory_cd8_t_cells = c("SELL", "CCR7", "KLRG1", "	GZMB", "TBX21", "Blimp1", "ID2", "STAT4", "ITGAE", "ITGA1", "CXCR3"),
    terminal_effector_cd8_t_cells = c("KLRG1", "IL7R", "B3GAT1", "CXCR3", "CD27", "CD28", "IFNG", "IL2", "CDKN1B")
  )
  
}


# 2. Custom Single-Cell Analysis & Plotting Functions

# 2.1 Downstream Analysis Utilities

#Functions for analysis
{
  
  
  get_singler_main <- function(seu){
    
    bped <- celldex::BlueprintEncodeData()
    
    DefaultAssay(seu) <- "SCT"
    pred_main <- as.SingleCellExperiment(seu) %>%
      SingleR(ref = bped, labels = bped$label.main)
    
    
    seu$SingleR.labels <- pred_main$labels
    cat(paste0("\033[31m","\nSummarize of SingleR annotation table\n"))
    print(table(seu$SingleR.labels))
    return(seu)
    
    
  }
  
  get_singler_fine <- function(seu){
    
    monaco <- celldex::MonacoImmuneData()
    
    DefaultAssay(seu) <- "SCT"
    pred_fine <- as.SingleCellExperiment(seu) %>%
      SingleR(ref = monaco, labels = monaco$label.fine)
    
    
    seu$SingleR_sub_labels <- pred_fine$labels
    cat(paste0("\033[31m","\nSummarize of SingleR annotation table\n"))
    # print(table(seu$SingleR_sub_labels))    
    return(seu)
    
    
  }
  
  get_four_umap_basic <- function(seu,size=0.8) {
    p_clusters <- dittoDimPlot(seu, "seurat_clusters", reduction.use = "umap", size = size) +
      ggplot2::ggtitle(label = "UMAP by Seurat Clusters")
    
    p_sample <- dittoDimPlot(seu, "sample", reduction.use = "umap", size = size) +
      ggplot2::ggtitle(label = "UMAP by sample")
    
    p_singler <- dittoDimPlot(seu, "SingleR.labels", reduction.use = "umap", size = size) +
      ggplot2::ggtitle(label = "UMAP by SingleR Labels")
    
    p_condition <- dittoDimPlot(seu, "condition", reduction.use = "umap", size = size, opacity = 0.8) +
      ggplot2::ggtitle(label = "UMAP by Sample")
    
    # Combine plots into one
    p_final <- ((p_clusters + p_condition) / (p_singler + p_sample))
    
    # Return the combined plot invisibly to avoid printing in the function
    invisible(p_final)
  }
  
  
  
  get_p_avg_expr_log2 <- function(seu) {
    Idents(seu) <- "condition"
    DefaultAssay(seu) <- "SCT"
    
    avg_exp <- AverageExpression(seu, assays = "SCT", group.by = "condition")
    avg_exp <- as.data.frame(avg_exp)
    
    # Making sure the columns are numeric and making it ready for log
    avg_exp$SCT.control <- as.numeric(as.character(avg_exp$SCT.control))
    avg_exp$SCT.patient <- as.numeric(as.character(avg_exp$SCT.patient))
    
    avg_exp$SCT.control <- avg_exp$SCT.control + 1
    avg_exp$SCT.patient <- avg_exp$SCT.patient + 1
    
    line_data <- data.frame(x = range(avg_exp$SCT.control, avg_exp$SCT.patient, na.rm = TRUE))
    line_data$y <- line_data$x  # displaying the line y = x
    avg_exp <- avg_exp %>% tibble::rownames_to_column("gene")
    
    gene_to_be_displayed <- avg_exp %>%
      filter(SCT.patient >= 1.5, SCT.patient >= SCT.control) %>%
      pull(gene)
    
    mt_genes <- (length(grep("^MT", gene_to_be_displayed)) / length(gene_to_be_displayed)) * 100
    ribo_genes <- (length(grep("^RP", gene_to_be_displayed)) / length(gene_to_be_displayed)) * 100
    
    p_avg_mon <- ggplot(avg_exp, aes(x = SCT.control, y = SCT.patient)) +
      geom_point(aes(color = ifelse(gene %in% gene_to_be_displayed, "significant", "non-significant"))) +
      geom_line(data = line_data, aes(x = x, y = y), color = "blue") +
      scale_x_continuous(trans = scales::log2_trans()) +  # Log2 scale for x-axis
      scale_y_continuous(trans = scales::log2_trans()) +  # Log2 scale for y-axis
      scale_color_manual(values = c("significant" = "red", "non-significant" = "black")) +
      labs(x = "Log2(Control)", y = "Log2(Patient)",
           title = "Control vs Patient Averaged Expression Levels") +
      theme_minimal() +
      guides(color = guide_legend(title = NULL)) +  # Remove legend title
      annotate("text", x = 0, y = Inf, label = sprintf("Percent of ribo_genes: %.2f%%", ribo_genes), 
               hjust = -0.1, vjust = 2.1, size = 5, color = "blue") +
      annotate("text", x = 0, y = Inf, label = sprintf("Percent of mt_genes: %.2f%%", mt_genes), 
               hjust = -0.1, vjust = 4.1, size = 5, color = "blue")
    
    # Return the plot invisibly to avoid printing
    return(p_avg_mon)
  }
  
  
  
  convert_to_msigdbr_format <- function(pathway_name,db) {
    paste0(db,"_", gsub(" ", "_", toupper(pathway_name)))
  }
  
  
  get_volcano <- function(de_results, key_genes, plottitle, group_a, group_b) {
    
    # Load necessary libraries and source the custom theme
    {
      suppressPackageStartupMessages(library(tidyverse))
      suppressPackageStartupMessages(library(gridExtra))
      suppressPackageStartupMessages(library(ggrepel))
      source('ggplot_theme_Publication-2.R')
    } 
    
    # Calculate the maximum absolute value on the x-axis to determine symmetric limits
    max_x_value <- max(abs(de_results$avg_log2FC))
    
    # Create the color mapping dynamically
    color_values <- c(
      group_a = "#00B0BA",
      group_b = "#f87f01",
      "Other" = "#cccccc"
    )
    names(color_values)[1] <- group_a
    names(color_values)[2] <- group_b
    
    # Create the volcano plot
    p_mon <- ggplot(de_results, aes(x = avg_log2FC, y = y_display, color = color_group)) +
      geom_point() +
      geom_text_repel(data = de_results %>% filter(gene %in% key_genes),
                      aes(label = gene, color = color_group),  
                      point.padding = 5,  #  point padding
                      nudge_y = 0.5,  # Slightly nudge the text labels up
                      max.overlaps = Inf,  # Allow infinite overlaps for key genes
                      size = 4  # Set a fixed size for all gene labels
      ) +
      scale_x_continuous(limits = c(-max_x_value - 1, max_x_value + 1)) +  # Expanding x-axis limits
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +  # Expanding y-axis on both ends
      geom_hline(yintercept = 1.30, linetype = "dashed", color = "blue", size = 1) +
      geom_vline(xintercept = -1, linetype = "dashed", color = "blue", size = 1) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "blue", size = 1) +
      labs(title = plottitle,
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value",
           size = "Percentage in Group 1",
           color = "Gene Group") +
      guides(size = guide_legend(title = "Percentage in Patient cells")) +
      scale_color_manual(values = color_values,
                         breaks = c(group_a, group_b, "Other")) +
      theme_dark_blue() +
      annotate("text", x = -Inf, y = -Inf, label = "The points were capped at 280 for display",
               hjust = 0, vjust = -1, color = "white", size = 3.5)
    
    return(p_mon)
  }
  
  
  
  get_four_umaps_pathway <- function (pathway_name, de_results ,seu, threshold,db){
    
    full_db <- paste0("CP:",db)
    specific_pathway  <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = full_db) %>% filter(gs_name == convert_to_msigdbr_format(pathway_name,db))
    
    key_genes <- list(
      key_genes = de_results %>%
        filter(gene %in% specific_pathway$human_gene_symbol) %>%
        filter(abs(avg_log2FC) > threshold) %>%
        pull(gene) %>%
        unique())
    
    seu <- AddModuleScore_UCell(seu, features = key_genes, name = NULL)      
    pathway_name <- paste('"',pathway_name,'"',"Pathway Key Genes Expression")
    
    #And finally, print it
    p_key_genes <- FeaturePlot_scCustom(seurat_object = seu, features = "key_genes")  + ggtitle(pathway_name) +
      theme(plot.title = element_text(face = "plain", size = 12))
    
    #Creating the other regular 3 umaps and printing them all togther
    {
      p_clusters <- dittoDimPlot(seu, "seurat_clusters", reduction.use = "umap") +
        ggplot2::ggtitle(label = "UMAP by Seurat Clusters")
      
      p_condition <- dittoDimPlot(seu, "condition", reduction.use = "umap", opacity = 0.8) +
        ggplot2::ggtitle(label = "UMAP by Sample")
      
      p_subcluster <- dittoDimPlot(seu,"adaptive_type", reduction.use = "umap", color.panel = brewer.pal(6, "Paired")) +
        ggplot2::ggtitle(label = "UMAP by subtype and condition")
      
      
      }
    
    p_umap_key_genes <- ((p_clusters | p_condition) / (p_key_genes | p_subcluster)) 
    #Boxplot
    {
      
      boxplot <- dittoBoxPlot(object = seu, var = "key_genes", group.by = "condition", main = pathway_name)
      
      
    }
    
    return(  list(p_umap_key_genes = p_umap_key_genes, key_genes = key_genes, boxplot = boxplot))
    
    
    
  }
  
  get_four_umaps_subclustring <- function(seu) {
    
    p_clusters <- dittoDimPlot(seu, "seurat_clusters", reduction.use = "umap") +
      ggplot2::ggtitle(label = "UMAP by Seurat Clusters")
    
    p_condition <- dittoDimPlot(seu, "condition", reduction.use = "umap", opacity = 0.8) +
      ggplot2::ggtitle(label = "UMAP by Sample")
    
    p_subcluster <- dittoDimPlot(seu,"SingleR_sub_labels", reduction.use = "umap", color.panel = brewer.pal(9, "Set1")) +
      ggplot2::ggtitle(label = "UMAP by SingleR sub cluster")
    
    p_phase <- dittoDimPlot(seu,"Phase", reduction.use = "umap", ) +
      ggplot2::ggtitle(label = "UMAP by Celly Cycle")
    
    p_4_umaps_sub_clustring <- ((p_clusters | p_condition) / (p_subcluster | p_phase)) 
    
    return(p_4_umaps_sub_clustring)
  }
  
  get_scpa <- function(seu, fgroup1, fgroup1_population) {
    #kegg
    {
      
      pathways_kegg <- msigdbr("Homo sapiens", "C2", "CP:KEGG") %>%
        format_pathways()
      
      results_kegg <- compare_seurat(seu,
                                     assay = "SCT",
                                     group1 = fgroup1, 
                                     group1_population = fgroup1_population,
                                     pathways = pathways_kegg,
                                     parallel = TRUE,
                                     cores = 4)   %>%
        mutate(Pathway = gsub("^KEGG_", "", Pathway),    # remove "KEGG_MEDICUS" 
               Pathway = tolower(Pathway),               # convert to lower case
               Pathway = gsub("_", " ", Pathway),        # replace "_" with ""
               Pathway = str_to_title(Pathway))          # Convert to title case
      
      
      
      results_kegg <- results_kegg %>%
        filter(FC > -0.3) %>%
        arrange(desc(qval)) %>%
        slice_head(n = 10)
      
      p_pathways_kegg <- ggdotchart(results_kegg, x = "Pathway", y = "FC",  title = "Pathway Analysis by KEGG",
                                    color = "adjPval",
                                    sorting = "descending",
                                    add = "segments",
                                    rotate = TRUE,
                                    dot.size = 10,
                                    gradient.cols = c("darkgreen", "grey"),
                                    font.label = list(color = "white", size = 9, vjust = 0.5),              
                                    ggtheme = theme_pubr(),
                                    ylab = FALSE) +
        theme(legend.key.height = unit(1.5, "cm"),
              legend.key.width = unit(1, "cm"),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              plot.margin = margin(1, 1, 1, 1, "cm"))
      
    }
    
    #reactom
    {
      pathways_reactom <- msigdbr("Homo sapiens", "C2", "CP:REACTOME") %>%
        format_pathways()
      
      results_reactom <- compare_seurat(seu,
                                        assay = "SCT",
                                        group1 = fgroup1, 
                                        group1_population = fgroup1_population,
                                        pathways = pathways_reactom,
                                        parallel = TRUE,
                                        cores = 4)   %>%
        mutate(Pathway = gsub("^REACTOME_", "", Pathway),    # remove "REACTOME_" 
               Pathway = tolower(Pathway),               # convert to lower case
               Pathway = gsub("_", " ", Pathway),        # replace "_" with ""
               Pathway = str_to_title(Pathway))          # convert to title case
      
      
      
      
      results_reactom <- results_reactom %>%
        filter(FC > -0.3) %>%
        arrange(desc(qval)) %>%
        slice_head(n = 10)
      
      
      p_pathways_reactom <- ggdotchart(results_reactom, x = "Pathway", y = "FC", title = "Pathway Analysis by REACTOME",
                                       color = "adjPval",
                                       sorting = "descending",
                                       add = "segments",
                                       rotate = TRUE,
                                       dot.size = 10,
                                       gradient.cols = c("darkblue", "grey"),
                                       font.label = list(color = "white", size = 9, vjust = 0.5),              
                                       ggtheme = theme_pubr(),
                                       ylab = FALSE) +
        theme(legend.key.height = unit(1.5, "cm"),
              legend.key.width = unit(1, "cm"),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              plot.margin = margin(1, 1, 1, 1, "cm"))
    }
    
    
    print("These are the kegg results:")
    print(results_kegg)
    print("These are the reactom results:")
    print(results_reactom)
    
    # Return a list with the two plots
    result_list <- list(p_pathways_kegg = p_pathways_kegg,
                        p_pathways_reactom = p_pathways_reactom,
                        results_kegg = results_kegg,
                        results_reactom = results_reactom)
    
    return(result_list)
    
    
    
  }
  
  get_nk_subtype <- function(seu) {
    
    library(BiocParallel)
    rucket_reference <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/data_for_get_nk_subtype/seu_rueckert_refrence.qs")
    
    DefaultAssay(seu) <- "RNA"
    
    # SingleR
    pred_fine <- as.SingleCellExperiment(seu) %>%
      SingleR(ref = rucket_reference, 
              labels = rucket_reference$type, 
              assay.type.ref = "counts", 
              assay.type.test = "counts",
              de.method = "wilcox", 
              BPPARAM = MulticoreParam(100))
    
    # Fine tuning + updating Seurat
    {
      adaptive_indices <- which(pred_fine$labels == "Adaptive" & pred_fine$delta.next < 0.11)
      
      # Reassign these cells to "Dim"
      pred_fine$labels[adaptive_indices] <- "Dim"
      
      # If you've already stored these labels in your Seurat object:
      seu$SingleR_NK_sub_labels <- pred_fine$labels
      }
    
    
    print(plotDeltaDistribution(pred_fine, show = "delta.next"))
    print(plotScoreHeatmap(pred_fine))
    
    
    table(seu$SingleR_NK_sub_labels)
    
    color_palette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
    print(dittoDimPlot(seu, var = "SingleR_NK_sub_labels", color.panel = color_palette))
    
    # My_dotplot
    solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')
    my_dotplot <- DotPlot(seu, 
                          features = c("GZMK", "XCL1", "IL7R", "TCF7", "GPR183",
                                       "GZMB", "PRF1", "CX3CR1", "CD7", "FCER1G", "KLRB1",
                                       "KLRC2", "CD3E", "PATL2", "ZBTB38"), 
                          group.by = "SingleR_NK_sub_labels", 
                          dot.min = 0.01) + 
      scale_color_gradientn(colours = solarExtra) + 
      theme(axis.text.x = element_text(angle = 45))
    
    print(my_dotplot)
    
    my_dotplot <- DotPlot(seu, 
                          features = c("KLRC2", "FCGR3A", "B3GAT1", "KIR3DL1", 
                                       "FCER1G", "NCR1", "IFNG"), 
                          group.by = "SingleR_NK_sub_labels", 
                          dot.min = 0.01) + 
      scale_color_gradientn(colours = solarExtra) + 
      theme(axis.text.x = element_text(angle = 45))
    print(my_dotplot)
    
    
    return(seu)
  }
  
  
  get_percentage_barplot <- function(seu, groupby, data_name){
    
    title =paste("Distribution of NK cells subtypes across conditions in the",data_name,"data")
    color_palette <- c("#E69F00", "#56B4E9", "#009E73",  "#CC79A7")
    
    
    # Create the dittoBarPlot and extract the plot and data
    plot_data <- dittoBarPlot(seu, "SingleR_NK_sub_labels", group.by = groupby, main = title,  color.panel = color_palette)
    
    data <- plot_data$data %>%
      group_by(grouping) %>%
      arrange(label) %>%
      mutate(mid_y = case_when(
        row_number() == 1 ~ 1 - (percent / 2),
        row_number() == 2 ~ (1 - lag(percent)) - (percent / 2),
        row_number() == 3 ~ (1 - lag(percent, 2) - lag(percent, 1)) - (percent / 2),
        row_number() == 4 ~ (1 - lag(percent, 3) - lag(percent, 2)- lag(percent, 1)) - (percent / 2)    
      ))
    
    
    # cumsum(percent) - 0.5 * percent
    plot_data <- plot_data + geom_text(data = data, aes(x = grouping, y = mid_y, label = sprintf("%.1f%%", percent * 100)), color = "white", size = 6, fontface = "bold")
    
    return (plot_data)
    
  }
  
  
  get_update_nk_prcnt <- function(nk_prcnt,seu, sample, data) {
    
    
    if (!(sample %in% unique(colnames(seu@meta.data)))) {
      print(paste("Column", sample, "does not exist in the data frame."))
      break
    } else if (!("condition" %in% unique(colnames(seu@meta.data)))) {
      
      print("Column condition does not exist in the data frame.")
      break
    }
    
    
    
    else {
      
      # Get unique conditions
      unique_samples <- unique(seu@meta.data[[sample]])
      seu$temp <- seu@meta.data[[sample]]
      
      i=0
      # Loop through each unique condition
      for (current_sample in unique_samples) {
        # Filter seu_mevorach for the current condition
        condition_data <- subset(seu, temp == current_sample)
        
        # Get dataset and disease (which is the same as condition in this case)
        dataset <- unique(condition_data$orig.ident)
        disease <- unique(condition_data$condition)
        
        # Get the percentages from the 'data' table
        dim_prcnt <- data %>% filter(label == "Dim" & grouping == current_sample)  %>% pull(percent)
        bright_prcnt <- data %>% filter(label == "Bright" & grouping == current_sample) %>% pull(percent)
        adaptive_prcnt <- data %>% filter(label == "Adaptive" & grouping == current_sample) %>% pull(percent)
        proliferating_prcnt <- data %>% filter(label == "Proliferating" & grouping == current_sample) %>% pull(percent)  
        
        
        # If no data is found, set percentages to 0
        if (length(dim_prcnt) == 0) dim_prcnt <- 0
        if (length(bright_prcnt) == 0) bright_prcnt <- 0
        if (length(adaptive_prcnt) == 0) adaptive_prcnt <- 0
        if (length(proliferating_prcnt) == 0) proliferating_prcnt <- 0
        
        # Create a new row
        new_row <- data.frame(
          sample = current_sample,
          dataset = dataset,
          disease = disease,
          dim_prcnt = dim_prcnt,
          bright_prcnt = bright_prcnt,
          adaptive_prcnt = adaptive_prcnt,
          prolifirating_prcnt = proliferating_prcnt,
          stringsAsFactors = FALSE
        )
        
        # Append the new row to nk_prcnt
        nk_prcnt <- rbind(nk_prcnt, new_row)
        i=i+1
      }
      
    }
    
    print(paste("nk_prcnt get updated",i,"times")) 
    return (nk_prcnt)
    
  }
  
  get_confusion_matrix <- function(seu_cd59, seu_cd55, method_column) {
    
    #Checking that the method_column indeed exsists in each of the seurat objects
    if (!(method_column %in% unique(colnames(seu_cd59@meta.data)))) {
      print(paste("Column", method_column, "does not exist in the seu_cd59 data."))
      break }
    else if (!(method_column %in% unique(colnames(seu_cd55@meta.data)))) {
      print(paste("Column", method_column, "does not exist in the seu_cd55 data."))
      break } 
    
    
    else {
      
      # Create a temp column that will hold the method_column variables
      seu_cd59$temp <- seu_cd59@meta.data[[method_column]]
      seu_cd55$temp <- seu_cd55@meta.data[[method_column]]
      
      # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), True Negatives (TN)
      tp <- sum(seu_cd59$temp == "NK cells" & seu_cd59$gold_standart == "NK cells") +
        sum(seu_cd55$temp == "NK cells" & seu_cd55$gold_standart == "NK cells")
      
      fp <- sum(seu_cd59$temp == "NK cells" & seu_cd59$gold_standart != "NK cells") +
        sum(seu_cd55$temp == "NK cells" & seu_cd55$gold_standart != "NK cells")
      
      fn <- sum(seu_cd59$temp != "NK cells" & seu_cd59$gold_standart == "NK cells") +
        sum(seu_cd55$temp != "NK cells" & seu_cd55$gold_standart == "NK cells")
      
      tn <- sum(seu_cd59$temp != "NK cells" & seu_cd59$gold_standart != "NK cells") +
        sum(seu_cd55$temp != "NK cells" & seu_cd55$gold_standart != "NK cells")
      
      precision <- tp / (tp + fp)
      sensitivity <- tp / (tp + fn)
      f1_score <- 2 * (precision * sensitivity) / (precision + sensitivity)
      accuracy <- (tp + tn) / (tp + tn + fp + fn)
      
      # Create a matrix for the confusion matrix
      confusion_matrix <- matrix(c(tp, fn, fp, tn), nrow = 2, byrow = TRUE,
                                 dimnames = list(
                                   "Actual" = c("Actual Positive", "Actual Negative"),
                                   "Predicted" = c("Predicted Positive", "Predicted Negative")
                                 ))
      
      # Create a data frame for the metrics
      metrics <- data.frame(
        Metric = c("Sensitivity", "Precision", "F1 Score", "Accuracy"),
        Value = c(sensitivity, precision, f1_score, accuracy)
      )
      
      # Combine the confusion matrix and metrics into a list
      output <- list(
        ConfusionMatrix = confusion_matrix,
        Metrics = metrics
      )
      
      # Return the output
      return(output)
    }
    
  }
  
  
  get_method_umap <- function (seu_cd55, seu_cd59, method_column){
    
    #Checking that the method_column indeed exsists in each of the seurat objects
    if (!(method_column %in% unique(colnames(seu_cd59@meta.data)))) {
      print(paste("Column", method_column, "does not exist in the seu_cd59 data."))
      break }
    else if (!(method_column %in% unique(colnames(seu_cd55@meta.data)))) {
      print(paste("Column", method_column, "does not exist in the seu_cd55 data."))
      break } 
    
    
    else {
      
      # Create a temp column that will hold the method_column variables
      # (It's the only way to access it)
      seu_cd59$temp <- seu_cd59@meta.data[[method_column]]
      seu_cd55$temp <- seu_cd55@meta.data[[method_column]]
      
      size <- 0.8
      
      cd55_method_title <- paste("CD55 UMAP by",method_column)
      cd59_method_title <- paste("CD59 UMAP by",method_column)
      
      
      
      p_sample_cd55 <- dittoDimPlot(seu_cd55, "gold_standart", reduction.use = "umap", size = size) +
        ggplot2::ggtitle(label = "CD55 UMAP by Gold standard")
      
      p_sample_cd59 <- dittoDimPlot(seu_cd59, "gold_standart", reduction.use = "umap", size = size) +
        ggplot2::ggtitle(label = "CD59 UMAP by Gold standard")
      
      p_temp_cd55 <- dittoDimPlot(seu_cd55, "temp", reduction.use = "umap", size = size) +
        ggplot2::ggtitle(label = cd55_method_title)
      
      p_temp_cd59 <- dittoDimPlot(seu_cd59, "temp", reduction.use = "umap", size = size) +
        ggplot2::ggtitle(label = cd59_method_title)
      
      # Combine plots into one
      p_final <- ((p_sample_cd55 + p_temp_cd55) / (p_sample_cd59 + p_temp_cd59))
      
    }
  }
  
  
  get_nk <- function(seu){
    
    
    
    seu <- FindClusters(seu, resolution = 80)
    
    #SingleR with Blueprint
    {
      singler_ref <- celldex::BlueprintEncodeData()
      
      
      singler_results <- SingleR(test = GetAssayData(seu, assay = "SCT"),
                                 ref = singler_ref,
                                 labels = singler_ref@colData@listData$label.main,
                                 clusters = seu@meta.data$SCT_snn_res.80
      )
      
      x <- data.frame(row.names = rownames(singler_results), labels = singler_results$labels)
      
      seu <- AddMetaData(seu, metadata = x[seu$SCT_snn_res.80, 1], col.name = "SingleR.labels.blueprint.cluster80")
      
    }
    
    #SingleR with Monaco
    
    {
      singler_ref <- celldex::MonacoImmuneData()
      
      
      singler_results <- SingleR(test = GetAssayData(seu, assay = "SCT"),
                                 ref = singler_ref,
                                 labels = singler_ref@colData@listData$label.main,
                                 clusters = seu@meta.data$SCT_snn_res.80
      )
      
      x <- data.frame(row.names = rownames(singler_results), labels = singler_results$labels)
      
      seu <- AddMetaData(seu, metadata = x[seu$SCT_snn_res.80, 1], col.name = "SingleR.labels.monaco.cluster80")
    }
    
    #Ucell with Viver's NK Signature
    
    {
      positive_nk_signature <- c("CD160+", "CD244+", "CHST12+", "CST7+", "GNLY+", "IL18RAP+", "IL2RB+", "KLRC1+", "KLRC3+", "KLRD1+", "KLRF1+", "PRF1+", "XCL2+")
      
      
      signatures <- list(
        positive_nk_signature = positive_nk_signature)
      
      seu <- AddModuleScore_UCell(seu, features = signatures, name = NULL)
      
      seu$nk_by_ucell <- ifelse(seu$positive_nk_signature >= 0.4, "NK cells", "other")
      
    }
    
    # Final idetification
    
    seu$is_nk <- dplyr::case_when(
      (seu$nk_by_ucell == "NK cells" & seu$SingleR.labels.blueprint.cluster80 == "NK cells") |
        (seu$nk_by_ucell == "NK cells" & seu$SingleR.labels.monaco.cluster80 == "NK cells") |
        (seu$SingleR.labels.blueprint.cluster80 == "NK cells" & seu$SingleR.labels.monaco.cluster80 == "NK cells") ~ "NK cells",
      TRUE ~ "other"
    )
    
    return(seu)
  }
  
  
}


# 2.2 Data Extraction & Matrix Merging Utilities
#Functions for extraction

{
  # Download function with quite mode and 3 attempts to download each file (it assumes network problems etc)
  download_file <- function(url, output_file) {
    retries <- 3
    for(i in 1:retries) {
      command <- paste("wget -q --show-progress -O", shQuote(output_file), shQuote(url))
      exit_status <- system(command)
      if (exit_status == 0) {
        break
      }
    }
  } 
  
  is_tar_gz <- function(file_name) {
    grepl("\\.tar\\.gz$", file_name)
  }
  
  one_verified_cm <- function(dir_name) {
    
    list_sample <- list.files(path = dir_name, recursive = FALSE)
    list_sample <- list_sample[!file.info(file.path(dir_name, list_sample))$isdir]  # Filter out directories
    list_sample <- sort(list_sample[grepl("_cm.qs", list_sample)])  # Keep only files with "_cm.qs" in their names
    
    
    # Load all the count matrixes to my environment 
    for (i in seq_along(list_sample)) {
      variable_name <- gsub(".qs", "", list_sample[i])  # Remove ".qs" from the filename
      assign(variable_name, qread(file.path(dir_name, list_sample[i])))  # Use full path
    }
    
    
    #Let's change list_sample so it will have the current names of the count matrixes loaded in the enviorment
    for (i in seq_along(list_sample)) {
      list_sample[i] <- gsub(".qs", "", list_sample[i]) # Remove ".qs" from the filename
    }
    #User verification and confermation  
    cat(paste0("\033[31m", "\nThese are all ", length(list_sample), " count matrices I'm going to inspect before merging:\n", paste(list_sample, collapse=", "), "\033[39m"))
    flush.console()
    response <- readline(prompt = paste0("\033[31m", "Do you wish to continue? (enter 'y' for yes or 'n' for no):\n"))
    
    if(tolower(response) != "y") {
      cat("Operation aborted.\n")
      return(invisible())
    }
    
    
    # Checking that all the rows (genes) are the same in all the count matrixes
    {
      cat("\nChecking that all the rows (genes) are the same in all the count matrixes\n")
      
      # Initialize a list to store pairs of matrices with different rownames
      mismatched_pairs <- list()
      
      # Nested loop to compare each matrix with every other matrix
      for (i in 1:(length(list_sample) - 1)) {
        for (j in (i + 1):length(list_sample)) {
          matrix1_name <- list_sample[i]
          matrix2_name <- list_sample[j]
          matrix1 <- get(matrix1_name)
          matrix2 <- get(matrix2_name)
          
          if (!all(rownames(matrix1) == rownames(matrix2))) {
            mismatched_pairs <- c(mismatched_pairs, paste(matrix1_name, matrix2_name, sep = " and "))
          }
        }
      }
      
      # Check if any mismatches were found
      if (length(mismatched_pairs) > 0) {
        cat("Mismatched rownames found between the following pairs of count matrices:\n")
        print(mismatched_pairs)
        return(invisible())} 
      else {
        cat("All count matrices have matching rownames.")
      }
    }
    
    # Creating a list of all the sample ids
    sample_ids <- gsub("_cm$", "", list_sample)
    
    # Making every column (cell) in each sample unique using the sample name
    # And printing the first 5 colnames of each sample to check that it worked
    {
      cat("\nMaking every column (cell) in each sample unique using the sample name\n")
      for (i in seq_along(list_sample)) {
        matrix1_name <- list_sample[i]
        matrix1 <- get(matrix1_name)
        
        # Prepend the unique ID to the column names
        colnames(matrix1) <- paste(sample_ids[i], colnames(matrix1), sep = "_")
        
        # Assign the modified matrix back to its original name in the environment
        assign(matrix1_name, matrix1)
        
        # Print the first five column names for verification
        if (ncol(matrix1) >= 5) {
          first_five_colnames <- colnames(matrix1)[1:5]
        } else {
          first_five_colnames <- colnames(matrix1)
        }
        cat("For sample",i,",here are the first five colnames\n",first_five_colnames,"\n")
      }
      
    }
    
    # Checking that each column in unique comparing to the rest of the data set
    # (So we will have no BALAGAN upon binding the count matrixes together)
    {
      cat("\nChecking that each column in unique comparing to the rest of the data set\n")
      
      # Initialize a list to keep track of all column names
      all_col_names <- vector("list", length(list_sample))
      names(all_col_names) <- list_sample
      
      # Loop to store column names of each matrix
      for (i in 1:length(list_sample)) {
        matrix_name <- list_sample[i]
        matrix <- get(matrix_name)
        all_col_names[[matrix_name]] <- colnames(matrix)
      }
      
      # Check for duplicate column names across matrices
      duplicate_found <- FALSE
      for (i in 1:(length(list_sample) - 1)) {
        for (j in (i + 1):length(list_sample)) {
          intersect_cols <- intersect(all_col_names[[list_sample[i]]], all_col_names[[list_sample[j]]])
          if (length(intersect_cols) > 0) {
            cat("Duplicate column names found between", list_sample[i], "and", list_sample[j], ":\n", paste(intersect_cols, collapse = ", "), "\n")
            duplicate_found <- TRUE
          }
        }
      }
      
      if (!duplicate_found) {
        cat("No duplicate column names found across matrices.")
      } else{
        cat("Operation aborted since duplicates were found.\n")
        return(invisible())
      }
      
      
    }
    
    #Finally, let's merge all these count matrixes into one big count matrix
    {
      
      cat(paste0("\033[31m","\nFinally, let's merge all these count matrixes into one big count matrix.\nPlease wait as it may take some time.\n"))
      
      # Assuming the first matrix in the list has the correct number of rows
      num_rows <- as.integer(nrow(get(list_sample[1])))
      
      # Initialize an empty sparse matrix with the correct dimensions
      combined_matrix <- Matrix(0, nrow = num_rows, ncol = 0, sparse = TRUE)
      
      # Iterate over each matrix name and bind them
      for (matrix_name in list_sample) {
        matrix <- get(matrix_name)
        combined_matrix <- cbind(combined_matrix, matrix)
      }
      
      # combined_matrix now contains all the matrices combined, let's save it!
      
      cat(paste0("\033[31m","\nReady to save the count matrix\n"))
      flush.console()
      response <- readline(prompt = "Please enter the name of the Seurat object to be saved:\n")
      cat(paste0("\033[31m", "\nCreating the cm\n"))
      
      # Check if a response was provided
      if (!is.na(response) && nzchar(response)) {
        # Save the file to the correct directory
        full_save_path <- file.path(dir_name, paste0("cm_", response, ".qs"))
        qsave(combined_matrix, full_save_path)
        
        # Inform the user where the file was saved
        cat(paste0("\033[31m","\n", response, "_cm.qs", " waits for you at ", dir_name, "\n"))
      } else {
        cat(paste0("\033[31m","\nNo valid name was entered. Saving aborted.\n"))
        return(invisible())
      }
      
      
      
    }
    
    
  }
  
  
  
}

# 3. Pre-Processed Healthy Cohort Databases (Fast Pass)
# 3.1 Loading Compiled Datasets
#Fast Pass to all the databases

{
# 4. Reference Data Ingestion & Processing (ETL)

# 4.1 Yang Cohort

#Yang
  seu_yang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
  
# 4.2 Witkowski Cohort

#Witkowski
  
  seu_witkowski_h6 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h6.qs")
  seu_witkowski_h8 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h8.qs")
  seu_witkowski_h21_h18_h2 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h_2_18_21.qs")
  
# 4.3 Smith Cohort

#Smith
  
  seu_smith_NRHV193_NRHV187 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith_NRHV193_NRHV187.qs")
  seu_smith_em <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith_em.qs")
  
# 4.4 Degagne Cohort
  
#Degagne
  seu_degagne <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_degagne.qs")
  
# 4.5 Zhang Cohort

#Zhang
  seu_zhang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
  
# 4.6 Gou Cohort

#Gou
  seu_gou <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou.qs")
  
}

# 4. Reference Data Ingestion & Processing (ETL)


# 4.1 Yang Cohort
#Yang

{

  untar("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/GSE130430_RAW.tar", exdir ="/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/")

  cm_ <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/yang/GSM3738542/")
  str(cm_)
  qsave(cm_,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/GSE130430/25y_F_2_PB_cm.qs")
  cm_ <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/yang/GSM3738543/")
  str(cm_)
  qsave(cm_,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/yang/26y_F_PB_cm.qs")

  one_verified_cm("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/yang/")
  
  cm_yang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cm_yang.qs")
 
  seu_yang <- CreateSeuratObject(cm_yang)
  qsave(seu_yang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
  
  seu_kang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
  seu <- seu_kang
  seu$sample <- seu$orig.ident
  seu$orig.ident <- "Yang"
  seu$condition <- "healthy"
  seu$sample <- ifelse(seu$sample == "25y" ,"25y_F_2_PB","26y_F_PB")
  qsave(seu,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
  
  #General pipeline
  {  
    seu_kang <- seu
    seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
    
    feats <- c("percent_mito")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = percent_mito < 20)
    
    
    
    options(future.globals.maxSize = 8 * 1024^3)  
    
    
    seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
    qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
    
    feats <- c("nFeature_SCT")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
    
    qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
    
    seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
    seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
    VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
    
    ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
    
    seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
    seu_kang <- FindClusters(seu_kang, resolution = 0.6)
    seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                        n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
    
    qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_yang.qs")
    
    
  }
}

# 4.2 Witkowski Cohort

#Witkowski
{
  # library(Seurat)
  # library(dplyr)
  # library(Matrix)
  # library(rhdf5)
  #H6
  {
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5584154&format=file&file=GSM5584154%5FNK%5Frun1%2Efiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584154_NK_run1.filtered_feature_bc_matrix.h5")
    
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5584154&format=file&file=GSM5584154%5FNK%5Frun1%2Etxt%2Egz",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584154_NK_run1.txt.gz")
    
    gunzip("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H6/GSM5584154_NK_run1.txt.gz", remove = FALSE)
    
    
    gene_expression_data <- Read10X_h5("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H6/GSM5584154_NK_run1.filtered_feature_bc_matrix.h5")
    
    
    # Load the metadata from the txt file
    metadata <- read.table("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H6/GSM5584154_NK_run1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Make sure the Barcode column is the row names of the metadata
    rownames(metadata) <- metadata$Barcode
    
    
    # Create the Seurat object
    seurat_obj <- CreateSeuratObject(counts = gene_expression_data$`Gene Expression`, project = "NK_cells", min.cells = 3, min.features = 200)
    
    # Add the metadata to the Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata)
    
    seurat_obj <- subset(seurat_obj,subset = Donor == "donor_H6")
    
    qsave(seurat_obj,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/seu_h6.qs")
    
    seu <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h6.qs")
    
    seu$orig.ident <- "Witkowski"
    seu$Barcode <- NULL
    seu$Timepoint <- NULL
    seu$Donor <- NULL
    seu$Status <- NULL
    seu$condition <- "healthy"
    
    seu_kang <- seu
    
    
    #General pipeline
    {  
      seu_kang <- seu
      seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
      
      feats <- c("percent_mito")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = percent_mito < 20)
      
      
      
      options(future.globals.maxSize = 8 * 1024^3)  
      
      
      seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
      qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h6.qs")
      
      feats <- c("nFeature_SCT")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
      
      qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h6.qs")
      
      seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                   s.features = cc.genes$s.genes)
      seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
      VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
      
      ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
      
      seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
      seu_kang <- FindClusters(seu_kang, resolution = 0.6)
      seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                          n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
      
      qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h6.qs")
      
      
    }
    
    
  }
  
  #H8
  {  
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5584155&format=file&file=GSM5584155%5FNK%5Frun2%2Efiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584155_NK_run2.filtered_feature_bc_matrix.h5")
    
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5584155&format=file&file=GSM5584155%5FNK%5Frun2%2Etxt%2Egz",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584155_NK_run2.txt.gz")
    
    gunzip("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584155_NK_run2.txt.gz", remove = FALSE)
    
    
    
    gene_expression_data <- Read10X_h5("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H8/GSM5584155_NK_run2.filtered_feature_bc_matrix.h5")
    
    
    # Load the metadata from the txt file
    metadata <- read.table("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H8/GSM5584155_NK_run2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Make sure the Barcode column is the row names of the metadata
    rownames(metadata) <- metadata$Barcode
    
    
    # Create the Seurat object
    seurat_obj <- CreateSeuratObject(counts = gene_expression_data$`Gene Expression`, project = "NK_cells", min.cells = 3, min.features = 200)
    
    # Add the metadata to the Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata)
    
    seurat_obj <- subset(seurat_obj,subset = Donor == "donor_H8")
    
    qsave(seurat_obj,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/seu_h8.qs")
    
    seu <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h8.qs")
    
    seu$orig.ident <- "Witkowski"
    seu$Barcode <- NULL
    seu$Timepoint <- NULL
    seu$Donor <- NULL
    seu$Status <- NULL
    seu$condition <- "healthy"
    seu$sample <- "H8"
    
    #General pipeline
    {  
      seu_kang <- seu
      seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
      
      feats <- c("percent_mito")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = percent_mito < 20)
      
      
      
      options(future.globals.maxSize = 8 * 1024^3)  
      
      
      seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
      qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h8.qs")
      
      feats <- c("nFeature_SCT")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
      
      qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h8.qs")
      
      seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                   s.features = cc.genes$s.genes)
      seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
      VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
      
      ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
      
      seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
      seu_kang <- FindClusters(seu_kang, resolution = 0.6)
      seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                          n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
      
      qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h8.qs")
      
      
    }
  }
  
  # H21 H18 H2
  {  
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5584156&format=file&file=GSM5584156%5FNK%5Frun3%2Efiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584156_NK_run3.filtered_feature_bc_matrix.h5")
    
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5584156&format=file&file=GSM5584156%5FNK%5Frun3%2Etxt%2Egz",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584156_NK_run3.txt.gz")
    
    gunzip("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/GSM5584156_NK_run3.txt.gz", remove = FALSE)
    
    
    
    gene_expression_data <- Read10X_h5("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H18,H21,H2/GSM5584156_NK_run3.filtered_feature_bc_matrix.h5")
    
    
    # Load the metadata from the txt file
    metadata <- read.table("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/H18,H21,H2/GSM5584156_NK_run3.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Make sure the Barcode column is the row names of the metadata
    rownames(metadata) <- metadata$Barcode
    
    
    # Create the Seurat object
    seurat_obj <- CreateSeuratObject(counts = gene_expression_data$`Gene Expression`, project = "NK_cells", min.cells = 3, min.features = 200)
    
    # Add the metadata to the Seurat object
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata)
    
    seurat_obj_h18 <- subset(seurat_obj,subset = Donor %in% c("donor_H18", "donor_H21", "donor_H2"))
    
    qsave(seurat_obj_h18,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Witkowski/seu_h_2_18_21.qs")
    
  seu <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h_2_18_21.qs")
    

  seu$orig.ident <- "Witkowski"
  seu$Barcode <- NULL
  seu$Timepoint <- NULL
  seu$Status <- NULL
  seu$condition <- "healthy"
  
  seu$sample <- case_when(seu$Donor == "donor_H18" ~ "H18",
                          seu$Donor == "donor_H21" ~ "H21",
                          seu$Donor == "donor_H2" ~ "H2")
  
  seu$Donor <- NULL
    
  seu_kang <- seu
  qsave (seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h_2_18_21.qs")
    
  #General pipeline
  {  
    seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
    
    feats <- c("percent_mito")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = percent_mito < 20)
    
    
    
    options(future.globals.maxSize = 8 * 1024^3)  
    
    
    seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
    qsave (seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h_2_18_21.qs")
    
    feats <- c("nFeature_SCT")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
    
    qsave (seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h_2_18_21.qs")
    
    seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
    seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
    VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
    
    ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
    
    seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
    seu_kang <- FindClusters(seu_kang, resolution = 0.6)
    seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                        n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
    
    qsave (seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_witkowski_h_2_18_21.qs")
    
    
  }
  }
  
}

# 4.3 Smith Cohort
#Smith
{
  #EM_HD1_NK_scRNAseq
  {
    
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144430&format=file",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/GSE144430_RAW.tar")

    untar(tarfile = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/GSE144430_RAW.tar", exdir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/")
    
    cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/GSE144430/")
    
    qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/smith_EM_cm.qs")
    
    seu_kang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith_em.qs")
    #General pipeline
    {  
      seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
      
      feats <- c("percent_mito")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = percent_mito < 20)
      
      
      
      options(future.globals.maxSize = 8 * 1024^3)  
      
      
      seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith_em.qs")
      
      feats <- c("nFeature_SCT")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
      
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith_em.qs")
      
      seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                   s.features = cc.genes$s.genes)
      seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
      VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
      
      ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
      
      seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
      seu_kang <- FindClusters(seu_kang, resolution = 0.6)
      seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                          n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
      
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith_em.qs")
      
      
    }   
  }
  
  # NRHV187_unstim and NRHV193_unstim
  {
    download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144191&format=file&file=GSE144191%5Fprocessed%5Fseurat%5Fobject%2Erds%2Egz",
                  "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/GSE144191_processed_seurat_object.rds.gz	")
    
    gunzip("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/GSE144191_processed_seurat_object.rds.gz", remove = FALSE)
    
    seurat_obj <- readRDS("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Smith/GSE144191_processed_seurat_object.rds")
    
    seu <- UpdateSeuratObject(seurat_obj)
    
    qsave(seu,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith.qs")
    
    seu <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith.qs")
    
    seu$percent.mito <- NULL
    seu$res.0.8 <- NULL
    seu$var.ratio.pca <- NULL
    seu$stim_cond <- NULL
    seu$cluster.stim <- NULL
    seu$donor <- NULL
    
    seu$orig.ident <- "Smith"
    seu$sample <- seu$donor_cond
    seu$donor_cond <- NULL
    seu$condition <- "healthy"
    
    seu_kang <- seu
    seu_kang@assays$RNA$scale.data  <- NULL
    seu_kang@reductions$cca <- NULL 
    seu_kang@reductions$cca.aligned <- NULL 
    seu_kang@reductions$tsne <- NULL 
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith.qs")
    
    #General pipeline
    {  
      seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
      
      feats <- c("percent_mito")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = percent_mito < 20)
      
      
      
      options(future.globals.maxSize = 8 * 1024^3)  
      
      
      seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith.qs")
      
      feats <- c("nFeature_SCT")
      VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
        NoLegend()
      
      seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
      
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith.qs")
      
      seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                   s.features = cc.genes$s.genes)
      seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
      VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
      
      ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
      
      seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
      seu_kang <- FindClusters(seu_kang, resolution = 0.6)
      seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                          n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
      
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_smith.qs")
      
    
  }
 
    
  }
}

# 4.4 Degagne Cohort

#Degagne
{
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7933938&format=file&file=GSM7933938%5FscCB000165%5Fbarcodes%2Etsv%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/GSM7933938_scCB000165_barcodes.tsv.gz")
  
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7933938&format=file&file=GSM7933938%5FscCB000165%5Ffeatures%2Etsv%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/GSM7933938_scCB000165_features.tsv.gz")
  
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM7933938&format=file&file=GSM7933938%5FscCB000165%5Fmatrix%2Emtx%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/GSM7933938_scCB000165_matrix.mtx.gz")
  
  
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE249316&format=file&file=GSE249316%5Fcell%5Fmetadata%2Etxt%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/GSE249316_cell_metadata.txt.gz")
  

  #  Read the gzipped metadata file into a variable -
  # Not usefull - weird format and has redundent information
  # cell_metadata <- read.table(gzfile("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/metadata.txt.gz"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/")
  qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/degagne_cm.qs")
  
  cm <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Degagne/degagne_cm.qs")
  
  seu_degagne <- CreateSeuratObject(cm)
  
  seu_degagne$orig.ident <- "Degagne"
  seu_degagne$condition <- "healthy"
  seu_degagne$sample <- "pNK_pre_infusion"
  
  qsave(seu_degagne, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_degagne.qs")
  #General pipeline
  {
    seu_kang <- seu_degagne
    seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
    
    feats <- c("percent_mito")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = percent_mito < 20)
    
    
    
    options(future.globals.maxSize = 8 * 1024^3)  
    
    
    seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_degagne.qs")
    
    feats <- c("nFeature_SCT")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
    
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_degagne.qs")
    
    seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
    seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
    VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
    
    ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
    
    seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
    seu_kang <- FindClusters(seu_kang, resolution = 0.6)
    seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                        n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
    
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_degagne.qs")
    
    
  }  
}

# 4.5 Zhang Cohort

#Zhang
{
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4798208&format=file&file=GSM4798208%5FD0%2D10%2Dcrl%5Fbarcodes%2Etsv%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Zhang/barcodes.tsv.gz")
  
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4798208&format=file&file=GSM4798208%5FD0%2D10%2Dcrl%5Ffeatures%2Etsv%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Zhang/features.tsv.gz")
  
  download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4798208&format=file&file=GSM4798208%5FD0%2D10%2Dcrl%5Fmatrix%2Emtx%2Egz",
                "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Zhang/matrix.mtx.gz")
  

  cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Zhang/")
  cm_rna <- cm$`Gene Expression`
  seu_zhang <- CreateSeuratObject(cm_rna)
  adt_assay <- CreateAssay5Object(counts = cm$`Antibody Capture`)
  
  seu_zhang[["ADT"]] <- adt_assay
  qsave(seu_zhang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
  qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Zhang/cm_zhang.qs")
  
  seu_kang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
  
  seu_kang$orig.ident <- "Zhang"
  seu_kang$condition <- "healthy"
  seu_kang$sample <- "D0-10-Crl"
  qsave(seu_kang, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
  
  #General pipeline
  {
    seu_kang <- PercentageFeatureSet(seu_kang, "^MT-", col.name = "percent_mito")
    
    feats <- c("percent_mito")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = percent_mito < 20)
    
    
    
    options(future.globals.maxSize = 8 * 1024^3)  
    
    
    seu_kang <- SCTransform(seu_kang, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = F)
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
    
    feats <- c("nFeature_SCT")
    VlnPlot(seu_kang, group.by = "orig.ident", features = feats, ncol = 1) +
      NoLegend()
    
    seu_kang <- subset(seu_kang, subset = nFeature_SCT > 500)
    
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
    
    seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
                                 s.features = cc.genes$s.genes)
    seu_kang <- RunPCA(seu_kang, npcs = 50, verbose = F)
    VizDimLoadings(seu_kang, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
    
    ElbowPlot(seu_kang, reduction = "pca", ndims = 50)
    
    seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
    seu_kang <- FindClusters(seu_kang, resolution = 0.6)
    seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                        n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
    
    qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_zhang.qs")
    
    
  }  
}

# 4.6 Gou Cohort

#Gou
{
  #Nk
  {
    #NK_elder
    {
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970633&format=file&file=GSM5970633%5FNK%5FElder%5Fbarcodes%2Etsv%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_elder/barcodes.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970633&format=file&file=GSM5970633%5FNK%5FElder%5Ffeatures%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_elder/features.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970633&format=file&file=GSM5970633%5FNK%5FElder%5Fmatrix%2Emtx%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_elder/matrix.mtx.gz")
      
      cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_elder/")
      qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/nk_elder_cm.qs")
      
    }
    
    #NK_young
    {
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970634&format=file&file=GSM5970634%5FNK%5FYoung%5Fbarcodes%2Etsv%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_young/barcodes.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970634&format=file&file=GSM5970634%5FNK%5FYoung%5Ffeatures%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_young/features.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970634&format=file&file=GSM5970634%5FNK%5FYoung%5Fmatrix%2Emtx%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_young/matrix.mtx.gz")
      
      cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/NK_young/")
      qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/nk_young_cm.qs")
      
    }
    one_verified_cm("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/nk/")
    cm <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/cm_gou_nk.qs")
    
    seu <- CreateSeuratObject(cm)
    seu$barcode <- colnames(seu@assays$RNA)
    seu$sample <- sapply(strsplit(seu$barcode, "_"), function(x) paste(x[1:2], collapse = "_"))
    seu$barcode <- NULL
    seu$orig.ident <- "Gou"
    seu$condition <- "healthy"
    
    qsave(seu,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/seu_gou.qs")
    
    seu_kang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou.qs")
    seu_kang$orig.ident <- "Gou"
    seu_kang$condition <- "healthy"
    qsave(seu_kang)
    
    # General Pipeline
    {
      #This countmatrix had already under went QC and normalization,
      # as written in the GEO sample pages GSM5970634 and GSM5970633:
      # "For the initial QC step,
      # we created Seurat objects for the young and elderly groups and filtered out
      # the cells that expressed <400 genes, >3000 genes or > 10% mitochondrial genes.
      # The gene expression count matrix for each cell was normalized using the “pp.normalize_total” 
      # function of Scanpy and log-transformed to the product (base e) after the addition of a pseudocount of 1."
      
      #So, we will go stright to the PCA (CellCycleScoring is not working for some reason)
      
      # seu_kang <- CellCycleScoring(object = seu_kang, g2m.features = cc.genes$g2m.genes,
      #                              s.features = cc.genes$s.genes,)
      
      
      LayerData(seu_kang, assay = "RNA", layer = "scale.data") <- LayerData(seu_kang, assay = "RNA", layer = "counts")
      seu_kang <- FindVariableFeatures(seu_kang, selection.method = "vst", nfeatures = 2000)
      
      
      # Access the scale.data layer from the RNA assay
      scale_data <- LayerData(seu_kang, assay = "RNA", layer = "scale.data")
      
      # Replace NA values with 0 in the scale.data matrix
      scale_data[is.na(scale_data)] <- 0
      
      # Save the modified scale.data back into the Seurat object
      LayerData(seu_kang, assay = "RNA", layer = "scale.data") <- scale_data
      
      
      #manual PCA to overpass some Seurat bullshit
      
      {
        # Extract the scale.data matrix
        scale_data_matrix <- LayerData(seu_kang, assay = "RNA", layer = "scale.data")
        
        # Ensure that the scale.data matrix has no NAs (this step might be redundant based on previous replacements)
        scale_data_matrix[is.na(scale_data_matrix)] <- 0
        
        # Manually run PCA using prcomp
        pca_result <- prcomp(t(scale_data_matrix), rank. = 50)
        
        # Store the PCA results back into the Seurat object as a dimensional reduction
        seu_kang[["pca"]] <- CreateDimReducObject(
          embeddings = pca_result$x, 
          loadings = pca_result$rotation, 
          assay = "RNA", 
          key = "PC_"
        )
        
        # Check the PCA result by plotting
        DimPlot(seu_kang, reduction = "pca")
        
      }
      
      
      
      seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
      seu_kang <- FindClusters(seu_kang, resolution = 0.6)
      seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                          n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou.qs")  
    }
  }
                         
# 4.7 Lym Cohort

  #Lym
  {
    #Lym_elder
    {
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970631&format=file&file=GSM5970631%5FLymphocytes%5FElder%5Fbarcodes%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_elder/barcodes.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970631&format=file&file=GSM5970631%5FLymphocytes%5FElder%5Ffeatures%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_elder/features.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970631&format=file&file=GSM5970631%5FLymphocytes%5FElder%5Fmatrix%2Emtx%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_elder/matrix.mtx.gz")
      
      cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_elder/")
      qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/lym_elder_cm.qs")
      
    }
    
    #	Lym_young
    {
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970632&format=file&file=GSM5970632%5FLymphocytes%5FYoung%5Fbarcodes%2Etsv%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_young/barcodes.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970632&format=file&file=GSM5970632%5FLymphocytes%5FYoung%5Ffeatures%2Etsv%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_young/features.tsv.gz")
      
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5970632&format=file&file=GSM5970632%5FLymphocytes%5FYoung%5Fmatrix%2Emtx%2Egz",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_young/matrix.mtx.gz")
      
      cm <- Read10X(data.dir = "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/Lym_young/")
      qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/lym_young_cm.qs")
      
    }
    
    one_verified_cm("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/lymphocytes/")
    cm <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/Gou/cm_gou_lym.qs")
    
    seu <- CreateSeuratObject(cm)
    seu$barcode <- colnames(seu@assays$RNA)
    seu$sample <- sapply(strsplit(seu$barcode, "_"), function(x) paste(x[1:2], collapse = "_"))
    seu$barcode <- NULL
    seu$orig.ident <- "Gou"
    seu$condition <- "healthy"
    
    qsave(seu,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou_lym.qs")
    
    seu_kang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou_lym.qs")

    
    # General Pipeline
    {
      #This countmatrix had already under went QC and normalization,
      # as written in the GEO sample pages GSM5970634 and GSM5970633:
      # "For the initial QC step,
      # we created Seurat objects for the young and elderly groups and filtered out
      # the cells that expressed <400 genes, >3000 genes or > 10% mitochondrial genes.
      # The gene expression count matrix for each cell was normalized using the “pp.normalize_total” 
      # function of Scanpy and log-transformed to the product (base e) after the addition of a pseudocount of 1."
      
      #So, we will go stright to the PCA (CellCycleScoring is not working for some reason)
      
      # Load the Seurat object
      seu_kang <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou_lym.qs")
      
      # Normalize the RNA assay data
      seu_kang <- NormalizeData(seu_kang, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
      
      # Find the 2000 most variable features
      seu_kang <- FindVariableFeatures(seu_kang, selection.method = "vst", nfeatures = 2000)
      
      # Scale the data to create the scale.data layer
      seu_kang <- ScaleData(seu_kang, assay = "RNA", features = VariableFeatures(seu_kang))
      
      # Access the scale.data layer from the RNA assay
      scale_data <- LayerData(seu_kang, assay = "RNA", layer = "scale.data")
      
      # Replace NA values with 0 in the scale.data matrix (if necessary)
      scale_data[is.na(scale_data)] <- 0
      
     
      
      #manual PCA to overpass some Seurat bullshit
      
      {
        
        # Subset the scale.data matrix to include only the variable genes
        variable_genes <- VariableFeatures(seu_kang)
        scale_data_subset <- scale_data[variable_genes, ]
        
        # Perform PCA on the variable genes
        pca_result <- prcomp(t(scale_data_subset), rank. = 50)
        
        # Save the PCA result
        qsave(pca_result, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/pca_result_lym.qs")
        
        # Store the PCA results back into the Seurat object as a dimensional reduction
        seu_kang[["pca"]] <- CreateDimReducObject(
          embeddings = pca_result$x, 
          loadings = pca_result$rotation, 
          assay = "RNA", 
          key = "PC_"
        )
        
        # Check the PCA result by plotting
        DimPlot(seu_kang, reduction = "pca")
        
        
      }
      
      
      
      seu_kang <- FindNeighbors(seu_kang, dims = 1:25)
      seu_kang <- FindClusters(seu_kang, resolution = 0.6)
      seu_kang <- RunUMAP(seu_kang, reduction = "pca", dims = 1:25, n.components = 2, n.neighbors = 30,
                          n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou_lym.qs")  
    }
    
    #Identifying the NK cells
    {
     nk_signature <- c("CD3G-", "CD3D-", "CD3E-", "CD14-", "CD19-", "MS4A1-", "CD7+")
      
      
      signatures <- list(
        nk_signature = nk_signature)
      
      seu_kang <- AddModuleScore_UCell(seu_kang, features = signatures, name = NULL)
      
      #Plot
      p_marker_gou_lym <- FeaturePlot_scCustom(seurat_object = seu_kang, features = "nk_signature")  + 
        ggtitle("Gou Lymphocytes UMAP by Nk  Signature") +
        theme(plot.title = element_text(face = "plain"))
      
      p_clusters <- dittoDimPlot(seu_kang, "seurat_clusters", reduction.use = "umap", size = size) +
        ggplot2::ggtitle(label = "Gou Lymphocytes UMAP by Seurat Clusters")
      
      p_final <- p_marker_gou_lym | p_clusters

      
      print(p_final)
      
      qsave(seu_kang,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/seu_gou_lym.qs")
      
    }
    
  }
 

}

#GEO research query - "Healthy NK" AND "Homo sapiens"[porgn:__txid9606] AND "scRNA-seq" AND "PBMCS" 

seu <- seu_kang
size = 0.8
p_clusters <- dittoDimPlot(seu, "seurat_clusters", reduction.use = "umap", size = size) +
  ggplot2::ggtitle(label = "UMAP by Seurat Clusters")

p_sample <- dittoDimPlot(seu, "sample", reduction.use = "umap", size = size) +
  ggplot2::ggtitle(label = "UMAP by sample")

p_singler <- dittoDimPlot(seu, "orig.ident", reduction.use = "umap", size = size) +
  ggplot2::ggtitle(label = "UMAP by orig.ident")

p_condition <- dittoDimPlot(seu, "condition", reduction.use = "umap", size = size, opacity = 0.8) +
  ggplot2::ggtitle(label = "UMAP by Sample")

# Combine plots into one
p_final <- ((p_clusters + p_condition) / (p_singler + p_sample))
print(p_final)


# 7. CITE-seq Benchmarking Datasets

{
  
  #Pont
  {
    #Download and untar
    {
      download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144434&format=file",
                    "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/GSE144434_RAW.tar")
      untar("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/GSE144434_RAW.tar", exdir ="/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/")
      
    }
    
    # create the seurat object
    {
      setwd("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/")
      
      # Read UMI gene expression counts and convert to sparase
      cm <- read.delim(gzfile("GSM4288255_Raw_UMI_genes.tsv.gz"), header = TRUE, row.names = 1)
      cm <- t(cm)
      seu_pont <- CreateSeuratObject(counts = cm, project = "Pont")
      
      cm_adt <- read.csv(gzfile("GSM4288255_Antibodies_count.csv.gz"), sep = "\t", row.names = 1)
      cm_adt <- as.matrix(cm_adt)
      cm_adt <- t(cm_adt)
      
      adt_assay <- CreateAssay5Object(counts = cm_adt)
      seu_pont[["ADT"]] <- adt_assay
      seu_pont[["ADT"]] <- CreateAssayObject(counts = cm_adt)
      
      qsave(seu_pont,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/seu_pont.qs")
      
 
      hto <- read.csv(gzfile("GSM4288255_HTO_identification.csv.gz"), sep = "\t", row.names = 1)
      seu_pont <- AddMetaData(
        object = seu_pont,
        metadata = hto,
        col.name = 'HTO'
      )
      
      umap_coordinate <- read.csv(gzfile("GSM4288255_UMAP_coordinate.csv.gz"), sep = "\t", row.names = 1)
      
      
      


      
      
    }

  }
  
  
  #Hao
  {
    #3P
    {
      #Download and untar
      {
        download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file",
                      "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/GSE164378_RAW.tar")
        
        untar("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/GSE164378_RAW.tar", exdir ="/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/")
        
        cm <- Read10X("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSM5008737_RNA_3P/")
        qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/cm_GSM5008737_RNA_3P.qs")
        
        cm_adt <- Read10X("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSM5008738_ADT_3P/")
        qsave(cm_adt,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/cm_adt_GSM5008738_ADT_3P.qs")
        
        cm_hto <- Read10X("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSM5008739_HTO_3P/")
        qsave(cm_hto,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/cm_hto_GSM5008739_HTO_3P.qs")
        
        seu_hao_3p <- CreateSeuratObject(counts=cm, assay = "RNA")
        seu_hao_3p[["ADT"]] <- CreateAssayObject(counts = cm_adt)
        seu_hao_3p[["HTO"]] <- CreateAssayObject(counts = cm_hto)
        
        qsave(seu_hao_3p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs")
        
        
        download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz",
                      "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/meta_data_3P.csv.gz")
        
        download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5FAntibody%5FHTO%5Finfo%2Exlsx",
                      "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSE164378_Antibody_HTO_info.xlsx")
        
        setwd("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/")
        meta_data <- read.csv(gzfile("meta_data_3P.csv.gz"))
        
        seu_hao_3p <- AddMetaData(object = seu_hao_3p, metadata = meta_data)
        
        seu_hao_3p$X <- NULL
        qsave(seu_hao_3p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs")
        
      }
      
      #Processing the seurat object:
      {
        
        seu_hao_3p <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs") 

        DefaultAssay(seu_hao_3p) <- 'RNA'
        seu_hao_3p <- PercentageFeatureSet(seu_hao_3p, pattern = "^MT-", col.name = "percent_mito")
        feats <- c("percent_mito")
        VlnPlot(seu_hao_3p, group.by = "orig.ident", features = feats, ncol = 1) +
          NoLegend()
        seu_hao_5p <- seu_hao_3p
        seu_hao_5p <- subset(seu_hao_5p, subset = percent_mito < 20)
        
        
        
        options(future.globals.maxSize = 32 * 1024^3)  # 32 GiB
        
        
        seu_hao_5p <- SCTransform(seu_hao_5p, vars.to.regress = "percent_mito", verbose = TRUE)
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs")
        
        seu_hao_5p <- RunPCA(seu_hao_5p, npcs = 50, verbose = T)
        VizDimLoadings(seu_hao_5p, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
        
        ElbowPlot(seu_hao_5p, reduction = "pca", ndims = 50)
        
        DefaultAssay(seu_hao_5p) <- 'ADT'
        # we will use all ADT features for dimensional reduction
        # we set a dimensional reduction name to avoid overwriting the 
        VariableFeatures(seu_hao_5p) <- rownames(seu_hao_5p[["ADT"]])
        seu_hao_5p <- NormalizeData(seu_hao_5p, normalization.method = 'CLR', margin = 2) %>% 
          ScaleData() %>% RunPCA(reduction.name = 'apca')
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs")
        ElbowPlot(seu_hao_5p, reduction = "apca", ndims = 50)
        
        
        seu_hao_5p <- FindMultiModalNeighbors(
          seu_hao_5p, reduction.list = list("pca", "apca"), 
          dims.list = list(1:25, 1:25), modality.weight.name = "RNA.weight"
        )
        
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs")
        
        seu_hao_5p <- RunUMAP(seu_hao_5p, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        seu_hao_5p <- FindClusters(seu_hao_5p, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
        
        plot_p5 <- DimPlot(seu_hao_5p, reduction = 'wnn.umap', group.by = 'celltype.l1', label = TRUE, repel = TRUE, label.size = 2.5)
        
        print(plot_p5)
        
        seu_hao_5p$gold_standart <- seu_hao_5p$celltype.l1
        
        seu_hao_3p$gold_standart <- ifelse(seu_hao_3p$gold_standart=="NK","NK cells",seu_hao_3p$gold_standart)
        qsave(seu_hao_3p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_3p.qs")
        
      }
    }
    
    #5P
    {
      #Download and untar
      {
        
        cm <- Read10X("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSM5008740_RNA_5P/")
        qsave(cm,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/cm_GSM5008740_RNA_5P.qs")
        
        cm_adt <- Read10X("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSM5008741_ADT_5P/")
        qsave(cm_adt,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/cm_adt_GSM5008741_ADT_5P.qs")
        
        cm_hto <- Read10X("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/GSM5008742_HTO_5P/")
        qsave(cm_hto,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/cm_hto_GSM5008742_HTO_5P.qs")
        
        seu_hao_5p <- CreateSeuratObject(counts=cm, assay = "RNA")
        seu_hao_5p[["ADT"]] <- CreateAssayObject(counts = cm_adt)
        seu_hao_5p[["HTO"]] <- CreateAssayObject(counts = cm_hto)
        
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_5p.qs")
        
        
        download_file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F5P%2Ecsv%2Egz",
                      "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/meta_data_5P.csv.gz")
        

        setwd("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/")
        meta_data <- read.csv(gzfile("meta_data_5P.csv.gz"))
        
        seu_hao_5p <- AddMetaData(object = seu_hao_5p, metadata = meta_data)
        
        seu_hao_5p$X <- NULL
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_5p.qs")
        
      }
      #Processing the seurat object:
      {
        DefaultAssay(seu_hao_5p) <- 'RNA'
        seu_hao_5p <- PercentageFeatureSet(seu_hao_5p, pattern = "^MT-", col.name = "percent_mito")
        feats <- c("percent_mito")
        VlnPlot(seu_hao_5p, group.by = "orig.ident", features = feats, ncol = 1) +
          NoLegend()
        
        seu_hao_5p <- subset(seu_hao_5p, subset = percent_mito < 20)
        
        
        
        options(future.globals.maxSize = 8 * 1024^3)  
        
        
        seu_hao_5p <- SCTransform(seu_hao_5p, vars.to.regress = "percent_mito", verbose = TRUE)
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_5p.qs")
        
        seu_hao_5p <- RunPCA(seu_hao_5p, npcs = 50, verbose = T)
        VizDimLoadings(seu_hao_5p, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
        
        ElbowPlot(seu_hao_5p, reduction = "pca", ndims = 50)
        
        DefaultAssay(seu_hao_5p) <- 'ADT'
        # we will use all ADT features for dimensional reduction
        # we set a dimensional reduction name to avoid overwriting the 
        VariableFeatures(seu_hao_5p) <- rownames(seu_hao_5p[["ADT"]])
        seu_hao_5p <- NormalizeData(seu_hao_5p, normalization.method = 'CLR', margin = 2) %>% 
          ScaleData() %>% RunPCA(reduction.name = 'apca')
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_5p.qs")
        ElbowPlot(seu_hao_5p, reduction = "apca", ndims = 50)
        
        
        seu_hao_5p <- FindMultiModalNeighbors(
          seu_hao_5p, reduction.list = list("pca", "apca"), 
          dims.list = list(1:25, 1:25), modality.weight.name = "RNA.weight"
        )
        
        qsave(seu_hao_5p,"/bigdata/bgutman/Collabs/Dror mevorach/data/healthy_nks/cite-seq_benchmark/Hao/seu_hao_5p.qs")
        
        seu_hao_5p <- RunUMAP(seu_hao_5p, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
        seu_hao_5p <- FindClusters(seu_hao_5p, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
        
        plot_p5 <- DimPlot(seu_hao_5p, reduction = 'wnn.umap', group.by = 'celltype.l1', label = TRUE, repel = TRUE, label.size = 2.5)
        
        print(plot_p5)
        
        seu_hao_5p$gold_standart <- seu_hao_5p$celltype.l1
        seu_hao_5p$gold_standart <- ifelse(seu_hao_5p$gold_standart=="NK","NK cells",seu_hao_5p$gold_standart)
      }
      
    }
 
    
    
  }
  
 
}













