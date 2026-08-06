
#libraries
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
    library(future.apply)
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
    library(scGate)
    library(scater)
    library(biomaRt)
    library(readxl)
    library(writexl)
    library(jpeg)
    library(cowplot)
    library(gridExtra)
    library(ggtext)
    library(tools)
    library(ReactomePA)
    # library(msigdb)
    library(GSEABase)
    library(stringr)
    library(AnnotationDbi)
    
    # library(edgeR)
    # library(zinbwave)
    # library(limma)
    # library(MAST)
    
  })
  
  
}



# Isolating Cell types-----------------------------------------------------------

get_all_cell_types <- function(seu, n_cores = 8) {
  ## ------------------------------------------------------------------ ##
  ## 0  Helpers                                                          ##
  ## ------------------------------------------------------------------ ##
  agree_2of3 <- function(seu,a, b, c, target,new_col) {
    seu@meta.data[[new_col]] <- (( seu@meta.data[[a]] == target) + ( seu@meta.data[[b]] == target) + ( seu@meta.data[[c]] == target))
    seu@meta.data[[new_col]] <- ifelse( seu@meta.data[[new_col]] >= 2, target, "other")
    seu
  }
  
  ## ------------------------------------------------------------------ ##
  ## 1  Pre-processing                       ##
  ## ------------------------------------------------------------------ ##
  plan(multicore, workers = n_cores)
  seu <- FindClusters(seu, resolution = 80)
  
  #––– cluster-level SingleR (Blueprint + Monaco) ––––––––––––––––––––––
  for (ref_fun in list(celldex::BlueprintEncodeData,
                       celldex::MonacoImmuneData)) {
    
    ref      <- ref_fun()
    ref_name <- if (identical(ref_fun, celldex::BlueprintEncodeData))
      "blueprint" else "monaco"
    
    singler  <- SingleR(
      test     = GetAssayData(seu, assay = "SCT"),
      ref      = ref,
      labels   = ref@colData@listData$label.main,
      clusters = seu@meta.data$SCT_snn_res.80
    )
    
    new_col  <- paste0("SingleR.labels.", ref_name, ".cluster80")
    x <- data.frame(row.names = rownames(singler), labels = singler$labels)
    seu <- AddMetaData(seu, metadata = x[seu$SCT_snn_res.80, 1], col.name = new_col)
    seu@meta.data[[new_col]] <- gsub("-", " ", seu@meta.data[[new_col]], fixed = TRUE)
    
    assign(paste0(ref_name, "_levels"),
           sort(unique(gsub("-", " ", ref$label.main))),
           envir = environment())
  }
  
  ## ------------------------------------------------------------------ ##
  ## 2  UCell scores                 ##
  ## ------------------------------------------------------------------ ##
  # 2.1 NK
  nk_sig <- c("CD160+","CD244+","CHST12+","CST7+","GNLY+","IL18RAP+",
              "IL2RB+","KLRC1+","KLRC3+","KLRD1+","KLRF1+","PRF1+","XCL2+")
  seu <- AddModuleScore_UCell(seu, list(nk_sig = nk_sig), name = NULL)
  seu$nk_by_ucell <- ifelse(seu$nk_sig >= .4, "NK cells", "other")
  
  # 2.2 Monocytes
  mono_sig <- c("CD14+","CD163+","FPR1+","S100A9+","S100A8+","CLEC9A-")
  seu <- AddModuleScore_UCell(seu, list(mono_sig = mono_sig), name = NULL)
  seu$monocyte_by_ucell <- ifelse(seu$mono_sig >= .4, "Monocytes", "other")
  
  # 2.3 B cells
  b_sig <- c("CD19+","MS4A1+","CD79A+","CD79B+","CD22+","CD24+")
  seu <- AddModuleScore_UCell(seu, list(b_sig = b_sig), name = NULL)
  seu$b_cells_by_ucell <- ifelse(seu$b_sig >= .4, "B cells", "other")
  
  ##  2.4 CD4 / CD8 signatures ----------
  cd4_sig  <- c("CD4+","IL7R+","TRAC+","CD3D+","CD3E+","CD2+","CXCR4+","LTB+")
  cd8_bur  <- c("CD8A+","CD8B+","CD8B2+")
  seu <- AddModuleScore_UCell(seu, list(cd4_sig = cd4_sig, cd8_bur = cd8_bur), name = NULL)
  seu$CD4_score_primary <- seu$cd4_sig - seu$cd8_bur
  seu$cd_4_signature_by_ucell <- ifelse(seu$CD4_score_primary >= .3, "CD4+ T cells","other")
  
  cd8_sig  <- c("CD8A+","CD8B+","CD8B2+","TRAC+","LCK+","RUNX3+","CD3E+","EOMES+")
  cd4_bur  <- c("CD4+","CCR7+","IL7R")
  nk_bur   <- c("NCR1","FCGR3A","TYROBP")
  seu <- AddModuleScore_UCell(seu, list(cd8_sig = cd8_sig, cd4_bur = cd4_bur, nk_bur = nk_bur), name = NULL)
  seu$CD8_to_CD4_score_primary <- seu$cd8_sig - seu$cd4_bur
  seu$CD8_to_NK_score_primary  <- seu$cd8_sig - seu$nk_bur
  seu$cd_8_signature_by_ucell  <- ifelse(
    seu$CD8_to_CD4_score_primary >= .3 | seu$CD8_to_NK_score_primary >= .3,
    "CD8+ T cells","other")
  
  ## ------------------------------------------------------------------ ##
  ## 3  Per-cell SingleR (Blueprint + Monaco)                            ##
  ## ------------------------------------------------------------------ ##
  for (ref_fun in list(celldex::BlueprintEncodeData,
                       celldex::MonacoImmuneData)) {
    ref      <- ref_fun()
    ref_name <- if (identical(ref_fun, celldex::BlueprintEncodeData)) "blueprint" else "monaco"
    singler  <- SingleR(test = GetAssayData(seu, assay = "SCT"),
                        ref = ref,
                        labels = ref@colData@listData$label.main)
    new_col  <- paste0("SingleR.labels.", ref_name)
    seu      <- AddMetaData(seu, metadata = setNames(singler$labels, rownames(singler)), col.name = new_col)
    seu@meta.data[[new_col]] <- gsub("-", " ", seu@meta.data[[new_col]], fixed = TRUE)
  }
  
  ## ------------------------------------------------------------------ ##
  ## 4  All flags → final_cell_type                                      ##
  ## ------------------------------------------------------------------ ##
  seu <- agree_2of3(seu, "nk_by_ucell",
                    "SingleR.labels.blueprint.cluster80",
                    "SingleR.labels.monaco.cluster80",
                    "NK cells", "nk_flag")
  
  seu <- agree_2of3(seu, "monocyte_by_ucell",
                    "SingleR.labels.blueprint.cluster80",
                    "SingleR.labels.monaco.cluster80",
                    "Monocytes", "mono_flag")
  
  seu <- agree_2of3(seu, "b_cells_by_ucell",
                    "SingleR.labels.blueprint.cluster80",
                    "SingleR.labels.monaco.cluster80",
                    "B cells", "b_flag")
  
  # CD4: pure majority vote (no overrides)
  seu <- agree_2of3(seu, "cd_4_signature_by_ucell",
                    "SingleR.labels.blueprint",
                    "SingleR.labels.monaco",
                    "CD4+ T cells", "cd4_flag")
  
  # CD8: pure majority vote (no overrides; use original SingleR labels)
  seu <- agree_2of3(seu, "cd_8_signature_by_ucell",
                    "SingleR.labels.blueprint",
                    "SingleR.labels.monaco",
                    "CD8+ T cells", "cd8_flag")
  
  ## 4.5 final priority order (NK > Mono > B > [CD4&CD8 both => Other] > CD4 > CD8 > Other)
  seu$final_cell_type = case_when(
    seu$nk_flag == "NK cells" ~ "NK cells",
    seu$mono_flag == "Monocytes" ~ "Monocytes",
    seu$b_flag == "B cells"   ~ "B cells",
    seu$cd4_flag == "CD4+ T cells" & seu$cd8_flag == "CD8+ T cells" ~ "Other",
    seu$cd4_flag == "CD4+ T cells" ~ "CD4+ T cells",
    seu$cd8_flag == "CD8+ T cells"  ~ "CD8+ T cells",
    TRUE      ~ "Other"
  )
  
  return(seu)
}


# get_healthy_list-----------------------------------------------------------

get_healthy_list <- function(seurat_list) {
  ## 1. sanity checks ----------------------------------------------------------
  required <- c("seu_cd55", "seu_cd59", "seu_hong", "seu_kang", "seu_ramesh")
  
  if (!all(required %in% names(seurat_list))) {
    missing <- setdiff(required, names(seurat_list))
    stop("get_healthy_list(): missing object(s): ", paste(missing, collapse = ", "))
  }
  if (length(seurat_list) != length(required)) {
    extra <- setdiff(names(seurat_list), required)
    stop("get_healthy_list(): list must contain **only** ",
         paste(required, collapse = ", "),
         if (length(extra)) paste0("; extra element(s): ", paste(extra, collapse = ", ")))
  }
  
  ## 2. pull objects -----------------------------------------------------------
  seu_cd55   <- seurat_list[["seu_cd55"]]
  seu_cd59   <- seurat_list[["seu_cd59"]]
  seu_hong   <- seurat_list[["seu_hong"]]
  seu_kang   <- seurat_list[["seu_kang"]]
  seu_ramesh <- seurat_list[["seu_ramesh"]]
  
  ## 3. subset helper ----------------------------------------------------------
  subset_one <- function(seu, smp) subset(seu, subset = sample == smp)
  
  ## 4. create all subsets -----------------------------------------------------
  hong_samples   <- paste0("HC-", 1:5)              # HC-1 … HC-5
  kang_samples   <- c("JSB", "KYO")
  ramesh_samples <- c("2", "3", "32")
  
  hong_sub   <- lapply(hong_samples,   subset_one, seu = seu_hong)
  kang_sub   <- lapply(kang_samples,   subset_one, seu = seu_kang)
  ramesh_sub <- lapply(ramesh_samples, subset_one, seu = seu_ramesh)
  cd55_sub   <- list(subset_one(seu_cd55, "CD55-control"))
  cd59_sub   <- list(subset_one(seu_cd59, "CD59-control"))
  
  ## 5. name each subset -------------------------------------------------------
  names(hong_sub)   <- sprintf("seu_hong_healthy_%s", hong_samples)
  names(kang_sub)   <- sprintf("seu_kang_healthy_%s", kang_samples)
  names(ramesh_sub) <- sprintf("seu_ramesh_healthy_%s", ramesh_samples)
  names(cd55_sub)   <- "seu_mevorach_healthy_CD55_Healthy"
  names(cd59_sub)   <- "seu_mevorach_healthy_CD59_Healthy"
  
  ## 6. bundle into the final healthy_list ------------------------------------
  make_entry <- function(seu, dataset, sample_label) {
    list(seu       = seu,
         condition = "Healthy",
         dataset   = dataset,
         sample    = sample_label)
  }
  
  healthy_list <- c(
    ## Hong
    mapply(make_entry, hong_sub, "Hong",   hong_samples,   SIMPLIFY = FALSE),
    ## Kang
    mapply(make_entry, kang_sub, "Kang",   kang_samples,   SIMPLIFY = FALSE),
    ## Ramesh
    mapply(make_entry, ramesh_sub, "Ramesh", ramesh_samples, SIMPLIFY = FALSE),
    ## Mevorach CD55/CD59
    list(
      seu_mevorach_healthy_CD55_Healthy = make_entry(cd55_sub[[1]], "Mevorach", "CD55-control"),
      seu_mevorach_healthy_CD59_Healthy = make_entry(cd59_sub[[1]], "Mevorach", "CD59-control")
    )
  )
  
  ## keep list element order identical to original script
  healthy_list <- healthy_list[
    c(names(hong_sub),
      names(kang_sub),
      names(ramesh_sub),
      names(cd55_sub),
      names(cd59_sub))
  ]
  
  return(healthy_list)
}



# get_confusion_matrix and identifaction plots--------------------------------------------------------------

get_confusion_matrix <- function(seu, 
                                 metric_column, 
                                 ground_truth_column = "celltype.l1",
                                 ground_truth_label) {
  # Extract logical vectors: predicted positive? and actual positive?
  
  preds  <- seu@meta.data[[metric_column]] != "other"
  truths <- seu@meta.data[[ground_truth_column]] == ground_truth_label
  
  # Compute confusion‐matrix entries
  TP <- sum( preds  &  truths,  na.rm = TRUE)
  FP <- sum( preds  & !truths,  na.rm = TRUE)
  TN <- sum(!preds  & !truths, na.rm = TRUE)
  FN <- sum(!preds  &  truths,  na.rm = TRUE)
  
  # Build a named 2×2 matrix
  conf_mat <- matrix(
    c(TP, FP,
      FN, TN),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      Prediction = c("Positive", "Negative"),
      Reference  = c("Positive", "Negative")
    )
  )
  
  # Calculate metrics (guarding against zero‐division)
  total      <- TP + TN + FP + FN
  accuracy   <- (TP + TN) / total
  precision  <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall     <- if ((TP + FN) > 0) TP / (TP + FN) else NA      # also called Sensitivity
  specificity<- if ((TN + FP) > 0) TN / (TN + FP) else NA
  f1_score   <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA
  
  list(
    confusion_matrix = conf_mat,
    accuracy         = accuracy,
    precision        = precision,
    recall           = recall,
    specificity      = specificity,
    f1_score         = f1_score
  )
}


get_cell_type_idenfication_4_umaps <- function(seu,metric_column, metric_label, ucell_column)
{
  seu$SingleR.labels.blueprint.cluster80 <- ifelse(seu$SingleR.labels.blueprint.cluster80==metric_label,metric_label,"other")
  
  seu$SingleR.labels.monaco.cluster80 <- ifelse(seu$SingleR.labels.monaco.cluster80 ==metric_label,metric_label,"other")
  
  #Checking reducation type
  {
    # Get reduction names in the object
    available_reds <- names(seu@reductions)          # character vector
    
    # Check which of the two we have
    has_umap      <- "umap"      %in% available_reds
    has_wnn_umap  <- "wnn.umap"  %in% available_reds
    
    # Decide, or error out
    if (has_umap && !has_wnn_umap) {
      reduction_type <- "umap"
    } else if (!has_umap && has_wnn_umap) {
      reduction_type <- "wnn.umap"
    } else {
      stop(
        "reduction_type error: the Seurat object must contain **exactly one** of ",
        "'umap' or 'wnn.umap'.\nFound: ",
        paste(intersect(available_reds, c('umap', 'wnn.umap')), collapse = ', ')
      )
    }
    
    colours <- setNames(
      c("#1f78b4", "#e31a1c"),          # the colour values
      c("other",  metric_label)         # the names – note metric_label is *not* quoted
    )
    
    
  }
 
  
  
  # Modify titles with specific words in bold
  
  p_blueprint_title <- paste(metric_label,"Identification by <b>SingleR</b> with <b>Blueprint</b> as Reference")
  
  p_blueprint <- dittoDimPlot(seu, "SingleR.labels.blueprint.cluster80", reduction.use = reduction_type, size = 0.4) +
    ggplot2::scale_color_manual(values = colours, drop = FALSE) +
    ggplot2::ggtitle(label = p_blueprint_title) +
    theme(plot.title = element_markdown())
  
  p_monaco_title <- paste(metric_label,"Identification by <b>SingleR</b> with <b>MonacoImmuneData</b> as Reference")
  
  p_monaco <- dittoDimPlot(seu, "SingleR.labels.monaco.cluster80", reduction.use = reduction_type, size = 0.4) +
    ggplot2::scale_color_manual(values = colours, drop = FALSE) +
    ggplot2::ggtitle(label = p_monaco_title) +
    theme(plot.title = element_markdown())
  
  p_signature_title <- paste(metric_label,"Identification with <b>Markers Signature</b> with a threshold of 0.4")
  
  
  p_signature <- dittoDimPlot(seu, ucell_column, reduction.use = reduction_type, size = 0.4) +
    ggplot2::scale_color_manual(values = colours, drop = FALSE) +
    ggplot2::ggtitle(label = p_signature_title) +
    theme(plot.title = element_markdown())
  
  p_signature_title <- paste("<b>Final</b>", metric_label ,"Identification")
  
  
  p_final_identification <- dittoDimPlot(seu, metric_column, reduction.use = reduction_type, size = 0.4) +
    ggplot2::scale_color_manual(values = colours, drop = FALSE) +
    ggplot2::ggtitle(label = p_signature_title) +
    theme(plot.title = element_markdown())
  
  
  
  
  
  # Add letters to each UMAP plot on the leftmost side
  p_blueprint <- p_blueprint + 
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
  
  p_monaco <- p_monaco + 
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
  
  p_signature <- p_signature + 
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
  
  p_final_identification <- p_final_identification + 
    ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
  
  # Combine plots into one
  p_cell_type_idenfication <- ((p_blueprint + p_monaco) / (p_signature + p_final_identification))
  
}

get_confusion_matrix_and_idenfication_4_umaps <- function(seu,metric_column, metric_label,                                  
                                                          ground_truth_label,  ground_truth_column = "celltype.l1", ucell_column){
  
  confusion_list<- get_confusion_matrix(seu, 
                          metric_column, 
                          ground_truth_column,
                          ground_truth_label)
  
  p_cell_type_idenfication <- get_cell_type_idenfication_4_umaps(seu,metric_column, metric_label, ucell_column)
  
  list(confusion_list = confusion_list,
       p_cell_type_idenfication =p_cell_type_idenfication)
  
  
}
# get_full_seurat--------------------------------------------------------------


get_full_seurat <- function(seu, nbin_for_cell_cycle = 23, compute_cell_cycle = TRUE) {
  seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
  seu <- SCTransform(seu, vars.to.regress = "percent_mito", verbose = FALSE)
  
  # QC subset
  seu <- subset(seu, subset = nFeature_SCT > 500 & percent_mito < 20)
  DefaultAssay(seu) <- "SCT"
  
  if (compute_cell_cycle) {
    # cell-cycle on SCT
    seu <- CellCycleScoring(
      seu,
      s.features   = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      nbin         = nbin_for_cell_cycle,
      assay        = "SCT"
    )
  }
  
  # PCA dims safe-guard
  max_pcs  <- max(2, min(50, ncol(seu) - 1, nrow(seu[["SCT"]]) - 1))
  dims_use <- 1:min(25, max_pcs)
  
  # compute a valid neighbour count
  n_cells <- ncol(seu)
  nn_use  <- max(2, min(30, n_cells - 1))
  
  # downstream
  seu <- RunPCA(seu, npcs = max_pcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = dims_use, k.param = nn_use, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
  seu <- RunUMAP(
    seu,
    reduction     = "pca",
    dims          = dims_use,
    n.components  = 2,
    n.neighbors   = nn_use,
    n.epochs      = 200,
    min.dist      = 0.3,
    learning.rate = 1,
    spread        = 1,
    verbose       = FALSE
  )
  
  seu
}

#  one_on_one and it's helper function----------------------------------------------------------

one_on_one <- function(seu_patient,seu_control,patient_type="patient",control_type="Healthy",
                       patient_orig_ident ="Mevorach") {
  # browser()
  #intial_check
  {
    
    #Checking that the method_column indeed exsists in each of the seurat objects
    if (!("condition" %in% unique(colnames(seu_patient@meta.data)))) {
      print(paste("The condition column does not exist in seu_patient"))
      break }
    else if (!("condition" %in% unique(colnames(seu_control@meta.data)))) {
      print(paste("The condition column does not exist in seu_control"))
      break }
    
    else if (!(control_type %in% unique(seu_control$condition))) {
      print(paste("The condition", control_type,"does not exist in the seu_control"))
      break }
    
    else if (!(patient_type %in% unique(seu_patient$condition))) {
      print(paste("The condition", patient_type,"does not exist in the seu_patient"))
      break }
    else if (!("sample" %in% unique(colnames(seu_control@meta.data)))) {
      print("The sample column does not exist in  seu_control")
      break }  
    else if (!("sample" %in% unique(colnames(seu_patient@meta.data)))) {
      print("The sample column does not exist in  seu_patient")
      break }  
    
  }  
  
  
  
  sample_output_list <- list()
  
  seu_control <- subset(seu_control,subset= condition == control_type)
  seu_patient <- subset(seu_patient,subset= condition == patient_type)
  
  #For each "healthy" sample in seu_control 
  for (control_sample in unique(seu_control$sample)) {
    
    #Subset to the that one healthy sample
    seu_control_sample <- subset(seu_control,subset= sample == control_sample)
    
    de_markers_list <- list()
    
    #Run downsample + SCT + FindMarkers 3 times
    #(Because in each downsample iteration, the count matrix changes up a little bit)
    for (i in 1:3) {
      
      #call the downsample function
      seurats_post_downsample <-  tryCatch({
        downsample_by_porp (seu_control_sample, seu_patient)
      }, error = function(e) {
        stop("The function one_on_one failed because the function downsample_by_porp encountered an error: ", e$message)
      })
      
      
      
      
      # creating and manipulating seu_merge
      seu_merge <- merge(seurats_post_downsample$seu_control_sample
                         ,seurats_post_downsample$seu_patient)
      
      
      seu_merge@active.assay <- "RNA"
      seu_merge[["SCT"]] <- NULL
      
      seu_merge[["RNA"]] <- JoinLayers(seu_merge[["RNA"]])
      
      seu_merge <- SCTransform(seu_merge, vars.to.regress = "percent_mito", verbose = TRUE, conserve.memory = TRUE)
      
      Idents(seu_merge) <- "condition"
      
      #maybe we should remove the "matching_item"
      patient_group <-   matching_item <- unique(seu_merge@active.ident)[grep(patient_orig_ident, unique(seu_merge@active.ident))]
      
      de_markers <- FindMarkers(seu_merge, ident.1 = patient_type, group.by = 'condition',logfc.threshold = 0, min.pct = 0,min.cells.feature = 0)
      head(de_markers)
      
      # Save the result of this iteration into the list
      de_markers_list[[i]] <- de_markers
      gc()
      
    }
    
    
    # Find the intersection of genes across all 3 tables
    {
      reoccuring_genes <- Reduce(intersect, lapply(de_markers_list, rownames))
      
      # Filter each table to retain only the common genes
      filtered_tables <- lapply(de_markers_list, function(tbl) tbl[reoccuring_genes, ])
      
      # Combine the filtered tables into a single array for averaging
      # This assumes all tables are filtered to have the same genes and structure
      combined_data <- abind::abind(filtered_tables, along = 3)
      
      # Compute the average across the 3rd dimension (iterations) for each gene and each column
      averaged_table <- apply(combined_data, c(1, 2), mean)
      
      # Convert the resulting matrix into a data frame for easier handling
      final_table <- as.data.frame(averaged_table)
      
      # Set row names to the gene names
      rownames(final_table) <- reoccuring_genes
      
    }
    
    
    #Filter out housekeeping, ribosomal, mitochondrial and hemoglobin genes
    
    original_table <- final_table
    final_table <- filter_housekeeping_genes(original_table)
    #get the log2 expression scatter plot
    
    plot_and_data <- tryCatch({
      get_p_avg_expr_log2 (seu_merge , patient_type,control_type)
    }, error = function(e) {
      stop("The function one_on_one failed because the function get_p_avg_expr_log2 encountered an error: ", e$message)
    })
    
    #violin plot
    {
      plot <- VlnPlot(seu_merge, features = c("nFeature_RNA", "nCount_RNA", "percent_mito","nCount_SCT","nFeature_SCT"), group.by = "condition" ,ncol = 5)
      
      violin_title <- sub(" Patient Averaged Expression Levels$", "", plot_and_data$plot_title)
      
      # Add centered overarching title and print the plot
      print(
        plot + plot_annotation(
          title = violin_title,
          theme = theme(plot.title = element_text(hjust = 0.5))  # hjust = 0.5 centers the title
        )
      )
      
    }
    

    # saving plot_and_data in sample_output_list  
    sample_output_list[[control_sample]] <- c(
      plot_and_data, 
      list(final_table = final_table, original_table = original_table)
    )
    
    
    
    
  }
  
  
  seurat_output_list <- list()
  
  seurat_output_list[[as.character(unique(seu_control$orig.ident)[1])]] <- sample_output_list
  
  
  
  
  return(seurat_output_list) 
  
}

filter_housekeeping_genes <- function(original_table) {
  
  
  #Filter out housekeeping, ribosomal, mitochondrial and hemoglobin genes
  
  #Read housekeeping genes names
  
  joshi_hk_genes <- read_excel("/bigdata/bgutman/Collabs/Dror mevorach/data/joshi_low_gini_coefficients_genes.xlsx", col_names = TRUE)
  
  levanon_hk_genes <- read_excel("/bigdata/bgutman/Collabs/Dror mevorach/data/levanon_housekeeping_genes.xlsx", col_names = TRUE)
  
  final_table <- original_table 
  # Add a temporary column with rownames
  final_table$temp_rownames <- rownames(final_table)
  
  # Combine all housekeeping genes from both files into a single vector
  all_hk_genes <- c(joshi_hk_genes[[2]], levanon_hk_genes[[1]])
  
  # Store the initial number of rows in final_table
  initial_row_count <- nrow(final_table)
  
  # Filter out rows whose names are in the housekeeping genes
  final_table <- final_table[!final_table$temp_rownames %in% all_hk_genes, ]
  
  
  
  # Calculate and print the number of rows removed
  rows_removed <- initial_row_count - nrow(final_table)
  print(paste("Number of housekeeping genes removed:", rows_removed))
  
  initial_row_count <- nrow(final_table)
  
  
  final_table <- final_table[!grepl("^MT-|^RP[SL]|^HB[^(P|E|S)]|^AC[0-9]+\\.[0-9]+$|^AL[0-9]+\\.[0-9]+$", final_table$temp_rownames), ]
  
  
  final_table$temp_rownames <- NULL
  
  rows_removed <- initial_row_count - nrow(final_table)
  print(paste("Number of additonal ribosomal, mitochondrial and hemoglobin genes and some lncRNAs removed:", rows_removed))
  
  return (final_table)  
}

get_p_avg_expr_log2 <- function(seu, patient_condition ="Mevorach", control_condition) {
  
  # browser()
  
  Idents(seu) <- "condition"
  DefaultAssay(seu) <- "SCT"
  
  avg_exp <- AverageExpression(seu, assays = "SCT", group.by = "condition")
  avg_exp <- as.data.frame(avg_exp)
  
  # Making sure the columns are numeric and making it ready for log
  # Making sure the columns are numeric and making it ready for log
  
  column_names <- colnames(avg_exp)
  
  #Creating new, eazy to work with columns
  {
    {if (!(length(column_names) == 2)) {
      stop("Error in the function 'get_p_avg_expr_log2': The number of unique original idents in the merged Seurat object isn't 2 as expected. Please change it accordingly")
    } else if (any(grepl(patient_condition, column_names[1]))) {
      avg_exp$SCT.patient <- avg_exp[[column_names[1]]]
      avg_exp$SCT.control <- avg_exp[[column_names[2]]]
    } else if (any(grepl(patient_condition, column_names[2]))) {
      avg_exp$SCT.patient <- avg_exp[[column_names[2]]]
      avg_exp$SCT.control <- avg_exp[[column_names[1]]]
    }
      else {
        # If patient identifier is not found in either column name
        stop("Error in the function 'get_p_avg_expr_log2': patient_condition was not found in any of the original idents. Please make sure to call the function with the correct patient_condition")
      }
    }
    
    #Remove the original columns
    avg_exp[[column_names[1]]] <- NULL
    avg_exp[[column_names[2]]] <- NULL
    
  }
  
  
  avg_exp$SCT.control <- as.numeric(as.character(avg_exp$SCT.control))
  avg_exp$SCT.patient <- as.numeric(as.character(avg_exp$SCT.patient))
  
  avg_exp$SCT.control <- avg_exp$SCT.control + 1
  avg_exp$SCT.patient <- avg_exp$SCT.patient + 1
  
  line_data <- data.frame(x = range(avg_exp$SCT.control, avg_exp$SCT.patient, na.rm = TRUE))
  line_data$y <- line_data$x  # displaying the line y = x
  avg_exp <- avg_exp %>% tibble::rownames_to_column("gene")
  
  gene_to_be_displayed <- avg_exp %>%
    filter(SCT.patient >= 1.5,
           SCT.patient >= SCT.control,
           SCT.patient > 2,
           SCT.control > 2) %>% pull(gene)
  
  
  weird_genes <- avg_exp %>%
    filter(SCT.patient >= SCT.control,
           SCT.control <1.1,
           SCT.patient > 1.35
    ) %>% pull(gene)  
  
  
  
  
  housekeeping_genes <- c("B2M", "H3F3A", "MALAT1", "YBX1", "UBA52", 
                          "PSMB9", "ATP5E", "ATP5L", "UQCR11", "TOMM7")
  
  # Filter genes for labeling
  
  label_genes <- avg_exp %>%
    filter(gene %in% gene_to_be_displayed,
           !(gene %in% housekeeping_genes),
           !grepl("^MT", gene),
           !grepl("^RP", gene),
           (SCT.patient - SCT.control) >= 0.5)
  
  
  
  
  
  mt_genes <- (length(grep("^MT", gene_to_be_displayed)) / length(gene_to_be_displayed)) * 100
  ribo_genes <- (length(grep("^RP", gene_to_be_displayed)) / length(gene_to_be_displayed)) * 100
  
  #setting the plot title  
  {
    
    #The control sample is the without the "-pt"|" or the "-ctrl", as these are the suffix of all the sample names in seu_cd59 and seu-cd55
    
    
    # browser()
    
    control_sample_name <-    
      unique(seu@meta.data$sample[seu@meta.data$condition == control_condition])       
    
    patient_sample_name <-    
      unique(seu@meta.data$sample[seu@meta.data$condition == patient_condition]) 
    
    orig_ident_control_sample <- unique(seu@meta.data$orig.ident[seu@meta.data$sample == control_sample_name])
    
    orig_ident_patient_sample <- unique(seu@meta.data$orig.ident[seu@meta.data$condition == patient_condition])
    
    
    
    # disease_condition_name <- 
    # unique(seu@meta.data$condition)[unique(seu@meta.data$condition) != control_condition]
    
    plot_title <- paste0(control_condition, " Control (",control_sample_name, " from ",orig_ident_control_sample,")"," vs ",patient_condition," Patient (", patient_sample_name, " from ", orig_ident_patient_sample,") Averaged Expression Levels")    
    
  }
  
  
  p_avg_mon <- ggplot(avg_exp, aes(x = SCT.control, y = SCT.patient)) +
    geom_point(aes(color = ifelse(gene %in% gene_to_be_displayed, "significant", "non-significant"))) +
    geom_line(data = line_data, aes(x = x, y = y), color = "blue") +
    scale_x_continuous(trans = scales::log2_trans()) +  # Log2 scale for x-axis
    scale_y_continuous(trans = scales::log2_trans()) +  # Log2 scale for y-axis
    scale_color_manual(values = c("significant" = "red", "non-significant" = "black")) +
    labs(x = "Log2(Control)", y = "Log2(Patient)",
         title = plot_title) +
    theme_minimal() +
    guides(color = guide_legend(title = NULL)) +  # Remove legend title
    annotate("text", x = 0, y = Inf, label = sprintf("Percent of ribo_genes: %.2f%%", ribo_genes), 
             hjust = -0.1, vjust = 2.1, size = 5, color = "blue") +
    annotate("text", x = 0, y = Inf, label = sprintf("Percent of mt_genes: %.2f%%", mt_genes), 
             hjust = -0.1, vjust = 4.1, size = 5, color = "blue") +
    geom_text_repel(data = label_genes, aes(x = SCT.control, y = SCT.patient, label = gene),
                    size = 3, max.overlaps = Inf)  # Adding gene labels without overlap
  
  
  print(p_avg_mon)
  
  output_list <- list()
  output_list$plot <- p_avg_mon
  output_list$de_genes <- gene_to_be_displayed
  output_list$artificat <- weird_genes
  output_list$plot_title <- plot_title
  
  return(output_list)
}

downsample_by_porp <- function(seu_control_sample, seu_patient) {
  # browser()
  
  #downsample
  {
    
    #find the intersect genes and keep only them in the cms
    
    intersect_genes = intersect(rownames(seu_control_sample@assays$RNA$counts),rownames(seu_patient@assays$RNA$counts))
    
    control_downsampled_cm = seu_control_sample@assays$RNA$counts[intersect_genes,]
    patient_downsampled_cm = seu_patient@assays$RNA$counts[intersect_genes,]
    
    #downsample
    downsampled_cms <- scuttle::downsampleBatches(control_downsampled_cm,patient_downsampled_cm)
    
    
    seu_control_sample <- CreateSeuratObject(counts = downsampled_cms[[1]],meta.data = seu_control_sample@meta.data, min.cells =3,  min.features=200)
    
    seu_patient <- CreateSeuratObject(counts = downsampled_cms[[2]],meta.data = seu_patient@meta.data
                                      , min.cells =3,  min.features=200)
    
    
    
    seurats_post_downsample <- list("seu_control_sample"= seu_control_sample,
                                    "seu_patient" = seu_patient)
    
    
  }
  
  
  
  
  #sanity check density plot
  {
    
    #setting the plot title  
    
    plot_title <- paste("Density Plot of Number of Counts Between",
                        unique(seu_control_sample$sample)[[1]], "from",unique(seu_control_sample$orig.ident)[[1]] , "(Control) and", unique(seu_patient$orig.ident)[[1]], "(Patient)")
    
    
    df =data.frame(number_of_counts=c(colSums(downsampled_cms[[1]]),colSums(downsampled_cms[[2]])),group=c(rep('control',ncol(downsampled_cms[[1]])),rep('patient',ncol(downsampled_cms[[2]]))))
    
    
    # Create the density plot
    plot <- ggplot(data = df, aes(x = number_of_counts, group = group, fill = group)) +
      geom_density(adjust = 1.5, alpha = 0.6) +  
      scale_fill_manual(values = c("steelblue", "coral")) +  #
      labs(
        title = plot_title,
        x = "Number of Counts",
        y = "Density",
        fill = "Group"
      ) +  
      theme_minimal(base_size = 14) + 
      theme(
        legend.position = "top",      
        legend.title = element_text(face = "plain"),  
        legend.text = element_text(size = 12),      
        plot.title = element_text(hjust = 0.5, face = "bold"),  
        axis.title = element_text(face = "bold")    
      )
    
    print(plot)
    
    seurats_post_downsample <- list("seu_control_sample"= seu_control_sample,
                                    "seu_patient" = seu_patient)
    return(seurats_post_downsample)
  }
  
  
  
}

# meta table-------------------------------------------------------------

create_meta_table <- function(output_list) {
  # 0) Identify the patient/dataset name from the object
  #    e.g. "CD59.Def_Mevorach_CD59_Patient_outputlist"
  list_name <- deparse(substitute(output_list))
  is_mevorach_patient <- grepl("_Mevorach_", list_name, ignore.case = TRUE)
  
  # 1) Flatten output_list into a simple list of "sample" objects.
  samples_list <- unlist(output_list, recursive = FALSE)
  
  # 2) Combine all sample$original_table data into one big data frame:
  samples_names <- names(samples_list)
  combined_df <- bind_rows(lapply(seq_along(samples_list), function(i) {
    df <- samples_list[[i]]$original_table
    df$gene    <- rownames(df)
    df$dataset <- sub("\\..*$", "", samples_names[i])
    df
  }))
  
  # 3) Precompute columns to support dsi_scores:
  combined_df <- combined_df %>%
    mutate(
      is_up   = ifelse(!is.na(avg_log2FC) & avg_log2FC >= 1 &
                         !is.na(p_val_adj)  & p_val_adj <= 0.05, 1, 0),
      is_down = ifelse(!is.na(avg_log2FC) & avg_log2FC <= -1 &
                         !is.na(p_val_adj) & p_val_adj <= 0.05, 1, 0)
    )
  
  # 4) Summarize per gene:
  meta_table <- combined_df %>%
    group_by(gene) %>%
    summarize(
      comparisons_num    = n(),
      fisher_score       = sum(-2 * log(p_val), na.rm = TRUE),
      valid_comparisons  = sum(!is.na(p_val)),
      avg_log2FC = if (is_mevorach_patient) {
        mean(avg_log2FC[ !(dataset == "Mevorach" & !is.na(p_val) & p_val > 0.05) ], na.rm = TRUE)
      } else {
        mean(avg_log2FC, na.rm = TRUE)
      },
      median_log2FC = if (is_mevorach_patient) {
        median(avg_log2FC[ !(dataset == "Mevorach" & !is.na(p_val) & p_val > 0.05) ], na.rm = TRUE)
      } else {
        median(avg_log2FC, na.rm = TRUE)
      },
      num_upregulated    = sum(is_up),
      num_downregulated  = sum(is_down),
      dsi                = sum(is_up) - sum(is_down),
      pct.1              = mean(pct.1, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # 5) Compute Fisher’s combined p-values, then BH adjust:
  meta_table <- meta_table %>%
    mutate(
      p_val_fisher = if_else(
        valid_comparisons == 0,
        NA_real_,
        1 - pchisq(fisher_score, df = 2 * valid_comparisons)
      ),
      q_val_fisher = p.adjust(p_val_fisher, method = "BH")
    )
  
  # 6) Round numeric columns where needed:
  meta_table <- meta_table %>%
    mutate(
      fisher     = round(fisher_score, 2),
      avg_log2FC = round(avg_log2FC, 2)
    )
  
  # 7) Reorder columns:
  meta_table <- meta_table %>%
    dplyr::select(
      gene,
      comparisons_num,
      dsi,
      num_upregulated,
      num_downregulated,
      avg_log2FC,
      median_log2FC,
      pct.1,
      fisher,
      valid_comparisons,
      p_val_fisher,
      q_val_fisher
    )
  meta_table <- as.data.frame(meta_table)
  
  # 8) Make gene names row names:
  rownames(meta_table) <- meta_table$gene
  
  return(meta_table)
}

# integration and helper functions -------------------------------------------------------------

get_integration <- function(seu_patient, healthy_list_for_integration, patient_condition) {

  
  seu_patient <- get_slim_seu(seu_patient)
  
  #merging loop
  for (seu_healthy in healthy_list_for_integration) {
    
    seu_healthy <- get_slim_seu(seu_healthy)
    seu_patient <- merge(seu_patient,seu_healthy)
    
  }
  
  # seu_patient[["ADT"]] <- NULL
  
  #change datasets names        
  seu_patient@meta.data$orig.ident <- case_when(
    seu_patient@meta.data$sample == "EM_HD1_NK_scRNAseq" ~ "Smith_em",
    seu_patient@meta.data$sample %in% c("NRHV187_unstim", "NRHV193_unstim")  ~ "Smith_unstim",
    seu_patient@meta.data$sample == "H6" ~ "Witkowski_H6",
    seu_patient@meta.data$sample == "H8" ~ "Witkowski_H8",
    seu_patient@meta.data$sample %in% c("H21", "H18", "H2")  ~ "Witkowski_H21_H18_H2",  
    TRUE ~ seu_patient@meta.data$orig.ident
  ) 
  
  seu_patient <- SCTransform(seu_patient, vars.to.regress = c("percent_mito", "S.Score", "G2M.Score"), verbose = F,   conserve.memory = FALSE)
  gc()
  seu_patient <- RunPCA(seu_patient, npcs = 50, verbose = F)


  seu_patient <- IntegrateLayers(object = seu_patient, method = HarmonyIntegration,
                                 orig.reduction = "pca", new.reduction = 'harmony',
                                 assay = "SCT", verbose = F)

  
  seu_patient <- RunUMAP(seu_patient, reduction = "harmony", dims = 1:30)

  seu_patient <- FindNeighbors(seu_patient, dims = 1:30, reduction = "harmony")

  seu_patient <- FindClusters (seu_patient)

  
  
  integration_4_umap <- get_four_integration_plots(seu_patient)
  
  
  seu_patient <- PrepSCTFindMarkers(seu_patient)   
  
  patient_condition <- "CD59.Def"
  
  ## ----- 1. Cell types & short names --------------------------------------------
  cell_types  <- c("NK cells", "B cells", "Monocytes", "CD4+ T cells", "CD8+ T cells")
  short_names <- c("NK", "B_cells", "Monocytes", "CD4_T_cells", "CD8_T_cells")
  
  ## ----- 2. Subset + prepare SCT models in one pass -----------------------------
  seu_by_ct <- lapply(cell_types, function(ct) {
    obj <- subset(seu_patient, subset = final_cell_type == ct)
    # PrepSCTFindMarkers(obj)           # fixes the unequal‑library‑size error
  })
  names(seu_by_ct) <- short_names
  
  
  ## ----- 3. Differential expression per subset ----------------------------------
  de_tables <- lapply(seu_by_ct, function(obj) {
    FindMarkers(
      obj,
      ident.1  = patient_condition,
      group.by = "condition",
      assay    = "SCT",                # explicit SCT assay
      recorrect_umi = FALSE 
    )
  })
  
  ## ----- 4. Assemble the five return lists --------------------------------------
  return_lists <- mapply(
    function(de_tab, nm) list(
      integrated_seu     = seu_patient,      # full integrated object
      integration_4_umap = integration_4_umap,
      de_table           = de_tab
    ),
    de_tables,
    short_names,
    SIMPLIFY = FALSE
  )
  names(return_lists) <- short_names
  
  ## If this code is inside a function, make this the last line:
  return(return_lists)
  
}

get_slim_seu <- function(seu) {
  
  seu_slim <- seu
  DefaultAssay(seu_slim) <- "RNA"
  seu_slim[["SCT"]] <- NULL
  seu_slim@reductions$pca <- NULL
  seu_slim@reductions$umap <- NULL
  
  # Keep only specific metadata columns: orig.ident, nCount_RNA, nFeature_RNA, sample
  metadata_colnames <- colnames(seu_slim@meta.data)
  
  # Exclude these columns (removed the empty string "")
  excluded_names <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA",
    "sample",     "condition",   "percent_mito",
    "S.Score",    "G2M.Score",   "Phase", "final_cell_type"
  )
  # Find indices of the excluded names in the metadata columns
  indices_to_remove <- match(excluded_names, metadata_colnames)
  
  # Check for any names not found and print a message, then return if any are missing
  not_found <- excluded_names[is.na(indices_to_remove)]
  if (length(not_found) > 0) {
    for (name in not_found) {
      cat(paste0("\nThe column '", name, "' wasn't found.\nPlease make sure to create such a column and then rerun this command.\n"))
    }
    return(invisible())
  }
  
  # Determine which metadata columns to keep (i.e. remove the excluded ones)
  metadata_to_keep <- setdiff(metadata_colnames, excluded_names)
  
  # Remove all other metadata columns
  for (column in metadata_to_keep) {
    seu_slim@meta.data[[column]] <- NULL
  }
  
  # seu_slim <- DietSeurat(object = seu_slim,layers = "counts")
  
  return(seu_slim)
}

get_four_integration_plots <- function(seu){
  
  p_dataset <- dittoDimPlot(seu,"orig.ident", reduction.use = "umap", size = 0.5) +
    ggplot2::ggtitle(label = "UMAP by DataSet of Origin")
  
  p_phase <- dittoDimPlot(seu,"Phase", reduction.use = "umap", size = 0.5) +
    ggplot2::ggtitle(label = "UMAP by Cell Cycle Phase")
  
  p_condition<- dittoDimPlot(seu,"condition", reduction.use = "umap", size = 0.5) +
    ggplot2::ggtitle(label = "UMAP by condition")
  
  
  #Creating a color_palate big enough for all the sampels
  # n <- length(unique(seu@meta.data$sample))
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  
  p_cell_type <- dittoDimPlot(seu,"final_cell_type", reduction.use = "umap", size = 0.5 #, color.panel=col_vector
  ) +
    ggplot2::ggtitle(label = "UMAP by Cell Type")
  
  
  #Making one pdf out of it
  
  p_final <- ((p_dataset + p_phase) / (p_condition + p_cell_type) ) 
  
  return(p_final)
  
}



# forset plots------------------------------------------------------------------

get_table_for_forest <- function(output_list, gene) {
  # Builds a “forest” summary table for one gene across many samples
  
  forest_table <- data.frame(
    dataset    = character(),
    sample     = character(),
    avg_log2FC = numeric(),
    pct_2      = numeric(),
    CI_lower   = numeric(),
    CI_upper   = numeric(),
    p_val_adj  = numeric(),
    stringsAsFactors = FALSE
  )
  
  pct1_vals <- numeric()  # ← NEW: collect pct.1 values
  
  # -------- main loop --------
  for (dataset_name in names(output_list)) {
    for (sample_name in names(output_list[[dataset_name]])) {
      
      final_tbl <- output_list[[dataset_name]][[sample_name]]$original_table
      
      if (gene %in% rownames(final_tbl)) {
        row <- final_tbl[ gene, , drop = FALSE ]
        
        # Wald-type CI from p_val_adj (or use row$p_val if preferred)
        SE <- abs(row$avg_log2FC) / qnorm(1 - row$p_val / 2)
        
        forest_table <- rbind(
          forest_table,
          data.frame(
            dataset    = dataset_name,
            sample     = sample_name,
            avg_log2FC = row$avg_log2FC,
            pct_2      = row$pct.2,
            CI_lower   = row$avg_log2FC - 1.96 * SE,
            CI_upper   = row$avg_log2FC + 1.96 * SE,
            p_val_adj  = row$p_val_adj,
            stringsAsFactors = FALSE
          )
        )
        
        pct1_vals <- c(pct1_vals, row$pct.1)  # ← NEW: store pct.1 for this comparison
      } else {
        forest_table <- rbind(
          forest_table,
          data.frame(
            dataset    = dataset_name,
            sample     = sample_name,
            avg_log2FC = NA,
            pct_2      = NA,
            CI_lower   = NA,
            CI_upper   = NA,
            p_val_adj  = NA,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  # -------- mean pct.1 across all comparisons --------
  pct_1 <- if (length(pct1_vals) > 0) mean(pct1_vals, na.rm = TRUE) else NA  # ← CHANGED
  
  forest_table$pct_1 <- pct_1
  return(forest_table)
}

create_forest_plot <- function(forest_table, gene) {
  forest_table <- forest_table[order(forest_table$avg_log2FC, decreasing = FALSE), ]
  forest_table$sample <- factor(forest_table$sample, levels = forest_table$sample)
  
  y_lines <- c(0)
  
  forest_plot <- ggplot(forest_table, aes(y = avg_log2FC, x = sample)) +
    geom_point(aes(size = pct_2, color = dataset), shape = 16, alpha = 0.8) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, color = "black") +
    geom_hline(yintercept = y_lines, linetype = "dashed", color = "black", size = 0.8) +
    scale_size_continuous(range = c(3, 8)) +
    scale_y_continuous(breaks = c(scales::pretty_breaks(n = 10)(forest_table$avg_log2FC), y_lines)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      plot.title   = element_markdown(hjust = 0.5, size = 14, color = "black")
    ) +
    ylab("Log2 Fold Change") +
    ggtitle(
      paste0(
        nrow(forest_table),
        " Comparisons of <b>", gene, "</b> Expression (PCT.1 = ",
        round(unique(forest_table$pct_1)[1], 3), ")"
      )
    ) +
    scale_color_manual(values = c("#ffc857", "#e9724c", "#c5283d", "#481d24",
                                  "#976391", "#255f85", "#2a9d8f", "#283618", "#606c38")) +
    guides(
      size  = guide_legend(title = "PCT.2"),
      color = guide_legend(title = "Dataset")
    )
  
  return(forest_plot)
}

get_4_forest_plots <- function(output_list, four_genes_vector) {
  # --------------------------- sanity checks --------------------------- #
  if (!is.vector(four_genes_vector) || length(four_genes_vector) != 4) {
    stop(sprintf(
      "four_genes_vector must contain exactly four genes – you supplied %d.",
      length(four_genes_vector)
    ))
  }
  
  # --------------------------- build plots ----------------------------- #
  plot_list <- setNames(vector("list", 4), paste0("plot_", LETTERS[1:4]))
  
  for (i in seq_along(four_genes_vector)) {
    gene          <- four_genes_vector[i]
    forest_table  <- get_table_for_forest(output_list, gene)
    forest_plot   <- create_forest_plot(forest_table, gene)
    
    # save to the correct slot (A, B, C, D)
    plot_list[[i]] <- forest_plot
  }
  
  # --------------------------- tag panels ------------------------------ #
  tag_vec <- LETTERS[1:4]
  plot_list <- Map(
    function(p, tag) {
      p + ggplot2::annotate(
        "text",
        x = -Inf, y = Inf,
        label = tag,
        hjust = -0.5, vjust = 1.5,
        size = 6, fontface = "bold"
      )
    },
    plot_list, tag_vec
  )
  names(plot_list) <- paste0("plot_", tag_vec)
  
  # --------------------------- assemble grid --------------------------- #
  # Requires library(patchwork)
  combined_plot <- (plot_list$plot_A + plot_list$plot_B) /
    (plot_list$plot_C + plot_list$plot_D)
  
  combined_plot   # returned invisibly if you call without assignment
}



#meta_scatter_plot------------------------------------------------------------------
get_dsi_table <- function(paths, genes) {
  
  # Description:
  #   Reads a set of “meta_table” files, parses metadata from each filename,
  #   extracts DSI values for specified genes, and merges with a static
  #   `meta_info_df` of sample annotations.
  # Arguments:
  #   paths        – character vector of file paths to `.qs` meta_table files.
  #   genes        – character vector of gene symbols for which to extract DSI.
  # Returns:
  #   A data.frame with one row per sample, containing:
  #     • parsed `condition`, `dataset`, `sample`
  #     • DSI(<gene>) columns for each gene in `genes`
  #     • joined metadata: age, sex, EDSS, treatment, number_of_cells, etc.
  
  
  # Prepare an empty results data frame
  result_df <- data.frame()
  
  parse_path_info <- function(path) {
    file_name <- basename(path)
    file_name <- sub("_meta_table\\.qs$", "", file_name)
    parts <- unlist(strsplit(file_name, "_"))
    
    condition_parsed <- NA
    dataset_parsed   <- NA
    sample_parsed    <- NA
    
    if (length(parts) >= 3) {
      condition_parsed <- parts[1]
      condition_parsed <- case_when(
        condition_parsed == "CIS" ~ "Clinically Isolated Syndrome",
        condition_parsed == "MOG" ~ "Myelin Oligodendrocyte Glycoprotein",
        condition_parsed == "NMO" ~ "Neuromyelitis Optica",
        condition_parsed == "MS"  ~ "Multiple Sclerosis",
        TRUE ~ condition_parsed
      )
      dataset_parsed <- parts[2]
      
      if (dataset_parsed == "Ramesh" && length(parts) >= 4) {
        sample_parsed <- parts[4]
      } else if (dataset_parsed == "Esaulova" && length(parts) >= 4) {
        sample_parsed <- paste(parts[3], parts[4], sep = "_")
      } else {
        sample_parsed <- parts[3]
      }
    }
    
    list(condition = condition_parsed, dataset = dataset_parsed, sample = sample_parsed)
  }
  
  # Loop over each path
  for (path in paths) {
    # Parse condition, dataset, and sample from path
    info <- parse_path_info(path)
    # Skip if missing info
    if (is.na(info$condition) || is.na(info$dataset) || is.na(info$sample)) {
      next
    }
    
    # Load meta_table
    meta_table <- qread(path)
    
    # Initialize a single row for this sample
    row_data <- data.frame(
      condition = info$condition,
      dataset   = info$dataset,
      sample    = info$sample,
      stringsAsFactors = FALSE
    )
    
    # For each gene in the genes vector, create only the dsi_* columns
    for (gene in genes) {
      dsi_col <- paste0("DSI(", gene,")")
      
      if (gene %in% rownames(meta_table)) {
        row_data[[dsi_col]] <- meta_table[gene, "dsi"]
      } else {
        row_data[[dsi_col]] <- NA
      }
    }
    
    # Append to the cumulative results
    result_df <- rbind(result_df, row_data)
  }
  
  ################################################################
  # Build meta_info_df with consistent column names
  ################################################################
  meta_info_df <- data.frame(
    sample = c("JYJ", "KJS", "YYW",
               "10", "22", "24", "25", "27", "28",
               "29", "31a", "31b", "5", "6", "7",
               "8", "9", "16", "18", "mog_blood",
               "ms1_blood", "ms2_blood", "MS19270", "MS71658", "MS49131",
               "MS60249", "MS74594", "cd59patient", "cd55patient","CYG","KHW","KSH","SYG","YYJ","YYS"),
    condition_accurate = c("Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis",
                           "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis",
                           "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis",
                           "Multiple Sclerosis", "Ab-mediated IDD", "Multiple Sclerosis", "Multiple Sclerosis",
                           "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis", "Ab-mediated IDD",
                           "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis", "Multiple Sclerosis",
                           "Multiple Sclerosis", "Multiple Sclerosis", "CD59 def", "CD55 def","Ab-mediated IDD","Ab-mediated IDD","Ab-mediated IDD","Ab-mediated IDD","Ab-mediated IDD","Ab-mediated IDD"),
    condition_subtype = c("MS-unknown", "MS-unknown", "MS-unknown",
                          "RRMS", "CIS", "RRMS", "RRMS", "RRMS", "RRMS",
                          "RRMS", "CIS", "CIS", "NMOSD", "RRMS", "RRMS",
                          "RRMS", "RRMS", "RRMS", "RRMS", "MOGAD",
                          "RRMS", "RRMS", "CDMS", "RRMS", "RRMS",
                          "RRMS", "RRMS", "CD59 def", "CD55 def", "Ab-mediated IDD-unknown", "Ab-mediated IDD-unknown", "Ab-mediated IDD-unknown", "Ab-mediated IDD-unknown", "Ab-mediated IDD-unknown", "Ab-mediated IDD-unknown"),
    EDSS = c(1.5, 3, 3.5,
             2, 4, 0, 1.5, 1.5, 4,
             2, 2.5, 2.5, NA, 2, 1.5,
             1.5, 3.5, 4, 2, NA,
             NA, NA, 0.5, 6, 2.5,
             1, NA, NA, NA, 3, 4, 5, 3, 3, 3),
    number_of_cells = c(1101, 1273, 697,
                        453, 307, 833, 414, 460, 220,
                        1107, 543, 601, 2890, 1683, 529,
                        376, 228, 1278, 472, 1479,
                        176, 279, 670, 631, 472,
                        330, 1236, 205, 266,1253, 2785, 1742, 1171, 165, 2774),
    age = c(47, 25, 31,
            45, 35, 45, 42, 37, 50,
            54, 53, 53, 72, 44, 37,
            32, 27, 22, 29, 56,
            36, 36, 35, 47, 47,
            28, 42, NA, NA, 63, 39, 50, 31, 47, 67),
    treatment = c("Non", "Non", "Non",
                  "Non", "Non", "Non", "Non", "Non", "Non",
                  "Non", "Non", "Non", "Non", "Non", "Non",
                  "Non", "Non", "Non", "Steroids", "Non",
                  "Non", "Non", "Non", "Non", "Non",
                  "Non", "Non", "Non", "Non", "Non", "Non",
                  "Non", "Non", "Non", "Non"),
    sex = c("Female", "Male", "Female",
            "Female", "Female", "Female", "Female", "Female", "Female",
            "Male", "Female", "Female", "Female", "Male", "Male",
            "Male", "Female", "Male", "Female", "Female",
            "Female", "Female", "Female", "Female", "Female",
            "Female", "Male", "Undefined", "Undefined","Male", "Male","Male", "Male","Female", "Female"),
    stringsAsFactors = FALSE
  )
  
  # Filter out sample 18
  meta_info_df <- meta_info_df %>% filter(sample != "18")
  
  # Validate that row counts match
  if (nrow(result_df) != nrow(meta_info_df)) {
    stop("The row counts differ between 'result_df' and 'meta_info_df'. Check sample names.")
  }
  
  # Join by sample
  result_df <- left_join(result_df, meta_info_df, by = "sample")
  
  # Double-check row count
  if (nrow(result_df) != nrow(meta_info_df)) {
    stop("Row count changed after join. Check sample name mismatches.")
  }
  
  # Reorder columns (and also rename & delete one)
  result_df <- result_df %>% dplyr::select(-condition) %>% 
    dplyr::select(c(
      "sample",
      "condition_accurate",
      "condition_subtype",
      "dataset",
      "age",
      "sex",
      "EDSS",
      "treatment",
      "number_of_cells"), everything()) %>% 
    dplyr::rename(condition = "condition_accurate", condition_subtype = "condition_subtype")
  
  
  # Return the final assembled table
  return(result_df)
}

get_meta_scatter_plot <- function(df, x_var, y_var, color_var, shape_var,
                                  color_type = c("discrete", "scale"), color_palette = NULL,
                                  scale_min = NULL, scale_max = NULL, legend_order = c("no", "yes")) {
  
  
  # Description:
  #   Creates a customizable scatter plot of two numeric variables across samples,
  #   colored and shaped by metadata columns, with options for discrete vs.
  #   continuous coloring.
  # Arguments:
  #   df           – data.frame of metadata and numeric variables.
  #   x_var        – character; name of column for x‐axis.
  #   y_var        – character; name of column for y‐axis.
  #   color_var    – character; name of column to map to color.
  #   shape_var    – character; name of column to map to point shape.
  #   color_type   – "discrete" or "scale"; whether to use categorical or
  #                  continuous color scale.
  #   color_palette– NULL or vector; palette name or named vector of colors.
  #   scale_min    – numeric; lower bound for continuous scale (if used).
  #   scale_max    – numeric; upper bound for continuous scale (if used).
  #   legend_order – "yes" or "no"; whether to reorder legend levels for discrete.
  # Returns:
  #   A ggplot object of x_var vs y_var with specified aesthetics.
  
  color_type   <- match.arg(color_type)
  legend_order <- match.arg(legend_order)
  df_plot      <- df
  
  # Symmetrical x-axis and y-axis limits
  x_max <- max(abs(df_plot[[x_var]]), na.rm = TRUE)
  y_max <- max(abs(df_plot[[y_var]]), na.rm = TRUE)
  
  # Convert legend names to title case
  color_var_title <- tools::toTitleCase(color_var)
  shape_var_title <- tools::toTitleCase(shape_var)
  
  # Reorder factor levels for discrete if requested
  if (
    color_type == "discrete" &&
    legend_order == "yes" &&
    !is.null(color_palette) &&
    length(color_palette) > 1 &&
    !is.null(names(color_palette))
  ) {
    df_plot[[color_var]] <- factor(df_plot[[color_var]], levels = names(color_palette))
  }
  
  # Initialize ggplot
  p <- ggplot(
    data = df_plot,
    aes(
      x = !!sym(x_var),
      y = !!sym(y_var),
      color = !!sym(color_var),
      shape = !!sym(shape_var)
    )
  ) +
    geom_point(size = 6, alpha = 0.9, stroke = 1.2) +
    theme_minimal(base_size = 16) +
    labs(
      title = paste("<b>", x_var, "</b>", "vs", "<b>", y_var, "</b>", "Across Neuroimmunological Samples Colored by", "<b>", color_var_title, "</b>"),
      x     = x_var,
      y     = y_var,
      color = color_var_title,
      shape = shape_var_title
    ) +
    theme(
      legend.position   = "right",
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.5, "cm"),
      plot.title        = ggtext::element_markdown(hjust = 0.5, size = 18),
      axis.title        = element_text(size = 16),
      legend.text       = element_text(size = 14),
      legend.title      = element_text(size = 16),
      panel.grid.major  = element_line(linewidth = 0.5, linetype = "dotted", color = "gray80")
    ) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    # Ensure symmetrical axes
    scale_x_continuous(limits = c(-x_max, x_max)) +
    scale_y_continuous(limits = c(-y_max, y_max)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8)
  
  # Discrete vs. continuous coloring
  if (color_type == "discrete") {
    if (is.null(color_palette)) {
      p <- p + scale_color_brewer(palette = "Dark2")
    } else if (length(color_palette) == 1) {
      p <- p + scale_color_brewer(palette = color_palette)
    } else {
      p <- p + scale_color_manual(values = color_palette)
    }
    # Ensure shape legend appears **above** color legend
    p <- p + guides(
      shape = guide_legend(order = 1, override.aes = list(size = 6, alpha = 1), title = shape_var_title),
      color = guide_legend(order = 2, override.aes = list(size = 6, alpha = 1), title = color_var_title)
    )
  } else {
    # color_type == "scale" => continuous color
    if (is.null(color_palette)) {
      color_palette <- c("#ffcccc", "#ff0000")
    }
    # Possibly handle NA color if user gave a 3-element palette with a named "NA"
    na_col <- "grey"
    if (length(color_palette) >= 3 && !is.null(names(color_palette)) && "NA" %in% names(color_palette)) {
      na_col <- color_palette["NA"]
    }
    
    if (!is.null(scale_min) && !is.null(scale_max)) {
      p <- p + scale_color_gradient(
        low     = color_palette[1],
        high    = color_palette[2],
        limits  = c(scale_min, scale_max),
        oob     = scales::squish,
        na.value = na_col
      )
    } else {
      p <- p + scale_color_gradient(
        low     = color_palette[1],
        high    = color_palette[2],
        na.value = na_col
      )
    }
    # Ensure shape legend appears **above** color legend
    p <- p + guides(
      shape = guide_legend(order = 1, override.aes = list(size = 6, alpha = 1), title = shape_var_title),
      color = guide_colorbar(order = 2, barheight = unit(4, "cm"), title = color_var_title)  # single vertical color bar
    )
  }
  
  return(p)
}

get_summarize_meta_scatter_plot <- function(df, x_var, y_var){
  
  
  # Description:
  #   Builds a composite set of meta‐scatter plots for two variables:
  #     1. condition vs subtype colored panels (A/B)
  #     2. sex, EDSS, age, and cell count panels (C/D)
  # Arguments:
  #   df           – data.frame with metadata and numeric variables.
  #   x_var        – character; name of column for x‐axis.
  #   y_var        – character; name of column for y‐axis.
  # Returns:
  #   A list with:
  #     • p_condition_and_sub_condition – combined plot A/B
  #     • p_4_meta_plots               – 2×2 grid of C–D plots
  
  
  # p_condition_and_sub_condition
  
  {
    
    my_color_platte <- c(
      "CD59 def" = "#2F5233",   # Deep Forest Green  
      "CD55 def" = "#8C9048",   # Muted Olive Green (More Yellowish, Distinct from MS)
      "Multiple Sclerosis" = "#74D3AE",  # Soft Teal-Green  
      "Ab-mediated IDD" = "#DD9787"      # Warm Muted Coral  
    )
    
    
    
    p_condition_color <-  get_meta_scatter_plot(df,x_var = x_var, y_var = y_var, color_var="condition", shape_var = "dataset",color_type = "discrete",color_palette =my_color_platte, legend_order = "yes")
    
    
    
    
    
    my_color_platte <- c(
      "CD59 def" = "#AE2012",   # Dark Red  
      "CD55 def" = "#7D1D6E",   # Deep Magenta (distinct from CD59 def)
      "NMOSD" = "#A35759",       # Deep Orange
      "MOGAD" = "#EE9B00",      # Bright Orange
      "Ab-mediated IDD-unknown" = "#FFD65A",
      "CDMS" = "#00A9FE",       # Soft Teal
      "CIS" = "#003092",        # Dark Teal
      "RRMS" = "#001219",       # Deep Navy
      "MS-unknown" = "#94D2BD" # Pale Mint
    )
    
    
    p_condition_subtype_color <-  get_meta_scatter_plot(df,x_var = x_var, y_var = y_var, color_var="condition_subtype", shape_var = "dataset",color_type = "discrete",color_palette =my_color_platte, legend_order = "yes")
    
    
    p_condition_color <- p_condition_color +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
    
    p_condition_subtype_color <- p_condition_subtype_color +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
    
    
    
    p_condition_and_sub_condition<- p_condition_color | p_condition_subtype_color  
    
    # print(p_condition_and_sub_condition)
  }
  
  
  
  
  #p_4_meta_plots: sex, EDSS, age and num_cells
  
  {
    
    my_color_platte <- c(
      "Female" = "#9C3848",   # Dark Red  
      "Male" = "#1E3888",   # Deep Magenta (distinct from CD59 def)
      "Undefined" = "#6C757D")
    
    
    
    p_sex_color <-  get_meta_scatter_plot(df,x_var = x_var, y_var = y_var, color_var="sex", shape_var = "dataset",color_type = "discrete",color_palette =my_color_platte, legend_order = "yes")
    
    
    
    
    my_color_platte <- c(
      "min" = "#F6CDCD",  # Very light pink
      "max" = "#E61717",  # Deep, dark red
      "NA"  = "black"   # Gray for missing
    )
    p_edss_color <-  get_meta_scatter_plot(df,x_var = x_var, y_var = y_var, color_var="EDSS", shape_var = "dataset",color_type = "scale",color_palette =my_color_platte, scale_min = 0, scale_max = 6, legend_order = "yes")
    
    
    
    
    my_color_platte <- c(
      "min" = "#ABE7FE",  # Very light pink
      "max" = "#082777",  # Deep, dark red
      "NA"  = "black"   # Gray for missing
    )
    
    
    p_age_color <-  get_meta_scatter_plot(df,x_var = x_var, y_var = y_var, color_var="age", shape_var = "dataset",color_type = "scale",color_palette =my_color_platte, scale_min = 20, scale_max = 75,legend_order = "yes")
    
    
    
    
    
    my_color_platte <- c(
      "min" = "#82D523",   # Very Light Pink
      "max" = "#094F29",   # Bold, Deep Purple (Stronger High-End Contrast)
      "NA"  = "black"      # Keeps NA values black
    )
    
    
    
    
    
    p_num_cells_color <-  get_meta_scatter_plot(df,x_var = x_var, y_var = y_var, color_var="number.of.cells", shape_var = "dataset",color_type = "scale",color_palette =my_color_platte, scale_min = 100, scale_max = 3000,legend_order = "yes")
    
    # print(p_num_cells_color)
    
    # Add letters to each UMAP plot on the leftmost side
    p_sex_color <- p_sex_color +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
    
    p_edss_color <- p_edss_color +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
    
    p_age_color <- p_age_color +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
    
    p_num_cells_color <- p_num_cells_color +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold")
    
    # Combine plots into one
    p_4_meta_plots <- ((p_sex_color + p_edss_color) / (p_age_color + p_num_cells_color))
    
    # print(p_4_meta_plots)
    
    
  }
  
  return_list <- list(p_condition_and_sub_condition = p_condition_and_sub_condition, 
                      p_4_meta_plots = p_4_meta_plots)
  
  return(return_list)
  
  
}


#get_meta_int_comparison_pds_for_all_the_samples and for one sample------------------------------------------------------------------

get_meta_int_comparison_pds_for_all_the_samples <- function(
    integrated_path,
    meta_path,
    files_to_exclude_int = NULL,
    files_to_exclude_meta = NULL,
    cell_type = "NKs",
    save_dir_gsea
) {
  
  # 1) Get paths for each patient
  patient_path_list <- get_patient_path_list(
    integrated_path = integrated_path,
    meta_path = meta_path,
    files_to_exclude_int = files_to_exclude_int,
    files_to_exclude_meta = files_to_exclude_meta
  )
  
  
  # 2) Initialize an empty list to store results
  patient_results_list <- list()
  
  # 3) Iterate over all patients
  for (patient_name in names(patient_path_list)) {
    
    ## Parse condition, dataset, sample
    {
      parts <- unlist(strsplit(patient_name, "_"))
      condition <- parts[1]
      
      if (parts[2] == "Ramesh") {
        dataset_name <- "Ramesh"
        # Skip parts[3] ("5" or "3"), combine from parts[4] onwards
        sample_name <- paste(parts[4:length(parts)], collapse = "_")
      } else {
        dataset_name <- parts[2]
        sample_name  <- paste(parts[3:length(parts)], collapse = "_")
      }
    }
    
    # 4) Read meta_table and integrated_list
    {
      patient_data <- patient_path_list[[patient_name]]
      meta_table <- qread(patient_data$meta_table)
      
      integrated_list <- qread(patient_data$integrated_object)
      # Ensure 'gene' column in integrated_list$de_table
      integrated_list$de_table$gene <- rownames(integrated_list$de_table)
    }
    
    # 5) Generatede gene plots 
    get_de_genes_comparison(
      meta_table       = meta_table,
      integrated_list  = integrated_list,
      highlight_genes  = NULL,
      sample_name      = sample_name,
      dataset_name     = dataset_name,
      cell_type        = cell_type
    )
    condition = names(which.min(table(integrated_list$integrated_seu@meta.data$condition)))
    
    # 6) Collect GSEA and ORA results
    gsea_list <- get_gsea(
      meta_table      = meta_table,
      integrated_list = integrated_list,
      cell_type       = cell_type,
      condition       = condition,
      dataset_name    = dataset_name, 
      sample_name     = sample_name)
    
    clean_gsea_list <- function(lst) {
      for (nm in names(lst)) {
        x <- lst[[nm]]
        is_empty <-
          is.null(x) ||                                   # already NULL
          (is.atomic(x) && length(x) == 1 && is.na(x)) || # single NA
          (inherits(x, "data.frame") && nrow(x) == 0)     # 0‑row data.frame / tibble
        if (is_empty) lst[[nm]] <- NULL                   # remove element
      }
      lst
    }
    
    ## ── remove empty pathway dfs
    gsea_list <- clean_gsea_list(gsea_list)
    
    full_gsea_dfs <- get_full_gsea_and_ora_dfs(gsea_return_list = gsea_list,
                                          ora_return_list  = NULL,
                                          save_dir_gsea    = save_dir_gsea,
                                          save_dir_ora     = NULL,
                                          apply_filter     = FALSE,
                                          condition = condition,
                                          dataset_name = dataset_name, 
                                          sample_name = sample_name) 
   
    
    # ora_list <- get_ora(
    #   meta_table      = meta_table,
    #   integrated_list = integrated_list,
    #   sample_name     = sample_name,
    #   dataset_name    = dataset_name,
    #   cell_type       = cell_type
    # )
    
    
    
    
    # 7) Store them in patient_results_list
    patient_results_list[[patient_name]] <- list(
      gsea_list = gsea_list,
      full_gsea_dfs = full_gsea_dfs
      # ora_list  = ora_list
    )
    
    print(paste("just finished", patient_name))
  }
  
  # 9) Return the named list of results
  return(patient_results_list)
}


#get_de_genes_comparison and its helper functions------------------------------------------------------------------

get_de_genes_comparison <- function(meta_table,integrated_list, highlight_genes = NULL, sample_name,dataset_name, cell_type = "NKs") {
  
  #Creates de_genes_plots PDF: the first page is the dsi plot and below the volcano.
  #The second PDF page is the Jaccard barplot and below the linear regression plot.
  
  # Description:
  #   Runs a side-by-side comparison of DE results from a meta‐analysis table and an integrated Seurat object,
  #   generating scatter, volcano, Jaccard, and linear regression plots, and saves them as a multipanel PDF.
  # Arguments:
  #   meta_table      – data.frame of meta‐analysis results, one row per gene.
  #   integrated_list – list with elements:
  #                      $integrated_seu (Seurat object with exactly two conditions)
  #                      $de_table      (data.frame of integrated DE results, one row per gene).
  #   highlight_genes – optional character vector of gene symbols to label prominently.
  #   sample_name     – character; name of the sample (for titles and filenames).
  #   dataset_name    – character; name of the dataset (for titles and filenames).
  #   cell_type       – character; cell type label for plot grouping (default "NKs").
  # Returns:
  #   NULL (invisibly); side‐effect is creation of PDF with all four comparison panels.
  
 
  # check of intewgrated seu has only 2 conditions
  if (!(length(unique(integrated_list$integrated_seu@meta.data$condition)) == 2)) {
    
    unique_conditions <- unique(integrated_list$integrated_seu@meta.data$condition)
    num_conditions <- length(unique_conditions)
    
    error_message <- paste(
      "Error: there are", num_conditions, "unique conditions in the integrated Seurat object, which are:",
      paste(unique_conditions, collapse = ", "),
      "\nPlease make sure there are exactly two conditions."
    )
    
    stop(error_message)
  }
  
  
  de_table <- integrated_list$de_table
  
  #adding gene column in both dataframes
  meta_table$gene <- rownames(meta_table)
  de_table$gene <- rownames(de_table)
  meta_table$avg_log2FC <- NULL
  
  group_a = names(which.min(table(integrated_list$integrated_seu@meta.data$condition)))
  
  group_b = names(which.max(table(integrated_list$integrated_seu@meta.data$condition)))
  
  title_suffix <- paste0(group_a, " ",cell_type," (", sample_name, " Sample, ", dataset_name, " Dataset)")
  
  # scatter_plot
  if (!is.null(highlight_genes)) {
    scatter_plot_list <- get_scatter_plot(meta_table, costume_de_genes = highlight_genes, title_suffix = title_suffix)
  } else {
    scatter_plot_list <- get_scatter_plot(meta_table, title_suffix = title_suffix)
  }
  
  
  
  #volcano_plot
  {
    
    
    key_genes <- scatter_plot_list$de_genes_table
    
    plottitle <- paste(group_a, " ",cell_type," (", sample_name, " Sample, ", dataset_name, " Dataset) vs. ", group_b, " ",cell_type," - Volcano Plot", sep = "")
    
    
    
    volcano_plot <- get_volcano(de_results = de_table,key_genes = key_genes, plottitle = plottitle, group_a = group_a, group_b = group_b)
  }
  
  
  #jaccard_plot
  
  jaccard_list <- get_jaccard_plot(meta_table, de_table, thresholds = c(20, 50, 100), title_suffix = title_suffix)
  
  
  #linear_regression_plot
  linear_regression <- get_linear_regression(meta_table, de_table, highlight_genes, title_suffix = title_suffix)
  
  
  #PDF
  {
    dir_name <- paste0("/bigdata/bgutman/Collabs/Dror mevorach/plots/", cell_type)
    
    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
    }
    
    dir_name <- paste0("/bigdata/bgutman/Collabs/Dror mevorach/plots/", cell_type,"/",group_a,"_",dataset_name, "_", sample_name)
    
    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
    }
    
    #scatter_and_volcano_plots
    
    scatter_and_volcano_plots <-  ((scatter_plot_list$p_scatter) /  (volcano_plot)) 
    jaccard_and_linear_regression_plots <-  ((jaccard_list$jaccard_plot) /  (linear_regression$linear_regression_plot))
    
    full_pdf_name <- paste0(dir_name, "/",group_a,"_",dataset_name, "_", sample_name, "_de_genes_plots.pdf")
    
    pdf(full_pdf_name, width = 26, height = 20)
    print(scatter_and_volcano_plots)
    print(jaccard_and_linear_regression_plots)
    dev.off()
    
    
    
    
  }
  
  
  
}

get_volcano <- function(de_results, key_genes, plottitle, group_a, group_b) {
  
  # Description:
  #   Constructs a publication‐style volcano plot from DE results, highlighting housekeeping vs non‐housekeeping genes,
  #   capping –log10 p‐values at 280, and labeling a set of key genes.
  # Arguments:
  #   de_results – data.frame of DE results with columns avg_log2FC, p_val_adj, rownames as gene.
  #   key_genes  – character vector or data.frame with column gene_name of genes to label.
  #   plottitle  – character; plot title.
  #   group_a    – character; name of condition A.
  #   group_b    – character; name of condition B.
  # Returns:
  #   ggplot object of the volcano plot.
  
  # Load necessary libraries and source the custom theme
  {
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(gridExtra))
    suppressPackageStartupMessages(library(ggrepel))
    source('/bigdata/bgutman/Collabs/Dror mevorach/scripts/ggplot_theme_Publication-2.R')
  }
  
  # Calculate the maximum absolute value on the x-axis
  max_x_value <- max(abs(de_results$avg_log2FC), na.rm = TRUE)
  
  # Read in housekeeping gene lists to label "Housekeeping Gene" vs. "Non-Housekeeping Gene"
  joshi_hk_genes <- read_excel("/bigdata/bgutman/Collabs/Dror mevorach/data/joshi_low_gini_coefficients_genes.xlsx")
  levanon_hk_genes <- read_excel("/bigdata/bgutman/Collabs/Dror mevorach/data/levanon_housekeeping_genes.xlsx")
  hk_genes <- unique(c(joshi_hk_genes$hgnc_symbol, levanon_hk_genes$gene_name))
  
  # Create a gene_type column
  # Make sure rownames(de_results) are present
  de_results$gene <- rownames(de_results)
  de_results$gene_type <- ifelse(de_results$gene %in% hk_genes, 
                                 "Housekeeping Gene", 
                                 "Non-Housekeeping Gene")
  
  # Cap the -log10(p_val_adj) at 280 (per your code)
  de_results$y_display <- -log10(de_results$p_val_adj)
  de_results$y_display <- ifelse(de_results$y_display < 280, 
                                 de_results$y_display, 280)
  
  # If 'key_genes' is a dataframe with a 'gene_name' column, extract its vector
  # Otherwise, if it's already a vector, leave as is
  if ("gene_name" %in% colnames(key_genes)) {
    key_genes <- key_genes$gene_name
  }
  
  # Now build the plot
  # Now build the plot with minimal changes from your code
  p <- ggplot(de_results, aes(x = avg_log2FC, y = y_display)) +
    
    # Points with same style as get_scatter_plot
    geom_point(
      aes(fill = gene_type),
      shape = 21, size = 3, color = "black", stroke = 0.5, alpha = 0.9
    ) +
    
    # Label key genes, same style
    geom_label_repel(
      data = dplyr::filter(de_results, gene %in% key_genes),
      aes(label = gene, fill = gene_type),
      color = "white",            
      fontface = "bold",
      size = 5,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray40",
      segment.size = 0.4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    # X-axis scale
    scale_x_continuous(
      limits = c(-max_x_value - 1, max_x_value + 1),
      breaks = c(-10, -5, -1, 0, 1, 5, 10)
    ) +
    
    # Y-axis scale
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
    
    # Reference lines
    geom_hline(yintercept = 1.30, linetype = "dashed", color = "#03045e", size = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#03045e", size = 1) +
    
    # Axis labels + legend title
    labs(
      title = plottitle,              # same variable you used before
      x = "Log2 Fold Change",         # same text as original
      y = "-Log10 Adjusted P-value",  # same text as original
      fill = "Gene type"              # same as in second plot
    ) +
    
    # Map fill to your two categories
    scale_fill_manual(values = c("Housekeeping Gene" = "black",
                                 "Non-Housekeeping Gene" = "firebrick")) +
    
    # REPLACE theme_Publication() WITH theme_minimal & extra theme lines
    theme_minimal(base_size = 20) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title       = element_text(size = 20),
      axis.text        = element_text(size = 16),
      panel.grid.major = element_blank(),  # remove the default major grid
      panel.grid.minor = element_blank(),  # remove minor grid
      legend.position  = "bottom"
    ) +
    
    # Corner annotation
    annotate(
      "text", x = -Inf, y = -Inf,
      label = "The points were capped at 280 for display",
      hjust = 0, vjust = -1, color = "white", size = 3.5
    )
  
  
  return(p)
}

get_scatter_plot <- function(meta_table,
                             x_variable = "dsi",
                             y_variable = "median_log2FC",
                             upper_threshold= 0.75,
                             below_threshold= -0.75,
                             upper_y_line =1,
                             below_y_line =-1,
                             full_x_variable_name = "Differential Stability Index (DSI)",
                             full_y_variable_name = "Median log2 Fold Change",
                             display_house_keeping_labels = TRUE,
                             costume_de_genes = NULL,
                             title_suffix) {
  
  
  # Description:
  #   Creates a labeled scatter plot of two metrics (e.g. DSI vs median_log2FC),
  #   coloring points by housekeeping status and optionally highlighting a custom gene list.
  # Arguments:
  #   meta_table                  – data.frame with gene metrics and rownames as gene_name.
  #   x_variable                  – character; column name for x‐axis (default "dsi").
  #   y_variable                  – character; column name for y‐axis (default "median_log2FC").
  #   upper_threshold             – numeric; fraction of max(x) above which to auto‐label.
  #   below_threshold             – numeric; fraction of min(x) below which to auto‐label.
  #   upper_y_line, below_y_line  – numeric; y‐values for dashed horizontal lines.
  #   full_x_variable_name        – character; x‐axis label text.
  #   full_y_variable_name        – character; y‐axis label text.
  #   display_house_keeping_labels – logical; whether to exclude HK genes from auto labels.
  #   costume_de_genes            – optional character vector of gene symbols to force‐label.
  #   title_suffix                – character; suffix for plot title.
  # Returns:
  #   List with elements:
  #     p_scatter      (ggplot object)
  #     de_genes_table (data.frame of labeled genes)
  
  # browser()
  
  if (!(x_variable %in% colnames(meta_table))) {
    stop("The x_variable column does not exist in the meta_table")
  }
  if (!(y_variable %in% colnames(meta_table))) {
    stop("The y_variable column does not exist in the meta_table")
  }
  
  # needed changes for the meta_table
  {
    meta_table$gene_name <- rownames(meta_table)
    meta_table[[y_variable]] <- as.numeric(as.character(meta_table[[y_variable]]))
    meta_table[[x_variable]] <- as.numeric(as.character(meta_table[[x_variable]]))
  }
  
  joshi_hk_genes <- read_excel("/bigdata/bgutman/Collabs/Dror mevorach/data/joshi_low_gini_coefficients_genes.xlsx")
  levanon_hk_genes <- read_excel("/bigdata/bgutman/Collabs/Dror mevorach/data/levanon_housekeeping_genes.xlsx")
  hk_genes <- unique(c(joshi_hk_genes$hgnc_symbol, levanon_hk_genes$gene_name))
  meta_table$gene_type <- ifelse(meta_table$gene_name %in% hk_genes, "Housekeeping Gene", "Non-Housekeeping Gene")
  
  # Define axis values
  x_breaks <- seq(-24, 24, by = 8)
  y_breaks <- seq(-10, 10, by = 5)
  
  # Define y-axis limits (make sure they cover your range + margin)
  y_axis_limit <- max(abs(meta_table[[y_variable]]), na.rm = TRUE)
  y_axis_limit <- ceiling(max(y_axis_limit, 10))
  custom_y_breaks <- sort(unique(c(seq(-y_axis_limit, y_axis_limit, by = 5), below_y_line, upper_y_line)))
  
  # Define x-axis limit
  x_axis_limit <- max(meta_table$comparisons_num)
  
  # Calculate threshold for gene labeling
  computed_upper_threshold <- max(meta_table$comparisons_num) * upper_threshold
  computed_lower_threshold <- max(meta_table$comparisons_num) * below_threshold
  
  # Subset for labeling using vectorized operator '|'
  label_data <- meta_table %>%
    filter(.data[[x_variable]] >= computed_upper_threshold |
             .data[[x_variable]] <= computed_lower_threshold)
  
  if (!display_house_keeping_labels) {
    label_data <- label_data %>%
      filter(gene_type == "Non-Housekeeping Gene")
  }
  
  if (!is.null(costume_de_genes)) {
    label_data <- meta_table %>%
      filter(gene_name %in% costume_de_genes)
  }
  
  # Begin plot
  p_scatter_avg_log2fc <- ggplot(meta_table, aes_string(x = x_variable, y = y_variable)) + 
    # Points colored by gene type
    geom_point(aes(fill = gene_type),
               shape = 21, size = 3, color = "black", stroke = 0.5, alpha = 0.9) +
    scale_fill_manual(values = c("Housekeeping Gene" = "black", 
                                 "Non-Housekeeping Gene" = "firebrick")) +
    
    # Label selected genes
    geom_label_repel(
      data = label_data,
      aes(label = gene_name, fill = gene_type),
      color = "white",
      fontface = "bold",
      size = 5,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray40",
      segment.size = 0.4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    # Horizontal dashed lines at y = ±1
    geom_hline(yintercept = c(1, -1), linetype = "dashed", color = "gray40") +
    
    # Solid horizontal grid lines at -10, -5, 0, 5, 10
    geom_hline(yintercept = y_breaks, linetype = "solid", color = "gray85", size = 0.3) +
    
    # Solid vertical grid lines at x_breaks
    geom_vline(xintercept = x_breaks, linetype = "solid", color = "gray85", size = 0.3) +
    
    # Vertical line at x = 0 (stronger emphasis)
    geom_vline(xintercept = 0, color = "darkgrey", size = 0.8) +
    
    # Axis settings
    scale_x_continuous(
      name = full_x_variable_name,
      limits = c(-x_axis_limit, x_axis_limit),
      breaks = x_breaks
    ) +
    scale_y_continuous(
      name =  full_y_variable_name,
      limits = c(-y_axis_limit, y_axis_limit),
      breaks = custom_y_breaks
    ) +
    
    # Title and theme
    labs(
      title = paste(full_x_variable_name, "vs.", full_y_variable_name,"-",title_suffix),
      fill = "Gene type"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      panel.grid.major = element_blank(),  # Disable default grid
      panel.grid.minor = element_blank(),
      legend.position = "bottom"           # <<< Moved legend to bottom
    )
  
  return_list <- list(
    "p_scatter" = p_scatter_avg_log2fc,
    "de_genes_table" = label_data
  )
  
  return(return_list)
}





get_spearman <- function(
    meta_table, 
    de_table, 
    join_column = "gene",
    meta_main_rank_column,
    meta_tie_breaker_column = NULL,
    int_main_rank_column,
    int_tie_breaker_column = NULL
){
  

  # Description:
  #   Computes Spearman’s rank correlation between two ranking columns
  #   (with optional tie‐breaker columns) after merging meta and integrated tables.
  # Arguments:
  #   meta_table               – data.frame containing join_column and ranking columns.
  #   de_table                 – data.frame containing join_column and ranking columns.
  #   join_column              – character; column name on which to merge (default "gene").
  #   meta_main_rank_column    – character; column name in meta_table to rank.
  #   meta_tie_breaker_column  – optional character; second column in meta_table for tie‐breaking.
  #   int_main_rank_column     – character; column name in de_table to rank.
  #   int_tie_breaker_column   – optional character; second column in de_table for tie‐breaking.
  # Returns:
  #   List with element spearman_correlation (numeric).
  
  # 1) Ensure join_column exist in both data frames
  if (!(join_column %in% colnames(meta_table)) || !(join_column %in% colnames(de_table))) {
    stop(
      "Both meta_table and de_table must contain a '", join_column, "' column."
    )
  }
  
  # 2) Merge with suffixes so if a column appears in both tables,
  #    we get e.g. 'col.meta' and 'col.int'
  merged <- merge(
    meta_table, de_table,
    by        = join_column,
    suffixes  = c(".meta", ".int")
  )
  
  # ———— edge‐case: no shared keys ——————————————
  if (nrow(merged) == 0) {
    warning("No overlapping '", join_column, "' values; returning 0 spearman.")
    return(list(spearman_correlation = 0))
  }
  

  
  # 3) Identify the columns in the merged data we actually need
  meta_col <- find_merged_col(meta_main_rank_column,  was_in_meta = TRUE, meta_table = meta_table, de_table = de_table)
  int_col  <- find_merged_col(int_main_rank_column,   was_in_meta = FALSE, meta_table = meta_table, de_table = de_table)
  
  meta_tie_breaker_col <- find_merged_col(meta_tie_breaker_column, was_in_meta = TRUE, meta_table = meta_table, de_table = de_table)
  int_tie_breaker_col  <- find_merged_col(int_tie_breaker_column,  was_in_meta = FALSE, meta_table = meta_table, de_table = de_table)
  
  # 4) Check that all needed columns actually exist in merged
  needed_cols <- c(meta_col, int_col, meta_tie_breaker_col, int_tie_breaker_col)
  needed_cols <- needed_cols[!is.null(needed_cols)]
  
  missing_cols <- setdiff(needed_cols, colnames(merged))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in merged data:\n",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Helper function to parse->eval->numeric if needed
  safe_to_numeric <- function(x) {
    sapply(x, function(s) eval(parse(text = s)))
  }
  
  # 5) Convert and rank meta_col
  merged[[meta_col]] <- safe_to_numeric(merged[[meta_col]])
  merged$meta_rank   <- rank(-merged[[meta_col]], ties.method = "average")
  
  # 6) If meta tie-breaker is provided, incorporate it
  if (!is.null(meta_tie_breaker_col)) {
    merged[[meta_tie_breaker_col]] <- safe_to_numeric(merged[[meta_tie_breaker_col]])
    
    # Combine the two ranks to break ties
    merged$meta_rank <- rank(
      merged$meta_rank + rank(-merged[[meta_tie_breaker_col]], ties.method = "average"),
      ties.method = "average"
    )
  }
  
  # 7) Convert and rank int_col
  merged[[int_col]] <- safe_to_numeric(merged[[int_col]])
  merged$de_rank    <- rank(-merged[[int_col]], ties.method = "average")
  
  # 8) If integrated tie-breaker is provided, incorporate it
  if (!is.null(int_tie_breaker_col)) {
    merged[[int_tie_breaker_col]] <- safe_to_numeric(merged[[int_tie_breaker_col]])
    
    merged$de_rank <- rank(
      merged$de_rank + rank(-merged[[int_tie_breaker_col]], ties.method = "average"),
      ties.method = "average"
    )
  }
  
  # 9) Compute Spearman correlation
  spearman_result <- cor(merged$meta_rank, merged$de_rank, method = "spearman")
  
  return(list(spearman_correlation = spearman_result))
}

#Helper:
# The function locates the appropriate column name in merged
# for a given "desired_col" based on whether that col
# also existed in de_table.  e.g. if "dsi" is in *both* data frames,
# we want "dsi.meta" for the meta_table version or "dsi.int"
# for de_table version. If it was only in one data frame,
# then it won't have a suffix after merging and we just return "dsi".

# If 'desired_col' is NULL or empty, return NULL
find_merged_col <- function(desired_col, was_in_meta, meta_table, de_table) {
  
  
  # Description:
  #   Helper to determine the correct column name in a merged data.frame
  #   when a column exists in one or both input tables (adds .meta or .int suffix).
  # Arguments:
  #   desired_col – character; original column name to locate.
  #   was_in_meta – logical; TRUE if selecting from meta suffix, FALSE for integrated suffix.
  #   meta_table  – data.frame before merge.
  #   de_table    – data.frame before merge.
  # Returns:
  #   Character; actual column name in merged data.frame or NULL if desired_col is NULL.
  
  
  if (is.null(desired_col)) {
    return(NULL)
  }
  
  # Did meta_table have 'desired_col'?
  meta_has  <- desired_col %in% colnames(meta_table)
  # Did de_table have 'desired_col'?
  de_has    <- desired_col %in% colnames(de_table)
  
  # If the column was in *both* data frames, after merging we should see
  #   'desired_col.meta' or 'desired_col.int'
  # If it was in *only one*, it remains 'desired_col' in merged.
  
  # We only need the "meta" version if was_in_meta=TRUE;
  # otherwise we want the "int" version if was_in_meta=FALSE.
  
  if (meta_has && de_has) {
    # The suffix will be ".meta" or ".int"
    suffix <- if (was_in_meta) ".meta" else ".int"
    return(paste0(desired_col, suffix))
  } else {
    # Column is only in one data frame => no suffix
    return(desired_col)
  }
}

get_effect_size_correlation <- function(
    meta_table, de_table,
    join_column  = "gene",
    meta_effect_size_column,
    int_effect_size_column
) {
  # ------------------------------------------------------------
  # 0) Sanity‑check the merge key
  # ------------------------------------------------------------
  if (!(join_column %in% colnames(meta_table)) ||
      !(join_column %in% colnames(de_table))) {
    stop(sprintf(
      "Both meta_table and de_table must contain a '%s' column.", join_column
    ))
  }
  
  # ------------------------------------------------------------
  # 1) Merge, keeping .meta / .int suffixes for duplicated names
  # ------------------------------------------------------------
  merged <- merge(
    meta_table, de_table,
    by       = join_column,
    suffixes = c(".meta", ".int")
  )
  
  # ———— edge‐case: no overlapping keys ——————————————— 
  if (nrow(merged) == 0) { 
    warning("No overlapping '", join_column, "' values; returning 0 correlation.") 
    return(list(effect_size_pearson_correlation = 0)) 
  }
  
  # ------------------------------------------------------------
  # 2) Locate the correct column names *via the helper* ---------
  #    (this mirrors get_spearman’s logic)
  # ------------------------------------------------------------
  meta_col <- find_merged_col(
    desired_col = meta_effect_size_column,
    was_in_meta = TRUE,          # we want the meta‑table version
    meta_table  = meta_table,
    de_table    = de_table
  )
  
  int_col  <- find_merged_col(
    desired_col = int_effect_size_column,
    was_in_meta = FALSE,         # we want the integrated‑table version
    meta_table  = meta_table,
    de_table    = de_table
  )
  
  # ------------------------------------------------------------
  # 3) Verify that the columns actually exist
  # ------------------------------------------------------------
  needed_cols   <- c(meta_col, int_col)
  missing_cols  <- setdiff(needed_cols, colnames(merged))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in merged data:\n",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------
  # 4) Convert to numeric — same helper used in get_spearman ----
  # ------------------------------------------------------------
  safe_to_numeric <- function(x) sapply(x, function(s) eval(parse(text = s)))
  
  merged[[meta_col]] <- safe_to_numeric(merged[[meta_col]])
  merged[[int_col]]  <- safe_to_numeric(merged[[int_col]])
  
  # ------------------------------------------------------------
  # 5) Pearson correlation of the two effect‑size columns
  # ------------------------------------------------------------
  pearson_corr <- cor(
    merged[[meta_col]],
    merged[[int_col]],
    use    = "complete.obs",
    method = "pearson"
  )
  
  list(effect_size_pearson_correlation = pearson_corr)
}

get_jaccard_index <- function(
    meta_table, 
    de_table, 
    threshold,
    main_meta_rank_column = "dsi",
    secondry_meta_rank_column = "median_log2FC",
    main_int_rank_column = "avg_log2FC",
    secondry_int_rank_column = NULL
) {
  
  # Description:
  #   Computes Jaccard indices for the top and bottom N genes between two ranked tables.
  # Arguments:
  #   meta_table                 – data.frame with 'gene' and ranking columns.
  #   de_table                   – data.frame with 'gene' and ranking columns.
  #   threshold                  – integer; number of genes in top/bottom sets.
  #   main_meta_rank_column      – character; primary ranking column in meta_table.
  #   secondry_meta_rank_column  – optional character; tie‐breaker in meta_table.
  #   main_int_rank_column       – character; primary ranking column in de_table.
  #   secondry_int_rank_column   – optional character; tie‐breaker in de_table.
  # Returns:
  #   Named list with top_<N>_jaccard_index and bottom_<N>_jaccard_index (numeric).
  
  # 1) Ensure 'gene' column is present in both tables
  if (!("gene" %in% colnames(meta_table)) || !("gene" %in% colnames(de_table))) {
    stop("Both meta_table and de_table must contain a 'gene' column.")
  }
  
  # 2) Rank the meta_table
  #    - If secondry_meta_rank_column is NOT NULL, use it as tie-breaker.
  #    - Otherwise, rank solely by main_meta_rank_column.
  if (!is.null(secondry_meta_rank_column)) {
    meta_ranked <- meta_table[order(
      -meta_table[[main_meta_rank_column]],
      -meta_table[[secondry_meta_rank_column]]
    ), ]
  } else {
    meta_ranked <- meta_table[order(
      -meta_table[[main_meta_rank_column]]
    ), ]
  }
  
  # 3) Rank the de_table
  #    - If secondry_int_rank_column is NOT NULL, use it as tie-breaker.
  #    - Otherwise, rank solely by main_int_rank_column.
  if (!is.null(secondry_int_rank_column)) {
    de_ranked <- de_table[order(
      -de_table[[main_int_rank_column]],
      -de_table[[secondry_int_rank_column]]
    ), ]
  } else {
    de_ranked <- de_table[order(
      -de_table[[main_int_rank_column]]
    ), ]
  }
  
  # 4) Identify top and bottom genes in each ranked table
  top_genes_meta <- head(meta_ranked$gene, threshold)
  bottom_genes_meta <- tail(meta_ranked$gene, threshold)
  
  top_genes_de <- head(de_ranked$gene, threshold)
  bottom_genes_de <- tail(de_ranked$gene, threshold)
  
  # 5) Jaccard calculation helper
  jaccard_index <- function(set1, set2) {
    intersect_len <- length(intersect(set1, set2))
    union_len <- length(union(set1, set2))
    return(intersect_len / union_len)
  }
  
  # 6) Compute Jaccard indices for top and bottom sets
  top_index <- jaccard_index(top_genes_meta, top_genes_de)
  bottom_index <- jaccard_index(bottom_genes_meta, bottom_genes_de)
  
  # 7) Return a named list
  result <- list()
  result[[paste0("top_", threshold, "_jaccard_index")]] <- top_index
  result[[paste0("bottom_", threshold, "_jaccard_index")]] <- bottom_index
  
  return(result)
}


get_jaccard_plot <- function(meta_table, de_table, thresholds = c(20, 50, 100), title_suffix) {
  
  # function get_jaccard_plot
  # Description:
  #   Builds a bar plot of Jaccard indices across multiple thresholds,
  #   annotating Spearman’s rank correlation.
  # Arguments:
  #   meta_table   – data.frame for meta analysis.
  #   de_table     – data.frame for integrated DE.
  #   thresholds   – integer vector of top/bottom set sizes.
  #   title_suffix – character; suffix for plot title.
  # Returns:
  #   List with elements:
  #     jaccard_plot (ggplot object),
  #     jaccard_df   (data.frame of scores),
  #     spearman     (numeric correlation).
  
  
  all_jaccard <- lapply(thresholds, function(t) {
    get_jaccard_index(meta_table, de_table, t)
  })
  names(all_jaccard) <- as.character(thresholds)
  
  jaccard_df <- do.call(rbind, lapply(names(all_jaccard), function(th) {
    vals <- all_jaccard[[th]]
    data.frame(
      threshold = th,
      category  = c("top", "bottom"),
      value     = unlist(vals)
    )
  }))
  
  jaccard_df$category <- tools::toTitleCase(as.character(jaccard_df$category))
  jaccard_df$stall <- paste(jaccard_df$category, jaccard_df$threshold)
  
  unique_thresholds <- sort(as.numeric(unique(jaccard_df$threshold)))  
  stall_levels <- unlist(lapply(unique_thresholds, function(t) {
    c(paste("Top", t), paste("Bottom", t))
  }))
  jaccard_df$stall <- factor(jaccard_df$stall, levels = stall_levels)
  
  # Get stall names
  unique_stalls <- levels(jaccard_df$stall)
  
  # Define your default 3-threshold palette
  default_colors <- c(
    "Top 20"    = "#48CAE4",
    "Bottom 20" = "#03045e",
    "Top 50"    = "#B2DF8A",
    "Bottom 50" = "#2C6E49",
    "Top 100"   = "#FF8FA3",
    "Bottom 100"= "firebrick" 
  )
  
  # Apply the default palette if thresholds == 3
  if (length(thresholds) == 3) {
    my_colors <- default_colors[unique_stalls]
  } else {
    my_colors <- setNames(rainbow(length(unique_stalls)), unique_stalls)
  }
  
  #Calculate spearman
  spearman <- get_spearman(meta_table = meta_table, de_table = de_table, meta_main_rank_column = "dsi", meta_tie_breaker_column = "median_log2FC", int_main_rank_column ="avg_log2FC")
  
  jaccard_plot <- ggplot(jaccard_df, aes(x = stall, y = value, fill = stall)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.4, size = 4) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_manual(values = my_colors) +
    # Replaced theme_bw(base_size = 14) with theme_minimal(base_size = 20)
    theme_minimal(base_size = 20) +
    theme(
      # Matching the style from second plot:
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      panel.grid.major = element_blank(),  # already present, keep
      panel.grid.minor = element_blank(),
      legend.position = "none",            # same as your original code
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(x = "Gene Set", y = "Jaccard Index") +
    ggtitle(paste("Jaccard Index for Top & Bottom Genes -", title_suffix)) +
    annotate("text", 
             x = -Inf, 
             y = Inf, 
             label = paste("Spearman's rank correlation:", sprintf("%.2f", spearman$spearman_correlation)),
             hjust = -0.1, 
             vjust = 1.1, 
             size = 6, 
             fontface = "bold")

  
  
  return_list <- list(
    "jaccard_plot" = jaccard_plot,
    "jaccard_df" = jaccard_df,
    "spearman" = spearman$spearman_correlation
  )
  
  return(return_list)
}

get_linear_regression <- function(meta_table, 
                                  de_table, 
                                  highlight_genes = NULL,
                                  title_suffix) {
  
  # Description:
  #   Fits and plots a linear regression of integrated avg_log2FC vs meta median_log2FC,
  #   highlights special gene groups, and annotates slope and Pearson effect‐size correlation.
  # Arguments:
  #   meta_table      – data.frame with 'gene' and median_log2FC.
  #   de_table        – data.frame with 'gene' and avg_log2FC.
  #   highlight_genes – optional character vector to label on plot.
  #   title_suffix    – character; suffix for plot title.
  # Returns:
  #   List with elements:
  #     linear_regression_plot  (ggplot object),
  #     linear_regression_slope (numeric),
  #     effect_size_correlation (numeric).
  
  # Ensure gene names exist
  if (!("gene" %in% colnames(meta_table)) || !("gene" %in% colnames(de_table))) {
    stop("Both meta_table and de_table must contain a 'gene' column.")
  }
  
  # Inner join on gene
  merged <- merge(meta_table, de_table, by = "gene")
  
  # Check that necessary columns exist
  if (!("median_log2FC" %in% colnames(merged)) || !("avg_log2FC" %in% colnames(merged))) {
    stop("meta_table must contain 'median_log2FC'; de_table must contain 'avg_log2FC'.")
  }
  
  # Flag genes that start with "AC[digits]", "AL[digits]", "AP[digits]", "MT", or "RP"
  merged$special <- ifelse(
    grepl("^(AC[0-9]+|AL[0-9]+|AP[0-9]+|MT|RP)", merged$gene),
    "Non-Coding, Mitochondrial and Ribosomal Genes",
    "Other Genes"
  )
  
  # Fit linear model: avg_log2FC ~ median_log2FC
  lm_model <- lm(avg_log2FC ~ median_log2FC, data = merged)
  linear_regression_slope <- coef(lm_model)["median_log2FC"]
  
  # Prepare the subset of genes to highlight (if any)
  if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
    label_data <- subset(merged, gene %in% highlight_genes)
  } else {
    # Same columns, 0 rows => no errors from inherited aesthetics
    label_data <- merged[FALSE, ]
  }
  
  # Calculate correlation
  effect_size_correlation <- get_effect_size_correlation(
    meta_table             = meta_table,
    de_table               = de_table,
    meta_effect_size_column= "median_log2FC",
    int_effect_size_column = "avg_log2FC"
  )
  
  # Create the scatter plot
  p <- ggplot(merged, aes(x = median_log2FC, y = avg_log2FC, fill = special)) +
    geom_point(shape = 21, color = "black", size = 3, stroke = 0.5, alpha = 0.9) +
    
    geom_label_repel(
      data           = label_data,
      aes(label      = gene),
      color          = "white",
      fontface       = "bold",
      size           = 5,
      box.padding    = 0.5,
      point.padding  = 0.3,
      segment.color  = "gray40",
      segment.size   = 0.4,
      max.overlaps   = Inf,
      show.legend    = FALSE
    ) +
    
    geom_smooth(
      data   = merged,
      aes(x = median_log2FC, y = avg_log2FC, group = 1),
      method = "lm",
      se     = FALSE,
      color  = "#000000"
    ) +
    
    # Reference lines
    geom_vline(xintercept = 0, color = "darkgrey", size = 0.8) +
    geom_hline(yintercept = 0, color = "darkgrey", size = 0.8) +
    
    # Annotate slope
    annotate("text",
             x = -Inf, y = Inf,
             label = paste(
               "Linear Regression Slope:", round(linear_regression_slope, 2),
               "\nEffect Size Correlation (Pearson):", 
               round(effect_size_correlation$effect_size_pearson_correlation, 2)
             ),
             hjust = 0,
             vjust = 1.1,
             size  = 6,
             fontface = "bold") +
    
    scale_fill_manual(
      values = c(
        "Non-Coding, Mitochondrial and Ribosomal Genes" = "#48CAE4",
        "Other Genes" = "#B2DF8A"
      ),
      name = "Gene Group",
      guide = guide_legend(
        override.aes = list(
          shape = 21,
          color = "black",
          stroke = 0.5,
          linetype = 0,
          size = 5
        )
      )
    ) +
    
    labs(
      title = paste("Linear Regression of Integrated log2FC vs. Median Log2FC - ", title_suffix),
      x = "Median Log2FC",
      y = "Integrated log2FC"
    ) +
    
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title = element_text(size = 20),
      axis.text  = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position   = "bottom"
    )
  
  # Return plot and slope
  return_list <- list(
    "linear_regression_plot"  = p,
    "linear_regression_slope" = linear_regression_slope,
    "effect_size_correlation" = effect_size_correlation$effect_size_pearson_correlation
  )
  
  return(return_list)
}



# get_gsea and helper functions--------------------------------------------------------------

get_gsea <- function(meta_table,
                     integrated_list = NULL,  # now optional
                     cell_type = "NKs",
                     condition,
                     dataset_name, 
                     sample_name,
                     filter = TRUE) {
  
  BPPARAM <- BiocParallel::MulticoreParam(workers = 6)   # or SnowParam
  BiocParallel::register(BPPARAM)
  
  # ---- helper flag -----------------------------------------------------------
  has_int <- !is.null(integrated_list)
  
  # =========================================================================== #
  # 1 · META GSEA -------------------------------------------------------------- #
  # =========================================================================== #
  ## a) build vectors ----------------------------------------------------------
  df_meta           <- data.frame(t(meta_table$median_log2FC))
  rownames(df_meta) <- "median_log2FC"
  colnames(df_meta) <- rownames(meta_table)
  
  gsea_vec_meta_sym <- make_symbol_vector(df_meta,  "median_log2FC")
  gsea_vec_meta_ent <- make_entrez_vector(df_meta, "median_log2FC")
  
  ## b) run enrichment (pvalueCutoff=1 disables internal filtering) ----------- #
  meta_reactome <- ReactomePA::gsePathway(
    gsea_vec_meta_ent, organism = "human", pvalueCutoff = 1, BPPARAM = BPPARAM
  )
  meta_kegg     <- clusterProfiler::gseKEGG(
    gsea_vec_meta_ent, pvalueCutoff = 1, BPPARAM = BPPARAM
  )
  hallmark_df <- msigdbr::msigdbr(
    species  = "Homo sapiens",
    category = "H"
  ) |>
    dplyr::select(gs_name, gene_symbol)
  
  meta_h <- clusterProfiler::GSEA(
    gsea_vec_meta_sym, TERM2GENE = hallmark_df, pvalueCutoff = 1, BPPARAM = BPPARAM
  )
  
  meta_gobp <- clusterProfiler::gseGO(
    gsea_vec_meta_ent, OrgDb = org.Hs.eg.db, ont = "BP",
    keyType = "ENTREZID", pvalueCutoff = 1, BPPARAM = BPPARAM
  )
  meta_gocc <- clusterProfiler::gseGO(
    gsea_vec_meta_ent, OrgDb = org.Hs.eg.db, ont = "CC",
    keyType = "ENTREZID", pvalueCutoff = 1, BPPARAM = BPPARAM
  )
  
  biocarta_df <- msigdbr::msigdbr(
    species     = "Homo sapiens",
    category    = "C2",
    subcategory = "CP:BIOCARTA"
  ) |>
    dplyr::select(gs_name, entrez_gene)
  meta_biocarta <- clusterProfiler::GSEA(
    gsea_vec_meta_ent, TERM2GENE = biocarta_df, pvalueCutoff = 1, BPPARAM = BPPARAM
  )
  
  ## keep originals for bar-plots --------------------------------------------- #
  meta_reactome_original <- meta_reactome
  meta_kegg_original     <- meta_kegg
  meta_h_original        <- meta_h
  meta_gobp_original     <- meta_gobp
  meta_gocc_original     <- meta_gocc
  meta_biocarta_original <- meta_biocarta
  
  ## optional post-filter for plotting ---------------------------------------- #
  if (filter) {
    meta_reactome@result <- dplyr::filter(meta_reactome@result, p.adjust <= 0.1)
    meta_kegg@result     <- dplyr::filter(meta_kegg@result,     p.adjust <= 0.1)
    meta_h@result        <- dplyr::filter(meta_h@result,        p.adjust <= 0.1)
    meta_gobp@result     <- dplyr::filter(meta_gobp@result,     p.adjust <= 0.1)
    meta_gocc@result     <- dplyr::filter(meta_gocc@result,     p.adjust <= 0.1)
    meta_biocarta@result <- dplyr::filter(meta_biocarta@result, p.adjust <= 0.1)
  }
  
  # =========================================================================== #
  # 2 · META PLOTS (unchanged) ------------------------------------------------ #
  # =========================================================================== #
  title_suffix <- paste(cell_type, "meta-analysis")
  
  p_meta_reactome <- plot_top10(meta_reactome,
                                paste("REACTOME GSEA for Meta Method:", title_suffix),
                                xcol = "Description", color1 = "firebrick")
  p_meta_kegg <- plot_top10(meta_kegg,
                            paste("KEGG GSEA for Meta Method:", title_suffix),
                            xcol = "Description", color1 = "#03045e")
  
  top10_meta_h <- meta_h@result %>%
    arrange(desc(enrichmentScore)) %>%
    slice_head(n = 10) %>%
    mutate(ID_clean = ID %>%
             str_remove("^HALLMARK_") %>%
             str_replace_all("_", " ") %>%
             str_to_title())
  
  p_meta_h <- ggdotchart(
    data = top10_meta_h,
    x = "ID_clean", y = "enrichmentScore",
    title = paste("MSigDB Hallmark GSEA for Meta Method:", title_suffix),
    color = "p.adjust",
    sorting = "descending",
    add = "segments",
    rotate = TRUE,
    dot.size = 10,
    gradient.cols = c("#fb8500", "grey"),
    ggtheme = theme_pubr()
  ) +
    theme(
      legend.key.height = unit(1.5, "cm"),
      legend.key.width  = unit(1, "cm"),
      legend.text  = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.margin  = margin(1, 1, 1, 1, "cm")
    )
  
  top10_meta_biocarta <- meta_biocarta@result %>%
    arrange(desc(enrichmentScore)) %>%
    slice_head(n = 10) %>%
    mutate(
      raw   = ID %>% str_remove("^BIOCARTA_") %>% str_replace_all("_", " "),
      first = word(raw, 1),
      rest  = str_remove(raw, paste0("^", first, " ")) %>% str_to_title(),
      ID_clean = if_else(rest == "", toupper(first), paste(toupper(first), rest))
    ) %>%
    dplyr::select(-raw, -first, -rest)
  
  p_meta_biocarta <- ggdotchart(
    data = top10_meta_biocarta,
    x = "ID_clean", y = "enrichmentScore",
    title = paste("MSigDB Biocarta GSEA for Meta Method:", title_suffix),
    color = "p.adjust",
    sorting = "descending",
    add = "segments",
    rotate = TRUE,
    dot.size = 10,
    gradient.cols = c("#2C6E49", "grey"),
    ggtheme = theme_pubr()
  ) +
    theme(
      legend.key.height = unit(1.5, "cm"),
      legend.key.width  = unit(1, "cm"),
      legend.text  = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.margin  = margin(1, 1, 1, 1, "cm")
    )
  
  p_meta_gobp <- plot_top10(meta_gobp,
                            paste("GO:BP GSEA for Meta Method:", title_suffix),
                            xcol = "Description", color1 = "#7C4585")
  p_meta_gocc <- plot_top10(meta_gocc,
                            paste("GO:CC GSEA for Meta Method:", title_suffix),
                            xcol = "Description", color1 = "#000000")
  
  # =========================================================================== #
  # 3 · INTEGRATED SECTION (executed only if supplied) ------------------------ #
  # =========================================================================== #
  if (has_int) {
    conds <- unique(integrated_list$integrated_seu@meta.data$condition)
    if (length(conds) != 2) stop("integrated_seu must contain exactly two conditions")
    
    # GSEA for the "integrated" method
    {
      df_int <- data.frame(t(integrated_list$de_table$avg_log2FC))
      rownames(df_int) <- "log2FC"
      colnames(df_int) <- rownames(integrated_list$de_table)
      
      gsea_vec_int_sym <- make_symbol_vector(df_int,  "log2FC")
      gsea_vec_int_ent <- make_entrez_vector(df_int, "log2FC")
      
      # Reactome / KEGG / Hallmark / GO / BioCarta with pvalueCutoff=1
      int_reactome <- clusterProfiler::gsePathway(
        geneList = gsea_vec_int_ent, organism = "human", pvalueCutoff = 1, BPPARAM = BPPARAM
      )
      int_reactome_original <- int_reactome
      if (filter) int_reactome@result <- dplyr::filter(int_reactome@result, p.adjust <= 0.1)
      
      int_kegg <- clusterProfiler::gseKEGG(
        geneList = gsea_vec_int_ent, pvalueCutoff = 1, BPPARAM = BPPARAM
      )
      int_kegg_original <- int_kegg
      if (filter) int_kegg@result <- dplyr::filter(int_kegg@result, p.adjust <= 0.1)
      
      int_h <- clusterProfiler::GSEA(
        geneList = gsea_vec_meta_sym, TERM2GENE = hallmark_df, pvalueCutoff = 1, BPPARAM = BPPARAM
      )
      int_h_original <- int_h
      if (filter) int_h@result <- dplyr::filter(int_h@result, p.adjust <= 0.1)
      
      int_gobp <- clusterProfiler::gseGO(
        geneList = gsea_vec_int_ent, OrgDb = org.Hs.eg.db, ont = "BP",
        keyType = "ENTREZID", pvalueCutoff = 1, BPPARAM = BPPARAM
      )
      int_gobp_original <- int_gobp
      if (filter) int_gobp@result <- dplyr::filter(int_gobp@result, p.adjust <= 0.1)
      
      int_gocc <- clusterProfiler::gseGO(
        geneList = gsea_vec_int_ent, OrgDb = org.Hs.eg.db, ont = "CC",
        keyType = "ENTREZID", pvalueCutoff = 1, BPPARAM = BPPARAM
      )
      int_gocc_original <- int_gocc
      if (filter) int_gocc@result <- dplyr::filter(int_gocc@result, p.adjust <= 0.1)
      
      int_biocarta <- clusterProfiler::GSEA(
        geneList = gsea_vec_int_ent, TERM2GENE = biocarta_df, pvalueCutoff = 1, BPPARAM = BPPARAM
      )
      int_biocarta_original <- int_biocarta
      if (filter) int_biocarta@result <- dplyr::filter(int_biocarta@result, p.adjust <= 0.1)
    }
    
    # plots for "integrated" method (unchanged)
    {
      title_suffix <- paste(cell_type, "Integrated Method")
      p_int_reactome <- plot_top10(int_reactome,
                                   paste("REACTOME GSEA for Integrated Method:", title_suffix),
                                   xcol = "Description", color1 = "#FF8FA3")
      p_int_kegg <- plot_top10(int_kegg,
                               paste("KEGG GSEA for Integrated Method:", title_suffix),
                               xcol = "Description", color1 = "#48CAE4")
      top10_int_h <- int_h@result %>%
        arrange(desc(enrichmentScore)) %>%
        slice_head(n = 10) %>%
        mutate(ID_clean = ID %>%
                 str_remove("^HALLMARK_") %>%
                 str_replace_all("_", " ") %>%
                 str_to_title())
      p_int_h <- ggdotchart(
        data = top10_int_h,
        x = "ID_clean", y = "enrichmentScore",
        title = paste("MSigDB Hallmark GSEA for Integrated Method:", title_suffix),
        color = "p.adjust",
        sorting = "descending",
        add = "segments",
        rotate = TRUE,
        dot.size = 10,
        gradient.cols = c("#ffb703", "grey"),
        ggtheme = theme_pubr()
      ) +
        theme(
          legend.key.height = unit(1.5, "cm"),
          legend.key.width  = unit(1, "cm"),
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 10),
          plot.margin  = margin(1, 1, 1, 1, "cm")
        )
      top10_int_biocarta <- int_biocarta@result %>%
        arrange(desc(enrichmentScore)) %>%
        slice_head(n = 10) %>%
        mutate(
          raw   = ID %>% str_remove("^BIOCARTA_") %>% str_replace_all("_", " "),
          first = word(raw, 1),
          rest  = str_remove(raw, paste0("^", first, " ")) %>% str_to_title(),
          ID_clean = if_else(rest == "", toupper(first), paste(toupper(first), rest))
        ) %>%
        dplyr::select(-raw, -first, -rest)
      p_int_biocarta <- ggdotchart(
        data = top10_int_biocarta,
        x = "ID_clean", y = "enrichmentScore",
        title = paste("MSigDB Biocarta GSEA for Integrated Method:", title_suffix),
        color = "p.adjust",
        sorting = "descending",
        add = "segments",
        rotate = TRUE,
        dot.size = 10,
        gradient.cols = c("#B2DF8A", "grey"),
        ggtheme = theme_pubr()
      ) +
        theme(
          legend.key.height = unit(1.5, "cm"),
          legend.key.width  = unit(1, "cm"),
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 10),
          plot.margin  = margin(1, 1, 1, 1, "cm")
        )
      p_int_gobp <- plot_top10(int_gobp,
                               paste("GO:BP GSEA for Integrated Method:", title_suffix),
                               xcol = "Description", color1 = "#E9A5F1")
      p_int_gocc <- plot_top10(int_gocc,
                               paste("GO:CC GSEA for Integrated Method:", title_suffix),
                               xcol = "Description", color1 = "#A76545")
    }
    
    # Create p_bar_plot_percentage (unchanged data source = *original* results)
    {
      meta_reactome <- meta_reactome_original@result$ID
      int_reactome  <- int_reactome_original@result$ID
      meta_kegg     <- meta_kegg_original@result$ID
      int_kegg      <- int_kegg_original@result$ID
      meta_h        <- meta_h_original@result$ID
      int_h         <- int_h_original@result$ID
      meta_biocarta <- meta_biocarta_original@result$ID
      int_biocarta  <- int_biocarta_original@result$ID
      meta_gobp     <- meta_gobp_original@result$ID
      int_gobp      <- int_gobp_original@result$ID
      meta_gocc     <- meta_gocc_original@result$ID
      int_gocc      <- int_gocc_original@result$ID
      
      df_reactome <- calc_fractions(meta_reactome, int_reactome, "REACTOME")
      df_kegg     <- calc_fractions(meta_kegg,     int_kegg,     "KEGG")
      df_h        <- calc_fractions(meta_h,        int_h,        "HALLMARK")
      df_biocarta <- calc_fractions(meta_biocarta, int_biocarta, "BIOCARTA")
      df_gobp     <- calc_fractions(meta_gobp,     int_gobp,     "GO:BP")
      df_gocc     <- calc_fractions(meta_gocc,     int_gocc,     "GO:CC")
      
      plot_df <- dplyr::bind_rows(df_reactome, df_kegg, df_h, df_biocarta, df_gobp, df_gocc)
      plot_df$group <- factor(plot_df$group, levels = c("Meta-exclusive", "Shared", "Int-exclusive"))
      
      p_bar_plot_percentage <- ggplot(plot_df, aes(x = dataset, y = fraction, fill = group)) +
        geom_bar(stat = "identity", position = "stack") +
        geom_text(
          aes(label = ifelse(fraction == 0, "", scales::percent(fraction, accuracy = 1))),
          position = position_stack(vjust = 0.5),
          color = "white", size = 5
        ) +
        scale_fill_manual(values = c(
          "Meta-exclusive" = "#03045e",
          "Shared"         = "grey70",
          "Int-exclusive"  = "#fb8500"
        )) +
        labs(
          title = "Comparison of Meta vs. Integrated pathways",
          x = "Pathway type",
          y = "Fraction of union of pathways"
        ) +
        theme_minimal(base_size = 14)
    }
  }
  
  # =========================================================================== #
  # 4 · PATCHWORK / PDF ------------------------------------------------------- #
  # =========================================================================== #
  if (has_int) {
    p_pathways_h        <- p_meta_h        / p_int_h
    p_pathways_reactome <- p_meta_reactome / p_int_reactome
    p_pathways_kegg     <- p_meta_kegg     / p_int_kegg
    p_pathways_biocarta <- p_meta_biocarta / p_int_biocarta
    p_pathways_gobp     <- p_meta_gobp     / p_int_gobp
    p_pathways_gocc     <- p_meta_gocc     / p_int_gocc
  } else {
    p_pathways_h        <- p_meta_h 
    p_pathways_reactome <- p_meta_reactome 
    p_pathways_kegg     <- p_meta_kegg    
    p_pathways_biocarta <- p_meta_biocarta 
    p_pathways_gobp     <- p_meta_gobp   
    p_pathways_gocc     <- p_meta_gocc  
  }
  
  #  PDF
  {
    dir_name <- paste0("/bigdata/bgutman/Collabs/Dror mevorach/plots/", cell_type)
    if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)
    
    dir_name <- paste0("/bigdata/bgutman/Collabs/Dror mevorach/plots/", cell_type, "/", condition, "_", dataset_name, "_", sample_name)
    if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)
    
    full_pdf_name <- paste0(dir_name, "/", condition, "_", dataset_name, "_", sample_name, "_gsea_plots.pdf")
    if (has_int) {
      pdf(full_pdf_name, width = 26, height = 20)
      print(p_pathways_h);        print(p_pathways_reactome)
      print(p_pathways_kegg);     print(p_pathways_biocarta)
      print(p_pathways_gobp);     print(p_pathways_gocc)
      print(p_bar_plot_percentage)
      dev.off()
    } else {
      pdf(full_pdf_name, width = 26, height = 20)
      print(p_pathways_h);        print(p_pathways_reactome)
      print(p_pathways_kegg);     print(p_pathways_biocarta)
      print(p_pathways_gobp);     print(p_pathways_gocc)
      dev.off()
    }
  }
  
  # =========================================================================== #
  # 5 · RETURN LIST ----------------------------------------------------------- #
  # =========================================================================== #
  return(c(
    list(
      meta_reactome = meta_reactome_original,
      meta_kegg     = meta_kegg_original,
      meta_h        = meta_h_original,
      meta_gobp     = meta_gobp_original,
      meta_gocc     = meta_gocc_original,
      meta_biocarta = meta_biocarta_original
    ),
    if (has_int) list(
      int_reactome  = int_reactome_original,
      int_kegg      = int_kegg_original,
      int_h         = int_h_original,
      int_gobp      = int_gobp_original,
      int_gocc      = int_gocc_original,
      int_biocarta  = int_biocarta_original
    )
  ))
}



# Helper 1: Make a symbol-based GSEA vector (for Hallmark sets)
make_symbol_vector <- function(df, row_name) {
  
  # Description:
  #   Converts a one‐row data.frame of gene‐rank values into a named numeric vector
  #   keyed by gene symbol, sorted descending.
  # Arguments:
  #   df       – data.frame with one row of numeric values and gene symbols as columns.
  #   row_name – character; name of the row to extract.
  # Returns:
  #   Named numeric vector (names = gene symbols, values = metrics).
  
  # df columns = gene symbols, row_name = "median_log2FC" or "log2FC"
  gvec <- as.numeric(df[row_name, ])
  names(gvec) <- colnames(df)  # gene symbols
  sort(gvec, decreasing = TRUE)
}

# Helper 2: Make an ENTREZ-based GSEA vector (for KEGG / Reactome)
make_entrez_vector <- function(df, row_name) {
  
  # Description:
  #   Maps gene symbols to ENTREZ IDs, then converts a one‐row data.frame into a named numeric vector
  #   keyed by ENTREZ ID, sorted descending.
  # Arguments:
  #   df       – data.frame with one row of numeric values and gene symbols as columns.
  #   row_name – character; name of the row to extract.
  # Returns:
  #   Named numeric vector (names = ENTREZ IDs, values = metrics).
  
  map <- mapIds(org.Hs.eg.db, keys = colnames(df),
                keytype = "SYMBOL", column = "ENTREZID") %>% na.omit()
  df <- df[, names(map), drop = FALSE]
  colnames(df) <- map  # rename to ENTREZ
  gvec <- as.numeric(df[row_name, ])
  names(gvec) <- colnames(df)
  sort(gvec, decreasing = TRUE)
}

# Helper 3: Convert Hallmark GeneSetCollection -> data.frame for GSEA
msigCollectionToDf <- function(gsc) {
  
  # Description:
  #   Flattens an MSigDB GeneSetCollection into a data.frame suitable for clusterProfiler,
  #   with TERM2GENE format (columns term, gene).
  # Arguments:
  #   gsc – GeneSetCollection object (e.g. from msigdb.hs subset).
  # Returns:
  #   data.frame with columns term (pathway ID/name) and gene (symbol).
  
  do.call(rbind, lapply(gsc, function(gs) {
    data.frame(term = setName(gs), gene = geneIds(gs), stringsAsFactors = FALSE)
  }))
}

# Helper 4: Dotchart of top 10 pathways
plot_top10 <- function(gsea_res, title, xcol = "Description",
                       ycol = "enrichmentScore", color1 = "grey40") {
  
  # Description:
  #   Creates a ggdotchart of the top 10 pathways by enrichment score or gene ratio
  #   from a clusterProfiler GSEA/enrichResult object.
  # Arguments:
  #   gsea_res – GSEA or enrichResult object with @result data.frame.
  #   title    – character; plot title.
  #   xcol     – character; column name for x‐axis (e.g. "Description").
  #   ycol     – character; column name for y‐axis (e.g. "enrichmentScore" or "GeneRatio").
  #   color1   – character; primary color for gradient (second is grey).
  # Returns:
  #   ggplot object of the dot chart.
  
  
  gsea_res@result[[ycol]] <- sapply(gsea_res@result[[ycol]], function(s) eval(parse(text = s)))
  
  
  top10 <- gsea_res@result %>%
    arrange(desc(.data[[ycol]])) %>%
    slice_head(n = 10)
  
  p <- ggdotchart(
    data = top10,
    x = xcol, y = ycol,
    title = title,
    color = "p.adjust",
    sorting = "descending",
    add = "segments",
    rotate = TRUE,
    dot.size = 10,
    gradient.cols = c(color1, "grey"),
    ggtheme = theme_pubr(),
    ylab = FALSE
  ) +
    theme(
      legend.key.height = unit(1.5, "cm"),
      legend.key.width  = unit(1, "cm"),
      legend.text  = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.margin  = margin(1, 1, 1, 1, "cm")
    )
  
  return(p)
}

# Helper 5: A helper function to compute fraction of union (for the barplot)

calc_fractions <- function(meta_ids, int_ids, dataset_label) {
  

  # Description:
  #   Computes fractions of Meta‐exclusive, Shared, and Int‐exclusive pathways
  #   given two sets of pathway IDs.
  # Arguments:
  #   meta_ids      – character vector of meta method pathway IDs.
  #   int_ids       – character vector of integrated method pathway IDs.
  #   dataset_label – character; label for the dataset (e.g. "REACTOME", "KEGG", "HALLMARK").
  # Returns:
  #   data.frame with columns dataset, group (Meta‑exclusive/Shared/Int‑exclusive), and fraction.
  
  if (length(meta_ids)==0 && length(int_ids)==0) {
    return(data.frame(
      dataset  = dataset_label,
      group    = c("Meta-exclusive","Shared","Int-exclusive"),
      fraction = c(0, 0, 0)
    ))
  }
  
  
  union_set     <- union(meta_ids, int_ids)
  intersect_set <- intersect(meta_ids, int_ids)
  meta_only     <- setdiff(meta_ids, int_ids)
  int_only      <- setdiff(int_ids, meta_ids)
  
  total         <- length(union_set)
  frac_meta     <- length(meta_only) / total
  frac_shared   <- length(intersect_set) / total
  frac_int      <- length(int_only) / total
  
  data.frame(
    dataset  = dataset_label,
    group    = c("Meta-exclusive", "Shared", "Int-exclusive"),
    fraction = c(frac_meta, frac_shared, frac_int)
  )
}


# get_ora and helper functions--------------------------------------------------------------
get_ora <- function(meta_table_list,
                    log2FC_threshold_for_gene = 1,
                    positive = TRUE,
                    q_value_threshold_for_gene = 0.05,
                    cols = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
                    BPPARAM = BiocParallel::MulticoreParam(workers = 6),
                    q_value_threshold_for_pathway = 0.05) {
  
  ## ── sanity checks (unchanged) ─────────────────────────────────────────────
  {
    expected_names <- c("NK_meta_table",
                        "monocytes_meta_table",
                        "B_cells_meta_table",
                        "CD4_T_cells_meta_table",
                        "CD8_T_cells_meta_table")
    if (!is.list(meta_table_list) || is.null(names(meta_table_list))) {
      stop("meta_table_list must be a *named* list holding the five cell-type meta tables.", call. = FALSE)
    }
    missing <- setdiff(expected_names, names(meta_table_list))
    extra   <- setdiff(names(meta_table_list), expected_names)
    if (length(missing) || length(extra)) {
      stop(paste0(
        "meta_table_list must contain exactly:\n",
        paste0(" • ", expected_names, collapse = "\n"), "\n\n",
        if (length(missing)) paste0("Missing: ", paste(missing, collapse = ", "), "\n") else "",
        if (length(extra))   paste0("Extra:   ", paste(extra,   collapse = ", "), "\n") else ""
      ), call. = FALSE)
    }
    if (!identical(names(meta_table_list), expected_names)) {
      stop("Elements present but not in required order: ",
           paste(expected_names, collapse = ", "),
           call. = FALSE)
    }
  }
  
  ## ── 1 · parallel filter per meta table ────────────────────────────────────
  meta_table_list_filtered <- BiocParallel::bplapply(
    meta_table_list,
    function(df) {
      if (positive) {
        dplyr::filter(df,
                      median_log2FC  >  log2FC_threshold_for_gene,
                      q_val_fisher   <= q_value_threshold_for_gene)
      } else {
        dplyr::filter(df,
                      median_log2FC  < -log2FC_threshold_for_gene,
                      q_val_fisher   <= q_value_threshold_for_gene)
      }
    },
    BPPARAM = BPPARAM
  )
  
  ## ── 2 · presence matrix ──────────────────────────────────────────────────
  gene_upset_df <- build_presence_df(meta_table_list_filtered, id_col = "gene")
  
  ## foreground
  {
    if (missing(cols) || !length(cols))
      stop("cols must be a non-empty character vector.", call. = FALSE)
    not_found <- setdiff(cols, names(gene_upset_df))
    if (length(not_found))
      stop("Missing columns in gene_upset_df: ", paste(not_found, collapse = ", "), call. = FALSE)
    
    gene_ORA_df <- dplyr::filter(gene_upset_df, dplyr::if_all(dplyr::all_of(cols), \(x) x))
    rownames(gene_ORA_df) <- gene_ORA_df$gene
    ora_vec_meta_sym <- unique(gene_ORA_df$gene)
    
    ora_vec_meta_ent <- {
      tmp <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys      = ora_vec_meta_sym,
        keytype   = "SYMBOL",
        column    = "ENTREZID",
        multiVals = "first"
      )
      unique(unname(stats::na.omit(tmp)))
    }
  }
  
  ## universe
  {
    genes_in_celltype <- function(meta_tbl) {
      meta_tbl %>%
        dplyr::filter(comparisons_num > 0) %>%
        dplyr::pull(gene) %>%
        unique()
    }
    gene_sets <- purrr::map(meta_table_list, genes_in_celltype)
    universe_genes_sym <- sort(purrr::reduce(gene_sets, intersect))
    
    universe_genes_entrez <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = universe_genes_sym,
      column    = "ENTREZID",
      keytype   = "SYMBOL",
      multiVals = "first"
    ) |>
      stats::na.omit() |>
      unname() |>
      unique() |>
      sort()
  }
  
  ## ── 3 · Run enrichments in parallel ──────────────────────────────────────
  enrichment_tasks <- list(
    reactome = list(
      fun = function() {
        x <- ReactomePA::enrichPathway(gene = ora_vec_meta_ent,
                                       organism = "human",
                                       universe = universe_genes_entrez)
        x@result <- dplyr::filter(x@result, p.adjust <= q_value_threshold_for_pathway)
        x
      }
    ),
    kegg = list(
      fun = function() {
        x <- clusterProfiler::enrichKEGG(
          gene     = ora_vec_meta_ent,
          organism = "hsa",            # ← FIX
          universe = universe_genes_entrez
        )
        if (is.null(x) || nrow(x@result) == 0) return(x)   # graceful
        x@result <- dplyr::filter(x@result,
                                  p.adjust <= q_value_threshold_for_pathway)
        x
      }
    ),
    hallmark = list(
      fun = function() {
        hallmark_df <- msigdbr::msigdbr(
          species  = "Homo sapiens",
          category = "H"
        ) %>% dplyr::select(gs_name, gene_symbol)
        x <- clusterProfiler::enricher(gene = ora_vec_meta_sym,
                                       TERM2GENE = hallmark_df,
                                       universe  = universe_genes_sym)
        x@result <- dplyr::filter(x@result, p.adjust <= q_value_threshold_for_pathway)
        x
      }
    ),
    gobp = list(
      fun = function() {
        x <- clusterProfiler::enrichGO(
          gene         = ora_vec_meta_ent,
          OrgDb      = org.Hs.eg.db,
          ont          = "BP",
          keyType      = "ENTREZID"
          # pvalueCutoff = 0.05
        )
        x@result <- dplyr::filter(x@result, p.adjust <= q_value_threshold_for_pathway)
        x
      }
    ),
    gocc = list(
      fun = function() {
        x <- clusterProfiler::enrichGO(
          gene         = ora_vec_meta_ent,
          OrgDb      = org.Hs.eg.db,
          ont          = "CC",
          keyType      = "ENTREZID"
          # pvalueCutoff = 0.05
        )
        x@result <- dplyr::filter(x@result, p.adjust <=q_value_threshold_for_pathway)
        x
      }
    ),
    biocarta = list(
      fun = function() {
        biocarta_df <- msigdbr::msigdbr(
          species     = "Homo sapiens",
          category    = "C2",
          subcategory = "CP:BIOCARTA"
        ) %>% dplyr::select(gs_name, entrez_gene)
        x <- clusterProfiler::enricher(gene = ora_vec_meta_ent,
                                       TERM2GENE = biocarta_df,
                                       universe  = universe_genes_entrez)
        x@result <- dplyr::filter(x@result, p.adjust <= q_value_threshold_for_pathway)
        x
      }
    )
  )
  
  safe_run <- function(f) {
    res <- f()
    if (is.null(res)) return(res)
    if (!methods::is(res, "enrichResult") && !methods::is(res, "gseaResult")) return(res)
    if (!"p.adjust" %in% names(res@result)) return(res)
    res
  }
  
  enrichment_results <- BiocParallel::bplapply(
    enrichment_tasks,
    function(task) safe_run(task$fun),  # pass the function itself
    BPPARAM = BPPARAM
  )
  
  meta_reactome <- enrichment_results$reactome
  meta_kegg     <- enrichment_results$kegg
  meta_h        <- enrichment_results$hallmark
  meta_gobp     <- enrichment_results$gobp
  meta_gocc     <- enrichment_results$gocc
  meta_biocarta <- enrichment_results$biocarta
  
  ## ── 4 · return list (unchanged structure) ─────────────────────────────────
  return(list(
    meta_reactome  = meta_reactome,
    meta_kegg      = meta_kegg,
    meta_h         = meta_h,
    meta_gobp      = meta_gobp,
    meta_gocc      = meta_gocc,
    meta_biocarta  = meta_biocarta
  ))
}


get_ora_scatter_plot <- function(ora_list, join_column = "ID", meta_main_rank_column = "GeneRatio", int_main_rank_column ="GeneRatio",
                                 meta_effect_size_column = "GeneRatio", int_effect_size_column = "GeneRatio", effect_size_threshold = 0.04) {

    # Description:
  #   Creates a scatter plot comparing GeneRatio from meta vs integrated ORA results
  #   for Reactome, KEGG, and Hallmark, with linear fit, Spearman & Pearson annotations,
  #   and selective pathway labels.
  # Arguments:
  #   ora_list                – list of enrichResult objects (meta/int for reactome, kegg, h).
  #   join_column             – character; column name to merge on (default "ID").
  #   meta_main_rank_column   – character; column in meta result for ranking (default "GeneRatio").
  #   int_main_rank_column    – character; column in int result for ranking (default "GeneRatio").
  #   meta_effect_size_column – character; column for effect size in meta (default "GeneRatio").
  #   int_effect_size_column  – character; column for effect size in int (default "GeneRatio").
  #   effect_size_threshold   – numeric; min GeneRatio to label points (default 0.04).
  # Returns:
  #   ggplot object of the ORA scatter comparison.
  
  # 1. Define your “collections” and their pretty labels
  collections <- c("reactome", "kegg", "h", "gobp", "gocc", "biocarta")
  labels      <- c("REACTOME", "KEGG", "H", "GO:BP", "GO:CC", "Biocarta")
  
  # 2. For each collection, pull out meta/int, compute both correlations
  results_tbl <- map2_dfr(
    .x = collections,
    .y = labels,
    .f = ~ {
      meta_df <- ora_list[[paste0("meta_", .x)]]@result
      int_df  <- ora_list[[paste0("int_",  .x)]]@result
      
      spearman_val <- get_spearman(
        meta_table            = meta_df,
        de_table              = int_df,
        join_column           = join_column,
        meta_main_rank_column = meta_main_rank_column,
        int_main_rank_column  = int_main_rank_column
      )$spearman_correlation
      
      pearson_val <- get_effect_size_correlation(
        meta_table                = meta_df,
        de_table                  = int_df,
        join_column               = join_column,
        meta_effect_size_column   = meta_effect_size_column,
        int_effect_size_column    = int_effect_size_column
      )$effect_size_pearson_correlation
      
      # return one row per collection
      tibble(
        collection = .x,
        label      = .y,
        spearman   = spearman_val,
        pearson    = pearson_val
      )
    }
  )
  
  # 3. Build your annotation strings in one go
  spearman_annotation <- sprintf(
    "Spearman's Rank Correlation (%s): %s",
    paste(results_tbl$label, collapse = ", "),
    paste(round(results_tbl$spearman,  2), collapse = ", ")
  )
  
  pearson_annotation <- sprintf(
    "Effect Size Correlation (%s): %s",
    paste(results_tbl$label, collapse = ", "),
    paste(round(results_tbl$pearson,  2), collapse = ", ")
  )
  
  
  
  # create ora_merged
  {
    #create merged_reactome, merged_kegg and merged_h
    
    for (ds in c("reactome", "kegg", "h", "gobp", "gocc", "biocarta")) {
      # Dynamically get the dataframes by name
      meta_df <- ora_list[[ paste0("meta_", ds) ]]@result
      int_df  <- ora_list[[ paste0("int_",  ds) ]]@result
      

      # Merge, select columns, and rename within a single pipe
      merged_obj <- merge(meta_df, int_df, by = join_column, suffixes = c(".meta", ".int")) %>% 
        dplyr::select(ID, Description.meta, GeneRatio.meta, GeneRatio.int, p.adjust.meta, p.adjust.int) %>% 
        dplyr::rename(Description = Description.meta)
      
      # For each specified column, evaluate the text as an expression
      for (col_name in c("GeneRatio.meta", "GeneRatio.int", "p.adjust.meta", "p.adjust.int")) {
        merged_obj[[col_name]] <- sapply(merged_obj[[col_name]], function(s) eval(parse(text = s)))
      }
      
      # Store the result with a dynamic name like "merged_reactome", "merged_kegg", etc.
      assign(paste0("merged_", ds), merged_obj)
    }
    
    
    
    # Vector of dataset names
    datasets <- c("reactome", "kegg", "h", "gobp", "gocc", "biocarta")
    
    # Create an empty list to store data frames
    df_list <- list()
    
    # Loop over each dataset name
    for (ds in datasets) {
      meta_df <- ora_list[[paste0("meta_", ds)]]@result
      int_df  <- ora_list[[paste0("int_",  ds)]]@result
      
      merged_obj <- merge(
        meta_df, int_df,
        by       = join_column,          # e.g. "ID"
        suffixes = c(".meta", ".int")
      ) |>
        dplyr::select(
          {{ join_column }},              # keep whichever column you joined on
          Description.meta,
          GeneRatio.meta, GeneRatio.int,
          p.adjust.meta,  p.adjust.int
        ) |>
        dplyr::rename(Description = Description.meta) |>
        dplyr::mutate(                         # safe even when nrow == 0
          across(
            c(GeneRatio.meta, GeneRatio.int,
              p.adjust.meta,  p.adjust.int),
            ~ sapply(., \(s) eval(parse(text = s)))
          ),
          dataset = ds
        )
      
      df_list[[ds]] <- merged_obj          # store for later bind_rows()
    }
    
    
    #improve the descriptions of hallmark and biocarta
    df_list$h <-  df_list$h %>% 
      mutate(Description = Description %>%
               str_remove("^HALLMARK_") %>%
               str_replace_all("_", " ") %>%
               str_to_title()
      )
    
    df_list$biocarta <- df_list$biocarta %>% 
      mutate(Description = str_remove(Description, "^BIOCARTA_") %>% 
               str_replace_all("_", " ") %>% 
               { 
                 h <- str_to_upper(word(., 1, -2))
                 t <- str_to_sentence(word(., -1))
                 paste(h, t) 
               } )
    
    
    # Combine into a single data frame
    ora_merged <- dplyr::bind_rows(df_list)
    
    ora_merged <- ora_merged %>%
      mutate(
        q_value_status = case_when(
          p.adjust.meta < 0.05 & p.adjust.int < 0.05 ~ "Meta & Int Significant",
          p.adjust.meta < 0.05 & p.adjust.int >= 0.05 ~ "Meta Only Significant",
          p.adjust.meta >= 0.05 & p.adjust.int < 0.05 ~ "Int Only Significant",
          TRUE ~ "None Significant"
        )
      ) %>% 
      mutate(
        display_label = case_when(
          (GeneRatio.meta >= effect_size_threshold & q_value_status != "None Significant") | 
            (GeneRatio.int >= effect_size_threshold & q_value_status != "None Significant") ~ TRUE,
          TRUE ~ FALSE
        )
      )
    
    
  }
  
  
  
  
  
  #p_ora_scatter
  {
    # 1. Find the maximum GeneRatio values
    x_max <- max(ora_merged$GeneRatio.meta, na.rm = TRUE)
    y_max <- max(ora_merged$GeneRatio.int, na.rm = TRUE)
    
    # 2. Create breaks from 0 in steps of 0.02 up to each max
    x_breaks <- seq(0, x_max, by = 0.02)
    y_breaks <- seq(0, y_max, by = 0.02)
    
    # Fit linear model: avg_log2FC as a function of median_log2FC
    lm_model <- lm(GeneRatio.int ~ GeneRatio.meta, data = ora_merged)
    # Extract the slope (coefficient for median_log2FC)
    linear_regression_slope <- coef(lm_model)["GeneRatio.meta"]
    
    
    
    
    # 3. Create scatter plot
    p_ora_scatter <- ggplot(
      data = ora_merged,
      aes(
        x = GeneRatio.meta,
        y = GeneRatio.int,
        fill = q_value_status,    # Interior color mapped to q_value_status
        shape = dataset           # Distinct shapes for each dataset
      )
    ) +
      # Points: shape glyphs that allow fill, black outline
      geom_point(
        color = "black",   # outline color for shapes
        size = 3,
        stroke = 0.5,
        alpha = 0.9
      ) +
      
      # Horizontal & vertical grid lines
      geom_hline(yintercept = y_breaks, linetype = "solid", color = "gray85", size = 0.3) +
      geom_vline(xintercept = x_breaks, linetype = "solid", color = "gray85", size = 0.3) +
      
      # Axis limits & breaks
      scale_x_continuous(limits = c(0, x_max), breaks = x_breaks) +
      scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
      
      # Assign specific point shapes for each dataset (optional; customize as desired)
      scale_shape_manual(
        values = c(
          "h"        = 21,  # circle       – Hallmark
          "kegg"     = 22,  # square
          "reactome" = 23,  # diamond
          "gobp"     = 24,  # ▲ triangle-up
          "gocc"     = 24,  # ▲ triangle-up (same shape)
          "biocarta" = 25   # ▼ triangle-down
        ),
        labels = c(
          "h"        = "H",
          "kegg"     = "KEGG",
          "reactome" = "REACTOME",
          "gobp"     = "GO:BP",
          "gocc"     = "GO:CC",
          "biocarta" = "Biocarta"
        )
      ) +
      
      # Manual fill colors for each q_value_status
      scale_fill_manual(
        values = c(
          "Int Only Significant"   = "#B2DF8A",
          "Meta & Int Significant" = "firebrick",
          "Meta Only Significant"  = "#48CAE4",
          "None Significant"       = "#000000"
        )
      ) +
      
      # Plot labels
      labs(
        title  = "Shared ORA Pathways: GeneRatio (Meta) vs. GeneRatio (Int)",
        x      = "GeneRatio (Meta)",
        y      = "GeneRatio (Int)",
        fill   = "Q-Value Status",  # Legend title for fill
        shape  = "Dataset"
      ) +
      
      # Label only rows with display_label == TRUE
      geom_label_repel(
        data = subset(ora_merged, display_label == TRUE),
        aes(label = Description, fill = q_value_status),
        color = "white",
        fontface = "bold",
        size = 5,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "gray40",
        segment.size = 0.4,
        max.overlaps = Inf,
        show.legend = FALSE  # Hides legend for label fill
      ) +
      geom_smooth(
        data = ora_merged,
        aes(x = GeneRatio.meta, y = GeneRatio.int, group = 1),
        method = "lm",
        se = FALSE,
        color = "#000000",
        show.legend = FALSE 
      ) + 
      # Minimal theme with white background
      theme_minimal(base_size = 20) +
      theme(
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.title       = element_text(size = 20),
        axis.text        = element_text(size = 16),
        panel.grid.major = element_blank(),  # Remove default major grid
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        legend.position  = "bottom"
      ) +
      # Annotate slope
      annotate("text",
               x = -Inf, 
               y = Inf,
               label = paste(
                 " Linear Regression Slope:", round(linear_regression_slope, 2),
                 "\n", pearson_annotation,
                 "\n", spearman_annotation),
               hjust = 0,        # left-align each line to x = -Inf
               vjust = 1.1,      # move it slightly down from the top
               size = 6, 
               fontface = "bold") + 
      # >>> ADDED: override legend aesthetics so fill colors show up
      guides(
        fill = guide_legend(
          override.aes = list(shape = 21, color = "black", alpha = 1, stroke = 0.5)
        )
      )
  }
  
  return(p_ora_scatter)
  
}


# get_meta_int_comparison and helper functions--------------------------------------------------------------

get_meta_int_comparison_boxplot <- function(
    integrated_path, meta_path,
    files_to_exclude_int = NULL,
    files_to_exclude_meta = NULL,
    cell_type = "NKs",
    general_condition = "Neuroimmunulogical"
) {
  
  
  # Description:
  #   Invokes get_comparison_results, reshapes the metrics to long form,
  #   and builds a combined ggplot boxplot + jitter of Spearman, Pearson,
  #   and Jaccard metrics across conditions and datasets.
  # Arguments:
  #   integrated_path      – character; directory of integrated objects.
  #   meta_path            – character; directory of meta_table files.
  #   files_to_exclude_int – optional vector; files to skip in integrated_path.
  #   files_to_exclude_meta– optional vector; files to skip in meta_path.
  #   cell_type            – character; label for cell type (default "NKs").
  #   general_condition    – character; overarching condition label (default "Neuroimmunulogical").
  # Returns:
  #   List with elements:
  #     • comparison_results (data.frame as from get_comparison_results)
  #     • comparison_plot    (ggplot object)
  
  comparison_results <- tryCatch({
    get_comparison_results(
      integrated_path = integrated_path,
      meta_path       = meta_path,
      files_to_exclude_int  = files_to_exclude_int,
      files_to_exclude_meta = files_to_exclude_meta
    )
  }, error = function(e) {
    stop(
      "The function get_meta_int_comparison_plot failed because get_comparison_results encountered an error: ",
      e$message
    )
  })
  
  title <- paste(
    "Comparison of Integrated and Meta Approaches Across",
    general_condition, cell_type, "Samples"
  )
  
  # Convert wide data to long
  df_long <- reshape2::melt(
    comparison_results,
    id.vars = c("condition", "sample", "dataset"),
    measure.vars = c(
      "spearman", "pearson",
      "jaccard_top_20", "jaccard_bottom_20",
      "jaccard_top_50", "jaccard_bottom_50",
      "jaccard_top_100", "jaccard_bottom_100"
    ),
    variable.name = "metric",
    value.name    = "value"
  )
  
  df_long$condition <- factor(
    df_long$condition,
    levels = c(
      "CD59.Def", "CD55.Def", "MS",
      "Ab.mediated.IDD", "MOG", "CIS", "NMO"
    )
  )
  
  # Custom colors for condition
  my_condition_colors <- c(
    "Ab.mediated.IDD" = "#2C6E49",
    "CD55.Def"        = "#FF8FA3",
    "CD59.Def"        = "firebrick",
    "CIS"             = "#ffb703",
    "MOG"             = "#B2DF8A",
    "MS"              = "#48CAE4",
    "NMO"             = "#fb8500"  
  )
  
  # Use shapes 21–25 (filled shapes)
  my_dataset_shapes <- c(
    "Esaulova"   = 25,  # filled triangle down
    "Mevorach"   = 22,  # filled square
    "Ramesh"     = 21,  # filled circle
    "Schafflick" = 24,  # filled triangle up
    "Kang"       = 23   # filled diamond
  )
  
  # Rename x-axis labels
  metric_labels <- c(
    "spearman"          = "Spearman Correlation",
    "pearson"           = "Pearson Correlation",
    "jaccard_top_20"    = "Jaccard Top 20",
    "jaccard_bottom_20" = "Jaccard Bottom 20",
    "jaccard_top_50"    = "Jaccard Top 50",
    "jaccard_bottom_50" = "Jaccard Bottom 50",
    "jaccard_top_100"   = "Jaccard Top 100",
    "jaccard_bottom_100"= "Jaccard Bottom 100"
  )
  
  p <- ggplot(df_long, aes(x = metric, y = value)) +
    geom_boxplot(
      aes(group = metric),
      fill = NA,
      color = "black",
      outlier.shape = NA
    ) +
    geom_jitter(
      aes(shape = dataset, fill = condition),
      color = "black",
      position = position_jitter(width = 0.15),
      size = 4,
      alpha = 0.7
    ) +
    # Manually set fill colors AND explicitly set breaks + labels
    scale_fill_manual(
      name   = "Condition",
      values = my_condition_colors,
      breaks = c("CD59.Def", "CD55.Def", "MS", "Ab.mediated.IDD", "MOG", "CIS", "NMO"),
      labels = c("CD59 Def.", "CD55 Def.", "MS", "Ab Mediated IDD", "MOG", "CIS", "NMO"),
      guide = guide_legend(
        override.aes = list(shape = 21, color = "black", size = 4)
      )
    ) +
    scale_shape_manual(
      name   = "Dataset",
      values = my_dataset_shapes
    ) +
    scale_x_discrete(labels = metric_labels) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5, face = "bold")  # keep it bold
    ) +
    labs(
      title = title,
      x     = "Metric",
      y     = "Value"
    )
  
  
  return_list <- list(
    comparison_results = comparison_results,
    comparison_plot    = p
  )
  
  return(return_list)
}



get_patient_path_list <- function(integrated_path, meta_path, 
                                  files_to_exclude_int = NULL,
                                  files_to_exclude_meta = NULL) {
  
  # Description:
  #   Scans two directories of serialized objects (integrated Seurat files and meta_table files),
  #   matches files by base name, and returns a named list of paired file paths per patient.
  # Arguments:
  #   integrated_path      – character; directory containing “_integrated.qs” files.
  #   meta_path            – character; directory containing “_meta_table.qs” files.
  #   files_to_exclude_int – optional character vector; full paths to skip in integrated_path.
  #   files_to_exclude_meta– optional character vector; full paths to skip in meta_path.
  # Returns:
  #   Named list (one element per patient base name) of lists with elements:
  #     • meta_table        (character path)
  #     • integrated_object (character path)
  
  # Get all file names in the integrated and meta paths
  {
    # Integrated files
    integrated_file_names <- list.files(path = integrated_path, full.names = TRUE)
    integrated_file_names <- normalizePath(integrated_file_names, winslash = "/")
    
    if (!is.null(files_to_exclude_int)) {
      integrated_file_names <- integrated_file_names[ !integrated_file_names %in% files_to_exclude_int ]
    }
    
    # Meta files
    meta_file_names <- list.files(path = meta_path, full.names = TRUE)
    meta_file_names <- normalizePath(meta_file_names, winslash = "/")
    
    if (!is.null(files_to_exclude_meta)) {
      meta_file_names <- meta_file_names[ !meta_file_names %in% files_to_exclude_meta ]
    }
    
    # Debug: print the resulting paths
    print(integrated_file_names)
    print(meta_file_names)
  }
  
  # Create patient_path_list
  {
    # Function to strip off a suffix from the base name
    strip_suffix <- function(path, suffix) {
      sub(suffix, "", basename(path))
    }
    
    # Extract base names
    meta_bases <- strip_suffix(meta_file_names, "_meta_table\\.qs$")
    integrated_bases <- strip_suffix(integrated_file_names, "_integrated\\.qs$")
    
    # Name the vectors for easier lookup
    names(meta_file_names) <- meta_bases
    names(integrated_file_names) <- integrated_bases
    
    # Check if the base names are identical
    missing_in_meta <- setdiff(integrated_bases, meta_bases)
    missing_in_integrated <- setdiff(meta_bases, integrated_bases)
    
    if (length(missing_in_meta) > 0 || length(missing_in_integrated) > 0) {
      error_message <- "Mismatch between meta and integrated base names:"
      
      if (length(missing_in_meta) > 0) {
        error_message <- paste0(error_message, "\nMissing in meta: ", paste(missing_in_meta, collapse = ", "))
      }
      
      if (length(missing_in_integrated) > 0) {
        error_message <- paste0(error_message, "\nMissing in integrated: ", paste(missing_in_integrated, collapse = ", "))
      }
      
      stop(error_message)
    }
    
    # Find the common base identifiers (should be all at this point)
    common_bases <- intersect(meta_bases, integrated_bases)
    
    # Build the patient_path_list
    patient_path_list <- lapply(common_bases, function(b) {
      list(
        meta_table        = meta_file_names[[b]],
        integrated_object = integrated_file_names[[b]]
      )
    })
    
    # Name the top-level list by the common base identifier
    names(patient_path_list) <- common_bases
  }
  
  return(patient_path_list)
}




get_comparison_results <- function(integrated_path, meta_path, files_to_exclude_int = NULL,
                                   files_to_exclude_meta =NULL){
  
  # Description:
  #   For each patient pairing from two directories, reads meta and integrated DE results,
  #   computes Spearman, Pearson, and Jaccard indices (top/bottom 20,50,100), and
  #   returns a summary data.frame.
  # Arguments:
  #   integrated_path      – character; directory of integrated objects.
  #   meta_path            – character; directory of meta_table files.
  #   files_to_exclude_int – optional character vector of paths to exclude from integrated_path.
  #   files_to_exclude_meta– optional character vector of paths to exclude from meta_path.
  # Returns:
  #   Data.frame with one row per patient, columns:
  #     condition, sample, dataset,
  #     spearman, pearson,
  #     jaccard_top_20, jaccard_bottom_20,
  #     jaccard_top_50, jaccard_bottom_50,
  #     jaccard_top_100, jaccard_bottom_100
  
  patient_path_list <-  tryCatch({
    get_patient_path_list(integrated_path = integrated_path, meta_path = meta_path, files_to_exclude_int = files_to_exclude_int,
                          files_to_exclude_meta =files_to_exclude_meta)
  }, error = function(e) {
    stop("The function get_comparison_results failed because the function get_patient_path_list encountered an error: ", e$message)
  })
  


  # Initialize results df
  comparison_results <- data.frame(
    condition = character(),
    sample = character(),
    dataset = character(),
    spearman = numeric(),
    pearson = numeric(),
    jaccard_top_20 = numeric(),
    jaccard_bottom_20 = numeric(),
    jaccard_top_50 = numeric(),
    jaccard_bottom_50 = numeric(),
    jaccard_top_100 = numeric(),
    jaccard_bottom_100 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through patients
  for (patient_name in names(patient_path_list)) {
    
    # get condition, sample and dataset
    {  
      # Split by "_"
      parts <- unlist(strsplit(patient_name, "_"))
      
      # 1. Extract condition (always first part)
      condition <- parts[1]
      
      # 2. Handle the "Ramesh" special case
      if (parts[2] == "Ramesh") {
        dataset <- "Ramesh"
        # Skip parts[3] (the "5" or "3"), and take everything from parts[4] onward
        sample <- paste(parts[4:length(parts)], collapse = "_")
      } else {
        # Normal case
        dataset <- parts[2]
        sample <- paste(parts[3:length(parts)], collapse = "_")
      }
    }
    
    patient_data <- patient_path_list[[patient_name]]
    
    
    #get meta_table and de_table
    {  meta_table <- qread(patient_data$meta_table)
      integrated_object <- qread(patient_data$integrated_object)
      de_table <- integrated_object$de_table
      de_table$gene <- rownames(de_table)
    }
    
    #get spearman, pearson,jaccard
    
    spearman <- get_spearman(
      meta_table = meta_table, 
      de_table = de_table, 
      join_column = "gene",
      meta_main_rank_column = "dsi",
      meta_tie_breaker_column = "median_log2FC",
      int_main_rank_column = "avg_log2FC"
    )
    
    pearson <- get_effect_size_correlation(meta_table = meta_table, de_table = de_table, meta_effect_size_column = "median_log2FC", int_effect_size_column = "avg_log2FC")
    
    all_jaccard <- lapply(c(20,50,100), function(t) {
      get_jaccard_index(meta_table, de_table, t)
    })
    names(all_jaccard) <- as.character(c(20,50,100))
    
    
    
    # Add row to results
    comparison_results <- rbind(
      comparison_results,
      data.frame(
        condition = condition,
        sample = sample,
        dataset = dataset,
        spearman = spearman$spearman_correlation,
        pearson = pearson$effect_size_pearson_correlation,
        jaccard_top_20 = all_jaccard$'20'$top_20_jaccard_index,
        jaccard_bottom_20 = all_jaccard$'20'$bottom_20_jaccard_index,
        jaccard_top_50 = all_jaccard$'50'$top_50_jaccard_index,
        jaccard_bottom_50 = all_jaccard$'50'$bottom_50_jaccard_index,
        jaccard_top_100 = all_jaccard$'100'$top_100_jaccard_index,
        jaccard_bottom_100 = all_jaccard$'100'$bottom_100_jaccard_index,
        stringsAsFactors = FALSE
      )
    )
    gc()
  }
  
  return(comparison_results)  
}












# gene heatmap and helper functions--------------------------------------------------------------

get_pdf_gene_heatmaps <- function(path_list, files_to_exclude_list = NULL, condition = "MS", full_condition_name = "Multiple Sclerosis",
                                    pdf_path, genes_to_display   = NULL){
  
  # get meta_list_of_condition_objects and int_list_of_condition_objects
  {
    meta_list_of_condition_objects <- get_list_of_condition_objects(
      path = path_list$meta_path, 
      files_to_exclude = files_to_exclude_list$meta, 
      condition = condition
    )
    
    int_list_of_condition_objects <- get_list_of_condition_objects(
      path = path_list$int_path, 
      files_to_exclude = files_to_exclude_list$int, 
      condition = condition
    )
  }
  

  # ---- meta heatmap ----
  ht_meta <- get_gene_heatmap(
    results_list       = meta_list_of_condition_objects$results_list,
    metric_col         = "median_log2FC",
    method_type        = "meta",
    metric_name        = "Median log2FC",
    condition          = "Multiple Sclerosis",
    min_abs_metric = 0,
    genes_to_display   = genes_to_display 
  )
  
  # ---- integrated heatmap----
  ht_int <- get_gene_heatmap(
    results_list       = int_list_of_condition_objects$results_list,
    metric_col         = "avg_log2FC",
    method_type        = "integrated",
    metric_name        = "avg_log2FC",
    condition          = full_condition_name,
    min_abs_metric = 0,
    genes_to_display   = genes_to_display
  )
  
  pdf(pdf_path, width = 26, height = 20)
  draw(ht_meta$heat_map, annotation_legend_list = attr(ht_meta$heat_map, "extra_legends"))
  draw(ht_int$heat_map, annotation_legend_list = attr(ht_int$heat_map, "extra_legends"))
  dev.off()
 
print(paste("The pdf awaits at", pdf_path))

}



get_gene_heatmap <- function(results_list,
                         metric_col,        # e.g. "median_log2FC" or "avg_log2FC"
                         method_type,       # "meta" | "integrated"
                         metric_name,
                         condition,
                         genes_to_display,
                         min_prop    = 0.75,
                         min_abs_metric = 1,
                         top_n       = 50) {
  


  # 1) long form: one row per gene-sample
  long_df <- imap_dfr(
    results_list,
    ~ .x %>%
      dplyr::select(gene, !!sym(metric_col)) %>%   # tidy-eval on column name
      mutate(sample = .y)
  ) %>% rename(metric = !!sym(metric_col))
  
  # 2) keep genes present in >= min_prop of samples
  min_required <- ceiling(min_prop * length(results_list))
  genes_keep <- long_df %>%
    dplyr::count(gene, name = "n_present") %>%
    dplyr::filter(n_present >= min_required) %>%
    dplyr::pull(gene)
  
  # 3) wide matrix for ComplexHeatmap
  wide_mat <- long_df %>%
    dplyr::filter(gene %in% genes_keep) %>%
    pivot_wider(names_from = sample, values_from = metric) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # 4) build the heatmap (returns the same list as before)
  draw_gene_heatmap(
    wide_mat,
    method_type      = method_type,
    metric_name      = metric_name,
    condition        = condition,
    min_abs_metric   = min_abs_metric,
    top_n            = top_n,
    genes_to_display = genes_to_display
  )
}




draw_gene_heatmap <- function(
    gene_mat,                          # genes × samples matrix
    method_type,                       # "meta" | "integrated"
    metric_name        = "Median log2FC",
    metric_name_plural = "Median log2FC",
    condition          = "Multiple Sclerosis",
    min_abs_metric     = 0,            # (>= 0). e.g. 0.5 keeps |log2FC| ≥ 0.5
    top_n              = 40,           # how many genes to label if auto-selecting
    show_all_gene_labels = FALSE,
    genes_to_display   = NULL          # <- NEW
) {
  

  stopifnot(is.matrix(gene_mat) || inherits(gene_mat, "dgCMatrix"))
  
  # ---------------------------------------------------------------------------
  # 0)  Metadata  (unchanged)
  # ---------------------------------------------------------------------------
  information_df <- qs::qread(
    "/bigdata/bgutman/Collabs/Dror mevorach/data/information_df_neuroimmunulogical_samples.qs"
  )
  
  extract_sample <- if (method_type == "meta") {
    function(x) sub("^[^_]+_[^_]+_", "", x)
  } else {
    function(x) sub("_integrated\\.qs$", "",
                    sub("^[^_]+_[^_]+_", "", x))
  }
  
  sample_key <- extract_sample(colnames(gene_mat))
  key_tbl    <- tibble::tibble(matrix_name = colnames(gene_mat),
                               sample      = sample_key)
  
  filt_meta <- information_df |>
    dplyr::filter(condition == condition) |>
    dplyr::mutate(sample = as.character(sample)) |>
    dplyr::inner_join(key_tbl, by = "sample")
  
  mat_keep <- filt_meta$matrix_name
  if (!length(mat_keep))
    stop("No overlapping samples between gene_mat and metadata for condition ",
         condition, ".")
  
  gene_mat <- gene_mat[, mat_keep, drop = FALSE]
  colnames(gene_mat) <- filt_meta$sample
  
  # ---------------------------------------------------------------------------
  # 1)  Optional DE filtering (unchanged)
  # ---------------------------------------------------------------------------
  if (min_abs_metric > 0) {
    gene_keep <- apply(abs(gene_mat), 1, max, na.rm = TRUE) >= min_abs_metric
    gene_mat  <- gene_mat[gene_keep, , drop = FALSE]
  }
  
  # ---------------------------------------------------------------------------
  # 2)  Decide which genes to label
  # ---------------------------------------------------------------------------
  row_median <- apply(gene_mat, 1, median, na.rm = TRUE)
  
  if (!is.null(genes_to_display)) {            # ---- NEW ----
    ## keep only genes present in the matrix
    genes_to_display <- intersect(genes_to_display, rownames(gene_mat))
    if (!length(genes_to_display))
      stop("None of the genes in 'genes_to_display' are found in 'gene_mat'.")
    
    pos_genes <- genes_to_display[row_median[genes_to_display] >  0]
    neg_genes <- genes_to_display[row_median[genes_to_display] <  0]
    
    label_genes <- genes_to_display
  } else {
    ## original automatic selection
    n_pos <- floor(top_n / 2)
    n_neg <- top_n - n_pos
    
    pos_genes <- names(sort(row_median, decreasing = TRUE))
    pos_genes <- pos_genes[row_median[pos_genes] > 0][seq_len(min(n_pos,
                                                                  sum(row_median > 0)))]
    
    neg_genes <- names(sort(row_median, decreasing = FALSE))
    neg_genes <- neg_genes[row_median[neg_genes] < 0][seq_len(min(n_neg,
                                                                  sum(row_median < 0)))]
    
    need_extra <- top_n - length(c(pos_genes, neg_genes))
    if (need_extra > 0) {
      remainder <- setdiff(names(sort(abs(row_median), decreasing = TRUE)),
                           c(pos_genes, neg_genes))
      extra_genes <- remainder[seq_len(min(need_extra, length(remainder)))]
      if (row_median[extra_genes[1]] >= 0) pos_genes <- c(pos_genes, extra_genes)
      else                                 neg_genes <- c(neg_genes, extra_genes)
    }
    
    label_genes <- c(pos_genes, neg_genes)
  }
  label_idx  <- match(label_genes, rownames(gene_mat))
  label_cols <- vapply(
    row_median[label_genes],
    \(x) {
      if      (x >  3) "firebrick"
      else if (x >  1) "#e95c4b"
      else if (x >=-1) "#000000"
      else if (x >=-3) "#2166ac"
      else             "#00296b"   # ← darker blue, clearly distinct from black
    },
    FUN.VALUE = character(1)
  )
  
  # ---------------------------------------------------------------------------
  # 3)  Column annotation (unchanged palettes but NA = grey85 everywhere)
  # ---------------------------------------------------------------------------
  subtype_pal <- c(CIS = "#bb9457", RRMS = "#432818",
                   "MS-unknown" = "#99582a", CDMS = "#ffe6a7")
  dataset_pal <- c(Esaulova = "#fb8500", Kang = "#8ecae6",
                   Ramesh   = "#ffb703", Schafflick = "#219ebc")
  sex_pal     <- c(Female = "#ff8daa", Male = "#32b7ff")
  edss_pal    <- colorRamp2(c(0, 6), c("#70e000", "#e01e37"))
  age_pal     <- colorRamp2(
    range(filt_meta$age, na.rm = TRUE),
    c("#48CAE4", "#03045e")
  )
  
  col_anno <- HeatmapAnnotation(
    Dataset = filt_meta$dataset,
    Sex     = filt_meta$sex,
    Age     = filt_meta$age,
    EDSS    = filt_meta$EDSS,
    Subtype = filt_meta$condition_subtype,
    col = list(
      Dataset = dataset_pal,
      Sex     = sex_pal,
      Age     = age_pal,
      EDSS    = edss_pal,
      Subtype = subtype_pal
    ),
    na_col = "grey85",
    annotation_legend_param = list(
      EDSS = list(
        title     = "EDSS",
        at        = seq(0, 8, 2),
        labels    = as.character(seq(0, 8, 2)),
        legend_gp = gpar(fill = edss_pal(seq(0, 8, 2)))
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 7)
  )
  # --------------------------------------------------------------------------
  
  # ---------------------------------------------------------------------------
  # 4)  Row annotation – place labels near their rows, avoid collisions
  # ---------------------------------------------------------------------------
  row_anno <- rowAnnotation(
    Gene = anno_mark(
      at        = label_idx,
      labels    = label_genes,
      labels_gp = gpar(fontsize = 7, col = label_cols),
      link_gp   = gpar(lwd = 0.4, col = label_cols),
      
      ## pull the text closer to the heat-map while keeping enough space to read
      extend    = unit(c(0.15, 0.55), "cm")   # ← was c(0.4, 1.0)
      
    ),
    annotation_name_side = "top"
  )
  
  # ---------------------------------------------------------------------------
  # 5)  Draw the heat-map (unchanged)
  # ---------------------------------------------------------------------------
  column_title <- paste(
    metric_name_plural, "Across", condition,
    "Samples in the", tools::toTitleCase(method_type), "Pipeline"
  )
  
  ht_body <- Heatmap(
    gene_mat,
    name                = metric_name,
    col                 = colorRamp2(c(-2, 0, 2),
                                     c("#2166ac", "white", "#b2182b")),
    na_col              = "grey85",
    top_annotation      = col_anno,
    show_row_names      = FALSE,
    show_column_names   = FALSE,
    cluster_rows        = TRUE,
    cluster_columns     = TRUE,
    heatmap_legend_param = list(
      at             = legend_breaks(gene_mat),
      title          = nice_two_line_title(metric_name),
      title_position = "leftcenter-rot"
    ),
    column_title    = column_title,
    column_title_gp = gpar(fontface = "bold", fontsize = 14),
    use_raster       = TRUE,          # keep raster for speed
    raster_by_magick = FALSE          # ← bypasses ImageMagick completely
  )
  
  na_leg <- legend_NA()
  full_ht <- row_anno + ht_body
  attr(full_ht, "extra_legends") <- list(na_leg)
  
  list(
    heat_map = full_ht,
    ht_body = ht_body,
    de_genes = list(pos_genes = pos_genes, neg_genes = neg_genes)
  )
  
  
}


# biological analysis functions and helper functions--------------------------------------------------------------

# helper 1 ── bind_results()  +  built-in checks
bind_results <- function(obj_names, data_list, prefix) {
  
  if (length(obj_names) == 0)
    stop(sprintf("No objects with prefix “%s” found.", prefix), call. = FALSE)
  
  dfs <- lapply(obj_names, function(nm) {
    obj <- data_list[[nm]]
    if (is.null(obj))
      stop(sprintf("Object “%s” is missing in the list.", nm), call. = FALSE)
    if (!"result" %in% slotNames(obj))
      stop(sprintf("Object “%s” lacks a @result slot.", nm), call. = FALSE)
    
    df <- obj@result
    if (!is.data.frame(df))
      stop(sprintf("@result of “%s” is not a data.frame.", nm), call. = FALSE)
    
    # ── column checks ────────────────────────────────────────────────
    if (!"p.adjust" %in% names(df))
      stop(sprintf("Result of “%s” is missing required column: p.adjust", nm),
           call. = FALSE)
    
    has_NES       <- "NES"       %in% names(df)
    has_GeneRatio <- "GeneRatio" %in% names(df)
    
    if (has_NES == has_GeneRatio) {      # both present or both absent → error
      stop(sprintf(
        "Result of “%s” must contain exactly one of 'NES' or 'GeneRatio' (found %s).",
        nm,
        if (has_NES) "both" else "neither"
      ), call. = FALSE)
    }
    
    # ── make GeneRatio numeric if it exists ──────────────────────────
    if (has_GeneRatio) {
      df$GeneRatio <- vapply(df$GeneRatio, function(x) {
        if (is.numeric(x)) return(x)               # already numeric
        if (grepl("^\\d+/\\d+$", x))               # "a/b" -> a/b
          eval(parse(text = x))
        else                                       
          as.numeric(x)                            # fallback (NA if not convertible)
      }, numeric(1))
    }
    

    
    # ── tag dataset and return ───────────────────────────────────────
    df$dataset <- sub(paste0("^", prefix), "", nm)
    df
  })
  
  out <- dplyr::bind_rows(dfs)
  out <- out %>% dplyr::mutate(ID = as.character(ID))
  out                    
}

# ------------------------------------------------------------------
# helper 2 ── get_top_filtered_pathways_per_dataset()
# ------------------------------------------------------------------
get_top_filtered_pathways_per_dataset <- function(full_df) {
  
  # ── 1) basic column checks ─────────────────────────────────────────
  # require both p.adjust and dataset
  required <- c("p.adjust", "dataset")
  miss     <- setdiff(required, names(full_df))
  if (length(miss) > 0) {
    stop("Input is missing required column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  # require exactly one of NES or GeneRatio
  
  has_NES       <- "NES"       %in% names(full_df)
  has_GeneRatio <- "GeneRatio" %in% names(full_df)
  # must have exactly one of NES or GeneRatio
  if (has_NES == has_GeneRatio) {
    stop("Input must contain exactly one of 'NES' or 'GeneRatio' (found ",
         if (has_NES) "both" else "neither", ").", call. = FALSE)
  }
  
  # if we're in the GeneRatio branch, also require Count
  if (has_GeneRatio && !("Count" %in% names(full_df))) {
    stop("ORA input has 'GeneRatio' but is missing required column 'Count'.", call. = FALSE)
  }
  
  # ── 2) branch on which metric we have ──────────────────────────────
  if (has_NES) {
    filtered <- full_df %>%
      dplyr::filter(
        abs(NES) > 1.5,
        p.adjust   < 0.05
      ) %>%
      dplyr::group_by(dataset) %>%
      dplyr::slice_max(order_by = abs(NES), n = 25) %>%
      dplyr::ungroup()
    
  } else {
    filtered <- full_df %>%
      dplyr::filter(
        p.adjust   < 0.05,
        GeneRatio  >= 0.05,
        Count      >= 5
      ) %>%
      dplyr::group_by(dataset) %>%
      dplyr::slice_max(order_by = GeneRatio, n = 25) %>%
      dplyr::ungroup()
  }
  
  # ── 3) return the result ───────────────────────────────────────────
  filtered
}

# ------------------------------------------------------------------
# main  ───────────────── get_full_gsea_and_ora_dfs()
# ------------------------------------------------------------------

get_full_gsea_and_ora_dfs <- function(gsea_return_list = NULL,
                                      ora_return_list  = NULL,
                                      save_dir_gsea    = NULL,
                                      save_dir_ora     = NULL,
                                      apply_filter     = FALSE,
                                      condition,
                                      dataset_name, 
                                      sample_name) {
  ## ── helper: prettify dataset names ----------------------------------------
  prettify_dataset <- function(.tbl) {
    if (is.null(.tbl) || nrow(.tbl) == 0) return(.tbl)
    dplyr::mutate(
      .tbl,
      dataset = dplyr::case_when(
        dataset == "h"    ~ "Hallmark",
        dataset == "gobp" ~ "GO:BP",
        dataset == "gocc" ~ "GO:CC",
        TRUE              ~ stringr::str_to_title(dataset)
      )
    )
  }
  
  ## ── helper: unwrap enrichResult / gseaResult ------------------------------
  get_df <- function(x) {
    if (is.data.frame(x))                  return(x)
    if (methods::is(x, "enrichResult") ||
        methods::is(x, "gseaResult"))      return(x@result)
    if (methods::isS4(x) && "result" %in% slotNames(x))
      return(x@result)
    NULL
  }
  
  ## ── helper: bind safely ---------------------------------------------------
  safe_bind <- function(lst, prefix) {
    sel <- names(lst)[startsWith(names(lst), prefix)]
    if (length(sel) == 0) return(NULL)
    
    sel <- sel[
      purrr::map_lgl(sel, \(nm) {
        df <- get_df(lst[[nm]])
        is.data.frame(df) && nrow(df) > 0
      })
    ]
    if (length(sel) == 0) return(NULL)
    
    bind_results(sel, lst, prefix)
  }
  
  ## output list --------------------------------------------------------------
  out <- list()
  
  #────────────────────────────────────────────────────────────────────────────
  # 1 · GSEA
  #────────────────────────────────────────────────────────────────────────────
  if (!is.null(gsea_return_list)) {
    
    meta_gsea_full <- safe_bind(gsea_return_list, "meta_") |> prettify_dataset()
    
    has_int_gsea   <- !is.null(names(gsea_return_list)) &&
      any(startsWith(names(gsea_return_list), "int_"))
    int_gsea_full  <- if (has_int_gsea)
      safe_bind(gsea_return_list, "int_") |> prettify_dataset()
    else
      NULL

    ## optional saving
    if (!is.null(save_dir_gsea) &&
        (!is.null(meta_gsea_full) || !is.null(int_gsea_full))) {
      
      if (!dir.exists(save_dir_gsea))
        dir.create(save_dir_gsea, recursive = TRUE, showWarnings = FALSE)
      
      if (!is.null(meta_gsea_full))
        full_path_name <- paste0(save_dir_gsea, "/",condition,"_",dataset_name, "_", sample_name, "_full_meta_gsea_df.qs")
      qs::qsave(meta_gsea_full, full_path_name)
      
      if (!is.null(int_gsea_full))
        full_path_name <- paste0(save_dir_gsea, "/",condition,"_",dataset_name, "_", sample_name, "_full_int_gsea_df.qs")
      qs::qsave(int_gsea_full,  full_path_name)
    }
    
    if (!is.null(meta_gsea_full)) {
      out$meta_gsea_top <- if (apply_filter)
        get_top_filtered_pathways_per_dataset(meta_gsea_full)
      else meta_gsea_full
    }
    
    if (!is.null(int_gsea_full)) {
      out$int_gsea_top <- if (apply_filter)
        get_top_filtered_pathways_per_dataset(int_gsea_full)
      else int_gsea_full
    }
  }
  
  #────────────────────────────────────────────────────────────────────────────
  # 2 · ORA
  #────────────────────────────────────────────────────────────────────────────
  if (!is.null(ora_return_list)) {
    
    meta_ora_full <- safe_bind(ora_return_list, "meta_") |> prettify_dataset()
    
    has_int_ora   <- !is.null(names(ora_return_list)) &&
      any(startsWith(names(ora_return_list), "int_"))
    int_ora_full  <- if (has_int_ora)
      safe_bind(ora_return_list, "int_") |> prettify_dataset()
    else
      NULL
    
    if (!is.null(save_dir_ora) &&
        (!is.null(meta_ora_full) || !is.null(int_ora_full))) {
      
      if (!dir.exists(save_dir_ora))
        dir.create(save_dir_ora, recursive = TRUE, showWarnings = FALSE)
      
      if (!is.null(meta_ora_full))
        full_path_name <- paste0(save_dir_ora, "/",condition,"_",dataset_name, "_", sample_name, "_full_meta_ora_df.qs")
      
      qs::qsave(meta_ora_full, full_path_name)
      if (!is.null(int_ora_full))
        full_path_name <- paste0(save_dir_ora, "/",condition,"_",dataset_name, "_", sample_name, "_full_int_ora_df.qs")
      
      qs::qsave(int_ora_full,  full_path_name)
    }
    
    if (!is.null(meta_ora_full)) {
      out$meta_ora_top <- if (apply_filter)
        get_top_filtered_pathways_per_dataset(meta_ora_full)
      else meta_ora_full
    }
    
    if (!is.null(int_ora_full)) {
      out$int_ora_top <- if (apply_filter)
        get_top_filtered_pathways_per_dataset(int_ora_full)
      else int_ora_full
    }
  }
  
  return(out)
}



get_pathway_lollipop <- function(
    df,
    list_of_highlight_pathways,
    manuall_changing_names = NULL,
    title,
    lollipop_color   = "#7C4585",
    axis_title_size  = 28,   # axis titles (x & y)
    legend_size      = 24,   # legend text size (title = legend_size + 2)
    title_size       = 26,   # plot title size
    axis_text_size   = 26    # NEW: axis values (tick labels) size
    # ,is_meta   = TRUE
) {
  
  # ------------------------------------------------------------------ #
  # 1) isolate the pathways of interest --------------------------------
  # ------------------------------------------------------------------ #
  highlight_df <- df %>%
    dplyr::filter(ID %in% list_of_highlight_pathways) %>%
    dplyr::mutate(Description_for_lollipop = Description)      # ← keep a copy
  
  # ————————————————————————————————————————————————
  # 2) manual renames (operate ONLY on the copy)
  # ————————————————————————————————————————————————
  if (length(manuall_changing_names) > 0) {
    for (old_name in names(manuall_changing_names)) {
      new_name <- manuall_changing_names[[old_name]]
      highlight_df$Description_for_lollipop[
        highlight_df$Description_for_lollipop == old_name
      ] <- new_name
    }
  }
  
  # after manual renames, decide which rows still need rule-based tidying
  descriptions_for_change <- setdiff(
    highlight_df$Description_for_lollipop,
    manuall_changing_names
  )
  
  # ————————————————————————————————————————————————
  # 3) rule-based tidying   (Hallmark, Biocarta, everything else)
  # ————————————————————————————————————————————————
  for (description in descriptions_for_change) {
    
    idx <- which(highlight_df$Description_for_lollipop == description)
    if (length(idx) == 0) next
    
    if (highlight_df$dataset[idx] == "Hallmark") {
      
      highlight_df$Description_for_lollipop[idx] <-
        description %>%
        stringr::str_remove("^HALLMARK_") %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_title()
      
    } else if (highlight_df$dataset[idx] == "Biocarta") {
      
      raw   <- description %>%
        stringr::str.remove("^BIOCARTA_") %>%
        stringr::str_replace_all("_", " ")
      first <- stringr::word(raw, 1)
      rest  <- stringr::str_remove(raw, paste0("^", first, " ")) %>% stringr::str_to_title()
      
      highlight_df$Description_for_lollipop[idx] <-
        ifelse(rest == "", toupper(first), paste(toupper(first), rest))
      
    } else {
      
      highlight_df$Description_for_lollipop[idx] <- stringr::str_to_title(description)
      
    }
  }
  
  # ------------------------------------------------------------------ #
  # 4) helper columns --------------------------------------------------
  # ------------------------------------------------------------------ #
  highlight_df <- highlight_df %>%
    dplyr::mutate(display_dataset = ifelse(dataset %in% c("GO:BP", "GO:CC"),
                                           "GO", dataset))
  
  metric <- intersect(c("NES", "GeneRatio"), names(highlight_df))
  if (length(metric) != 1)
    stop("highlight_df must contain **exactly one** of 'NES' or 'GeneRatio'.")
  metric <- metric[[1]]
  
  # ------------------------------------------------------------------ #
  # 4.5) plot title + helper columns ------------------------------------
  # ------------------------------------------------------------------ #
  # title <- if (is_meta) {
  #   paste("PDE-Based GSEA: Key Pathways Enriched in CD59-Deficient", cell_type)
  # } else {
  #   paste("Integration-Based GSEA: Key Pathways Enriched in CD59-Deficient", cell_type)
  # }
  
  # ------------------------------------------------------------------ #
  # 5) lollipop plot  (uses Description_for_lollipop) -----------------
  # ------------------------------------------------------------------ #
  p_pathway_lollipop <- ggplot2::ggplot(
    highlight_df,
    ggplot2::aes(x = .data[[metric]],
                 y = forcats::fct_reorder(Description_for_lollipop, .data[[metric]]),
                 shape = display_dataset,
                 fill  = p.adjust)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = 0,
                                       y    = Description_for_lollipop,
                                       yend = Description_for_lollipop),
                          colour = "grey80", size = 0.5) +
    ggplot2::geom_point(size = 5, colour = "black", stroke = 0.8) +
    ggplot2::scale_shape_manual(
      name   = "Dataset",
      values = c(Biocarta = 25, `GO` = 24, Kegg = 22,
                 Hallmark = 21, Reactome = 23),
      guide  = ggplot2::guide_legend(order = 1)
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c(lollipop_color, "grey"),
      limits  = c(0, 0.1),
      breaks  = c(0, 0.05, 0.1),                 # ticks at 0, 0.05, 0.10
      labels  = function(x) sprintf("%.2f", x),  # 2 decimal places
      name    = "p.adjust",
      guide   = ggplot2::guide_colorbar(order = 2)
    ) +
    ggplot2::labs(
      y     = "Pathway",
      x     = metric,
      title = title
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(
      legend.key.height = grid::unit(1.5, "cm"),
      legend.key.width  = grid::unit(1,   "cm"),
      legend.text       = ggplot2::element_text(size = legend_size),
      legend.title      = ggplot2::element_text(size = legend_size + 2),
      axis.text.y       = ggplot2::element_text(size = axis_text_size),
      axis.text.x       = ggplot2::element_text(size = axis_text_size),
      axis.title.y      = ggplot2::element_text(size = axis_title_size),
      axis.title.x      = ggplot2::element_text(size = axis_title_size),
      plot.title        = ggplot2::element_text(size = title_size,
                                                face = "bold",
                                                hjust = 0.5),
      plot.margin       = ggplot2::margin(1, 1, 1, 1, "cm")
    ) +
    ggplot2::guides()
  
  # ------------------------------------------------------------------ #
  # 6) return --------------------------------------------------------- #
  # ------------------------------------------------------------------ #
  list(
    highlight_df       = highlight_df,
    p_pathway_lollipop = p_pathway_lollipop
  )
}




#  Helper function: Converters one '/'-separated string → vector of symbols ────
convert_to_symbol <- function(core_str) {
  
  ids <- str_split(core_str, "/", simplify = TRUE) |> str_trim()
  
  # A) numeric? → treat as Entrez and map
  is_num <- str_detect(ids, "^[0-9]+$")
  if (any(is_num)) {
    sym_num <- select(org.Hs.eg.db,
                      keys      = ids[is_num],
                      keytype   = "ENTREZID",
                      columns   = "SYMBOL") |>
      deframe()                 # named vector: ENTREZ → SYMBOL
    ids[is_num] <- sym_num[ids[is_num]]
  }
  
  # B) non-numeric but not official? map aliases → SYMBOL
  is_alias <- !is_num & !ids %in% keys(org.Hs.eg.db, "SYMBOL")
  if (any(is_alias)) {
    sym_alias <- select(org.Hs.eg.db,
                        keys      = ids[is_alias],
                        keytype   = "ALIAS",
                        columns   = "SYMBOL") |>
      deframe()
    ids[is_alias] <- coalesce(sym_alias[ids[is_alias]], ids[is_alias])
  }
  
  unique(ids[!is.na(ids)])          # drop NA & duplicates
}

#  Helper function: get gene overlaps

get_overlap <- function(sets, gene_sets){
  ints <- Reduce(intersect, gene_sets[sets])         # use gene_sets
  tibble(pathways = paste(sets, collapse = " ∩ "),
         n_genes  = length(ints),
         genes    = list(ints))
}


get_overlap_genes_meta_and_int_result_df <- function(
    overlap_tbl,   # data.frame with columns `pathways`, `genes`, `n_genes`
    meta_table,    # must have columns `gene`, `dsi`, `median_log2FC`, `q_val_fisher`
    de_table       # must have columns `gene`, `avg_log2FC`, `p_val_adj`
) {
  # 1) nothing to do?
  if (is.null(overlap_tbl) || !any(overlap_tbl$n_genes != 0)) {
    return(list())
  }
  
  # 2) prepare the two look-up tables
  mt <- meta_table %>%
    { if (!"gene" %in% names(.)) rownames_to_column(., "gene") else . } %>%
    dplyr::select(gene, dsi, median_log2FC, q_val_fisher)
  
  dt <- de_table %>%
    { if (!"gene" %in% names(.)) rownames_to_column(., "gene") else . } %>%
    rename(log2FC = avg_log2FC) %>%
    dplyr::select(gene, log2FC, p_val_adj)
  
  # 3) keep only non-zero overlaps
  ov <- overlap_tbl %>% filter(n_genes != 0)
  
  # 4) map over each row, splitting the comma-list, joining metadata
  result <- imap(
    ov$genes,
    ~ {
      pathway_name <- ov$pathways[.y]
      gene_vec     <- str_split(.x, ",", simplify = TRUE)[1, ]
      
      tibble(gene = gene_vec) %>%
        left_join(mt, by = "gene") %>%
        left_join(dt, by = "gene") %>%
        mutate(pathway = pathway_name) %>%
        dplyr::select(pathway, gene, dsi, median_log2FC, q_val_fisher, log2FC, p_val_adj)
    }
  )
  
  # 5) name each list element with a safe variable name
  names(result) <- make.names(ov$pathways)
  
  # convert each tibble → plain data.frame
  result <- lapply(result, as.data.frame)
  
  result
}

get_overlap_genes_between_pathways <- function(
    highlight_df,
    meta_table,     
    de_table,       
    targets_pair   = NULL,
    targets_triple = NULL
) {
  # 1. convert to HUGO & build gene_sets -------------------------------
  highlight_df <- highlight_df |>
    rowwise() |>
    mutate(core_enrichment_hugo = paste(convert_to_symbol(core_enrichment),
                                        collapse = "/")) |>
    ungroup()
  
  gene_sets <- highlight_df %>%
    mutate(gene_vec = str_split(core_enrichment_hugo, "/")) %>%
    { setNames(.$gene_vec, .$Description) }
  
  # 2. if BOTH are NULL, auto-fill targets_pair with every unique pair --
  if (is.null(targets_pair) && is.null(targets_triple)) {
    all_paths    <- names(gene_sets)
    targets_pair <- combn(all_paths, 2, simplify = FALSE)
  }
  
  # 3. compute pair overlaps (only if targets_pair is non-NULL) --------
  if (!is.null(targets_pair)) {
    pairs_tbl <- map_dfr(targets_pair,  ~ get_overlap(.x, gene_sets)) %>%
      arrange(desc(n_genes)) %>%
      mutate(genes = sapply(genes, paste, collapse = ",")) %>%
      as.data.frame()
  } else {
    pairs_tbl <- NULL
  }
  
  # 4. compute triple overlaps (only if targets_triple is non-NULL) ------
  if (!is.null(targets_triple)) {
    triples_tbl <- map_dfr(targets_triple,  ~ get_overlap(.x, gene_sets)) %>%
      arrange(desc(n_genes)) %>%
      mutate(genes = sapply(genes, paste, collapse = ",")) %>%
      as.data.frame()
  } else {
    triples_tbl <- NULL
  }
  
  pairs_list <- list()
  triples_list <- list()
  
  if (!"gene" %in% names(meta_table)) {
    meta_table <- meta_table %>% rownames_to_column("gene")
  }
  if (!"gene" %in% names(de_table)) {
    de_table <- de_table %>% rownames_to_column("gene")
  }
  
  
  
  if (!is.null(pairs_tbl) && any(pairs_tbl$n_genes != 0)) {
    
    # 1) keep only the non-zero overlaps
    pairs_nonzero <- subset(pairs_tbl, n_genes != 0)
    
    pairs_list   <- get_overlap_genes_meta_and_int_result_df(pairs_nonzero,   meta_table, de_table)
    
  }
  
  
  if (!is.null(triples_tbl) && any(triples_tbl$n_genes != 0)) {
    
    # 1) keep only the non-zero overlaps
    triples_nonzero <- subset(triples_tbl, n_genes != 0)
    
    triples_list   <- get_overlap_genes_meta_and_int_result_df(triples_nonzero,   meta_table, de_table)
    
  }
  
  
  # 5. return both tables -----------------------------------------------
  list(
    pairs_tbl   = pairs_tbl,
    triples_tbl = triples_tbl,
    pairs_list = pairs_list,
    triples_list = triples_list
  )
}  


# Unified cell type analysis-----------------------------------------------------------


# remove pathways whose 90% of thier genes are shared with another pathway
dedup <- function(df, overlap = 0.9) {
  keep <- rep(TRUE, nrow(df))
  for (i in seq_len(nrow(df) - 1)) if (keep[i]) {
    gi <- strsplit(df$core_enrichment[i], "/")[[1]]
    dup <- map_lgl((i+1):nrow(df), function(j) {
      gj <- strsplit(df$core_enrichment[j], "/")[[1]]
      min(length(intersect(gi, gj))/length(gi),
          length(intersect(gi, gj))/length(gj)) >= overlap
    })
    keep[(i+1):nrow(df)] <- keep[(i+1):nrow(df)] & !dup
  }
  df[keep, ]
}

#' Merge multiple group-map objects into one - for grouping of the pathways
#'
#' @param ...   Any number of group_map objects (named lists).
#'              You may also supply a single list that *contains* group_maps.
#' @return      A merged group_map (named list). If the same pathway is found
#'              under different group names, a warning is emitted and the
#'              pathway is kept in **all** groups so you can decide later.
merge_group_maps <- function(...) {
  
  # collect inputs ------------------------------------------------------------
  maps <- list(...)
  if (length(maps) == 1 && is.list(maps[[1]]) && !is.character(maps[[1]])) {
    # allow merge_group_maps(list(map1, map2, map3))
    maps <- maps[[1]]
  }
  if (!all(vapply(maps, is.list, logical(1))))
    stop("All inputs must be group_map lists.")
  
  # initialise the merged container ------------------------------------------
  merged <- list()
  
  # helper to add unique values into merged[[grp]]
  add_to_group <- function(grp, desc_vec) {
    if (is.null(merged[[grp]])) {
      merged[[grp]] <<- unique(desc_vec)
    } else {
      merged[[grp]] <<- unique(c(merged[[grp]], desc_vec))
    }
  }
  
  # build the merged map ------------------------------------------------------
  for (gm in maps) {
    for (grp in names(gm)) {
      add_to_group(grp, gm[[grp]])
    }
  }
  
  # optional: detect pathways assigned to >1 group ---------------------------
  desc_to_groups <- list()
  for (grp in names(merged)) {
    for (desc in merged[[grp]]) {
      desc_to_groups[[desc]] <- c(desc_to_groups[[desc]], grp)
    }
  }
  multi_group <- vapply(desc_to_groups, length, integer(1)) > 1
  if (any(multi_group)) {
    dup_msgs <- vapply(
      names(desc_to_groups)[multi_group],
      function(d) sprintf("  • %s : %s",
                          d, paste(desc_to_groups[[d]], collapse = ", ")),
      character(1)
    )
    warning(
      "Some pathways appear in more than one group:\n",
      paste(dup_msgs, collapse = "\n"),
      call. = FALSE
    )
  }
  
  merged
}


# -------------------------------------------------------------------------
# 1)  FUNCTION:  copy the `group` values from any number of reference GSEA
#                tables into a primary table (here B_cells_gsea_clean).
#                If conflicting non-NA group labels are found for the same
#                pathway ID, a warning is issued and the conflicting rows
#                are printed.
# -------------------------------------------------------------------------
copy_groups_from_refs <- function(primary_df, ..., 
                                  id_col    = "ID",
                                  group_col = "group",
                                  desc_col  = "Description") {
  
  refs <- list(...)
  if (length(refs) == 1 && is.list(refs[[1]]) && !is.data.frame(refs[[1]])) {
    refs <- refs[[1]]                           # allow a list-of-refs argument
  }
  ref_names <- names(refs)
  if (is.null(ref_names)) ref_names <- paste0("ref", seq_along(refs))
  
  # 1) gather reference ID–group pairs --------------------------------------
  ref_groups <- bind_rows(lapply(seq_along(refs), function(i) {
    refs[[i]] %>% 
      dplyr::select(
        !!sym(id_col),
        !!sym(group_col),
        !!sym(desc_col)
      ) %>% 
      mutate(source = ref_names[i])
  }))
  
  # 2) detect conflicting non-NA group labels --------------------------------
  conflicts <- ref_groups %>% 
    filter(!is.na(.data[[group_col]])) %>% 
    group_by(.data[[id_col]]) %>% 
    filter(n_distinct(.data[[group_col]]) > 1) %>% 
    ungroup()
  
  if (nrow(conflicts) > 0) {
    warning("Conflicting group assignments detected:")
    print(conflicts %>% arrange(.data[[id_col]]))
  }
  
  # 3) build ID → group lookup (first non-NA wins) ---------------------------
  lookup <- ref_groups %>% 
    group_by(.data[[id_col]]) %>% 
    summarise(
      !!group_col := dplyr::first(.data[[group_col]]),
      .groups = "drop"
    )
  
  # 4) join lookup back onto primary table -----------------------------------
  updated_primary <- primary_df %>% 
    dplyr::select(-dplyr::any_of(group_col)) %>%   # safe even if 'group' absent
    dplyr::left_join(lookup, by = id_col)
  
  updated_primary
}

#subset the pathways into positive and negative
make_subsets <- function(df, sign = c("positive", "negative")) {
  sign <- match.arg(sign)
  if (sign == "positive") {
    subset(df, NES > 0 & !is.na(group))
  } else {
    subset(df, NES < 0 & !is.na(group))
  }
}



# ════════════════════════════════════════════════════════════════════
#  helper: get_intersection_ylim()
# --------------------------------------------------------------------
#  • figures out the highest bar that will appear in the “Intersection
#    size” panel (given df, sets, intersections & mode) and rounds it
#    **up** to the next multiple of 100
#  • returns a single numeric value → feed it straight into
#       scale_y_continuous(limits = c(0, y_lim))
# ════════════════════════════════════════════════════════════════════
get_intersection_ylim <- function(df, sets, intersections, mode) {
  # ── 1. decide which intersections we will actually plot ────────────
  if (is.null(intersections)) {
    # create every combination of ≥ 3 sets
    intersections <- unlist(
      lapply(3:length(sets), \(k) combn(sets, k, simplify = FALSE)),
      recursive = FALSE
    )
    # and later we will keep only those with ≥ 3 members anyway
  }
  
  # ── 2. helper to count one intersection’s size ─────────────────────
  count_int <- function(selected_sets) {
    sel <- df |>
      dplyr::select(dplyr::all_of(sets))
    
    if (mode == "inclusive_intersection") {
      # in ALL selected sets
      dplyr::filter(sel, dplyr::if_all(dplyr::all_of(selected_sets), \(x) x)) |>
        nrow()
    } else {                        # "exclusive_intersection" & everything else
      other <- setdiff(sets, selected_sets)
      dplyr::filter(
        sel,
        dplyr::if_all(dplyr::all_of(selected_sets), \(x) x),   # in the chosen
        dplyr::if_all(dplyr::all_of(other),          \(x) !x)  # in no others
      ) |>
        nrow()
    }
  }
  
  sizes <- vapply(intersections, count_int, integer(1))
  
  ## round UP to the next multiple of 100
  ceiling(max(sizes, 1) / 100) * 100
}


# ───────────────────────────────────────────────────────────────────────────────
#  DOT-PLOT  ▸  NES → colour   •   q-value → size
#             size breaks: 0.20, 0.15, 0.10, 0.05, ≤0.01
#             pathways ordered top-to-bottom by **decreasing** mean NES
# ───────────────────────────────────────────────────────────────────────────────
get_pathway_dotplot <- function(
    df,
    manuall_changing_names = NULL,          # named list: old → new labels
    title                  = NULL,          # plot title
    row_order              = NULL,          # explicit pathway ordering (top → bottom)
    celltype_order         = NULL,          # explicit x-axis order
    colour_low   = "#56B4E9",               # low  NES  (depleted)
    colour_mid   = "beige",                 # neutral
    colour_high  = "#B62055",               # high NES  (enriched)
    size_range   = c(2, 14),                # dot-size range (pts)
    colourbar_width_cm = 5                  # width of NES legend bar
) {
  
  
  ## 1)  PRETTIFY PATHWAY NAMES  ─────────────────────────────────────────── ##
  df <- df %>%
    mutate(Pathway_label = Description)      # keep editable copy
  
  # 1a) manual renames ----------------------------------------------------- #
  if (length(manuall_changing_names) > 0) {
    for (old_name in names(manuall_changing_names)) {
      new_name <- manuall_changing_names[[old_name]]
      df$Pathway_label[df$Pathway_label == old_name] <- new_name
    }
  }
  
  # 1b) rule-based tidy ---------------------------------------------------- #
  to_change <- setdiff(df$Pathway_label, unlist(manuall_changing_names))
  
  for (desc in to_change) {
    
    idx <- which(df$Pathway_label == desc)
    if (length(idx) == 0) next
    
    dataset_tag <- df$dataset[idx[1]]           # take first occurrence
    
    tidy_name <- if (dataset_tag == "Hallmark") {
      
      desc %>%
        str_remove("^HALLMARK_") %>%
        str_replace_all("_", " ") %>%
        str_to_title()
      
    } else if (dataset_tag == "Biocarta") {
      
      raw   <- desc %>% str_remove("^BIOCARTA_") %>% str_replace_all("_", " ")
      first <- word(raw, 1)
      rest  <- str_remove(raw, paste0("^", first, " ")) %>% str_to_title()
      ifelse(rest == "", toupper(first), paste(toupper(first), rest))
      
    } else {
      
      str_to_title(desc)
      
    }
    
    df$Pathway_label[idx] <- tidy_name
  }
  
  ## 2)  FACTOR LEVELS  (ordering)  ──────────────────────────────────────── ##
  if (!is.null(row_order)) {
    df$Pathway_label <- factor(df$Pathway_label, levels = rev(row_order))
  } else {
    lvl <- df %>%
      group_by(Pathway_label) %>%
      summarise(meanNES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(meanNES)) %>%
      pull(Pathway_label)
    df$Pathway_label <- factor(df$Pathway_label, levels = rev(lvl))
  }
  
  if (!is.null(celltype_order)) {
    df$Cell_type <- factor(df$Cell_type, levels = celltype_order)
  }
  
  ## 3)  SIZE MAPPING  (q-value)  ───────────────────────────────────────── ##
  df <- df %>%
    mutate(
      q_val_cap = pmin(pmax(p.adjust, 0.01), 0.20),
      size_map  = 1 - q_val_cap              # smaller q → larger size
    )
  
  breaks_q  <- c(0.20, 0.15, 0.10, 0.05, 0.01)
  breaks_sz <- 1 - breaks_q
  labels_q  <- c("0.20", "0.15", "0.10", "0.05", "≤0.01")
  
  ## 4)  TITLE  ─────────────────────────────────────────────────────────── ##
  if (is.null(title)) title <- "Pathway Enrichment Across Cell Types"
  
  ## 5)  PLOT  ───────────────────────────────────────────────────────────── ##
  p <- ggplot(
    df,
    aes(x    = Cell_type,
        y    = Pathway_label,
        fill = NES,
        size = size_map)) +
    geom_point(shape = 21, colour = "black", stroke = 0.3) +
    
    # 5a) colour scale (NES) -------------------------------------------- #
    scale_fill_gradient2(
      name     = "Enrichment Score (NES)",
      low      = colour_low,
      mid      = colour_mid,
      high     = colour_high,
      midpoint = 0,
      guide    = guide_colourbar(
        barwidth  = grid::unit(colourbar_width_cm, "cm"),
        barheight = grid::unit(0.6, "cm"),
        label.hjust = 0.5)
    ) +
    
    # 5b) size scale (q-value) ------------------------------------------ #
    scale_size_continuous(
      name    = "q-value",
      breaks  = breaks_sz,
      labels  = labels_q,
      limits  = c(min(breaks_sz), max(breaks_sz)),
      range   = size_range,
      guide   = guide_legend(override.aes = list(shape = 21), reverse = TRUE)
    ) +
    
    labs(x = NULL, y = NULL, title = title) +
    
    theme_pubr() +
    theme(
      plot.title        = element_text(size = 28, face = "bold", hjust = 0.5),
      axis.text.y       = element_text(size = 22),
      axis.text.x       = element_text(size = 22, angle = 45, hjust = 1),
      legend.title      = element_text(size = 22),
      legend.text       = element_text(size = 20),
      legend.key.height = unit(1.3, "cm"),
      legend.key.width  = unit(0.6, "cm"),
      panel.grid.major  = element_blank(),
      panel.border      = element_rect(colour = "black", fill = NA)
    )
  
  return(p)
}






# ════════════════════════════════════════════════════════════════════════
#  get_upset_plot()
#  ----------------------------------------------------------------------
#  A convenience wrapper around ComplexUpset that:
#    1. renames `group → pathway_group`
#    2. optionally re-codes selected genes to a user-supplied group
#    3. standardises the “Stress” label to “Replication-Stress Response”
#    4. draws a publication-ready UpSet plot with sensible defaults
#       (font sizes, palettes, ordering, etc.)
#  ----------------------------------------------------------------------
#  Arguments
#  • output_df      – data-frame with at least ‘gene’, ‘group’, and the
#                     boolean set columns (Monocytes, NK, …)
#  • recode_genes   – *named* character vector, e.g.
#                     c("RFC2" = "Proliferation", "NSD2" = "Proliferation")
#                     Default = NULL (no recoding)
#  • sets           – character vector naming the set columns
#  • group_palette  – named named vector colour palette
#  • legend_breaks  – factor order for the legend.  NULL ⇒ all palette
#                     names except “Migration / Chemotaxis”
#  • intersections  – list of intersections to keep.  NULL ⇒ the 5-way,
#                     4-way, 3-way combo used in your example
#  • big_title / big_axis / big_txt – font sizes
#  • my_themes      – list passed to ComplexUpset::upset_modify_themes().
#                     NULL ⇒ sensible defaults built from the size args.
#  • mode           – "inclusive_intersection" | "exclusive_intersection"
#                     (or any other valid ComplexUpset mode)
#  • plot_title     – character string for ggtitle
#  • ...            – extra arguments forwarded to ComplexUpset::upset()
#  ----------------------------------------------------------------------
#  Returns a ggplot object.
# ════════════════════════════════════════════════════════════════════════
get_upset_plot <- function(
    df,
    plot_title ,
    sets           = c("Monocytes", "NKs", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
    group_palette  = c(
      "Other / No Pathway"           = "#393E46",
      "Metabolism"                   = "#B2DF8A",
      "Immune Processes"             = "#ff4d6d",
      "Stress and Genome Regulation" = "#E9A5F1"
    ),
    legend_breaks  = NULL,
    legend_order   = NULL,
    intersections  = NULL,
    big_title  = 28,
    big_axis   = 26,
    big_txt    = 22,
    big_legend = 22,
    my_themes  = NULL,
    mode       = "inclusive_intersection",
    legend_title = "Pathway Group",
    group_column = "group",
    ...
) {
  
  ## 0. rename & coerce factor so FIRST color in group_palette is at the BOTTOM
  if (group_column %in% names(df)) {
    df <- dplyr::rename(df, pathway_group = !!rlang::sym(group_column))
  } else {
    stop(sprintf("`%s` not found in `df`.", group_column), call. = FALSE)
  }
  
  palette_order <- names(group_palette)       # user-given order for legend
  stack_levels  <- rev(palette_order)         # reversed → first color goes to bottom
  
  df$pathway_group <- factor(df$pathway_group, levels = stack_levels)
  
  ## 0.5 reorder `sets` by number of TRUEs (shortest first)
  if (!all(sets %in% names(df))) {
    missing_cols <- setdiff(sets, names(df))
    stop("These `sets` columns are missing from df: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  true_counts <- vapply(sets, function(s) sum(isTRUE(df[[s]]) | df[[s]] == TRUE, na.rm = TRUE), numeric(1))
  sets <- sets[order(true_counts, decreasing = FALSE)]
  
  ## 1. legend order/breaks
  if (is.null(legend_order)) legend_order <- palette_order
  
  groups_in_df <- levels(df$pathway_group)
  if (!setequal(groups_in_df, legend_order)) {
    missing <- setdiff(groups_in_df, legend_order)
    extra   <- setdiff(legend_order, groups_in_df)
    stop(
      paste0(
        "legend_order must contain exactly the labels in df$pathway_group.\n",
        if (length(missing)) paste0(" • Missing: ", paste(missing, collapse = ", "), "\n") else "",
        if (length(extra))   paste0(" • Extra:   ", paste(extra,   collapse = ", "), "\n") else ""
      ),
      call. = FALSE
    )
  }
  
  base_breaks    <- legend_order
  breaks_to_show <- if (is.null(legend_breaks)) base_breaks else base_breaks[!base_breaks %in% legend_breaks]
  
  ## 2. intersections/min-degree
  if (is.null(intersections)) {
    intersections_param <- "all"
    min_degree_param    <- 0
  } else {
    intersections_param <- intersections
    min_degree_param    <- 0
  }
  
  ## 3. themes
  if (is.null(my_themes)) {
    my_themes <- ComplexUpset::upset_modify_themes(list(
      'intersections_matrix' = ggplot2::theme(
        axis.text.y  = ggplot2::element_text(size = big_txt),
        axis.text.x  = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = big_axis),
        axis.title.y = ggplot2::element_blank()
      ),
      'Intersection size' = ggplot2::theme(
        axis.text.y  = ggplot2::element_text(size = big_txt, family = "sans"),
        axis.title.y = ggplot2::element_text(size = big_axis, family = "sans"),
        legend.title = ggplot2::element_text(size = big_legend, family = "sans"),
        legend.text  = ggplot2::element_text(size = big_legend, family = "sans")
      ),
      'overall_sizes' = ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = big_axis),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_text(size = big_txt),
        axis.text.y  = ggplot2::element_blank()
      ),
      'default' = ggplot2::theme(
        text = ggplot2::element_text(size = big_txt, family = "sans")
      )
    ))
  }
  
  ## 4. y-axis limit helper
  y_lim <- get_intersection_ylim(
    df            = df,
    sets          = sets,
    intersections = if (is.null(intersections)) NULL else intersections,
    mode          = mode
  )
  
  ## 5. plot
  ComplexUpset::upset(
    df,
    intersect             = tidyselect::all_of(sets),
    name                  = "Cell Type",
    encode_sets           = FALSE,
    intersections         = intersections_param,
    mode                  = mode,
    sort_intersections_by = "cardinality",
    sort_intersections    = "ascending",
    sort_sets             = FALSE,
    keep_empty_groups     = TRUE,
    width_ratio           = 0.12,
    wrap                  = TRUE,
    themes                = my_themes,
    min_degree            = min_degree_param,
    base_annotations = list(
      "Intersection size" =
        ComplexUpset::intersection_size(
          counts  = TRUE,
          mapping = ggplot2::aes(fill = pathway_group),
          text    = list(size = 6),
          mode    = mode
        ) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_fill_manual(
          values = group_palette,     # keep palette as given
          breaks = breaks_to_show,    # legend order
          limits = stack_levels,      # stacking order (reversed)
          name   = legend_title
        ) +
        ggplot2::theme(
          legend.title    = ggplot2::element_text(size = big_axis, hjust = .5),
          legend.text     = ggplot2::element_text(size = big_legend),
          legend.key.size = grid::unit(0.8, "cm")
        )
    ),
    ...
  ) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme_minimal(base_size = big_txt + 4) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size  = big_title,
        face  = "bold",
        hjust = 0.5
      )
    )
}
















# ──────────────────────────────────────────────────────────────────────────────
#  get_de_pathways_for_core_enrichment_genes()
#  ---------------------------------------------------------------------------
#  Given:
#    • pathway_list – the same named list of GSEA data-frames used before
#    • genes_vec    – character vector of HGNC symbols to inspect
#
#  Returns a named list.  Each entry is a tibble with the pathways in which
#  the given gene is part of the *core_enrichment* set, along with:
#      cell_type   | ID | Description | NES | group
# ──────────────────────────────────────────────────────────────────────────────
get_de_pathways_for_core_enrichment_genes <- function(pathway_list,genes_vec) {
  
  ## 1 ─ map list-element names → canonical cell-type strings -----------------
  cell_lookup <- c(
    NK          = "NK",
    monocytes   = "Monocytes",
    B_cells     = "B_cells",
    CD4_T_cells = "CD4_T_cells",
    CD8_T_cells = "CD8_T_cells"
  )
  
  ## 2 ─ build Entrez→Symbol map (reuse logic from build_gene_presence_df) ----
  entrez_pool <- unlist(lapply(pathway_list, function(df) {
    df %>%
      dplyr::filter(!stringr::str_detect(dataset,
                                         stringr::regex("^hallmark$", TRUE))) %>%
      dplyr::pull(core_enrichment) %>%
      stringr::str_split("/") %>%
      unlist()
  }), use.names = FALSE)
  
  if (length(entrez_pool)) {
    mart <- biomaRt::useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl",
      mirror  = "useast"
    )
    map_tbl <- biomaRt::getBM(
      attributes = c("entrezgene_id", "hgnc_symbol"),
      filters    = "entrezgene_id",
      values     = unique(entrez_pool),
      mart       = mart
    )
    entrez2symbol <- setNames(map_tbl$hgnc_symbol,
                              as.character(map_tbl$entrezgene_id))
  } else {
    entrez2symbol <- character()
  }
  
  ## 3 ─ helper to convert the core_enrichment string into HGNC symbols -------
  split_to_symbols <- function(core_str, ds_is_hallmark) {
    raw <- stringr::str_split(core_str, "/")[[1]]
    if (ds_is_hallmark) {
      raw                                   # already symbols
    } else {
      sy <- unname(entrez2symbol[raw])
      sy[!is.na(sy) & nzchar(sy)]           # drop NAs / empty
    }
  }
  
  ## 4 ─ construct result list ------------------------------------------------
  results <- vector("list", length(genes_vec))
  names(results) <- genes_vec
  
  # pre-allocate empty tibble for genes that never appear
  empty_tbl <- tibble::tibble(
    cell_type   = character(),
    ID          = character(),
    Description = character(),
    NES         = numeric(),
    group       = character()
  )
  
  for (gene in genes_vec) {
    holder <- list()
    
    for (elt in names(pathway_list)) {
      
      generic   <- sub("_gsea_clean_.*$", "", elt)
      cell_type <- cell_lookup[[generic]]
      if (is.null(cell_type))
        next                                           # unknown / skip
      
      df <- pathway_list[[elt]]
      if (!all(c("ID", "Description", "NES", "group", "dataset",
                 "core_enrichment") %in% names(df))) {
        stop("Missing expected columns in list element ‘", elt, "’.")
      }
      
      # iterate over pathways, but vectorised for speed
      gene_in_pathway <- mapply(function(core_str, ds) {
        gene %in% split_to_symbols(core_str,
                                   stringr::str_detect(ds,
                                                       stringr::regex("^hallmark$", TRUE)))
      },
      core_str = df$core_enrichment,
      ds       = df$dataset,
      USE.NAMES = FALSE)
      
      if (any(gene_in_pathway)) {
        hit_tbl <- df[gene_in_pathway, ] %>%
          dplyr::transmute(
            cell_type   = cell_type,
            ID,
            Description,
            NES,
            group
          )
        holder[[length(holder) + 1L]] <- hit_tbl
      }
    }
    
    results[[gene]] <- if (length(holder))
      dplyr::bind_rows(holder)
    else
      empty_tbl
  }
  
  results
}






# ──────────────────────────────────────────────────────────────────────────────
# General presence-matrix builder
# ──────────────────────────────────────────────────────────────────────────────
build_presence_df <- function(
    obj_list,
    id_col      = NULL,
    meta_cols   = NULL,
    base_map    = c(
      "NK"          = "NKs",
      "monocytes"   = "Monocytes",
      "B_cells"     = "B Cells",
      "CD4_T_cells" = "CD4+ T Cells",
      "CD8_T_cells" = "CD8+ T Cells"
    ),
    cell_order  = c("NKs", "Monocytes", "B Cells",
                    "CD4+ T Cells", "CD8+ T Cells"),
    filter_nocoding = FALSE        # NEW
) {
  
  ## ── 0 · infer id column ───────────────────────────────────────────────────
  if (is.null(id_col)) {
    sample_names <- names(obj_list[[1L]])
    if      ("ID"   %in% sample_names) id_col <- "ID"
    else if ("gene" %in% sample_names) id_col <- "gene"
    else stop("Couldn’t infer `id_col`. Please supply it explicitly.")
  }
  
  ## ── 1 · long table ────────────────────────────────────────────────────────
  long_tbl <- purrr::imap_dfr(obj_list, function(df, nm) {
    base <- nm %>%
      stringr::str_remove("_gsea_(positive|negative)$") %>%
      stringr::str_remove("_meta_table$")
    if (!base %in% names(base_map)) return(NULL)
    
    keep_cols <- c(id_col, meta_cols)
    keep_cols <- keep_cols[keep_cols %in% names(df)]
    
    df %>%
      dplyr::select(dplyr::all_of(keep_cols)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(cell = base_map[base],
                    presence = TRUE)
  })
  
  ## ── 2 · widen ─────────────────────────────────────────────────────────────
  wide_tbl <- long_tbl %>%
    tidyr::pivot_wider(
      names_from  = cell,
      values_from = presence,
      values_fill = FALSE
    )
  
  ## ── 3 · column ordering ───────────────────────────────────────────────────
  wide_tbl <- wide_tbl %>%
    dplyr::relocate(dplyr::all_of(id_col))
  
  if (!is.null(meta_cols)) {
    wide_tbl <- wide_tbl %>%
      dplyr::relocate(dplyr::all_of(meta_cols),
                      .after = !!rlang::sym(id_col))
  }
  
  wide_tbl <- wide_tbl %>%
    dplyr::relocate(dplyr::all_of(cell_order[cell_order %in% names(wide_tbl)]),
                    .after = dplyr::last_col())
  
  ## ── 4 · optional: filter to protein‑coding genes ──────────────────────────
  if (isTRUE(filter_nocoding) && identical(id_col, "gene") && "gene" %in% names(wide_tbl)) {
    # safely attempt biomaRt lookup
    coding_syms <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = unique(wide_tbl$gene),
        mart       = mart
      )
      annot %>%
        dplyr::filter(.data$gene_biotype == "protein_coding") %>%
        dplyr::pull(.data$hgnc_symbol) %>%
        unique()
    }, error = function(e) {
      warning("Protein-coding filter skipped (biomaRt error: ", conditionMessage(e), ")")
      NULL
    })
    
    if (!is.null(coding_syms)) {
      wide_tbl <- dplyr::filter(wide_tbl, .data$gene %in% coding_syms)
    }
  }
  
  wide_tbl
}






get_core_enrichment_genes_for_de_pathways <- function(path_ids, pathway_list) {
  
  ## ── 1. map list-element names → canonical cell-type strings ──────────
  cell_lookup <- c(
    NK          = "NKs",
    monocytes   = "Monocytes",
    B_cells     = "B Cells",
    CD4_T_cells = "CD4+ T Cells",
    CD8_T_cells = "CD8+ T Cells"
  )
  cell_cols <- unname(cell_lookup)      # keeps output-column order
  
  ## ── 2. Entrez → HGNC map (offline, no Ensembl call) ──────────────────
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
  
  # collect every core_enrichment string, split on “/”, flatten
  entrez_pool <- unique(unlist(
    strsplit(
      unlist(lapply(pathway_list, `[[`, "core_enrichment"), use.names = FALSE),
      "/", fixed = TRUE
    ), use.names = FALSE
  ))
  
  if (length(entrez_pool)) {
    entrez2symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys     = entrez_pool,
      column   = "SYMBOL",
      keytype  = "ENTREZID",
      multiVals = "first"
    )
  } else {
    entrez2symbol <- character()
  }
  
  split_to_symbols <- function(core_str, ds_is_hallmark) {
    raw <- strsplit(core_str, "/", fixed = TRUE)[[1]]
    if (ds_is_hallmark) raw
    else {
      sy <- unname(entrez2symbol[raw])
      sy[!is.na(sy) & nzchar(sy)]
    }
  }
  
  ## ── 3. iterate over requested pathway IDs ────────────────────────────
  result_list <- vector("list", length(path_ids))
  names(result_list) <- path_ids
  
  empty_df <- tibble::tibble(
    gene   = character(),
    !!!setNames(rep(list(logical()), length(cell_cols)), cell_cols),
    count  = integer()
  )
  
  for (pid in path_ids) {
    gene_sets <- vector("list", length(cell_cols))
    names(gene_sets) <- cell_cols
    
    for (elt in names(pathway_list)) {
      generic <- sub("_gsea_.*$", "", elt)
      cell    <- cell_lookup[[generic]]
      if (is.null(cell)) next
      
      df <- pathway_list[[elt]] |> dplyr::filter(ID == pid)
      if (nrow(df) == 0) next
      
      genes <- unlist(mapply(
        split_to_symbols,
        core_str       = df$core_enrichment,
        ds_is_hallmark = stringr::str_detect(df$dataset, 
                                             stringr::regex("^hallmark$", TRUE)),
        USE.NAMES = FALSE))
      
      gene_sets[[cell]] <- unique(c(gene_sets[[cell]], genes))
      if (!exists("path_desc", inherits = FALSE) && nrow(df))
        path_desc <- df$Description[1]
    }
    
    if (!length(unlist(gene_sets))) {
      result_list[[pid]] <- empty_df
      next
    }
    
    all_genes <- sort(unique(unlist(gene_sets)))
    bool_mat  <- sapply(cell_cols, \(cl) all_genes %in% gene_sets[[cl]])
    count_vec <- rowSums(bool_mat)
    
    out_df <- tibble::tibble(
      gene  = all_genes,
      !!!as.data.frame(bool_mat),
      count = count_vec
    ) |> 
      dplyr::arrange(dplyr::desc(count), gene)
    
    result_list[[path_desc %||% pid]] <- out_df
    rm(path_desc, inherits = FALSE)
  }
  
  result_list
}




add_pathway_group_reactome_and_goslim <- function(
    gene_upset_df,
    id_col      = "gene",
    other_label = "Other / No Pathway",
    ann         = NULL,
    BPPARAM     = NULL          # optional BiocParallel param
){
  ## ---- 1 · SYMBOL → ENTREZ --------------------------------------------------
  gene_syms <- gene_upset_df[[id_col]]
  gene_entrezid <- {
    tmp <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = gene_syms,
      keytype   = "SYMBOL",
      column    = "ENTREZID",
      multiVals = "first")
    unique(unname(na.omit(tmp)))
  }
  if (!length(gene_entrezid)) return(
    dplyr::mutate(gene_upset_df,
                  group_reactome = other_label,
                  group_go       = other_label)
  )
  
  ## ---- 2 · Reactome mapping -------------------------------------------------
  ext2path <- get("reactomeEXTID2PATHID", asNamespace("reactome.db"))
  mapL     <- mget(gene_entrezid, ext2path, ifnotfound = NA)
  
  lp_fun <- if (!is.null(BPPARAM)) {
    function(X, FUN) BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM)
  } else lapply
  
  gene2path <- data.table::rbindlist(
    lp_fun(names(mapL), function(ez){
      pid <- mapL[[ez]]
      if (all(is.na(pid))) return(NULL)
      data.table::data.table(ENTREZID = ez, PathwayID = pid)
    }),
    use.names = TRUE
  )
  if (!nrow(gene2path)) {
    return(dplyr::mutate(gene_upset_df,
                         group_reactome = other_label,
                         group_go       = other_label))
  }
  
  ## ---- 2b · Reactome hierarchy helpers -------------------------------------
  get_top_ids <- function(species = "Homo sapiens",
                          cache_dir = tools::R_user_dir("reactomeTop","cache")){
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    rel_file <- file.path(cache_dir, "ReactomePathwaysRelation.txt")
    pth_file <- file.path(cache_dir, "ReactomePathways.txt")
    if (!file.exists(rel_file))
      download.file("https://reactome.org/download/current/ReactomePathwaysRelation.txt",
                    rel_file, quiet = TRUE)
    if (!file.exists(pth_file))
      download.file("https://reactome.org/download/current/ReactomePathways.txt",
                    pth_file, quiet = TRUE)
    rel <- data.table::fread(rel_file, col.names = c("Parent","Child"),
                             showProgress = FALSE)
    pth <- data.table::fread(pth_file, col.names = c("PathwayID","PathwayName","Species"),
                             showProgress = FALSE)
    list(
      top_ids = intersect(setdiff(rel$Parent, rel$Child),
                          pth[Species == species, PathwayID]),
      rel = rel
    )
  }
  cache   <- get_top_ids()
  top_ids <- cache$top_ids
  rel     <- cache$rel
  
  data.table::setorder(rel, Child, Parent)
  child_parent <- rel$Parent
  names(child_parent) <- rel$Child
  child_parent <- child_parent[!duplicated(names(child_parent))]
  
  ## ---- Vectorised root finding ----------------------------------------------
  ids <- gene2path$PathwayID
  unresolved <- !is.na(ids) & !(ids %in% top_ids)
  iter <- 0L; max_iter <- 50L
  while (any(unresolved) && iter < max_iter) {
    ids[unresolved] <- unname(child_parent[ids[unresolved]])
    ids[lengths(ids) == 0] <- NA_character_
    unresolved <- !is.na(ids) & !(ids %in% top_ids)
    iter <- iter + 1L
  }
  gene2path[, TopID := ids]
  gene2path  <- gene2path[!is.na(TopID)]
  
  ## ---- 2c · Names & collapse per gene --------------------------------------
  all_ids <- unique(c(gene2path$TopID, gene2path$PathwayID))
  id2name <- unlist(mget(na.omit(all_ids), reactomePATHID2NAME, ifnotfound = NA))
  
  gene2path[, Cluster := data.table::fifelse(!is.na(TopID),
                                             id2name[TopID],
                                             id2name[PathwayID])]
  gene2path[, Cluster := sub("^Homo sapiens:\\s*", "", Cluster)]
  
  gene2path_final <- gene2path[
    , .(Cluster = {
      tab <- table(na.omit(Cluster))
      if (length(tab)) names(tab)[which.max(tab)] else NA_character_
    }),
    by = ENTREZID
  ]
  id2sym <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                  keys      = gene2path_final$ENTREZID,
                                  keytype   = "ENTREZID",
                                  column    = "SYMBOL",
                                  multiVals = "first")
  gene2path_final[, HUGO := id2sym[ENTREZID]]
  
  ## ---- 3 · GO + Slim mapping (ALL ontologies) -------------------------------
  if (is.null(ann)) {
    ann <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys    = gene_entrezid,
                                 keytype = "ENTREZID",
                                 columns = c("GO","ONTOLOGY"))
    ann <- ann[!is.na(ann$GO), c("ENTREZID","GO","ONTOLOGY")]
    ann <- ann[ann$ONTOLOGY %in% c("BP", "MF"), ]   # keep only BP & MF
    
  } else {
    ann <- ann[ann$ENTREZID %in% gene_entrezid, ]
  }
  
  # Load slim & ancestor environments
  slim_obo   <- "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"
  goslim     <- getOBOCollection(slim_obo)
  goslim_ids <- unlist(ids(goslim))
  
  data("GOBPANCESTOR", package = "GO.db")
  data("GOCCANCESTOR", package = "GO.db")
  data("GOMFANCESTOR", package = "GO.db")
  
  # Worker to map a single GO term (with its ontology) to a Slim
  worker_fun <- function(idx){
    go  <- ann$GO[idx]
    ont <- ann$ONTOLOGY[idx]
    anc <- switch(ont,
                  "BP" = GOBPANCESTOR[[go]],
                  "CC" = GOCCANCESTOR[[go]],
                  "MF" = GOMFANCESTOR[[go]],
                  NULL)
    # Include the term itself (some slim terms are direct children)
    chain <- unique(c(go, anc))
    hit <- chain[chain %in% goslim_ids]
    if (length(hit)) hit[1] else NA_character_
  }
  
  if (nrow(ann)) {
    if (!is.null(BPPARAM)) {
      ann$Slim <- unlist(BiocParallel::bplapply(seq_len(nrow(ann)), worker_fun,
                                                BPPARAM = BPPARAM),
                         use.names = FALSE)
    } else {
      ann$Slim <- vapply(seq_len(nrow(ann)), worker_fun, character(1))
    }
  } else {
    ann$Slim <- character(0)
  }
  
  ## Collapse to one Slim per ENTREZID (mode)
  gene2bucket <- data.table::as.data.table(ann)[,
                                                .(Slim = {
                                                  tab <- table(na.omit(Slim))
                                                  if (length(tab)) names(tab)[which.max(tab)] else NA_character_
                                                }),
                                                by = ENTREZID
  ]
  
  gene2bucket[, Symbol := AnnotationDbi::mapIds(org.Hs.eg.db,
                                                keys      = ENTREZID,
                                                column    = "SYMBOL",
                                                keytype   = "ENTREZID",
                                                multiVals = "first")]
  
  # Lookup Slim names
  if (any(!is.na(gene2bucket$Slim))) {
    go_lookup <- AnnotationDbi::select(GO.db,
                                       keys    = na.omit(unique(gene2bucket$Slim)),
                                       keytype = "GOID",
                                       columns = "TERM")
  } else {
    go_lookup <- data.frame(GOID = character(), TERM = character())
  }
  
  final_named <- gene2bucket |>
    dplyr::left_join(go_lookup, by = c("Slim" = "GOID")) |>
    dplyr::mutate(Bucket = TERM) |>
    dplyr::select(Symbol, Bucket)
  
  ## ---- 4 · Merge & output ---------------------------------------------------
  gene_upset_df |>
    dplyr::left_join(gene2path_final |>
                       dplyr::select(HUGO, Cluster),
                     by = setNames("HUGO", id_col)) |>
    dplyr::left_join(final_named,
                     by = setNames("Symbol", id_col)) |>
    dplyr::mutate(
      group_reactome = dplyr::if_else(is.na(Cluster), other_label, Cluster),
      group_go       = dplyr::if_else(is.na(Bucket),  other_label, Bucket),
      .keep = "unused"
    )
}




# ======================================================================
#  get_scatter_ORA()
# ----------------------------------------------------------------------
#  positive_ora_df  data.frame   ORA table for positively enriched genes
#  negative_ora_df  data.frame   ORA table for negatively enriched genes
#  title            character    plot title
#  big_title        numeric      title font size     (default 28)
#  big_axis         numeric      axis-label font     (default 26)
#  big_txt          numeric      axis-tick & legend  (default 22)
#  big_legend       numeric      legend-title font   (default 22)
#  label_map        named chr    IDs → label text
#  fill_pal         named chr    category colours
#  shape_pal        named int    dataset shapes
# ----------------------------------------------------------------------
#  RETURNS          ggplot       scatter plot object
# ======================================================================
get_scatter_ORA <- function(positive_ora_df,
                            negative_ora_df,
                            title      = "Pathway Gene Ratios: Positive vs. Negative DEGs",
                            big_title  = 28,   # plot title
                            big_axis   = 26,   # axis titles
                            big_txt    = 22,   # axis ticks, legend text **and label text**
                            big_legend = 22,   # legend titles
                            label_map = c(
                              "R-HSA-71291"                     = "AA Metabolism",
                              # "R-HSA-376176"                    = "ROBO Signaling",
                              "R-HSA-9633012"                   = "GCN2 – AA deficiency",
                              "HALLMARK_TNFA_SIGNALING_VIA_NFKB"= "TNF-α / NF-κB",
                              "HALLMARK_INFLAMMATORY_RESPONSE"  = "Inflammatory Response",
                              "R-HSA-2408557"                   = "Selenocysteine Synthesis",
                              "GO:0097191" = "Extrinsic Apoptosis"),
                            fill_pal = c(
                              "Positive only" = "#219ebc",
                              "Negative only" = "#B62055",
                              "Both ORAs"     = "#fb8500"
                            ),
                            shape_pal = c(
                              Biocarta = 25, GO = 24, Kegg = 22,
                              Hallmark = 21, Reactome = 23
                            )) {
  
  ## ---- 1 · Merge ORA tables -----------------------------------------
  ora_df_for_scatter <- dplyr::full_join(
    dplyr::select(positive_ora_df, ID, Description,
                  GeneRatio_p = GeneRatio, dataset_p = dataset),
    dplyr::select(negative_ora_df, ID, Description,
                  GeneRatio_n = GeneRatio, dataset_n = dataset),
    by = c("ID", "Description")
  ) |>
    dplyr::mutate(
      GeneRatio_p = tidyr::replace_na(GeneRatio_p, 0),
      GeneRatio_n = tidyr::replace_na(GeneRatio_n, 0),
      dataset     = dplyr::coalesce(dataset_p, dataset_n)
    ) |>
    dplyr::transmute(
      ID, Description, dataset,
      `GeneRatio(pDEGs)` = GeneRatio_p,
      `GeneRatio(nDEGs)` = GeneRatio_n
    )
  
  ## ---- 2 · Extra columns -------------------------------------------
  ora_df_for_scatter <- ora_df_for_scatter |>
    dplyr::mutate(
      display_dataset = ifelse(dataset %in% c("GO:BP", "GO:CC"), "GO", dataset),
      category = dplyr::case_when(
        `GeneRatio(pDEGs)` > 0 & `GeneRatio(nDEGs)` > 0 ~ "Both ORAs",
        `GeneRatio(pDEGs)` > 0                          ~ "Positive only",
        `GeneRatio(nDEGs)` > 0                          ~ "Negative only",
        TRUE                                            ~ "None"
      )
    )
  
  ## ---- 3 · Labels ---------------------------------------------------
  label_data <- ora_df_for_scatter |>
    dplyr::filter(ID %in% names(label_map)) |>
    dplyr::mutate(label_txt = label_map[ID])
  
  ## ---- 4 · Axis helpers --------------------------------------------
  x_max    <- max(ora_df_for_scatter$`GeneRatio(pDEGs)`)
  y_max    <- max(ora_df_for_scatter$`GeneRatio(nDEGs)`)
  axis_pad <- 0.02
  x_breaks <- pretty(c(0, x_max), n = 5)
  y_breaks <- pretty(c(0, y_max), n = 5)
  
  ## ---- 5 · Convert text size (theme uses pt, label size uses mm) ----
  mm_per_pt <- 1 / ggplot2::.pt                 # ≈ 0.3528
  label_size_mm <- big_txt * mm_per_pt          # pathway label size
  
  ## ---- 6 · Plot -----------------------------------------------------
  p <- ggplot2::ggplot(
    ora_df_for_scatter,
    ggplot2::aes(x = `GeneRatio(pDEGs)`, y = `GeneRatio(nDEGs)`)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = category, shape = display_dataset),
      size   = 3.5,
      colour = "black",
      stroke = 0.5,
      alpha  = 0.9
    ) +
    
    ## colour legend
    ggplot2::scale_fill_manual(
      values = fill_pal,
      name   = "Enrichment source",
      guide  = ggplot2::guide_legend(
        order = 1,
        override.aes = list(shape = 21, colour = "black",
                            size  = 4, stroke = 0.5)
      )
    ) +
    
    ## shape legend
    ggplot2::scale_shape_manual(
      values = shape_pal,
      name   = "Dataset",
      guide  = ggplot2::guide_legend(order = 2)
    ) +
    
    ## labels
    ggrepel::geom_label_repel(
      data          = label_data,
      ggplot2::aes(label = label_txt, fill = category, shape = NULL),
      colour        = "white",
      fontface      = "bold",
      size          = label_size_mm,
      box.padding   = 0.6,
      point.padding = 0.3,
      segment.color = "gray40",
      segment.size  = 0.4,
      show.legend   = FALSE,
      max.overlaps  = Inf
    ) +
    
    ## grid
    ggplot2::geom_hline(yintercept = y_breaks, linetype = "solid",
                        colour = "grey85", size = 0.3) +
    ggplot2::geom_vline(xintercept = x_breaks, linetype = "solid",
                        colour = "grey85", size = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "darkgrey", size = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "darkgrey", size = 0.8) +
    
    ## axes
    ggplot2::scale_x_continuous(
      name   = "Gene ratio in positive-DEG ORA",
      limits = c(0, x_max + axis_pad),
      breaks = x_breaks,
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(
      name   = "Gene ratio in negative-DEG ORA",
      limits = c(0, y_max + axis_pad),
      breaks = y_breaks,
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    ## titles & theme
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal(base_size = big_txt) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold",
                                              size = big_title),
      axis.title      = ggplot2::element_text(size = big_axis),
      axis.text       = ggplot2::element_text(size = big_txt),
      legend.title    = ggplot2::element_text(size = big_legend),
      legend.text     = ggplot2::element_text(size = big_txt),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box      = "horizontal"
    )
  
  return(p)
}





annotate_with_positive_ORA <- function(always_true_df, ora_tbl) {
  # ora_tbl: e.g. positive_ora_df$meta_ora_top
  
  stopifnot(all(c("geneID","Description","GeneRatio") %in% names(ora_tbl)))
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    BiocManager::install("org.Hs.eg.db")
  
  # 1. Split all Entrez IDs appearing in ORA
  entrez_vec <- unique(unlist(strsplit(ora_tbl$geneID, "/", fixed = TRUE), use.names = FALSE))
  entrez_vec <- entrez_vec[nzchar(entrez_vec)]
  
  # 2. Map Entrez → SYMBOL (reuse logic style from provided functions)
  entrez2symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys     = entrez_vec,
    column   = "SYMBOL",
    keytype  = "ENTREZID",
    multiVals = "first"
  )
  
  # 3. Long table: one row per (SYMBOL, pathway)
  ora_long <- ora_tbl %>%
    dplyr::transmute(
      ID,
      Description,
      GeneRatio,
      p.adjust = .data[["p.adjust"]],
      geneID
    ) %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::mutate(
      gene = unname(entrez2symbol[geneID])
    ) %>%
    dplyr::filter(!is.na(gene) & nzchar(gene))
  
  # 4. For each gene pick pathway with highest GeneRatio (tie‑break: lowest p.adjust, then alphabetical)
  best_path_per_gene <- ora_long %>%
    dplyr::arrange(dplyr::desc(GeneRatio), p.adjust, Description) %>%
    dplyr::group_by(gene) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(gene,
                  ORA_positive_name = Description,
                  ORA_positive_ID   = ID,
                  ORA_positive_GeneRatio = GeneRatio)
  
  # 5. Join onto original
  always_true_df %>%
    dplyr::left_join(best_path_per_gene, by = "gene") %>%
    dplyr::mutate(positive_ORA = !is.na(ORA_positive_name))
}







build_metric_df_for_fold_change <- function(
    obj_list,
    id_col        = NULL,
    metric_col,
    base_map      = c(
      "NK"          = "NKs",
      "monocytes"   = "Monocytes",
      "B_cells"     = "B Cells",
      "CD4_T_cells" = "CD4+ T Cells",
      "CD8_T_cells" = "CD8+ T Cells"
    ),
    cell_order    = c("NKs", "Monocytes", "B Cells",
                      "CD4+ T Cells", "CD8+ T Cells"),
    filter_nocoding = FALSE
) {
  ## ── 0 · infer id column ───────────────────────────────────────────────────
  if (is.null(id_col)) {
    sample_names <- names(obj_list[[1L]])
    if      ("ID"   %in% sample_names) id_col <- "ID"
    else if ("gene" %in% sample_names) id_col <- "gene"
    else stop("Couldn’t infer `id_col`. Please supply it explicitly.")
  }
  
  ## ── 1 · long table (keep metric instead of presence) ──────────────────────
  long_tbl <- purrr::imap_dfr(obj_list, function(df, nm) {
    base <- nm %>%
      stringr::str_remove("_gsea_(positive|negative)$") %>%
      stringr::str_remove("_meta_table$")
    if (!base %in% names(base_map)) return(NULL)
    
    keep_cols <- c(id_col, metric_col)
    keep_cols <- keep_cols[keep_cols %in% names(df)]
    
    df %>%
      dplyr::select(dplyr::all_of(keep_cols)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(cell = base_map[base])
  })
  
  ## ── 2 · widen: values are metric_col ──────────────────────────────────────
  wide_tbl <- long_tbl %>%
    tidyr::pivot_wider(
      names_from  = cell,
      values_from = !!rlang::sym(metric_col)
    )
  
  ## ── 3 · keep only genes present in *all* lists (all cell cols non‑NA) ─────
  cell_cols <- intersect(cell_order, names(wide_tbl))
  wide_tbl  <- wide_tbl %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(cell_cols), ~ !is.na(.)))
  
  ## ── 4 · column ordering ───────────────────────────────────────────────────
  wide_tbl <- wide_tbl %>%
    dplyr::relocate(dplyr::all_of(id_col))
  
  wide_tbl <- wide_tbl %>%
    dplyr::relocate(dplyr::all_of(cell_cols),
                    .after = dplyr::last_col())
  
  ## ── 5 · optional: filter to protein‑coding genes ──────────────────────────
  if (isTRUE(filter_nocoding) && identical(id_col, "gene") && "gene" %in% names(wide_tbl)) {
    coding_syms <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = unique(wide_tbl$gene),
        mart       = mart
      )
      annot %>%
        dplyr::filter(.data$gene_biotype == "protein_coding") %>%
        dplyr::pull(.data$hgnc_symbol) %>%
        unique()
    }, error = function(e) {
      warning("Protein-coding filter skipped (biomaRt error: ", conditionMessage(e), ")")
      NULL
    })
    
    if (!is.null(coding_syms)) {
      wide_tbl <- dplyr::filter(wide_tbl, .data$gene %in% coding_syms)
    }
  }
  
  wide_tbl
}







get_overlap_genes_significance <- function(
    obj_list          = meta_table_list,
    id_col            = "gene",
    metric_col        = "median_log2FC",
    base_map          = c(
      "NK"          = "NKs",
      "monocytes"   = "Monocytes",
      "B_cells"     = "B Cells",
      "CD4_T_cells" = "CD4+ T Cells",
      "CD8_T_cells" = "CD8+ T Cells"
    ),
    cell_order        = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
    filter_nocoding   = TRUE,
    list_of_overlaps  = list(
      "All Five"             = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude NKs"          = c("Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude Monocytes"    = c("NKs", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude B Cells"      = c("NKs", "Monocytes", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude CD4+ T Cells" = c("NKs", "Monocytes", "B Cells", "CD8+ T Cells"),
      "Exclude CD8+ T Cells" = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells")
    ),
    thresholds_list   = list(
      "Positive DE" = list(median_log2FC = 1.25, q_val_fisher = 0.05, Positive = TRUE),
      "Negative DE" = list(median_log2FC = 1.25, q_val_fisher = 0.05, Positive = FALSE)
    )
){
  ## ---- packages --------------------------------------------------------------
  for (pkg in c("dplyr","purrr","stringr","tidyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required.")
  }
  
  ## ---- sanity checks ---------------------------------------------------------
  stopifnot(is.list(obj_list), length(obj_list) > 0)
  needed_cols <- c("median_log2FC","q_val_fisher")
  if (!all(purrr::map_lgl(obj_list, ~all(needed_cols %in% names(.x))))) {
    stop("Every data frame in obj_list must contain: ", paste(needed_cols, collapse=", "))
  }
  all_cells_used <- unique(unlist(list_of_overlaps))
  if (!all(all_cells_used %in% cell_order)) {
    stop("Some cells in list_of_overlaps are not in cell_order/base_map target names.")
  }
  check_thr <- function(thr){
    stopifnot(
      is.list(thr),
      all(c("median_log2FC","q_val_fisher","Positive") %in% names(thr)),
      is.numeric(thr$median_log2FC), length(thr$median_log2FC)==1, thr$median_log2FC > 0,
      is.numeric(thr$q_val_fisher),  length(thr$q_val_fisher)==1,  thr$q_val_fisher  > 0,
      is.logical(thr$Positive),      length(thr$Positive)==1
    )
  }
  invisible(purrr::walk(thresholds_list, check_thr))
  
  ## ---- 1. rho calculation ----------------------------------------------------
  rho_metric <- build_metric_df_for_fold_change(
    obj_list        = obj_list,
    id_col          = id_col,
    metric_col      = metric_col,
    base_map        = base_map,
    cell_order      = cell_order,
    filter_nocoding = filter_nocoding
  )
  if (!id_col %in% names(rho_metric)) stop("rho_metric missing id_col '", id_col,"'.")
  S <- rho_metric %>% dplyr::select(-!!id_col) %>% as.data.frame()
  if (ncol(S) < 2) stop("Need at least 2 cell-type columns to compute correlations.")
  
  pairs <- combn(colnames(S), 2, simplify = FALSE)
  rho_df <- purrr::map_dfr(
    pairs,
    ~{
      a <- .x[1]; b <- .x[2]
      r_s <- suppressWarnings(stats::cor(S[[a]], S[[b]],
                                         method = "spearman",
                                         use = "pairwise.complete.obs"))
      dplyr::tibble(cell1 = a, cell2 = b, r_s = r_s, rho = 1 + r_s)
    }
  )
  geom_mean_rho <- function(rho_vec, na.rm = TRUE) {
    if (any(rho_vec <= 0, na.rm = na.rm))
      stop("All rho values must be > 0 for a geometric mean.")
    exp(mean(log(rho_vec), na.rm = na.rm))
  }
  rho_list <- purrr::imap(
    list_of_overlaps,
    ~{
      subset_cells <- .x
      rho_vec <- rho_df %>%
        dplyr::filter(cell1 %in% subset_cells, cell2 %in% subset_cells) %>%
        dplyr::pull(rho)
      geom_mean_rho(rho_vec)
    }
  )
  
  ## ---- 2. Filtered obj_list & intersections ----------------------------------
  filter_one_df <- function(df, thr){
    fc_ok <- if (thr$Positive) df$median_log2FC >=  thr$median_log2FC
    else                       df$median_log2FC <= -thr$median_log2FC
    q_ok  <- df$q_val_fisher <= thr$q_val_fisher
    dplyr::filter(df, fc_ok & q_ok)
  }
  filtered_obj_list <- purrr::imap(thresholds_list, ~{
    thr <- .x
    purrr::map(obj_list, filter_one_df, thr = thr)
  })
  
  filtered_intersection_list <- purrr::imap(
    filtered_obj_list,
    ~ build_presence_df(
      obj_list        = .x,
      id_col          = id_col,
      filter_nocoding = filter_nocoding,
      base_map        = base_map,
      cell_order      = cell_order
    )
  )
  
  rename_inner <- function(x){
    old_core <- stringr::str_remove(names(x), "_meta_table$")
    new_nm   <- base_map[old_core]
    if (any(is.na(new_nm)))
      stop("Some inner names not found in base_map.")
    stats::setNames(x, new_nm)
  }
  filtered_obj_list <- lapply(filtered_obj_list, rename_inner)
  
  ## ---- 4. E0_list ------------------------------------------------------------
  N <- nrow(rho_metric)  # universe size
  
  compute_E0 <- function(cell_subset, cell_df_list, N) {
    stopifnot(all(cell_subset %in% names(cell_df_list)))
    k <- vapply(cell_subset, function(nm) nrow(cell_df_list[[nm]]), integer(1))
    prod(k) / N^(length(k) - 1)
  }
  
  E0_list <- purrr::imap(filtered_obj_list, function(cell_df_list, thr_name) {
    purrr::imap(list_of_overlaps, function(cells, overlap_name) {
      compute_E0(cells, cell_df_list, N)
    })
  })
  
  ## ---- 5. x_list (observed overlaps) -----------------------------------------
  count_overlap <- function(df, cols) {
    stopifnot(all(cols %in% names(df)))
    sum(rowSums(df[cols]) == length(cols), na.rm = TRUE)
  }
  
  x_list <- purrr::imap(filtered_intersection_list, function(df_top, top_name) {
    purrr::imap(list_of_overlaps, function(cols, sub_name) count_overlap(df_top, cols))
  })
  
  ## ---- 6. EF_list with stats -------------------------------------------------
  stopifnot(all(names(x_list) == names(E0_list)))
  stopifnot(all(sort(names(rho_list)) == sort(names(list_of_overlaps))))
  
  all_pvals <- numeric(0)
  idx_map   <- list()
  counter   <- 0
  
  EF_list <- purrr::imap(x_list, function(x_sub, top_nm) {
    purrr::imap(x_sub, function(x_val, ov_nm) {
      mu  <- E0_list[[top_nm]][[ov_nm]] * rho_list[[ov_nm]]
      if (mu <= 0) stop("Non-positive expected count.")
      ef_val <- x_val / mu
      
      # exact Poisson test on rate ratio (EF)
      pt <- stats::poisson.test(x_val, T = mu, r = 1, alternative = "greater")
      
      counter <<- counter + 1
      all_pvals[counter] <<- pt$p.value
      idx_map[[counter]] <<- c(top_nm, ov_nm)
      
      list(
        ef     = ef_val,
        p_val  = pt$p.value,
        q_val  = NA_real_,
        CI_bot = pt$conf.int[1],
        CI_top = pt$conf.int[2]
      )
    })
  })
  
  # BH adjust
  q_adj <- stats::p.adjust(all_pvals, method = "BH")
  for (i in seq_along(q_adj)) {
    nm1 <- idx_map[[i]][1]
    nm2 <- idx_map[[i]][2]
    EF_list[[nm1]][[nm2]]$q_val <- q_adj[i]
  }
  
  ## ---- return ----------------------------------------------------------------
  list(
    EF_list                    = EF_list,
    x_list                     = x_list,
    E0_list                    = E0_list,
    rho_list                   = rho_list,
    rho_df                     = rho_df,
    filtered_obj_list          = filtered_obj_list,
    filtered_intersection_list = filtered_intersection_list,
    rho_metric                 = rho_metric,
    universe_size              = N
  )
}




# ────────────────────────────────────
# 1. helper  ─ get_intersection_matrix
# ────────────────────────────────────
get_intersection_matrix <- function(
    ord,
    list_of_overlaps,
    cell_order_mat,
    base_size         = 12,
    mat_row_text_size = 12
) {
  ## ── sanity checks ────────────────────────────────────────────────────────
  req_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2")
  miss <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss))
    stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  
  if (!all(names(list_of_overlaps) %in% ord))
    stop("All names(list_of_overlaps) must be present in 'ord'.", call. = FALSE)
  
  if (!all(cell_order_mat %in% unlist(list_of_overlaps)))
    stop("All entries in cell_order_mat must appear in list_of_overlaps.", call. = FALSE)
  
  ## ── build long matrix ────────────────────────────────────────────────────
  cells  <- cell_order_mat
  y_map  <- tibble::tibble(cell = cells,
                           y_num = length(cells):1)
  
  mat_long <- tibble::tibble(name = names(list_of_overlaps)) %>%
    dplyr::mutate(cells_in = list_of_overlaps[name]) %>%
    tidyr::unnest_longer(cells_in, values_to = "cell") %>%
    dplyr::mutate(present = TRUE) %>%
    dplyr::right_join(
      tidyr::expand_grid(name = names(list_of_overlaps), cell = cells),
      by = c("name", "cell")
    ) %>%
    dplyr::mutate(
      present = !is.na(present),
      name    = factor(name, levels = ord)
    ) %>%
    dplyr::left_join(y_map, by = "cell")
  
  seg_df <- mat_long %>%
    dplyr::filter(present) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      ymin = min(y_num),
      ymax = max(y_num),
      .groups = "drop"
    )
  
  ## ── plot  ────────────────────────────────────────────────────────────────
  ggplot2::ggplot(mat_long, ggplot2::aes(x = name, y = y_num)) +
    ggplot2::geom_segment(
      data        = seg_df,
      ggplot2::aes(x = name, xend = name, y = ymin, yend = ymax),
      inherit.aes = FALSE,
      size        = 0.4,
      colour      = "black"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = present),
      shape   = 21,
      size    = 2.8,
      colour  = "black",
      stroke  = 0.25
    ) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "black", `FALSE` = "grey85"),
      guide  = "none"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_map$y_num,
      labels = y_map$cell
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = mat_row_text_size),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 0)
    )
}

# ──────────────────────────────────────────────────────────────
# 2. main  ─ get_overlap_genes_significance_barplot (updated)
# ──────────────────────────────────────────────────────────────
get_overlap_genes_significance_barplot <- function(
    barplot_list,
    list_of_overlaps = list(
      "All Five"             = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude NKs"          = c("Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude Monocytes"    = c("NKs", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude B Cells"      = c("NKs", "Monocytes", "CD4+ T Cells", "CD8+ T Cells"),
      "Exclude CD4+ T Cells" = c("NKs", "Monocytes", "B Cells", "CD8+ T Cells"),
      "Exclude CD8+ T Cells" = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells")
    ),
    cell_order_mat   = c("Monocytes", "CD4+ T Cells", "CD8+ T Cells", "B Cells", "NKs"),
    my_dir_cols      = c("Positive DEGs" = "#219ebc", "Negative DEGs" = "#fb8500"),
    title            = "Dependence‑Adjusted Effect Sizes for Overlapping DE Protein-Coding Genes",
    title_size        = 12,
    y_title_size      = 12,
    y_text_size       = 12,
    mat_row_text_size = 12,
    legend_title_size = 12,
    legend_text_size  = 12,
    base_size         = 12,
    heights           = c(3, 1.3)
) {
  ## ── sanity checks (unchanged) ─────────────────────────────────────────────
  req_pkgs <- c("dplyr", "forcats", "tidyr", "ggplot2", "patchwork")
  miss <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  
  if (!all(c("Positive DE","Negative DE") %in% names(barplot_list)))
    stop("barplot_list must contain elements named 'Positive DE' and 'Negative DE'.", call. = FALSE)
  
  needed_cols <- c("name","log2ef")
  if (!all(needed_cols %in% names(barplot_list$`Positive DE`)))
    stop("barplot_list$`Positive DE` must have columns: ", paste(needed_cols, collapse=", "), call. = FALSE)
  if (!all(needed_cols %in% names(barplot_list$`Negative DE`)))
    stop("barplot_list$`Negative DE` must have columns: ", paste(needed_cols, collapse=", "), call. = FALSE)
  
  if (!all(names(list_of_overlaps) %in% unique(barplot_list$`Positive DE`$name)))
    stop("Names of list_of_overlaps must match the 'name' column of barplot_list dfs.", call. = FALSE)
  
  if (!all(cell_order_mat %in% unlist(list_of_overlaps)))
    stop("All entries in cell_order_mat must appear in list_of_overlaps.", call. = FALSE)
  
  if (length(my_dir_cols) != 2 || any(names(my_dir_cols) == ""))
    stop("my_dir_cols must be a named vector of length 2.", call. = FALSE)
  
  ## ── BAR DATA ──────────────────────────────────────────────────────────────
  pos <- barplot_list$`Positive DE` %>%
    dplyr::select(name, log2ef) %>%
    dplyr::mutate(direction = "Positive DEGs")
  
  neg <- barplot_list$`Negative DE` %>%
    dplyr::select(name, log2ef) %>%
    dplyr::mutate(direction = "Negative DEGs")
  
  ord <- neg %>% dplyr::arrange(dplyr::desc(log2ef)) %>% dplyr::pull(name)
  
  pdat <- dplyr::bind_rows(pos, neg) %>%
    dplyr::mutate(name = factor(name, levels = ord))
  
  ## ── BAR PLOT ──────────────────────────────────────────────────────────────
  pdodge <- ggplot2::position_dodge(width = 0.65)
  
  p_bar <- ggplot2::ggplot(pdat, ggplot2::aes(x = name, y = log2ef, fill = direction)) +
    ggplot2::geom_col(width = 0.6, position = pdodge) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(values = my_dir_cols, name = "Overlap Group") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      x = NULL,
      y = expression(log[2] * " EF"),
      title = title
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.title.y       = ggplot2::element_text(size = y_title_size),
      axis.text.y        = ggplot2::element_text(size = y_text_size),
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      legend.title       = ggplot2::element_text(size = legend_title_size),
      legend.text        = ggplot2::element_text(size = legend_text_size),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  
  ## ── INTERSECTION MATRIX (outsourced) ──────────────────────────────────────
  p_mat <- get_intersection_matrix(
    ord              = ord,
    list_of_overlaps = list_of_overlaps,
    cell_order_mat   = cell_order_mat,
    base_size        = base_size,
    mat_row_text_size = mat_row_text_size
  )
  
  ## ── COMBINE & RETURN ──────────────────────────────────────────────────────
  p_bar / p_mat + patchwork::plot_layout(heights = heights)
}







get_intersection_composition_df <- function(
    df,
    cell_cols   = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
    group_col   = "group_go",
    groups_map  = c("Immune Processes"             = "Immune_percent",
                    "Metabolism"                   = "Metabolism_percent",
                    "Stress and Genome Regulation" = "Stress_percent",
                    "Other / No Pathway"           = "Other_percent"),
    drop_empty_intersections = TRUE
) {
  ## ── sanity checks ──────────────────────────────────────────────────────────
  stopifnot(is.data.frame(df))
  miss <- setdiff(cell_cols, names(df))
  if (length(miss)) stop("Missing cell columns: ", paste(miss, collapse = ", "))
  if (!group_col %in% names(df)) stop("group_col not in df.")
  if (!all(vapply(df[cell_cols], is.logical, logical(1))))
    stop("All cell_cols must be logical (TRUE/FALSE).")
  
  ## ── build ALL combinations (inclusive intersections) ──────────────────────
  combo_list <- unlist(lapply(seq_along(cell_cols), function(k) {
    combn(cell_cols, k, simplify = FALSE)
  }), recursive = FALSE)
  
  all_groups <- names(groups_map)
  
  ## ── for each combo: count group_col, compute % ────────────────────────────
  combo_tabs <- lapply(combo_list, function(combo) {
    idx <- rowSums(df[, combo, drop = FALSE]) == length(combo)   # inclusive: all combo TRUE
    sub_df <- df[idx, , drop = FALSE]
    
    # If empty, either skip or create zero rows
    if (nrow(sub_df) == 0) {
      if (drop_empty_intersections) return(NULL)
      tmp <- tibble(
        .intersection = paste(combo, collapse = "|"),
        !!group_col   := all_groups,
        n             = 0
      )
    } else {
      tmp <- sub_df %>%
        dplyr::count(!!rlang::sym(group_col), name = "n") %>%
        tidyr::complete(!!rlang::sym(group_col) := all_groups, fill = list(n = 0)) %>%
        dplyr::mutate(.intersection = paste(combo, collapse = "|"))
    }
    
    tmp %>%
      dplyr::group_by(.intersection) %>%
      dplyr::mutate(total = sum(n),
                    percent = ifelse(total == 0, 0, 100 * n / total)) %>%
      dplyr::ungroup()
  })
  
  combo_tabs <- combo_tabs[!vapply(combo_tabs, is.null, logical(1))]
  if (length(combo_tabs) == 0) return(tibble())
  
  sum_tab <- dplyr::bind_rows(combo_tabs)
  
  perc_tab <- sum_tab %>%
    dplyr::select(.intersection, !!rlang::sym(group_col), percent) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(group_col),
                       values_from = percent,
                       values_fill = 0)
  
  ## ── rename percentage columns ─────────────────────────────────────────────
  for (g in names(groups_map)) {
    if (g %in% names(perc_tab)) {
      names(perc_tab)[names(perc_tab) == g] <- groups_map[[g]]
    }
  }
  
  ## ── boolean flags per cell type ───────────────────────────────────────────
  parts   <- strsplit(perc_tab$.intersection, "\\|")
  bool_df <- as.data.frame(do.call(rbind, lapply(parts, function(x) cell_cols %in% x)))
  names(bool_df) <- cell_cols
  
  ## ── bind, add intersection_num, order ─────────────────────────────────────
  out <- dplyr::bind_cols(bool_df, perc_tab[, setdiff(names(perc_tab), ".intersection")])
  
  out <- out %>%
    dplyr::mutate(intersection_num = rowSums(dplyr::across(all_of(cell_cols)))) %>%
    dplyr::arrange(dplyr::desc(intersection_num), dplyr::across(all_of(cell_cols)))
  
  out
}




get_intersection_composition_stack_barplot <- function(
    df,
    cell_cols    = c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells"),
    percent_cols = c("Immune_percent","Metabolism_percent",
                     "Proteins_Metabolism_percent","RNA_Metabolism_percent",
                     "Other_percent"),
    group_name_map   = NULL,    # named vector: old_name = "Pretty name"
    group_palette    = NULL,    # named vector: "Pretty name" = "#hex"
    legend_title     = "Category",
    plot_title       = NULL,    # title for the stacked barplot
    x_title          = "Intersection",
    y_title          = "Percent",
    x_angle          = 45,
    x_hjust          = 1,
    x_order          = NULL,    # optional vector to set the order of the bars
    group_order      = NULL,    # character vector: stacking order of bars
    legend_order     = NULL,    # character vector: legend order ONLY
    ## ── matrix args ──────────────────────────────────────────────────────────
    list_of_overlaps = NULL,    # named list: names = column titles, values = character vectors of cells
    cell_order_mat   = c("Monocytes","CD4+ T Cells","CD8+ T Cells","B Cells","NKs"),
    heights          = c(3, 1.3),  # relative heights barplot / matrix
    base_size        = 12,
    mat_row_text_size = 12,
    ## ── reference line ───────────────────────────────────────────────────────
    y_line           = NULL,    # if numeric, draw dashed line at this y and include it in y-axis breaks
    ## ── style sizes to match get_overlap_genes_significance_barplot ──────────
    title_size        = 16,
    y_title_size      = 14,
    y_text_size       = 14,
    legend_title_size = 14,
    legend_text_size  = 14,
    ## manual vertical shift for the legend (0 = top)
    legend_vshift     = 0.25
){
  ## ── sanity checks ──────────────────────────────────────────────────────────
  stopifnot(is.data.frame(df))
  miss_cells   <- setdiff(cell_cols, names(df))
  miss_percent <- setdiff(percent_cols, names(df))
  if (length(miss_cells))   stop("Missing cell columns: ",   paste(miss_cells, collapse = ", "))
  if (length(miss_percent)) stop("Missing percent columns: ",paste(miss_percent, collapse = ", "))
  
  ## ── readable intersection labels ───────────────────────────────────────────
  intersection_lbl <- apply(df[cell_cols], 1, function(x){
    paste(cell_cols[as.logical(x)], collapse = " ∩ ")
  })
  df <- df %>% dplyr::mutate(.intersection = intersection_lbl)
  
  ## ── long format ────────────────────────────────────────────────────────────
  long_df <- df %>% 
    dplyr::select(.intersection, dplyr::all_of(percent_cols)) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(percent_cols),
      names_to  = "group",
      values_to = "percent"
    )
  
  ## pretty names
  if (!is.null(group_name_map)) {
    stopifnot(all(names(group_name_map) %in% unique(long_df$group)))
    long_df$group <- dplyr::recode(long_df$group, !!!group_name_map)
  }
  
  ## stacking order for bars only
  if (!is.null(group_order)) {
    stopifnot(setequal(group_order, unique(long_df$group)))
    long_df$group <- factor(long_df$group, levels = group_order)
  } else {
    long_df$group <- factor(long_df$group)
  }
  groups_chr <- as.character(long_df$group)
  groups_unique <- unique(groups_chr)
  
  ## legend order validation (independent)
  if (!is.null(legend_order)) {
    stopifnot(setequal(legend_order, groups_unique))
  }
  
  ## x‑axis order
  if (is.null(x_order)) {
    x_order <- long_df %>% 
      dplyr::mutate(set_count = stringr::str_count(.intersection, "∩") + 1) %>% 
      dplyr::distinct(.intersection, set_count) %>% 
      dplyr::arrange(dplyr::desc(set_count), .intersection) %>% 
      dplyr::pull(.intersection)
  }
  long_df$.intersection <- factor(long_df$.intersection, levels = x_order)
  
  ## list_of_overlaps
  if (is.null(list_of_overlaps)) {
    list_of_overlaps <- lapply(split(df[cell_cols], df$.intersection), function(sub){
      cell_cols[as.logical(sub[1, ])]
    })
  }
  list_of_overlaps <- list_of_overlaps[x_order]
  
  ## y‑axis breaks — always include 100 and append '%'
  max_y <- max(
    long_df$percent,
    100,
    if (!is.null(y_line) && is.numeric(y_line) && is.finite(y_line)) y_line else -Inf,
    na.rm = TRUE
  )
  base_breaks <- seq(0, ceiling(max_y / 10) * 10, by = 10)
  if (!is.null(y_line) && is.numeric(y_line) && is.finite(y_line)) {
    y_breaks <- sort(unique(c(base_breaks, y_line)))
  } else {
    y_breaks <- base_breaks
  }
  
  ## ── stacked bar plot -------------------------------------------------------
  p_bar_core <- ggplot2::ggplot(long_df, ggplot2::aes(x = .intersection, y = percent, fill = group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_y_continuous(breaks = y_breaks, labels = function(b) paste0(b, "%")) +
    ggplot2::labs(
      title = plot_title,
      x     = x_title,
      y     = y_title,
      fill  = legend_title
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "panel",
      plot.title          = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                  size = title_size,
                                                  margin = ggplot2::margin(t = 10, b = 4)),
      axis.title.x        = ggplot2::element_text(size = y_title_size),
      axis.title.y        = ggplot2::element_text(size = y_title_size),
      axis.text.y         = ggplot2::element_text(size = y_text_size),
      axis.text.x         = ggplot2::element_blank(),
      axis.ticks.x        = ggplot2::element_blank(),
      legend.title        = ggplot2::element_text(size = legend_title_size),
      legend.text         = ggplot2::element_text(size = legend_text_size),
      panel.grid.major.x  = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      plot.margin         = ggplot2::margin(t = 6, r = 6, b = 4, l = 6)
    )
  if (!is.null(y_line) && is.numeric(y_line) && length(y_line) == 1 && is.finite(y_line)) {
    p_bar_core <- p_bar_core + ggplot2::geom_hline(yintercept = y_line, linetype = "dashed")
  }
  
  ## palette + legend order (independent from stacking)
  if (!is.null(group_palette)) {
    stopifnot(all(groups_unique %in% names(group_palette)))
    # legend breaks: legend_order if given, else previous default (reversed stacking)
    legend_breaks <- if (!is.null(legend_order)) legend_order else rev(levels(long_df$group))
    p_bar_core <- p_bar_core + ggplot2::scale_fill_manual(
      values = group_palette,
      breaks = legend_breaks
    )
  } else {
    if (!is.null(legend_order)) {
      p_bar_core <- p_bar_core + ggplot2::scale_fill_discrete(breaks = legend_order)
    } else {
      p_bar_core <- p_bar_core + ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))
    }
  }
  
  ## legend as separate grob
  if (!requireNamespace("cowplot", quietly = TRUE))
    stop("Package 'cowplot' is required for assembling the plots.", call. = FALSE)
  legend_g <- cowplot::get_legend(
    p_bar_core +
      ggplot2::theme(
        legend.position      = "right",
        legend.justification = c(0, 1),
        legend.box.margin    = ggplot2::margin(t = 0, r = 10, b = 0, l = 10),
        plot.margin          = ggplot2::margin(0, 0, 0, 0)
      )
  )
  p_bar <- p_bar_core + ggplot2::theme(legend.position = "none")
  
  ## ── intersection matrix ----------------------------------------------------
  matrix_plot <- get_intersection_matrix(
    ord               = x_order,
    list_of_overlaps  = list_of_overlaps,
    cell_order_mat    = cell_order_mat,
    base_size         = base_size,
    mat_row_text_size = mat_row_text_size
  ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 2, r = 6, b = 2, l = 6)
    )
  
  ## stack bar + matrix (same width)
  stack_vm <- cowplot::plot_grid(
    p_bar,
    matrix_plot,
    ncol        = 1,
    align       = "v",
    rel_heights = heights
  )
  
  ## vertical legend shift
  vshift <- max(0, min(1, legend_vshift))
  y_pos  <- 1 - vshift
  legend_box <- cowplot::ggdraw() +
    cowplot::draw_plot(legend_g, x = 0, y = y_pos, width = 1, hjust = 0, vjust = 1)
  
  ## combine
  combined_plot <- cowplot::plot_grid(
    stack_vm,
    legend_box,
    ncol        = 2,
    rel_widths  = c(1, 0.28),
    align       = "h",
    axis        = "t"
  )
  
  list(
    plot      = combined_plot,
    bar_plot  = p_bar,
    mat_plot  = matrix_plot,
    legend    = legend_g,
    data_long = long_df
  )
}







#───────────────────────────────────────────────────────────────────────────────
#  get_output_signature_score()
#  --------------------------------
#  Compute a pair‑wise UCell‑style score for one FindMarkers data‑frame
#───────────────────────────────────────────────────────────────────────────────

get_output_signature_score <- function(output_df,
                                       positive_genes,
                                       negative_genes = NULL,
                                       pct.threshold = 0,
                                       w_neg         = 1,
                                       max_rank      = 1500,
                                       filter_nocoding = TRUE) {
  ## 0 ·──────── sanity checks ────────────────────────────────────────────────
  if (!is.data.frame(output_df))
    stop("`output_df` must be a data‑frame produced by Seurat::FindMarkers().")
  req_cols <- c("avg_log2FC", "pct.1", "pct.2")
  if (!all(req_cols %in% colnames(output_df)))
    stop("`output_df` is missing required columns: ",
         paste(setdiff(req_cols, colnames(output_df)), collapse = ", "))
  if (is.null(rownames(output_df)))
    stop("`output_df` must have gene symbols as row names.")
  if (missing(positive_genes) || length(positive_genes) == 0)
    stop("`positive_genes` must be a non‑empty character vector.")
  if (!is.null(negative_genes) && !is.character(negative_genes))
    stop("`negative_genes` must be NULL or a character vector.")
  if (!is.numeric(w_neg) || w_neg < 0)
    stop("`w_neg` must be a non‑negative numeric value.")
  if (!is.numeric(pct.threshold) ||
      pct.threshold < 0 ||          # allow zero
      pct.threshold >= 1)
    stop("`pct.threshold` must be in [0, 1).")
  if (!is.numeric(max_rank) || max_rank < 1)
    stop("`max_rank` must be a positive integer.")

    ## 1 ·──────── optional: filter to protein‑coding genes ─────────────────────
  if (isTRUE(filter_nocoding)) {
    coding_syms <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      mart  <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = unique(rownames(output_df)),
        mart       = mart
      )
      unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
    }, error = function(e) {
      warning("Protein‑coding filter skipped (biomaRt error: ",
              conditionMessage(e), ").")
      NULL
    })
    
    if (!is.null(coding_syms))
      output_df <- output_df[rownames(output_df) %in% coding_syms, , drop = FALSE]
  }
  
  # browser()
  
  ## 2 ·──────── assemble the gene universe ───────────────────────────────────
  sig_genes <- unique(c(positive_genes, negative_genes))
  universe  <- unique(c(rownames(output_df), sig_genes))
  
  # add missing genes as zero‑rows (so they get worst rank later)
  missing_genes <- setdiff(universe, rownames(output_df))
  if (length(missing_genes)) {
    pad <- data.frame(
      p_val      = NA_real_,
      avg_log2FC = 0,
      pct.1      = 0,
      pct.2      = 0,
      p_val_adj  = NA_real_,
      stringsAsFactors = FALSE
    )
    pad <- pad[rep(1, length(missing_genes)), , drop = FALSE]  # replicate row
    rownames(pad) <- missing_genes                              # assign names
    output_df <- rbind(output_df, pad)
  }
  
  ## 3 ·──────── flag low‑expression genes (< pct.threshold in BOTH groups) ───
  low_expr <- with(output_df,
                   pct.1 < pct.threshold & pct.2 < pct.threshold)
  
  ## 4 ·──────── build the ranked vector ──────────────────────────────────────
  logfc_vec <- output_df$avg_log2FC
  names(logfc_vec) <- rownames(output_df)
  
  # base rank: descending log2FC (highest up‑reg gets rank 1)
  ranks      <- rank(-logfc_vec, ties.method = "average")
  ranks[ranks > max_rank] <- max_rank            # clip long tail
  
  # override low‑expression genes with worst‑rank
  ranks[low_expr] <- max_rank
  
  ## 5 ·──────── UCell helper to score a gene set ─────────────────────────────
  u_score <- function(gset) {
    if (is.null(gset) || !length(gset)) return(0)
    r <- ranks[gset]
    L <- length(r)
    U <- sum(r) - L * (L + 1) / 2
    Umax <- L * max_rank - L * (L + 1) / 2
    1 - U / Umax
  }
  
  s_pos <- u_score(intersect(positive_genes, universe))
  s_neg <- u_score(intersect(negative_genes, universe))
  
  ## 6 ·──────── final mixed score (canonical formula) ────────────────────────
  score <- if (is.null(negative_genes) || length(negative_genes) == 0) {
    s_pos
  } else {
    max(0, s_pos - w_neg * s_neg)
  }
  
  return(score)
}


get_meta_table_signature_score <- function(meta_table,
                                           positive_genes,
                                           negative_genes = NULL,
                                           pct.threshold  = 0,
                                           w_neg          = 1,
                                           max_rank       = 1500,
                                           filter_nocoding = TRUE) {
  
  
  ## 0 ·──────── sanity checks ────────────────────────────────────────────────
  if (!is.data.frame(meta_table))
    stop("`meta_table` must be a data-frame.")
  
  if (is.null(rownames(meta_table))) {
    if (!"gene" %in% colnames(meta_table))
      stop("`meta_table` must have gene symbols as row names or a 'gene' column.")
    rownames(meta_table) <- as.character(meta_table$gene)
  }
  
  if (!"median_log2FC" %in% colnames(meta_table))
    stop("`meta_table` is missing required column: 'median_log2FC'.")
  if (!"pct.1" %in% colnames(meta_table))
    stop("`meta_table` is missing required column: 'pct.1'.")
  
  if (missing(positive_genes) || length(positive_genes) == 0)
    stop("`positive_genes` must be a non-empty character vector.")
  if (!is.null(negative_genes) && !is.character(negative_genes))
    stop("`negative_genes` must be NULL or a character vector.")
  if (!is.numeric(w_neg) || w_neg < 0)
    stop("`w_neg` must be a non-negative numeric value.")
  if (!is.numeric(pct.threshold) || pct.threshold < 0 || pct.threshold >= 1)
    stop("`pct.threshold` must be in [0, 1).")
  if (!is.numeric(max_rank) || max_rank < 1)
    stop("`max_rank` must be a positive integer.")
  
  ## 1 ·──────── optional: filter to protein-coding genes ─────────────────────
  if (isTRUE(filter_nocoding)) {
    coding_syms <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      mart  <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = unique(rownames(meta_table)),
        mart       = mart
      )
      unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
    }, error = function(e) {
      warning("Protein-coding filter skipped (biomaRt error: ",
              conditionMessage(e), ").")
      NULL
    })
    if (!is.null(coding_syms))
      meta_table <- meta_table[rownames(meta_table) %in% coding_syms, , drop = FALSE]
  }
  
  ## 2 ·──────── assemble the gene universe ───────────────────────────────────
  sig_genes <- unique(c(positive_genes, negative_genes))
  universe  <- unique(c(rownames(meta_table), sig_genes))
  
  # add missing genes as zero-rows (so they get worst rank later)
  missing_genes <- setdiff(universe, rownames(meta_table))
  if (length(missing_genes)) {
    pad <- as.data.frame(matrix(0, nrow = length(missing_genes), ncol = ncol(meta_table)))
    colnames(pad) <- colnames(meta_table)
    rownames(pad) <- missing_genes
    if ("gene" %in% colnames(pad)) pad$gene <- missing_genes
    meta_table <- rbind(meta_table, pad)
  }
  
  ## 3 ·──────── flag low-expression genes (< pct.threshold) ──────────────────
  low_expr <- meta_table$pct.1 < pct.threshold
  
  ## 4 ·──────── build the ranked vector (by median_log2FC) ───────────────────
  logfc_vec <- meta_table$median_log2FC
  names(logfc_vec) <- rownames(meta_table)
  
  ranks <- rank(-logfc_vec, ties.method = "average")
  ranks[ranks > max_rank] <- max_rank
  ranks[low_expr] <- max_rank
  
  ## 5 ·──────── UCell helper to score a gene set ─────────────────────────────
  u_score <- function(gset) {
    if (is.null(gset) || !length(gset)) return(0)
    r <- ranks[gset]
    L <- length(r)
    U <- sum(r) - L * (L + 1) / 2
    Umax <- L * max_rank - L * (L + 1) / 2
    1 - U / Umax
  }
  
  s_pos <- u_score(intersect(positive_genes, universe))
  s_neg <- u_score(intersect(negative_genes, universe))
  
  ## 6 ·──────── final mixed score ───────────────────────────────────────────
  if (is.null(negative_genes) || length(negative_genes) == 0) {
    score <- s_pos
  } else {
    score <- max(0, s_pos - w_neg * s_neg)
  }
  
  return(score)
}




#───────────────────────────────────────────────────────────────────────────────
#  get_signature_score_boxplot()
#  --------------------------------
#  Traverse a 4‑level output_lists object (cell_type → dataset → sample →
#  original_table) and compute the UCell‑style oxidative‑stress score for every
#  sample using get_output_signature_score().  Returns a tidy data.frame ready
#  for ggplot2::geom_boxplot().
#───────────────────────────────────────────────────────────────────────────────
#───────────────────────────────────────────────────────────────────────────────
#  get_signature_score_boxplot()
#  ---------------------------------------------------------------------------
#  Walk through a nested output_lists object:
#    • level-1  = cell type   (e.g. "NKs")
#    • level-2  = dataset     (e.g. "Hong")
#    • level-3  = sample      (e.g. "HC-1")
#    • slot     = $original_table  (Seurat::FindMarkers data-frame)
#  For every original_table it calls get_output_signature_score() and returns
#  a tidy data-frame ready for ggplot2::geom_boxplot().
#───────────────────────────────────────────────────────────────────────────────

get_pathway_genes_by_id <- function(path_id, species = "hsa") {
  # ── 0) deps & helpers ──────────────────────────────────────────────────
  needs <- c("AnnotationDbi", "org.Hs.eg.db")
  for (pkg in needs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
  to_symbols <- function(entrez_vec) {
    entrez_vec <- unique(as.character(entrez_vec))
    entrez_vec <- entrez_vec[nzchar(entrez_vec)]
    if (!length(entrez_vec)) return(character())
    sy <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys     = entrez_vec,
      column   = "SYMBOL",
      keytype  = "ENTREZID",
      multiVals = "first"
    )
    sy <- unname(sy[!is.na(sy) & nzchar(sy)])
    sort(unique(sy))
  }
  
  # ── 1) KEGG pathway (e.g., "hsa04668") ─────────────────────────────────
  if (grepl(sprintf("^%s\\d{5}$", species), path_id, ignore.case = TRUE)) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("KEGGREST", ask = FALSE, update = FALSE)
    }
    kg <- tryCatch(KEGGREST::keggGet(path_id)[[1]], error = function(e) NULL)
    if (is.null(kg) || is.null(kg$GENE)) return(character())
    raw <- kg$GENE
    entrez_ids <- raw[seq(1, length(raw), by = 2)]
    entrez_ids <- sub(";.*$", "", entrez_ids)
    return(to_symbols(entrez_ids))
  }
  
  # ── 2) Reactome ID (e.g., "R-HSA-168256") ──────────────────────────────
  if (grepl("^R-HSA-\\d+$", path_id)) {
    if (!requireNamespace("reactome.db", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("reactome.db", ask = FALSE, update = FALSE)
    }
    map <- reactome.db::reactomePATHID2EXTID
    if (!is.null(map[[path_id]])) {
      entrez_ids <- as.character(map[[path_id]])
      return(to_symbols(entrez_ids))
    } else {
      return(character())
    }
  }
  
  # ── 3) MSigDB set by gs_name (works for all categories) ────────────────
  #     e.g., "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_APOPTOSIS",
  #     "KEGG_TGF_BETA_SIGNALING_PATHWAY", "GO_APOPTOTIC_PROCESS", etc.
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    install.packages("msigdbr")
  }
  gs <- toupper(path_id)
  df <- msigdbr::msigdbr(species = "Homo sapiens")
  
  df <- df %>%
    dplyr::filter(gs_name == !!gs) %>%
    dplyr::pull(gene_symbol) %>%
    unique() 
  
  if (length(df) > 0) {
    return(df)
  }
  
  # ── 4) Unknown ID type ─────────────────────────────────────────────────
  character()
}



get_signature_score_boxplot_df <- function(output_lists,
                                        positive_genes,
                                        negative_genes   = NULL,
                                        pct.threshold    = 0,
                                        w_neg            = 1,
                                        max_rank         = 1500,
                                        filter_nocoding  = TRUE)
{
  ## sanity ------------------------------------------------------------
  if (!is.list(output_lists))
    stop("`output_lists` must be a nested list as described.")
  if (missing(positive_genes) || length(positive_genes) == 0)
    stop("`positive_genes` must be a non-empty character vector.")
  
  ## helper to recurse -------------------------------------------------
  collect_scores <- function(cell_type, dataset, sample, df) {
    
    score <- get_output_signature_score(
      output_df       = df,
      positive_genes  = positive_genes,
      negative_genes  = negative_genes,
      pct.threshold   = pct.threshold,
      w_neg           = w_neg,
      max_rank        = max_rank,
      filter_nocoding = filter_nocoding
    )
    data.frame(cell_type = cell_type,
               dataset    = dataset,
               sample     = sample,
               oxidative_stress_score = score,
               stringsAsFactors = FALSE)
  }
  
  ## main traversal ----------------------------------------------------
  res_list <- list()
  
  for (ct in names(output_lists)) {
    print(ct)
    cell_lvl <- output_lists[[ct]]
    if (!is.list(cell_lvl)) next
    for (ds in names(cell_lvl)) {
      print(ds)
      ds_lvl <- cell_lvl[[ds]]
      if (!is.list(ds_lvl)) next
      for (smpl in names(ds_lvl)) {
        print(smpl)
        smpl_lvl <- ds_lvl[[smpl]]
        if (!is.list(smpl_lvl) || is.null(smpl_lvl$original_table)) {
          warning(sprintf("Skipping %s / %s / %s (no original_table)", ct, ds, smpl))
          next
        }
        res_list[[length(res_list) + 1]] <-
          collect_scores(ct, ds, smpl, smpl_lvl$original_table)
      }
    }
  }
  
  ## bind & return -----------------------------------------------------
  if (length(res_list) == 0)
    stop("No valid original_table objects found in `output_lists`.")
  
  scores_df <- do.call(rbind, res_list)
  return(scores_df)
}




get_signature_score_boxplot <- function(
    output_lists,
    positive_genes,
    negative_genes   = NULL,
    pct.threshold    = 0,
    w_neg            = 1,
    max_rank         = 1500,
    filter_nocoding  = TRUE,
    score_col   = "oxidative_stress_score",
    cell_col    = "cell_type",
    dataset_col = "dataset",
    box_colors = c(
      "NKs"          = "#9EB7E0",
      "Monocytes"    = "#219ebc",
      "B cells"      = "#48CAE4",
      "CD4+ T cells" = "#ffb703",
      "CD8+ T cells" = "#fb8500"
    ),
    shape_vals = c(
      Hong     = 24,
      Kang     = 22,
      Ramesh   = 21,
      Mevorach = 23
    ),
    title              = "Oxidative Stress Signatures in Pairwise Comparisons",
    base_size          = 24,
    title_size         = 28,
    y_title_size       = 26,
    y_text_size        = 22,
    legend_title_size  = 22,
    legend_text_size   = 19
) {
  
  df <- get_signature_score_boxplot_df (output_lists,
                                        positive_genes,
                                        negative_genes,
                                        pct.threshold,
                                        w_neg ,
                                        max_rank,
                                        filter_nocoding)
  stopifnot(
    is.data.frame(df),
    all(c(cell_col, dataset_col, score_col) %in% names(df)),
    all(unique(df[[cell_col]]) %in% names(box_colors)),
    all(unique(df[[dataset_col]]) %in% names(shape_vals))
  )
  
  cell_order <- names(box_colors)
  plot_df <- df |>
    dplyr::filter(.data[[cell_col]] %in% cell_order) |>
    dplyr::mutate(
      !!cell_col := factor(.data[[cell_col]], levels = cell_order)
    )
  
  ## separate zero‑score rows and give them a tiny positive offset so they remain visible
  zero_df <- plot_df |>
    dplyr::filter(.data[[score_col]] == 0) |>
    dplyr::mutate(score_jitter = 0.02)   # baseline just above zero
  
  nonzero_df <- dplyr::filter(plot_df, .data[[score_col]] != 0)
  
  signature_boxplot <- ggplot2::ggplot(plot_df,
                                       ggplot2::aes(x = .data[[cell_col]],
                                                    y = .data[[score_col]],
                                                    group = .data[[cell_col]])) +
    ## 1) boxplots
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data[[cell_col]]),
      width         = 0.6,
      outlier.shape = NA,
      colour        = "black",
      lwd           = 0.4,
      alpha         = 0.85
    ) +
    ## 2a) non‑zero points (regular position jitter)
    ggplot2::geom_jitter(
      data   = nonzero_df,
      ggplot2::aes(shape = .data[[dataset_col]]),
      width  = 0.15,
      size   = 2.4,
      stroke = 0.25,
      alpha  = 1,
      colour = "black",
      fill   = "black"
    ) +
    ## 2b) zero‑score points (nudged upward, then lightly jittered)
    ggplot2::geom_jitter(
      data   = zero_df,
      ggplot2::aes(y = score_jitter, shape = .data[[dataset_col]]),
      width  = 0.15,
      height = 0.01,               # tiny vertical jitter *above* zero
      size   = 2.4,
      stroke = 0.25,
      alpha  = 1,
      colour = "black",
      fill   = "black"
    ) +
    ggplot2::scale_fill_manual(values = box_colors, guide = "none") +
    ggplot2::scale_shape_manual(values = shape_vals, name = "Dataset") +
    ggplot2::scale_y_continuous(
      name   = "Oxidative‑stress score",
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(x = NULL, title = title) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(
        hjust  = 0.5, face = "bold", size = title_size
      ),
      axis.title.y       = ggplot2::element_text(size = y_title_size),
      axis.text.y        = ggplot2::element_text(size = y_text_size),
      axis.text.x        = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1, size = y_text_size
      ),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.direction   = "horizontal",
      legend.justification = "center",
      legend.box.just    = "center",
      legend.title       = ggplot2::element_text(size = legend_title_size),
      legend.text        = ggplot2::element_text(size = legend_text_size),
      plot.margin        = ggplot2::margin(t = 6, r = 20, b = 12, l = 10)
    )
  
  list(signature_df = df,
       signature_boxplot = signature_boxplot)
}

get_signature_score_boxplot_by_pathway <- function(
    output_lists,
    path_id,
    # positive_genes,
    # negative_genes   = NULL,
    pct.threshold    = 0,
    w_neg            = 1,
    max_rank         = 1500,
    filter_nocoding  = TRUE,
    score_col   = "oxidative_stress_score",
    cell_col    = "cell_type",
    dataset_col = "dataset",
    box_colors = c(
      "NKs"          = "#9EB7E0",
      "Monocytes"    = "#219ebc",
      "B cells"      = "#48CAE4",
      "CD4+ T cells" = "#ffb703",
      "CD8+ T cells" = "#fb8500"
    ),
    shape_vals = c(
      Hong     = 24,
      Kang     = 22,
      Ramesh   = 21,
      Mevorach = 23
    ),
    title              = "Oxidative Stress Signatures in Pairwise Comparisons",
    base_size          = 24,
    title_size         = 28,
    y_title_size       = 26,
    y_text_size        = 22,
    legend_title_size  = 22,
    legend_text_size   = 19
) {
  
  positive_genes <- get_pathway_genes_by_id(path_id = path_id)
  
  return_list <- get_signature_score_boxplot(
    output_lists = output_lists,
    positive_genes = positive_genes,
    pct.threshold    = pct.threshold,
    w_neg            = w_neg,
    max_rank         = max_rank,
    filter_nocoding  = filter_nocoding,
    score_col   = score_col,
    cell_col    = cell_col,
    dataset_col = dataset_col,
    box_colors = box_colors,
    shape_vals = shape_vals,
    title              = title,
    base_size          = base_size,
    title_size         = title_size,
    y_title_size       = y_title_size,
    y_text_size        = y_text_size,
    legend_title_size  = legend_title_size,
    legend_text_size   = legend_text_size)
    
  return_list
  
  
  
}





















get_cellular_level_signature <- function(seu,
                                         positive_top_genes,
                                         negative_top_genes = NULL,
                                         filter_nocoding    = TRUE,
                                         assay              = "SCT",
                                         signature_name     = "os_signature"){
  # ── sanity ───────────────────────────────────────────────────────────────
  if (!inherits(seu, "Seurat")) stop("`seu` must be a Seurat object.")
  if (missing(positive_top_genes) || length(positive_top_genes) == 0)
    stop("`positive_top_genes` must be a non-empty character vector.")
  if (!assay %in% names(seu@assays))
    stop("Assay '", assay, "' not found in object.")
  if (!is.character(signature_name) || length(signature_name) != 1)
    stop("`signature_name` must be a single string.")
  
  # get feature universe (works for Seurat v5/v4)
  counts_mat <- tryCatch(
    GetAssayData(seu, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(seu, assay = assay, slot = "counts")
  )
  feat_universe <- rownames(counts_mat)
  
  # ── build signature (UCell negative genes use trailing '-') ───────────────
  gene_sig <- unique(c(positive_top_genes,
                       if (!is.null(negative_top_genes)) paste0(negative_top_genes, "-")))
  
  # optional: keep only protein-coding genes (filter the signature, not the assay)
  if (isTRUE(filter_nocoding)) {
    coding_syms <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      mart  <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = unique(sub("-$", "", gene_sig)),
        mart       = mart
      )
      unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
    }, error = function(e) {
      warning("Protein-coding filter skipped (biomaRt error: ",
              conditionMessage(e), ").")
      NULL
    })
    
    if (!is.null(coding_syms)) {
      base <- sub("-$", "", gene_sig)
      neg  <- grepl("-$", gene_sig)
      keep <- base %in% coding_syms
      dropped_nc <- setdiff(unique(base[!keep]), character(0))
      if (length(dropped_nc))
        message("Dropped non-coding genes from signature: ",
                paste(dropped_nc, collapse = ", "))
      gene_sig <- paste0(base[keep], ifelse(neg[keep], "-", ""))
    }
  }
  
  # intersect with features present in the object
  base <- sub("-$", "", gene_sig)
  neg  <- grepl("-$", gene_sig)
  keep_present <- base %in% feat_universe
  if (!any(keep_present))
    stop("None of the signature genes are present in assay '", assay, "'.")
  dropped_absent <- setdiff(unique(base[!keep_present]), character(0))
  if (length(dropped_absent))
    message("Dropped genes not found in object: ",
            paste(dropped_absent, collapse = ", "))
  gene_sig <- paste0(base[keep_present], ifelse(neg[keep_present], "-", ""))
  gene_sig <- unique(gene_sig)
  
  # wrap in a named list
  signatures <- setNames(list(gene_sig), signature_name)
  
  # score (column will be '<signature_name>_UCell' in meta.data)
  DefaultAssay(seu) <- assay
  if (!requireNamespace("UCell", quietly = TRUE))
    stop("Package 'UCell' is required for AddModuleScore_UCell().")
  seu <- UCell::AddModuleScore_UCell(
    seu,
    features = signatures,
    name     = NULL
  )
  
  return(seu)
}







get_pathway_genes_hugo <- function(pathway_id) {
  # Load needed packages
  if (!requireNamespace("KEGGREST", quietly = TRUE))
    stop("Package 'KEGGREST' is required for KEGG IDs.")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Package 'org.Hs.eg.db' is required for ID mapping.")
  
  # Helper to map Entrez → HUGO
  entrez_to_hugo <- function(entrez_ids) {
    suppressMessages({
      mapped <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                      keys = entrez_ids,
                                      keytype = "ENTREZID",
                                      columns = "SYMBOL")
    })
    unique(na.omit(mapped$SYMBOL))
  }
  
  # KEGG case
  if (grepl("^hsa\\d+$", pathway_id)) {
    kegg_data  <- KEGGREST::keggGet(pathway_id)[[1]]
    genes_field <- kegg_data$GENE
    
    entrez_ids <- character(0)
    
    if (!is.null(genes_field) && length(genes_field) > 0) {
      # Pull only the id tokens (lines that actually start with hsa:)
      id_tokens <- genes_field[grepl("^hsa:", genes_field)]
      if (length(id_tokens) > 0) {
        entrez_ids <- base::sub("^hsa:", "", as.character(id_tokens))
      }
    }
    
    # Hard fallback when $GENE is empty or had no id tokens
    if (length(entrez_ids) == 0) {
      links <- KEGGREST::keggLink("hsa", pathway_id)
      if (!is.null(links) && length(links) > 0) {
        entrez_ids <- base::sub("^hsa:", "", as.character(unname(links)))
      }
    }
    
    return(entrez_to_hugo(entrez_ids))
  }
  
  
  # Reactome case
  if (grepl("^R-HSA-\\d+$", pathway_id)) {
    if (!requireNamespace("ReactomeContentService4R", quietly = TRUE))
      stop("Package 'ReactomeContentService4R' is required for Reactome IDs.")
    reactome_genes <- ReactomeContentService4R::reactome_entities(pathway_id, species = "Homo sapiens")
    entrez_ids <- unique(reactome_genes$identifier[reactome_genes$databaseName == "EntrezGene"])
    return(entrez_to_hugo(entrez_ids))
  }
  
  stop("Unrecognized pathway ID format. Use KEGG (hsaXXXXX) or Reactome (R-HSA-XXXXX).")
}

get_pathway_genes <- function(pathway_id) {
  
  # browser()
  if (!requireNamespace("KEGGREST", quietly = TRUE))
    stop("Package 'KEGGREST' is required for KEGG IDs.")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Package 'org.Hs.eg.db' is required for ID mapping.")
  
  # Helper to map Entrez → HUGO
  entrez_to_hugo <- function(entrez_ids) {
    suppressMessages({
      mapped <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                      keys = entrez_ids,
                                      keytype = "ENTREZID",
                                      columns = "SYMBOL")
    })
    unique(na.omit(mapped$SYMBOL))
  }
  
  # KEGG case
  if (grepl("^hsa\\d+$", pathway_id)) {
    kegg_data <- KEGGREST::keggGet(pathway_id)[[1]]
    entrez_ids <- sub("hsa:", "", kegg_data$GENE[seq(1, length(kegg_data$GENE), 2)])
    return(entrez_to_hugo(entrez_ids))
  }
  
  # Reactome case
  if (grepl("^R-HSA-\\d+$", pathway_id)) {
    if (!requireNamespace("ReactomeContentService4R", quietly = TRUE))
      stop("Package 'ReactomeContentService4R' is required for Reactome IDs.")
    reactome_genes <- ReactomeContentService4R::reactome_entities(pathway_id, species = "Homo sapiens")
    entrez_ids <- unique(reactome_genes$identifier[reactome_genes$databaseName == "EntrezGene"])
    return(entrez_to_hugo(entrez_ids))
  }
  
  stop("Unrecognized pathway ID format. Use KEGG (hsaXXXXX) or Reactome (R-HSA-XXXXX).")
}

get_SBGNview_df <- function(pathway_id, meta_table_list) {
  
  # browser()
  
  ## ── 0 · sanity checks ────────────────────────────────────────────────────
  if (!is.list(meta_table_list) || length(meta_table_list) == 0)
    stop("`meta_table_list` must be a non-empty list of data frames.")
  
  not_df <- vapply(meta_table_list, function(x) !is.data.frame(x), logical(1))
  if (any(not_df))
    stop("All elements of `meta_table_list` must be data frames. ",
         "Problems at positions: ", paste(which(not_df), collapse = ", "))
  
  # Each df must have either a 'gene' column or rownames, and a 'median_log2FC' column
  has_gene_ref <- vapply(
    meta_table_list,
    function(df) "gene" %in% colnames(df) || (!is.null(rownames(df)) && any(nchar(rownames(df)) > 0)),
    logical(1)
  )
  if (any(!has_gene_ref))
    stop("Each table must have a 'gene' column or non-empty rownames. ",
         "Issues at positions: ", paste(which(!has_gene_ref), collapse = ", "))
  
  has_med_col <- vapply(meta_table_list, function(df) "median_log2FC" %in% colnames(df), logical(1))
  if (any(!has_med_col))
    stop("Missing 'median_log2FC' column in tables at positions: ",
         paste(which(!has_med_col), collapse = ", "))
  
  # Give tables sensible names if missing
  if (is.null(names(meta_table_list)) || any(names(meta_table_list) == "")) {
    names(meta_table_list) <- paste0("table_", seq_along(meta_table_list))
  }
  
  ## ── 1 · fetch pathway genes (HUGO) ───────────────────────────────────────
  if (exists("get_pathway_genes_hugo")) {
    pathway_genes_vec <- get_pathway_genes_hugo(pathway_id)
  } else if (exists("get_pathway_genes")) {
    pathway_genes_vec <- get_pathway_genes(pathway_id)
  } else {
    stop("Define either get_pathway_genes_hugo() or get_pathway_genes() first.")
  }
  
  pathway_genes_vec <- unique(toupper(stats::na.omit(as.character(pathway_genes_vec))))
  if (length(pathway_genes_vec) == 0)
    stop("No genes returned for pathway_id: ", pathway_id)
  
  ## ── 2 · stack all median_log2FC per gene across tables ───────────────────
  long_list <- lapply(names(meta_table_list), function(nm) {
    df <- meta_table_list[[nm]]
    
    genes <- if ("gene" %in% colnames(df)) df[["gene"]] else rownames(df)
    if (is.null(genes)) return(NULL)
    
    genes <- toupper(as.character(genes))
    vals  <- suppressWarnings(as.numeric(df[["median_log2FC"]]))
    
    keep <- genes %in% pathway_genes_vec
    if (!any(keep)) return(NULL)
    
    data.frame(gene = genes[keep], score = vals[keep], src = nm, stringsAsFactors = FALSE)
  })
  
  long_df <- do.call(rbind, long_list)
  
  ## ── 3 · summarize per gene ───────────────────────────────────────────────
  if (is.null(long_df) || nrow(long_df) == 0) {
    # none of the pathway genes appear in any meta table
    SBGNview_df <- data.frame(
      gene         = pathway_genes_vec,
      median_score = NA_real_,
      avg_score    = NA_real_,
      appearances  = 0L,
      stringsAsFactors = FALSE
    )
    return(SBGNview_df[order(SBGNview_df$gene), ])
  }
  
  split_by_gene <- split(long_df, long_df$gene)
  
  summary_list <- lapply(names(split_by_gene), function(g) {
    dd <- split_by_gene[[g]]
    sc <- dd$score
    data.frame(
      gene          = g,
      median_score  = if (all(is.na(sc))) NA_real_ else stats::median(sc, na.rm = TRUE),
      avg_score     = if (all(is.na(sc))) NA_real_ else mean(sc, na.rm = TRUE),
      appearances   = as.integer(length(unique(dd$src))),   # in how many meta tables it appears
      stringsAsFactors = FALSE
    )
  })
  
  summary_df <- do.call(rbind, summary_list)
  
  # ensure all pathway genes are present (fill missing with NA/0)
  missing_genes <- setdiff(pathway_genes_vec, summary_df$gene)
  if (length(missing_genes) > 0) {
    filler <- data.frame(
      gene         = missing_genes,
      median_score = NA_real_,
      avg_score    = NA_real_,
      appearances  = 0L,
      stringsAsFactors = FALSE
    )
    summary_df <- rbind(summary_df, filler)
  }
  
  # order: more appearances first, then alphabetical
  summary_df <- summary_df[order(-summary_df$appearances, summary_df$gene), ]
  rownames(summary_df) <- NULL
  
  summary_df <- summary_df %>%
    filter(appearances >= length(meta_table_list)) %>%
    dplyr::select(meta_score = median_score, gene)
  
  # Prepare for SBGNview
  gene.data <- summary_df
  
  
  # Convert to matrix for SBGNview (gene names as rownames)
  gene.data <- summary_df %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  gene.data
}

get_pathview_full <- function(meta_table_list, pathway_id,out.suffix,
                              breaks = c(-Inf, -0.5, -0.25, 0.25, 0.5, Inf),
                              labels = c(-2L, -1L, 0L, 1L, 2L)){
  
  #get data
  {
  
    
    gene_data <- get_SBGNview_df(pathway_id = pathway_id, meta_table_list = meta_table_list)
    
    gene.vec <- setNames(as.numeric(gene_data[, "meta_score"]),
                         rownames(gene_data))
    
  }
  
  #get pathview
  {
    gene_disc5 <- cut(
      gene.vec,
      breaks = breaks,
      labels = labels
    )
    
    # Convert factor to numeric and keep gene names
    gene_disc5 <- setNames(as.numeric(as.character(gene_disc5)), names(gene.vec))
    
    # Run Pathview with 5 discrete bins
    pv <- pathview::pathview(
      gene.data   = gene_disc5,
      pathway.id  = pathway_id,
      species     = "hsa",
      gene.idtype = "SYMBOL",
      kegg.native = TRUE,
      node.sum    = "mean",      # <- mean across isoforms in each node
      out.suffix  = out.suffix,
      discrete    = list(gene = TRUE),
      bins        = list(gene = 5),
      limit       = list(gene = 2),
      low         = "red",
      midlow      = "#FF4646",
      mid         = "gray",
      midhigh     = "#B2DF8A",
      high        = "green",
      kegg.native.width  = 2000,   # increase width in pixels
      kegg.native.height = 1500    # increase height in pixels
    )
    
    
  }
  

  
}


draw_legend <- function(file = NULL, type = c("png","pdf"),
                        width_in = 6.0, height_in = 2.0, res = 600,
                        xlabel = "Average(Median (Log2FC))") {
  type <- match.arg(type)
  if (type == "png") {
    # Convert inches to pixels for PNG
    png(file %||% "legend_discrete_5bins.png", 
        width = width_in * res, height = height_in * res, res = res)
  } else {
    pdf(file %||% "legend_discrete_5bins.pdf", 
        width = width_in, height = height_in, useDingbats = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  
  par(mar = c(4, 0.5, 1, 0.5), oma = c(1, 1, 1, 1), bg = "white")
  
  plot(NA, xlim = c(0, 5), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  
  for (i in 1:5) {
    rect(i - 1, 0, i, 1, col = cols[i], border = "black", lwd = 1.5)
  }
  
  axis(1, at = 0:5, labels = labs, lwd = 0, lwd.ticks = 0, tck = 0, cex.axis = 0.9)
  mtext(xlabel, side = 1, line = 2.2, cex = 0.9)
  box(col = NA)
}

get_pathview_per_cell_type <- function(meta_table_list,
                                       pathway_id,
                                       out.suffix,
                                       breaks = c(-Inf, -0.5, -0.25, 0.25, 0.5, Inf),
                                       labels = c(-2L, -1L, 0L, 1L, 2L)) {
  # basic checks
  if (!is.list(meta_table_list) || length(meta_table_list) == 0) {
    stop("`meta_table_list` must be a non-empty list of data frames.")
  }
  
  # run once per cell type
  results <- lapply(names(meta_table_list), function(nm) {
    df <- meta_table_list[[nm]]
    
    # derive cell type from name like "meta_table_cd8_t_cells" -> "cd8_t_cells"
    cell_type <- sub("^meta_table_", "", nm)
    
    # build a single-table list so get_SBGNview_df uses only this cell type
    one_list <- list(df)
    names(one_list) <- nm
    
    # delegate to your existing function (handles SBGN prep + pathview call)
    get_pathview_full(
      meta_table_list = one_list,
      pathway_id      = pathway_id,
      out.suffix      = paste0(out.suffix, "_", cell_type),
      breaks          = breaks,
      labels          = labels
    )
  })
  
  # name results by the cleaned cell type for convenience
  names(results) <- sub("^meta_table_", "", names(meta_table_list))
  
  # pathview writes files; return invisibly in case you want the objects too
  invisible(results)
}




#───────────────────────────────────────────────────────────────────────────────
#  AUC functions ()
#───────────────────────────────────────────────────────────────────────────────




get_meta_table_auc_score <- function(meta_table,
                                     positive_genes,
                                     filter_nocoding = TRUE) {
  ## 0 ·──────── sanity checks ────────────────────────────────────────────────
  if (!is.data.frame(meta_table))
    stop("`meta_table` must be a data-frame.")
  
  if (is.null(rownames(meta_table))) {
    if (!"gene" %in% colnames(meta_table))
      stop("`meta_table` must have gene symbols as row names or a 'gene' column.")
    rownames(meta_table) <- as.character(meta_table$gene)
  }
  
  if (!"median_log2FC" %in% colnames(meta_table))
    stop("`meta_table` is missing required column: 'median_log2FC'.")
  
  if (missing(positive_genes) || length(positive_genes) == 0)
    stop("`positive_genes` must be a non-empty character vector.")
  
  ## 1 ·──────── optional: filter to protein-coding genes ─────────────────────
  if (isTRUE(filter_nocoding)) {
    coding_syms <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      mart  <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = unique(rownames(meta_table)),
        mart       = mart
      )
      unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
    }, error = function(e) {
      warning("Protein-coding filter skipped (biomaRt error: ",
              conditionMessage(e), ").")
      NULL
    })
    if (!is.null(coding_syms))
      meta_table <- meta_table[rownames(meta_table) %in% coding_syms, , drop = FALSE]
  }
  
  ## 2 ·──────── NO padding: use only genes that are found ────────────────────
  positive_genes <- unique(positive_genes)
  
  ## 3 ·──────── prepare the score vector (median_log2FC) ─────────────────────
  logfc_vec <- meta_table$median_log2FC
  names(logfc_vec) <- rownames(meta_table)
  logfc_vec <- logfc_vec[!is.na(logfc_vec)]
  
  pos <- intersect(positive_genes, names(logfc_vec))
  if (!length(pos)) {
    warning("None of the `positive_genes` are present after filtering; returning NA.")
    return(NA_real_)
  }
  neg <- setdiff(names(logfc_vec), pos)
  if (!length(neg)) {
    warning("Background set is empty; returning NA.")
    return(NA_real_)
  }
  
  ## 4 ·──────── compute AUC via Mann–Whitney/Wilcoxon equivalence ────────────
  r <- rank(logfc_vec, ties.method = "average")
  L <- length(pos)
  sum_r_pos <- sum(r[pos])
  U <- sum_r_pos - L * (L + 1) / 2
  auc <- U / (L * length(neg))
  
  ## clamp to [0,1]
  auc <- max(0, min(1, auc))
  return(auc)
}




get_meta_table_auc_score_for_many_cell_types_and_pathways <- function(
    meta_table_list,
    vec_path_ids,
    filter_nocoding = TRUE,
    is_msigdbr = TRUE
) {
  # 0) sanity checks
  if (!is.list(meta_table_list) || length(meta_table_list) == 0)
    stop("`meta_table_list` must be a non-empty list of meta tables.")
  needed <- c("meta_table_nk", "meta_table_monocytes", "meta_table_b_cells",
              "meta_table_cd4_t_cells", "meta_table_cd8_t_cells")
  # ── allow partial presence ────────────────────────────────────────────────
  missing_ct <- setdiff(needed, names(meta_table_list))
  if (length(missing_ct) > 0)
    warning("Missing cell types in `meta_table_list`: ",
            paste(missing_ct, collapse = ", "),
            ". Will fill NA for them.")
  # ──────────────────────────────────────────────────────────────────────────
  if (isTRUE(is_msigdbr)) {
    if (!is.character(vec_path_ids) || !length(vec_path_ids))
      stop("With is_msigdbr=TRUE, `vec_path_ids` must be a non-empty character vector of pathway IDs.")
  } else {
    if (!is.list(vec_path_ids) || length(vec_path_ids) == 0)
      stop("With is_msigdbr=FALSE, `vec_path_ids` must be a non-empty *named* list of character vectors (gene sets).")
    if (is.null(names(vec_path_ids)) || any(!nzchar(names(vec_path_ids))))
      stop("With is_msigdbr=FALSE, `vec_path_ids` must be a *named* list; all elements must have non-empty names.")
    if (!all(vapply(vec_path_ids, function(x) is.character(x) && length(x) > 0, logical(1))))
      stop("With is_msigdbr=FALSE, every list element must be a non-empty character vector of gene symbols.")
  }
  
  # 1) list-name -> display-name mapping & desired column order
  cell_name_map <- c(
    meta_table_nk          = "NKs",
    meta_table_monocytes   = "Monocytes",
    meta_table_b_cells     = "B Cells",
    meta_table_cd4_t_cells = "CD4+ T Cells",
    meta_table_cd8_t_cells = "CD8+ T Cells"
  )
  desired_cols <- unname(cell_name_map)
  
  # 2) helper: resolve a pretty pathway name from ID (used only for msigdbr mode)
  get_pathway_name_by_id <- function(path_id, species = "hsa") {
    if (grepl(sprintf("^%s\\d{5}$", species), path_id, ignore.case = TRUE)) {
      if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("KEGGREST", ask = FALSE, update = FALSE)
      }
      kg <- tryCatch(KEGGREST::keggGet(path_id)[[1]], error = function(e) NULL)
      nm <- if (!is.null(kg) && !is.null(kg$NAME)) kg$NAME else path_id
      return(as.character(nm))
    }
    if (grepl("^R-HSA-\\d+$", path_id)) {
      if (!requireNamespace("reactome.db", quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install("reactome.db", ask = FALSE, update = FALSE)
      }
      nm <- tryCatch(reactome.db::reactomePATHID2NAME[[path_id]], error = function(e) NULL)
      return(if (is.null(nm)) path_id else as.character(nm))
    }
    nm <- toupper(path_id)
    nm <- gsub("_", " ", nm, fixed = TRUE)
    nm <- tools::toTitleCase(tolower(nm))
    nm
  }
  
  # 3) prepare result container: rows = pathways, cols = cell types + found_genes
  if (isTRUE(is_msigdbr)) {
    path_ids  <- vec_path_ids
    path_names <- vapply(path_ids, get_pathway_name_by_id, character(1))
    n_paths <- length(path_ids)
  } else {
    # list-of-genes mode
    path_ids  <- names(vec_path_ids)  # use names as IDs/rownames
    path_names <- path_ids
    n_paths <- length(vec_path_ids)
  }
  
  res <- matrix(NA_real_, nrow = n_paths, ncol = length(desired_cols),
                dimnames = list(path_names, desired_cols))
  found_str_vec <- rep(NA_character_, n_paths)  # <── replaces num_genes
  
  # (optional) precompute protein-coding symbols once across all tables
  coding_syms_all <- NULL
  if (isTRUE(filter_nocoding)) {
    coding_syms_all <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      all_syms <- unique(unlist(lapply(meta_table_list, function(df) {
        rn <- rownames(df)
        if (is.null(rn) && "gene" %in% colnames(df)) rn <- as.character(df$gene)
        rn
      })))
      mart  <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters    = "hgnc_symbol",
        values     = all_syms,
        mart       = mart
      )
      unique(annot$hgnc_symbol[annot$gene_biotype == "protein_coding"])
    }, error = function(e) NULL)
  }
  
  # helper to get available genes in one meta table (respecting filter_nocoding)
  get_available_genes <- function(df) {
    rn <- rownames(df)
    if (is.null(rn) && "gene" %in% colnames(df)) rn <- as.character(df$gene)
    rn <- unique(rn)
    if (!is.null(coding_syms_all)) rn <- intersect(rn, coding_syms_all)
    rn
  }
  
  # Precompute union once
  available_union <- unique(unlist(lapply(meta_table_list, get_available_genes)))
  
  # 4) loop over pathways
  for (i in seq_len(n_paths)) {
    if (isTRUE(is_msigdbr)) {
      pid <- path_ids[i]
      positive_genes <- tryCatch(get_pathway_genes_by_id(pid), error = function(e) character())
    } else {
      pid <- path_ids[i]
      positive_genes <- vec_path_ids[[pid]]
      if (!is.character(positive_genes)) positive_genes <- character()
    }
    
    if (!length(positive_genes)) next
    
    # x/y: x = #found across union of available genes in all present cell types; y = pathway size
    x_found <- length(intersect(unique(positive_genes), available_union))
    y_all   <- length(unique(positive_genes))
    found_str_vec[i] <- sprintf("%d/%d", x_found, y_all)
    
    for (nm in names(cell_name_map)) {
      display <- cell_name_map[[nm]]
      mt <- meta_table_list[[nm]]
      if (is.null(mt)) {                 # missing cell type → NA
        res[i, display] <- NA_real_
        next
      }
      score <- tryCatch(
        get_meta_table_auc_score(
          meta_table      = mt,
          positive_genes  = positive_genes,
          filter_nocoding = filter_nocoding
        ),
        error = function(e) NA_real_
      )
      res[i, display] <- score
    }
  }
  
  # 5) return as data.frame with found_genes column
  df <- as.data.frame(res, check.names = FALSE, stringsAsFactors = FALSE)
  df$found_genes <- found_str_vec
  df
}








# genes: character vector of HGNC symbols
# categories: optional filter like c("H","C2","C5","CP:REACTOME","CP:BIOCARTA")
find_msig_sets_with_genes <- function(genes,
                                      species = "Homo sapiens",
                                      categories = NULL,
                                      require_all = TRUE,
                                      ignore_case = TRUE) {
  stopifnot(is.character(genes), length(genes) > 0)
  
  # normalize case if desired
  if (ignore_case) genes <- toupper(genes)
  
  df <- msigdbr(species = species, category = categories)
  if (ignore_case) df <- df %>% mutate(gene_symbol = toupper(gene_symbol))
  
  df %>%
    filter(gene_symbol %in% genes) %>%
    group_by(gs_cat, gs_subcat, gs_name) %>%
    summarise(
      set_size   = n_distinct(gene_symbol),            # how many of your query genes are in
      overlap    = paste(sort(unique(gene_symbol)), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(
      query_size  = length(unique(genes)),
      hits_needed = if (require_all) query_size else 1
    ) %>%
    filter(set_size >= hits_needed) %>%
    arrange(desc(set_size), gs_cat, gs_name)
}






# get_auc_heatmap and its helper functions--------------------------------------------------------------

# get_metric_heatmap helper functions:

## A) five equally spaced legend ticks, based on real min/max
legend_breaks <- function(mat) {
  rng <- range(mat, na.rm = TRUE)
  
  adjust_border <- function(x) {
    int <- round(x)
    if (abs(x - int) <= 0.45) int else x
  }
  rng <- vapply(rng, adjust_border, numeric(1))
  if (diff(rng) == 0) rng <- rng + c(-1, 1)
  
  round(seq(rng[1], rng[2], length.out = 5), 2)
}

## B) split any phrase roughly half‑way between words
nice_two_line_title <- function(x) {
  words <- strsplit(trimws(x), "\\s+")[[1]]
  if (length(words) <= 1) return(x)
  cut <- ceiling(length(words) / 2)
  paste(paste(words[1:cut], collapse = " "),
        paste(words[-(1:cut)], collapse = " "),
        sep = "\n")
}

## C) grey square legend for NA
legend_NA <- function(title = "NA", fill = "grey85") {
  Legend(labels    = title,
         title     = NULL,
         type      = "points",
         pch       = 15,
         size      = unit(4, "mm"),
         legend_gp = gpar(col = fill, fill = fill))
}




get_auc_heatmap <- function(df, title,pathway_names,row_order,  meta_name_fontsize = 11) {
  # --- Sanity checks ----------------------------------------------------------
  required_cols <- c("sample", "pathway", "AUC", "dataset",
                     "condition", "condition_subtype", "age", "sex", "EDSS")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input data.frame is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(df$AUC))  stop("Column 'AUC' must be numeric, but got: ", class(df$AUC)[1])
  if (!is.numeric(df$age))  stop("Column 'age' must be numeric, but got: ", class(df$age)[1])
  if (!is.numeric(df$EDSS)) stop("Column 'EDSS' must be numeric, but got: ", class(df$EDSS)[1])
  if (!is.character(df$sample) && !is.factor(df$sample))
    stop("Column 'sample' must be character or factor, but got: ", class(df$sample)[1])
  if (!is.character(df$pathway) && !is.factor(df$pathway))
    stop("Column 'pathway' must be character or factor, but got: ", class(df$pathway)[1])
  
  # --- 1. Prepare meta --------------------------------------------------------
  df <- df |>
    dplyr::mutate(sample  = as.character(.data$sample),
                  pathway = as.character(.data$pathway))
  
  ## ── 2. Rename pathways (shorter labels) & set desired row order ───────────

  
  # apply renaming (unmapped names stay as-is)
  df$pathway <- ifelse(df$pathway %in% names(pathway_names),
                       unname(pathway_names[df$pathway]),
                       df$pathway)
  
  # Per-sample metadata (one row per sample)
  filt_meta <- df |>
    dplyr::arrange(.data$sample) |>
    dplyr::distinct(.data$sample, .keep_all = TRUE)
  
  # --- 2. Build AUC matrix: pathways (rows) × samples (cols) -----------------
  auc_wide <- df |>
    dplyr::select(.data$pathway, .data$sample, .data$AUC) |>
    dplyr::group_by(.data$pathway, .data$sample) |>
    dplyr::summarise(AUC = mean(.data$AUC, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from  = .data$sample,
      values_from = .data$AUC,
      values_fill = NA_real_
    ) |>
    dplyr::arrange(.data$pathway)
  
  auc_mat <- as.matrix(auc_wide[ , -1, drop = FALSE])
  rownames(auc_mat) <- auc_wide$pathway
  
  # Reorder rows to desired order (keep any extras after)
  desired_rows <- intersect(row_order, rownames(auc_mat))
  other_rows   <- setdiff(rownames(auc_mat), desired_rows)
  auc_mat      <- auc_mat[c(desired_rows, other_rows), , drop = FALSE]
  
  # --- Allocate space for long left-side row labels --------------------------
  row_lab_gp <- grid::gpar(fontsize = 11)
  row_lab_w  <- ComplexHeatmap::max_text_width(rownames(auc_mat), gp = row_lab_gp) * 1.1
  
  # --- Custom column ordering -------------------------------------------------
  condition_order <- c("CD59 Def.", "Multiple Sclerosis", "Ab-mediated IDD")
  filt_meta$condition <- factor(filt_meta$condition, levels = condition_order)
  filt_meta <- filt_meta |>
    dplyr::arrange(.data$condition, .data$dataset, .data$sample)
  
  common_cols <- filt_meta$sample
  auc_mat     <- auc_mat[ , common_cols, drop = FALSE]
  filt_meta   <- filt_meta[match(common_cols, filt_meta$sample), , drop = FALSE]
  
  # --- 3. Colour palettes -----------------------------------------------------
  dataset_pal <- c(
    Esaulova = "#E9A5F1", Kang = "#FFB4B4",
    Ramesh = "#e76f51", Schafflick = "#000000",
    Mevorach ="#ff4d6d"
  )
  subtype_pal <- c(CIS = "#bb9457", RRMS = "#432818", "CD59 Def." = "#fb8500",
                   "MOGAD" = "#2C6E49", "NMOSD" = "#B2DF8A")
  condition_pal <- c("Multiple Sclerosis" = "#023047", "Ab-mediated IDD" = "#8ecae6", "CD59 Def." = "#fb8500")
  sex_pal     <- c(Female = "#ff8daa", Male = "#32b7ff")
  edss_pal <- circlize::colorRamp2(c(0, 6),  c("#70e000", "#e01e37"))
  age_pal  <- circlize::colorRamp2(c(10, 80), c("#48CAE4", "#03045e"))
  
  # --- 4. Column annotation (order + labels) ---------------------------------
  anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition       = filt_meta$condition,
    `Sub-Condition` = filt_meta$condition_subtype,
    Dataset         = filt_meta$dataset,
    Sex             = filt_meta$sex,
    Age             = filt_meta$age,
    EDSS            = filt_meta$EDSS,
    col = list(
      Condition       = condition_pal,
      `Sub-Condition` = subtype_pal,
      Dataset         = dataset_pal,
      Sex             = sex_pal,
      Age             = age_pal,
      EDSS            = edss_pal
    ),
    na_col = "grey85",
    annotation_legend_param = list(
      EDSS = list(
        title     = "EDSS",
        at        = seq(0, 6, 2),
        labels    = as.character(seq(0, 6, 2)),
        legend_gp = grid::gpar(fill = edss_pal(seq(0, 6, 2)))
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp   = grid::gpar(fontsize = meta_name_fontsize)  # <— control metadata label font
  )
  
  # --- 5. Heat-map ------------------------------------------------------------
  ht <- ComplexHeatmap::Heatmap(
    auc_mat,
    name                        = "AUC",
    top_annotation              = anno,
    show_row_names              = TRUE,
    row_names_side              = "left",
    row_names_gp                = row_lab_gp,
    row_names_max_width         = row_lab_w,
    show_column_names           = FALSE,
    cluster_rows                = FALSE,
    row_order                   = rownames(auc_mat),
    cluster_columns             = FALSE,
    column_order                = common_cols,
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    na_col                      = "grey85",
    heatmap_legend_param        = list(
      at    = legend_breaks(auc_mat),
      title = "AUC",
      title_position = "leftcenter-rot"
    ),
    column_title    = title,
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 14)
  )
  
  return(ht)
}


get_auc_boxplot <- function(
    df, title,
    point_size         = 4,
    title_size         = 24,
    axis_title_size    = 22,
    axis_text_size     = 20,
    legend_title_size  = 22,
    legend_label_size  = 20,
    plot_margin_left   = 130
) {
  # --- Sanity checks ----------------------------------------------------------
  required_cols <- c("pathway", "AUC", "dataset", "condition")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Input data.frame is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(df$AUC)) {
    stop("Column 'AUC' must be numeric, but got: ", class(df$AUC)[1])
  }
  if (!is.character(df$pathway) && !is.factor(df$pathway)) {
    stop("Column 'pathway' must be character or factor, but got: ",
         class(df$pathway)[1])
  }
  if (!is.character(df$dataset) && !is.factor(df$dataset)) {
    stop("Column 'dataset' must be character or factor, but got: ",
         class(df$dataset)[1])
  }
  if (!is.character(df$condition) && !is.factor(df$condition)) {
    stop("Column 'condition' must be character or factor, but got: ",
         class(df$condition)[1])
  }
  # Optional: warn if AUC outside [0,1]
  if (any(!is.na(df$AUC) & (df$AUC < 0 | df$AUC > 1))) {
    warning("Some 'AUC' values are outside [0, 1]; they may be clipped by y-axis limits.")
  }
  
  
  
  condition_pal <- c(
    "Multiple Sclerosis" = "#012A3A",
    "Ab-mediated IDD"    = "#2C7FB8",
    "CD59 Def."          = "#E66100"   # deep orange
  )
  
  # Palette coverage check
  missing_cond <- setdiff(unique(as.character(df$condition)), names(condition_pal))
  if (length(missing_cond) > 0) {
    warning("Conditions missing from color palette: ", paste(missing_cond, collapse = ", "))
  }
  
  # --- Ensure desired condition order for grouped boxes -----------------------
  df$condition <- factor(as.character(df$condition),
                         levels = c("CD59 Def.", "Multiple Sclerosis", "Ab-mediated IDD"))
  
  # --- Pathway renaming & ordering -------------------------------------------
  pathway_names <- c(
    "Kegg Medicus Reference Electron Transfer in Complex i"         = "KEGG-MR: ETC in Complex I",
    "Reactome Complex i Biogenesis"                                 = "Reactome: Complex I Biogenesis",
    "Wp Mitochondrial Complex Ii Assembly"                          = "WP: Complex II Assembly",
    "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "GO-BP:ETC in Complex II",
    "Wp Mitochondrial Complex Iii Assembly"                         = "WP: Complex III Assembly",
    "Kegg Medicus Reference Electron Transfer in Complex Iii"       = "KEGG-MR: ETC in Complex III",
    "Wp Mitochondrial Complex Iv Assembly"                          = "WP: Complex IV Assembly",
    "Kegg Medicus Reference Electron Transfer in Complex Iv"        = "KEGG-MR: ETC in Complex IV",
    "Gobp Mitochondrial Proton Transporting Atp Synthase Complex Assembly" = "GO-BP: Complex V Assembly",
    "Gobp Proton Motive Force Driven Atp Synthesis"                 = "GO-BP: PMF-Driven ATP Synthesis"
  )
  row_order <- c(
    "Reactome: Complex I Biogenesis",
    "KEGG-MR: ETC in Complex I",
    "WP: Complex II Assembly",
    "GO-BP:ETC in Complex II",
    "WP: Complex III Assembly",
    "KEGG-MR: ETC in Complex III",
    "WP: Complex IV Assembly",
    "KEGG-MR: ETC in Complex IV",
    "GO-BP: Complex V Assembly",
    "GO-BP: PMF-Driven ATP Synthesis"
  )
  
  # case-insensitive rename; keep originals if not in map
  map_lc <- stats::setNames(unname(pathway_names), tolower(names(pathway_names)))
  cur    <- as.character(df$pathway)
  renamed <- unname(map_lc[tolower(cur)])
  df$pathway <- ifelse(is.na(renamed), cur, renamed)
  
  # order: requested order first, then any extras (if present)
  level_order <- c(row_order, setdiff(unique(df$pathway), row_order))
  df$pathway  <- factor(df$pathway, levels = level_order)
  
  # --- Plot: grouped boxplots with darker colored outlines, no fill ----------
  dodge_w <- 0.8
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pathway, y = AUC, colour = condition)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(pathway, condition)),
      outlier.shape = NA,
      width = 0.6,
      fill  = NA,     # no interior fill
      size  = 0.8,    # slightly thicker outlines for visibility
      position = ggplot2::position_dodge(width = dodge_w)
    ) +
    ggplot2::scale_color_manual(
      name   = "Condition",
      values = condition_pal,
      drop   = FALSE
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "AUC score"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "#efefef"),
      axis.text.x  = ggplot2::element_text(angle = 25, hjust = 1, vjust = 1, size = axis_text_size),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      plot.title   = ggplot2::element_text(face = "bold", size = title_size),
      legend.position = "right",
      legend.box      = "vertical",
      legend.title    = ggplot2::element_text(size = legend_title_size),
      legend.text     = ggplot2::element_text(size = legend_label_size),
      plot.margin     = ggplot2::margin(t = 10, r = 10, b = 10, l = plot_margin_left, unit = "pt")
    )
  
  return(p)
}


get_auc_boxplot_for_pathway <- function(
    auc_list, pathway_name, title,
    title_size         = 24,
    axis_title_size    = 22,
    axis_text_size     = 20,
    legend_title_size  = 22,
    legend_label_size  = 20,
    plot_margin_left   = 130
) {
  # --- Sanity checks ----------------------------------------------------------
  if (!is.list(auc_list) || length(auc_list) == 0)
    stop("'auc_list' must be a non-empty list of data.frames.")
  
  needed_cols <- c("sample", "pathway", "AUC", "condition")
  has_cols <- vapply(auc_list, function(df)
    all(needed_cols %in% colnames(df)), logical(1))
  if (!all(has_cols)) {
    bad <- names(auc_list)[!has_cols]
    stop("These elements of 'auc_list' are missing required columns: ",
         paste(bad, collapse = ", "),
         ". Each must contain: ", paste(needed_cols, collapse = ", "))
  }
  
  if (!is.character(pathway_name) || length(pathway_name) != 1 || nchar(pathway_name) == 0)
    stop("'pathway_name' must be a single non-empty character string.")
  
  # --- Palette ---------------------------------------------------------------
  condition_pal <- c(
    "Multiple Sclerosis" = "#012A3A",
    "Ab-mediated IDD"    = "#2C7FB8",
    "CD59 Def."          = "#E66100"   # deep orange
  )
  
  # --- Desired order of x-axis (cell types) ----------------------------------
  desired_cell_types <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
  
  # --- Pathway renaming & ordering (reuse of your mapping) -------------------
  pathway_names <- c(
    "Kegg Medicus Reference Electron Transfer in Complex i"         = "KEGG-MR: ETC in Complex I",
    "Reactome Complex i Biogenesis"                                 = "Reactome: Complex I Biogenesis",
    "Wp Mitochondrial Complex Ii Assembly"                          = "WP: Complex II Assembly",
    "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "GO-BP:ETC in Complex II",
    "Wp Mitochondrial Complex Iii Assembly"                         = "WP: Complex III Assembly",
    "Kegg Medicus Reference Electron Transfer in Complex Iii"       = "KEGG-MR: ETC in Complex III",
    "Wp Mitochondrial Complex Iv Assembly"                          = "WP: Complex IV Assembly",
    "Kegg Medicus Reference Electron Transfer in Complex Iv"        = "KEGG-MR: ETC in Complex IV",
    "Gobp Mitochondrial Proton Transporting Atp Synthase Complex Assembly" = "GO-BP: Complex V Assembly",
    "Gobp Proton Motive Force Driven Atp Synthesis"                 = "GO-BP: PMF-Driven ATP Synthesis"
  )
  map_lc <- stats::setNames(unname(pathway_names), tolower(names(pathway_names)))
  
  prettify_path <- function(x) {
    cur <- as.character(x)
    renamed <- unname(map_lc[tolower(cur)])
    ifelse(is.na(renamed), cur, renamed)
  }
  
  pretty_target <- prettify_path(pathway_name)[1]
  
  # --- Build combined df across requested cell types --------------------------
  present_cells <- intersect(desired_cell_types, names(auc_list))
  if (length(present_cells) == 0)
    stop("None of the expected cell types were found in 'auc_list': ",
         paste(desired_cell_types, collapse = ", "))
  
  combined <- do.call(
    rbind,
    lapply(present_cells, function(ct) {
      df <- auc_list[[ct]]
      df$cell_type <- ct
      df
    })
  )
  
  if (!is.numeric(combined$AUC)) {
    stop("Column 'AUC' must be numeric, but got: ", class(combined$AUC)[1])
  }
  
  if (any(!is.na(combined$AUC) & (combined$AUC < 0 | combined$AUC > 1))) {
    warning("Some 'AUC' values are outside [0, 1]; they may be clipped by y-axis limits.")
  }
  
  combined$pathway <- prettify_path(combined$pathway)
  df <- combined[combined$pathway == pretty_target, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No rows found for pathway '", pathway_name,
         "'. (After normalization matched: '", pretty_target, "').")
  }
  
  # --- Ordering ---------------------------------------------------------------
  # For plotting/dodging (left->right): AbID, MS, CD55, CD59
  # This yields the requested right->left visual order: CD59, CD55, MS, AbID
  plot_level_order <- c("CD59 Def.", "Multiple Sclerosis", "Ab-mediated IDD")
  df$condition <- factor(as.character(df$condition), levels = plot_level_order)
  df$cell_type <- factor(df$cell_type, levels = desired_cell_types)
  
  # Palette coverage check
  missing_cond <- setdiff(unique(as.character(df$condition)), names(condition_pal))
  if (length(missing_cond) > 0) {
    warning("Conditions missing from color palette: ", paste(missing_cond, collapse = ", "))
  }
  
  # --- Plot -------------------------------------------------------------------
  dodge_w <- 0.8
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cell_type, y = AUC, colour = condition)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(group = interaction(cell_type, condition)),
      outlier.shape = NA,
      width = 0.6,
      fill  = NA,
      size  = 0.8,
      position = ggplot2::position_dodge(width = dodge_w)
    ) +
    ggplot2::scale_color_manual(
      name   = "Condition",
      values = condition_pal,
      drop   = FALSE,
      # Legend order (top->bottom) to match your requested right->left: CD59, CD55, MS, AbID
      breaks = c("CD59 Def.", "CD55 Def.", "Multiple Sclerosis", "Ab-mediated IDD")
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "AUC score"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "#efefef"),
      axis.text.x  = ggplot2::element_text(angle = 25, hjust = 1, vjust = 1, size = axis_text_size),
      axis.text.y  = ggplot2::element_text(size = axis_text_size),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      plot.title   = ggplot2::element_text(face = "bold", size = title_size),
      legend.position = "right",
      legend.box      = "vertical",
      legend.title    = ggplot2::element_text(size = legend_title_size),
      legend.text     = ggplot2::element_text(size = legend_label_size),
      plot.margin     = ggplot2::margin(t = 10, r = 10, b = 10, l = plot_margin_left, unit = "pt")
    )
  
  return(p)
}



# cell_type_list: named list of 5 data.frames (NKs, Monocytes, B Cells, CD4+ T Cells, CD8+ T Cells)
# strict: if TRUE, stop() when any required pathway is missing in a cell type; if FALSE, compute with available ones (score = mean of observed)

get_metabolic_shifits_list <- function(cell_type_list, strict = TRUE, verbose = TRUE) {
  # 0) Sanity checks on the outer structure -----------------------------------
  required_ct <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
  if (!is.list(cell_type_list) || length(cell_type_list) == 0) {
    stop("cell_type_list must be a non-empty list.")
  }
  if (!all(required_ct %in% names(cell_type_list))) {
    warning(sprintf(
      "cell_type_list is missing these cell types: %s",
      paste(setdiff(required_ct, names(cell_type_list)), collapse = ", ")
    ))
  }
  extra_ct <- setdiff(names(cell_type_list), required_ct)
  if (length(extra_ct) > 0 && verbose) {
    message("Note: extra cell types present and will be processed: ",
            paste(extra_ct, collapse = ", "))
  }
  
  # 1) Required columns for each data.frame ------------------------------------
  req_cols <- c("sample", "pathway", "AUC", "dataset",
                "condition", "condition_subtype", "age", "sex", "EDSS")
  
  check_df_cols <- function(df, ct_name) {
    if (!is.data.frame(df)) stop(sprintf("Element '%s' is not a data.frame.", ct_name))
    missing <- setdiff(req_cols, colnames(df))
    if (length(missing) > 0) {
      stop(sprintf("Cell type '%s' is missing required columns: %s",
                   ct_name, paste(missing, collapse = ", ")))
    }
    if (!is.numeric(df$AUC)) stop(sprintf("Column AUC must be numeric in '%s'.", ct_name))
    invisible(TRUE)
  }
  
  # 2) Pathway → metric mapping (EXACT STRINGS as provided) --------------------
  ci_label <- "CI-linked NADH:"  # keep the colon, as requested
  
  metrics_map <- list()
  metrics_map[[ci_label]] <- c(
    "Wp Thiamine Metabolic Pathways",
    "Gobp Pyruvate Metabolic Process",
    "Gomf Pyruvate Dehydrogenase Acetyl Transferring Phosphatase Activity",
    "Gobp Glutamate Metabolic Process",
    "Reactome Malate Aspartate Shuttle"
  )
  metrics_map[["FAO"]] <- c(
    "Gobp Fatty Acid Beta Oxidation using Acyl Coa Dehydrogenase",  # ACAD
    "Reactome Mitochondrial Fatty Acid Beta Oxidation",             # Mito FAO
    "Reactome Carnitine Shuttle",
    "Gomf Acetyl Coa c Acyltransferase Activity"                    # Thiolase
  )
  metrics_map[["Succinate"]] <- c(
    "Kegg Valine Leucine and Isoleucine Degradation",               # BCAA degradation
    "Kegg Propanoate Metabolism",                                   # Propionate
    "Gobp Succinyl Coa Catabolic Process"                           # Succinyl-CoA catabolism
  )
  
  # Prepare mapping table (metric-pathway rows)
  mapping_df <- purrr::imap_dfr(metrics_map, function(paths, metric_name) {
    tibble::tibble(metric = metric_name, pathway = paths)
  })
  
  # Helper to compute metrics for one cell type --------------------------------
  compute_one_ct <- function(df_ct, ct_name) {
    check_df_cols(df_ct, ct_name)
    
    # Normalize whitespace defensively
    df_ct <- df_ct %>%
      mutate(pathway = gsub("\\s+", " ", trimws(pathway)))
    
    # Check presence of required pathways (global)
    present <- unique(df_ct$pathway)
    missing_paths <- setdiff(mapping_df$pathway, present)
    if (length(missing_paths) > 0) {
      msg <- paste0(
        "Cell type '", ct_name, "' is missing ",
        length(missing_paths), " required pathway(s): ",
        paste(missing_paths, collapse = " | ")
      )
      if (strict) stop(msg) else if (verbose) warning(msg, call. = FALSE)
    }
    
    # Compute base metrics (three means) per sample
    base_long <- df_ct %>%
      inner_join(mapping_df, by = "pathway") %>%
      group_by(sample, dataset, condition, condition_subtype, age, sex, EDSS, metric) %>%
      summarise(score = mean(AUC, na.rm = TRUE), .groups = "drop")
    
    # Ensure all (sample x metric) combos exist; fill missing with NA
    sample_keys <- df_ct %>%
      distinct(sample, dataset, condition, condition_subtype, age, sex, EDSS)
    
    all_metric_grid <- tidyr::expand_grid(
      sample_keys,
      metric = names(metrics_map)
    )
    
    base_long_complete <- all_metric_grid %>%
      left_join(base_long, by = c(colnames(sample_keys), "metric"))
    
    # Compute shifts from wide form
    wide <- base_long_complete %>%
      tidyr::pivot_wider(names_from = metric, values_from = score)
    
    # Guarantee columns exist
    for (nm in c("FAO", "Succinate", ci_label)) {
      if (!nm %in% colnames(wide)) wide[[nm]] <- NA_real_
    }
    
    # Compute shifts without rlang/data.table complications
    wide$FAO_shift  <- wide$FAO - wide[[ci_label]]
    wide$Succ_shift <- wide$Succinate - wide[[ci_label]]
    
    # Back to long; order metrics
    out_long <- wide %>%
      tidyr::pivot_longer(
        cols = c(names(metrics_map), "FAO_shift", "Succ_shift"),
        names_to = "metric",
        values_to = "score"
      ) %>%
      arrange(sample, metric)
    
    # Final column order
    out_long %>%
      dplyr::select(sample, metric, score, dataset, condition, condition_subtype, age, sex, EDSS)
  }
  
  # 3) Apply to each cell type present -----------------------------------------
  output_list <- purrr::imap(cell_type_list, compute_one_ct)
  
  # Order result as canonical cell type order first, then extras
  ordered <- output_list[intersect(required_ct, names(output_list))]
  extras  <- output_list[setdiff(names(output_list), required_ct)]
  output_list <- c(ordered, extras)
  
  if (verbose) {
    message("Done. Returned a list with ", length(output_list), " item(s): ",
            paste(names(output_list), collapse = ", "))
  }
  
  return(output_list)
}






get_metabolic_shifit_heatmap <- function(metabolic_shifits_list, title, row_order = NULL, meta_name_fontsize = 11) {
  # --- Sanity: list & per-DF columns -----------------------------------------
  if (!is.list(metabolic_shifits_list) || length(metabolic_shifits_list) == 0) {
    stop("metabolic_shifits_list must be a non-empty list of data frames.")
  }
  req_cols <- c("sample","metric","score","dataset","condition","condition_subtype","age","sex","EDSS")
  lapply(names(metabolic_shifits_list), function(ct) {
    df <- metabolic_shifits_list[[ct]]
    if (!is.data.frame(df)) stop("Element '", ct, "' is not a data.frame.")
    miss <- setdiff(req_cols, colnames(df))
    if (length(miss) > 0) stop("Cell type '", ct, "' is missing columns: ", paste(miss, collapse=", "))
    if (!is.numeric(df$score)) stop("Column 'score' must be numeric in '", ct, "'.")
    if (!is.numeric(df$age))  stop("Column 'age' must be numeric in '", ct, "'.")
    if (!is.numeric(df$EDSS)) stop("Column 'EDSS' must be numeric in '", ct, "'.")
    if (!is.character(df$sample) && !is.factor(df$sample)) stop("'sample' must be character/factor in '", ct, "'.")
  })
  
  # --- Map cell-type labels for row names ------------------------------------
  ct_label_map <- c(
    "NKs"          = "NK",
    "Monocytes"    = "Monocytes",
    "B Cells"      = "B_Cells",
    "CD4+ T Cells" = "CD4_T_Cells",
    "CD8+ T Cells" = "CD8_T_Cells"
  )
  
  # --- 1) Bind list → long df with 'row' names (two shifts per cell type) ----
  wanted_metrics <- c("FAO_shift","Succ_shift")
  
  long_df <- do.call(dplyr::bind_rows, lapply(names(metabolic_shifits_list), function(ct) {
    df <- metabolic_shifits_list[[ct]]
    lbl <- ifelse(ct %in% names(ct_label_map), ct_label_map[[ct]], gsub("\\s+", "_", ct))
    df |>
      dplyr::filter(.data$metric %in% wanted_metrics) |>
      dplyr::mutate(
        sample = as.character(.data$sample),
        row    = paste0(.data$metric, "_", lbl)
      ) |>
      dplyr::select(.data$row, dplyr::all_of(req_cols))
  }))
  
  # Per-sample metadata (one row per sample)
  filt_meta <- long_df |>
    dplyr::arrange(.data$sample) |>
    dplyr::distinct(.data$sample, .keep_all = TRUE)
  
  # --- 2) Build SHIFT matrix: rows (row) × columns (sample) ------------------
  shift_wide <- long_df |>
    dplyr::select(.data$row, .data$sample, .data$score) |>
    dplyr::group_by(.data$row, .data$sample) |>
    dplyr::summarise(score = mean(.data$score, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from  = .data$sample,
      values_from = .data$score,
      values_fill = NA_real_
    ) |>
    dplyr::arrange(.data$row)
  
  shift_mat <- as.matrix(shift_wide[, -1, drop = FALSE])
  rownames(shift_mat) <- shift_wide$row
  
  # Desired row order (if not supplied, canonical FAO/Succ per CT order)
  if (is.null(row_order)) {
    base_rows <- unname(ct_label_map[names(ct_label_map)])
    row_order <- as.vector(rbind(paste0("FAO_shift_",  base_rows),
                                 paste0("Succ_shift_", base_rows)))
  }
  desired_rows <- intersect(row_order, rownames(shift_mat))
  other_rows   <- setdiff(rownames(shift_mat), desired_rows)
  shift_mat    <- shift_mat[c(desired_rows, other_rows), , drop = FALSE]
  
  # --- Allocate space for long left-side row labels --------------------------
  row_lab_gp <- grid::gpar(fontsize = 11)
  row_lab_w  <- ComplexHeatmap::max_text_width(rownames(shift_mat), gp = row_lab_gp) * 1.1
  
  # --- 3) Custom column ordering (same logic as your AUC heatmap) ------------
  condition_order <- c("CD59 Def.", "CD55 Def.", "Multiple Sclerosis", "Ab-mediated IDD")
  filt_meta$condition <- factor(filt_meta$condition, levels = condition_order)
  filt_meta <- filt_meta |>
    dplyr::arrange(.data$condition, .data$dataset, .data$sample)
  
  common_cols <- filt_meta$sample
  shift_mat   <- shift_mat[, common_cols, drop = FALSE]
  filt_meta   <- filt_meta[match(common_cols, filt_meta$sample), , drop = FALSE]
  
  # --- 4) Palettes (unchanged) + FIXED SHIFT SCALE ---------------------------
  dataset_pal <- c(
    Esaulova = "#E9A5F1", Kang = "#FFB4B4",
    Ramesh = "#e76f51", Schafflick = "#000000",
    Mevorach = "#ff4d6d"
  )
  subtype_pal <- c(
    CIS = "#bb9457", RRMS = "#432818",
    "CD55 Def." = "#ffb703", "CD59 Def." = "#fb8500",
    "MOGAD" = "#2C6E49", "NMOSD" = "#B2DF8A"
  )
  condition_pal <- c(
    "Multiple Sclerosis" = "#023047", "Ab-mediated IDD" = "#8ecae6",
    "CD55 Def." = "#ffb703", "CD59 Def." = "#fb8500"
  )
  sex_pal     <- c(Female = "#ff8daa", Male = "#32b7ff")
  edss_pal    <- circlize::colorRamp2(c(0, 6),  c("#70e000", "#e01e37"))
  age_pal     <- circlize::colorRamp2(c(20, 80), c("#48CAE4", "#03045e"))
  
  # >>> NEW: fixed color scale for shifts [-0.25, 0.25]
  shift_col_fun <- circlize::colorRamp2(
    c(-0.25, 0, 0.25),
    c("#4575b4", "white", "#d73027")  # choose your diverging colors
  )
  # >>> OPTIONAL: clip values outside the range so legend matches tiles
  shift_mat[shift_mat < -0.25] <- -0.25
  shift_mat[shift_mat >  0.25] <-  0.25
  
  # --- 5) Column annotation ---------------------------------------------------
  anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition       = filt_meta$condition,
    `Sub-Condition` = filt_meta$condition_subtype,
    Dataset         = filt_meta$dataset,
    Sex             = filt_meta$sex,
    Age             = filt_meta$age,
    EDSS            = filt_meta$EDSS,
    col = list(
      Condition       = condition_pal,
      `Sub-Condition` = subtype_pal,
      Dataset         = dataset_pal,
      Sex             = sex_pal,
      Age             = age_pal,
      EDSS            = edss_pal
    ),
    na_col = "grey85",
    annotation_legend_param = list(
      EDSS = list(
        title     = "EDSS",
        at        = seq(0, 6, 2),
        labels    = as.character(seq(0, 6, 2)),
        legend_gp = grid::gpar(fill = edss_pal(seq(0, 6, 2)))
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp   = grid::gpar(fontsize = meta_name_fontsize)
  )
  
  # --- 6) Heatmap -------------------------------------------------------------
  ht <- ComplexHeatmap::Heatmap(
    shift_mat,
    name                        = "Shift",
    top_annotation              = anno,
    show_row_names              = TRUE,
    row_names_side              = "left",
    row_names_gp                = row_lab_gp,
    row_names_max_width         = row_lab_w,
    show_column_names           = FALSE,
    cluster_rows                = FALSE,
    row_order                   = rownames(shift_mat),
    cluster_columns             = FALSE,
    column_order                = common_cols,
    clustering_distance_rows    = "euclidean",
    clustering_distance_columns = "euclidean",
    na_col                      = "grey85",
    # >>> CHANGE: fix legend ticks and pass fixed color scale
    col                         = shift_col_fun,
    heatmap_legend_param        = list(
      at    = seq(-0.25, 0.25, by = 0.05),
      title = "Shift",
      title_position = "leftcenter-rot"
    ),
    column_title    = title,
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 14)
  )
  
  return(ht)
}



#───────────────────────────────────────────────────────────────────────────────
#  signatures from shared DEGs functions ()
#───────────────────────────────────────────────────────────────────────────────




# Helper to map Entrez → HUGO
entrez_to_hugo <- function(entrez_ids) {
  suppressMessages({
    mapped <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys     = entrez_ids,
      keytype  = "ENTREZID",
      columns  = "SYMBOL"
    )
  })
  unique(na.omit(mapped$SYMBOL))
}

# Requires: extract_signature_genes() from earlier + your entrez_to_hugo() helper
# Also assumes 'negative_ora_df' and 'signatures' exist.

# --- map per-signature to HUGO (dedupe within signature) ---------------------
to_hugo_by_signature <- function(ora_df, signatures) {
  # 1) get raw tokens per signature (split on '/')
  raw_res <- extract_signature_genes(ora_df, signatures)
  raw_list <- raw_res$genes_by_signature
  
  is_entrez <- function(x) grepl("^[0-9]+$", x)  # purely digits
  
  hugo_list <- lapply(raw_list, function(tokens) {
    tokens <- unique(trimws(tokens))
    if (length(tokens) == 0) return(character(0))
    
    entrez_tokens <- tokens[is_entrez(tokens)]
    symbol_like   <- tokens[!is_entrez(tokens)]   # already looks like a symbol
    
    mapped_symbols <- character(0)
    if (length(entrez_tokens) > 0) {
      mapped_symbols <- entrez_to_hugo(entrez_tokens)
      mapped_symbols <- unique(trimws(mapped_symbols))
    }
    
    # union within signature, then dedupe
    unique(c(symbol_like, mapped_symbols))
  })
  
  # tidy long
  tidy_hugo <- do.call(rbind, lapply(names(hugo_list), function(sig) {
    if (length(hugo_list[[sig]]) == 0) return(NULL)
    data.frame(signature = sig, gene = hugo_list[[sig]], stringsAsFactors = FALSE)
  }))
  
  list(genes_by_signature_hugo = hugo_list, tidy_hugo = tidy_hugo)
}

# --- Function: extract unique gene tokens per signature ----------------------
extract_signature_genes <- function(ora_df, signatures, id_col = "ID", gene_col = "geneID") {
  stopifnot(all(c(id_col, gene_col) %in% colnames(ora_df)))
  
  # helper: split "a/b/c", trim, drop empties, unique
  split_genes <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(character(0))
    parts <- unlist(strsplit(x, "/", fixed = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts[nzchar(parts)]
  }
  
  out_list <- lapply(names(signatures), function(sig_name) {
    ids <- signatures[[sig_name]]
    subdf <- ora_df[ora_df[[id_col]] %in% ids, , drop = FALSE]
    if (nrow(subdf) == 0) {
      # no matching IDs in this df
      return(character(0))
    }
    uniq_genes <- unique(split_genes(subdf[[gene_col]]))
    uniq_genes
  })
  names(out_list) <- names(signatures)
  
  # also produce a tidy long data frame (signature, gene)
  tidy_df <- do.call(rbind, lapply(names(out_list), function(sig) {
    if (length(out_list[[sig]]) == 0) return(NULL)
    data.frame(signature = sig, gene = out_list[[sig]], stringsAsFactors = FALSE)
  }))
  
  list(genes_by_signature = out_list, tidy = tidy_df)
}

overlaps_for_category <- function(category = NULL, genes, species = "Homo sapiens") {
  # choose msigdbr source depending on category
  msig <- if (is.null(category)) {
    msigdbr(species = species)
  } else {
    msigdbr(species = species, category = category)
  }
  
  msig |>
    mutate(gene_symbol = toupper(gene_symbol)) |>
    group_by(gs_exact_source, gs_name) |>
    summarise(
      overlap       = n_distinct(intersect(gene_symbol, genes)),
      size          = n(),
      genes_matched = paste(sort(intersect(gene_symbol, genes)), collapse = ", "),
      .groups = "drop"
    ) |>
    arrange(desc(overlap))
}


compute_overlaps_for_sets <- function(genes, set_keys, species = "Homo sapiens") {
  set_keysU <- toupper(trimws(unique(set_keys)))
  
  msig <- msigdbr(species = species) |>
    mutate(
      GS_NAME_U   = toupper(gs_name),
      GS_ID_U     = toupper(gs_exact_source),
      gene_symbol = toupper(gene_symbol)
    ) |>
    # keep only what we need; tolerate missing cols across msigdbr versions
    dplyr::select(any_of(c(
      "gs_cat","gs_subcat","gs_name","gs_exact_source",
      "gene_symbol","GS_NAME_U","GS_ID_U"
    )))
  
  msig_sel <- msig |>
    filter(GS_NAME_U %in% set_keysU | GS_ID_U %in% set_keysU)
  
  # helpful notice for typos/missing keys
  found_keys <- unique(c(msig_sel$GS_NAME_U, msig_sel$GS_ID_U))
  missing <- setdiff(set_keysU, found_keys)
  if (length(missing)) message("Not found: ", paste(missing, collapse = ", "))
  
  if (nrow(msig_sel) == 0) {
    return(tibble(
      gs_cat = character(), gs_subcat = character(),
      gs_name = character(), gs_exact_source = character(),
      set_size = integer(), hits = integer(),
      genes_matched = character(),
      query_size = integer(), pct_of_set = numeric(), pct_of_query = numeric()
    ))
  }
  
  msig_sel |>
    group_by(across(any_of(c("gs_cat","gs_subcat","gs_name","gs_exact_source")))) |>
    summarise(
      set_size      = n_distinct(gene_symbol),
      hits          = n_distinct(intersect(gene_symbol, genes)),
      genes_matched = paste(sort(intersect(gene_symbol, genes)), collapse = ", "),
      .groups = "drop"
    ) |>
    mutate(
      query_size   = length(genes),
      pct_of_set   = ifelse(set_size > 0, hits / set_size, NA_real_),
      pct_of_query = ifelse(query_size > 0, hits / query_size, NA_real_)
    ) |>
    arrange(desc(hits), desc(pct_of_set), gs_name)
}

add_mappings <- function(overlaps_df, mapping_df, sig_tbl) {
  new_rows <- overlaps_df |>
    inner_join(mapping_df, by = "gs_name") |>
    filter(!is.na(genes_matched), genes_matched != "") |>
    separate_rows(genes_matched, sep = "\\s*,\\s*") |>
    transmute(signature, gene = trimws(genes_matched)) |>
    distinct(signature, gene)
  
  bind_rows(
    sig_tbl,
    anti_join(new_rows, sig_tbl, by = c("signature","gene"))
  )
}

# Power-match genes to a single signature and return a list of gene vectors per signature
power_match_sig_list <- function(sig_gene_hugo_long) {
  stopifnot(all(c("signature", "gene") %in% colnames(sig_gene_hugo_long)))
  
  # 0) Clean & coerce
  x <- sig_gene_hugo_long[, c("signature", "gene")]
  x$signature <- as.character(x$signature)
  x$gene      <- as.character(x$gene)
  x <- x[!is.na(x$signature) & !is.na(x$gene) & nzchar(x$signature) & nzchar(x$gene), , drop = FALSE]
  
  # 1) De-duplicate pairs
  df <- unique(x)
  
  # 2) Original signature sizes (unique genes per signature in the *original* table)
  tmp <- aggregate(gene ~ signature, data = x, FUN = function(v) length(unique(v)))
  orig_counts <- setNames(as.integer(tmp$gene), tmp$signature)
  
  # 3) Prepare accumulators
  genes <- sort(unique(df$gene))
  assigned_counts <- setNames(integer(length(orig_counts)), names(orig_counts))
  assigned_sig <- character(length(genes))  # parallel to genes
  
  # 4) Assign each gene to exactly one signature by "power matching"
  for (i in seq_along(genes)) {
    g <- genes[i]
    sigs_g <- unique(df$signature[df$gene == g])
    
    if (length(sigs_g) == 1L) {
      pick <- sigs_g
    } else {
      # choose the signature with the smallest original size
      sizes <- orig_counts[sigs_g]
      min_size <- min(sizes, na.rm = TRUE)
      candidates <- sigs_g[sizes == min_size]
      
      if (length(candidates) > 1L) {
        # tie-break 1: currently smallest assigned count
        cur <- assigned_counts[candidates]
        min_cur <- min(cur, na.rm = TRUE)
        candidates <- candidates[cur == min_cur]
        
        # tie-break 2: alphabetical (deterministic)
        pick <- sort(candidates)[1]
      } else {
        pick <- candidates[1]
      }
    }
    
    assigned_sig[i] <- pick
    assigned_counts[pick] <- assigned_counts[pick] + 1L
  }
  
  # 5) Build unique, power-matched output and return list-of-vectors
  out_df <- data.frame(
    signature = assigned_sig,
    gene = genes,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  split(out_df$gene, out_df$signature)
}




#───────────────────────────────────────────────────────────────────────────────
#  shared DEGs composition, boxplots and scatter
#───────────────────────────────────────────────────────────────────────────────


get_deg_signature_composition_plot <- function(
    sig_list,
    col_vec = c(
      "Amino_Starvation_ISR"          = "#FF5A78",
      "Apoptosis_Resistance"          = "#9C77B3",
      "MitoROS_Module"                = "#FF9BB3", 
      "Endomembrane_Trafficking"      = "#fb8500",
      "Mito_Function_Maintenance"     = "#ffb703",
      "NFkB_Immune_Activation"        = "#219ebc",
      "Extrinsic_Apoptosis"           = "#023047",
      "Other"                         = "grey"
    ),
    total_gene_num = list(n_positive = 75, n_negative = 175),
    title = "Signature Breakdown of Shared Multi-lineage DEGs in CD59 Deficiency",
    base_size = 12,
    # font styling (match your other plot)
    title_size        = 16,
    y_title_size      = 14,
    y_text_size       = 14,
    legend_title_size = 14,
    legend_text_size  = 14,
    # explicit orders (bottom → top)
    ndegs_order  = c(
      "NFkB_Immune_Activation",      # bottom
      "Mito_Function_Maintenance",
      "Endomembrane_Trafficking",
      "Extrinsic_Apoptosis"          # top
    ), 
    pdegs_order  = c(
      "MitoROS_Module",
      "Apoptosis_Resistance",
      "Amino_Starvation_ISR"         # top
    ),
    legend_order = c("Extrinsic_Apoptosis", "Endomembrane_Trafficking",
                     "Mito_Function_Maintenance", "NFkB_Immune_Activation",
                     "Apoptosis_Resistance", "Other"),
    # NEW: pretty legend display names
    legend_names = c(
      "Amino_Starvation_ISR"     = "Metabolic Stress Response",
      "Apoptosis_Resistance"     = "Anti-Apoptotic Signaling",
      "MitoROS_Module"           = "Mitochondrial ROS pathway",
      "Endomembrane_Trafficking" = "Vesicular Transport",
      "Mito_Function_Maintenance"= "Mitochondrial Homeostasis",
      "NFkB_Immune_Activation"   = "NF-κB–Inflammatory Response",
      "Extrinsic_Apoptosis"      = "Extrinsic Apoptosis",
      "Other"                    = "Other"
    )
) {
  # ── Dependencies ───────────────────────────────────────────────────────────
  req_pkgs <- c("ggplot2", "dplyr", "tidyr", "forcats", "tibble", "scales")
  miss <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  
  # ── Sanity checks ──────────────────────────────────────────────────────────
  if (!is.list(sig_list))
    stop("sig_list must be a list containing 'positive_sigs' and 'negative_sigs'.", call. = FALSE)
  if (!all(c("positive_sigs", "negative_sigs") %in% names(sig_list)))
    stop("sig_list must contain both 'positive_sigs' and 'negative_sigs'.", call. = FALSE)
  if (!is.list(sig_list$positive_sigs) || !is.list(sig_list$negative_sigs))
    stop("sig_list$positive_sigs and sig_list$negative_sigs must each be lists of signatures.", call. = FALSE)
  if (!is.list(total_gene_num) ||
      !all(c("n_positive", "n_negative") %in% names(total_gene_num)))
    stop("total_gene_num must be a list with named elements 'n_positive' and 'n_negative'.", call. = FALSE)
  if (!is.numeric(total_gene_num$n_positive) || !is.numeric(total_gene_num$n_negative))
    stop("n_positive and n_negative must be numeric values.", call. = FALSE)
  if (!is.character(col_vec) || is.null(names(col_vec)))
    stop("col_vec must be a named character vector mapping signatures to colors.", call. = FALSE)
  
  # ── Helper: build data frame from sig_list ─────────────────────────────────
  get_sig_df <- function(sig_sublist, direction) {
    tibble::tibble(
      signature = rep(names(sig_sublist), times = vapply(sig_sublist, length, integer(1))),
      gene      = unlist(sig_sublist),
      direction = direction
    )
  }
  
  df_pos <- get_sig_df(sig_list$positive_sigs, "Shared pDEGs")
  df_neg <- get_sig_df(sig_list$negative_sigs, "Shared nDEGs")
  df <- dplyr::bind_rows(df_pos, df_neg)
  
  # totals map
  totals_map <- c(
    "Shared pDEGs" = total_gene_num$n_positive,
    "Shared nDEGs" = total_gene_num$n_negative
  )
  
  # ── Calculate composition (counts) ─────────────────────────────────────────
  comp_df <- df %>%
    dplyr::group_by(direction, signature) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(total = as.numeric(totals_map[direction]))
  
  # ── Add 'Other' category for uncovered genes ───────────────────────────────
  comp_df <- comp_df %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(
      covered = sum(n),
      other_n = pmax(0, total[1] - covered[1])
    ) %>%
    dplyr::ungroup()
  
  others <- comp_df %>%
    dplyr::distinct(direction, total, other_n) %>%
    dplyr::filter(other_n > 0) %>%
    dplyr::mutate(signature = "Other", n = other_n) %>%
    dplyr::select(direction, signature, n, total)
  
  comp_df <- dplyr::bind_rows(comp_df, others)
  
  # ── Build per-bar stacking order (respect provided bottom→top order) ───────
  ndegs_present <- comp_df %>%
    dplyr::filter(direction == "Shared nDEGs", signature != "Other", n > 0) %>%
    dplyr::arrange(dplyr::desc(n)) %>% dplyr::pull(signature) %>% unique()
  if (!is.null(ndegs_order)) {
    ndegs_present <- c(ndegs_order[ndegs_order %in% ndegs_present],
                       setdiff(ndegs_present, ndegs_order))
  }
  levels_n <- c("Other", ndegs_present)
  
  pdegs_present <- comp_df %>%
    dplyr::filter(direction == "Shared pDEGs", signature != "Other", n > 0) %>%
    dplyr::arrange(dplyr::desc(n)) %>% dplyr::pull(signature) %>% unique()
  if (!is.null(pdegs_order)) {
    pdegs_present <- c(pdegs_order[pdegs_order %in% pdegs_present],
                       setdiff(pdegs_present, pdegs_order))
  }
  levels_p <- c("Other", pdegs_present)
  
  # Per-bar data with their own factor orders
  comp_n <- comp_df %>% dplyr::filter(direction == "Shared nDEGs")
  comp_p <- comp_df %>% dplyr::filter(direction == "Shared pDEGs")
  comp_n$signature <- factor(comp_n$signature, levels = rev(levels_n))
  comp_p$signature <- factor(comp_p$signature, levels = rev(levels_p))
  # explicit numeric stack order (bottom→top)
  comp_n$stack_ord <- match(as.character(comp_n$signature), rev(levels_n))
  comp_p$stack_ord <- match(as.character(comp_p$signature), rev(levels_p))
  
  # Direction order on x
  comp_df$direction <- factor(comp_df$direction, levels = c("Shared nDEGs", "Shared pDEGs"))
  
  # Legend order
  if (is.null(legend_order)) {
    n_legend <- comp_df %>%
      dplyr::filter(direction == "Shared nDEGs", signature != "Other", n > 0) %>%
      dplyr::arrange(n) %>% dplyr::pull(signature) %>% unique()
    p_legend <- pdegs_present
    legend_levels <- c(n_legend, setdiff(p_legend, n_legend), "Other")
  } else {
    legend_levels <- legend_order
    legend_levels <- legend_levels[legend_levels %in% unique(as.character(comp_df$signature))]
    if ("Other" %in% legend_levels) {
      legend_levels <- c(setdiff(legend_levels, "Other"), "Other")
    }
  }
  # Legend labels (use pretty names if provided; fall back to raw names)
  legend_labels <- {
    labs <- legend_names[legend_levels]
    labs[is.na(labs)] <- legend_levels[is.na(labs)]
    unname(labs)
  }
  
  # totals for top labels
  totals_df <- comp_df %>% dplyr::distinct(direction, total)
  
  # y-axis breaks
  y_breaks <- sort(unique(c(0, 50, 100, 150, total_gene_num$n_negative)))
  
  # ── Plot ───────────────────────────────────────────────────────────────────
  ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = comp_n,
      ggplot2::aes(x = direction, y = n, fill = signature, order = stack_ord),
      stat = "identity", width = 0.70, color = "white", linewidth = 0.35
    ) +
    ggplot2::geom_bar(
      data = comp_p,
      ggplot2::aes(x = direction, y = n, fill = signature, order = stack_ord),
      stat = "identity", width = 0.70, color = "white", linewidth = 0.35
    ) +
    ggplot2::scale_fill_manual(
      values = col_vec,
      breaks = legend_levels,
      labels = legend_labels,
      name   = "Signature"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      minor_breaks = NULL,
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::geom_text(
      data = totals_df,
      mapping = ggplot2::aes(x = direction, y = total, label = total),
      vjust = -0.6, fontface = "bold", size = 4,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Number of genes"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "panel",
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = title_size,
        margin = ggplot2::margin(t = 10, b = 4)
      ),
      axis.title.y        = ggplot2::element_text(size = y_title_size),
      axis.text.y         = ggplot2::element_text(size = y_text_size),
      axis.text.x         = ggplot2::element_text(size = y_text_size),
      legend.title        = ggplot2::element_text(size = legend_title_size),
      legend.text         = ggplot2::element_text(size = legend_text_size),
      legend.position     = "right",
      legend.title.align  = 0.5,
      panel.grid.major.x  = ggplot2::element_blank(),
      panel.grid.major.y  = ggplot2::element_line(color = "#e5e7eb", linewidth = 0.4),
      panel.grid.minor    = ggplot2::element_blank(),
      axis.line.y         = ggplot2::element_line(linewidth = 0.4, colour = "#2b2b2b"),
      plot.margin         = ggplot2::margin(t = 6, r = 10, b = 6, l = 8)
    )
}



plot_shared_degs_signatures_boxplots <- function(
    df,
    cell_col        = "cell_type",
    sig_col         = "sig_name",
    condition_col   = "condition",
    value_col       = "sig_score",
    
    # legacy controls (still supported)
    cell_order      = NULL,
    sig_order       = NULL,
    
    # explicit ordering controls
    row_order       = c("Apoptosis_Resistance",
                        "Extrinsic_Apoptosis", "NFkB_Immune_Activation", "Mito_Function_Maintenance",
                        "Endomembrane_Trafficking"),   # order of signatures (rows)
    col_order       = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),   # order of cell types (columns)
    row_decreasing  = FALSE,  # reverse row order
    col_decreasing  = FALSE,  # reverse column order
    
    # rename rows (facet strip labels); named character vector old -> new
    row_renames     = c(
      "Apoptosis_Resistance"     = "Anti-apoptotic Signaling",
      "Endomembrane_Trafficking" = "Vesicular Transport",
      "Mito_Function_Maintenance"= "Mitochondrial Homeostasis",
      "NFkB_Immune_Activation"   = "NF-κB–Inflammatory Response",
      "Extrinsic_Apoptosis"      = "Extrinsic Apoptosis"
    ),
    
    condition_order = c("CD59.Def","MS","Ab.mediated.IDD"),
    palette         = c("MS" = "#012A3A",
                        "Ab.mediated.IDD" = "#2C7FB8",
                        "CD59.Def" = "#E66100"),
    
    # NEW: legend condition label renaming (old level -> new label)
    condition_renames =  c(
      "CD59.Def"        = "CD59 Def.",
      "MS"              = "Multiple Sclerosis",
      "Ab.mediated.IDD" = "Ab-mediated IDD"
    ),
    
    share_y       = TRUE,
    show_points   = FALSE,
    y_limits      = NULL,
    
    # ---- Plot title & legend title controls ----------------------------------
    title               = "Signature AUC Scores from Shared DEGs Across Conditions", 
    legend_title        = "Condition",
    title_size          = 20,
    cell_type_size      = 17,
    sig_size            = 17,
    y_text_size         = 17,
    legend_title_size   = 17,
    legend_text_size    = 17,
    title_face          = "bold",
    title_hjust         = 0.5,
    
    row_gap       = 2.5,        # gap between rows
    row_gap_unit  = "lines",    # gap unit
    
    # horizontal guide at y = 0.5 in every panel
    add_row_y_line      = TRUE,
    row_y_line          = 0.5,
    row_y_line_color    = "grey60",
    row_y_line_size     = 0.9,
    row_y_line_linetype = "dashed"
) {
  stopifnot(is.data.frame(df))
  req_cols <- c(cell_col, sig_col, condition_col, value_col)
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))
  if (!is.numeric(df[[value_col]]))
    stop(sprintf("`%s` must be numeric.", value_col))
  
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(forcats)
  })
  
  # ----- Resolve desired orders ----------------------------------------------
  desired_row <- if (!is.null(row_order)) row_order else sig_order
  desired_col <- if (!is.null(col_order)) col_order else cell_order
  
  if (isTRUE(row_decreasing) && !is.null(desired_row)) desired_row <- rev(desired_row)
  if (isTRUE(col_decreasing) && !is.null(desired_col)) desired_col <- rev(desired_col)
  
  if (!is.null(desired_row)) {
    not_in_rows <- setdiff(desired_row, unique(df[[sig_col]]))
    if (length(not_in_rows))
      warning("These 'row_order' entries were not found in ", sig_col, ": ",
              paste(not_in_rows, collapse = ", "))
  }
  if (!is.null(desired_col)) {
    not_in_cols <- setdiff(desired_col, unique(df[[cell_col]]))
    if (length(not_in_cols))
      warning("These 'col_order' entries were not found in ", cell_col, ": ",
              paste(not_in_cols, collapse = ", "))
  }
  
  # ----- Factorize with requested orders -------------------------------------
  df <- df %>%
    mutate(
      !!sig_col       := if (is.null(desired_row)) .data[[sig_col]] else
        fct_relevel(.data[[sig_col]], desired_row),
      !!cell_col      := if (is.null(desired_col)) .data[[cell_col]] else
        fct_relevel(.data[[cell_col]], desired_col),
      !!condition_col := if (is.null(condition_order)) .data[[condition_col]] else
        fct_relevel(.data[[condition_col]], condition_order)
    )
  
  # ----- Optional row renaming (facet strip labels) ---------------------------
  if (!is.null(row_renames)) {
    if (!is.character(row_renames) || is.null(names(row_renames)) || any(names(row_renames) == "")) {
      stop("'row_renames' must be a NAMED character vector like c('old'='New Label', ...).")
    }
    sig_levels <- levels(factor(df[[sig_col]]))
    if (!setequal(sig_levels, names(row_renames))) {
      stop("`row_renames` must have names for ALL current signature names. ",
           "Expected names: ", paste(sig_levels, collapse = ", "))
    }
    renames_in_order <- unname(row_renames[sig_levels])
    df[[sig_col]] <- forcats::fct_relabel(factor(df[[sig_col]], levels = sig_levels),
                                          function(x) renames_in_order[match(x, sig_levels)])
  }
  
  # ----- Build plot -----------------------------------------------------------
  p <- ggplot(
    df,
    aes(
      x      = .data[[condition_col]],
      y      = .data[[value_col]],
      colour = .data[[condition_col]]
    )
  )
  
  if (isTRUE(add_row_y_line)) {
    p <- p +
      geom_hline(
        yintercept = row_y_line,
        colour     = row_y_line_color,
        linewidth  = row_y_line_size,
        linetype   = row_y_line_linetype,
        inherit.aes = FALSE
      )
  }
  
  p <- p +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = NA, size = 0.8)
  
  if (show_points) {
    p <- p + geom_point(
      position = position_jitter(width = 0.07, height = 0, seed = 1),
      size = 0.6, alpha = 0.6, color = "black"
    )
  }
  
  # Determine legend breaks in use
  legend_breaks <- if (!is.null(condition_order)) condition_order else levels(df[[condition_col]])
  
  # Build legend labels (apply renames if provided)
  legend_labels <- legend_breaks
  if (!is.null(condition_renames)) {
    if (!is.character(condition_renames) || is.null(names(condition_renames)) ||
        any(names(condition_renames) == "")) {
      stop("`condition_renames` must be a NAMED character vector like c('old'='New Label', ...).")
    }
    missing_names <- setdiff(legend_breaks, names(condition_renames))
    if (length(missing_names)) {
      stop("`condition_renames` must provide labels for ALL legend entries. Missing: ",
           paste(missing_names, collapse = ", "))
    }
    legend_labels <- unname(condition_renames[legend_breaks])
  }
  
  p <- p +
    facet_grid(
      rows   = vars(.data[[sig_col]]),
      cols   = vars(.data[[cell_col]]),
      scales = if (share_y) "fixed" else "free_y",
      switch = "y",
      drop   = FALSE,
      labeller = labeller(.rows = label_wrap_gen(width = 25))
    ) +
    scale_color_manual(values = palette,
                       breaks = legend_breaks,
                       labels = legend_labels,
                       drop   = FALSE,
                       name   = legend_title)
  
  # y scale / ensure baseline visible
  if (!is.null(y_limits)) {
    y_limits <- range(c(y_limits, row_y_line))
    p <- p + scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.02, 0.05)))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))
    if (isTRUE(add_row_y_line)) p <- p + expand_limits(y = row_y_line)
  }
  
  p <- p +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_line(color = "#efefef"),
      panel.grid.minor   = element_blank(),
      strip.placement    = "outside",
      strip.background   = element_blank(),
      strip.text.y.left  = element_text(size = sig_size, face = "bold", angle = 0, hjust = 1),
      strip.text.x       = element_text(size = cell_type_size, face = "bold"),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = y_text_size),
      panel.spacing.y    = grid::unit(row_gap, row_gap_unit),
      panel.spacing.x    = grid::unit(0.8, "lines"),
      legend.position    = "right",
      legend.box         = "vertical",
      legend.title       = element_text(size = legend_title_size, face = "bold"),
      legend.text        = element_text(size = legend_text_size),
      plot.title         = element_text(size = title_size, face = title_face, hjust = title_hjust)
    )
  
  return(p)
}




auc_scores_scatter_across_condtions <- function(
    df,
    
    # which signatures to plot (x vs y)
    x_sig,
    y_sig,
    # NEW: axis label override (vector). Examples:
    # c("X axis title","Y axis title")  OR  c(x="X axis", y="Y axis")
    axis_labels = NULL, 
    
    # column names (keep defaults for your sig_score_long_df)
    cell_col      = "cell_type",
    condition_col = "condition",
    sig_col       = "sig_name",
    score_col     = "sig_score",
    sample_col    = "sample",
    dataset_col   = "dataset",
    
    # visuals
    origin        = 0.5,
    condition_palette = c("MS" = "#012A3A",
                          "Ab.mediated.IDD" = "#2C7FB8",
                          "CD59.Def" = "#E66100"),
    shape_map = c("NKs" = 21, "Monocytes" = 22, "B Cells" = 23,
                  "CD4+ T Cells" = 25, "CD8+ T Cells" = 24),
    
    # size controls (match your previous plot’s feel)
    base_size         = 20,
    title             = "Mito_Function_Maintenance vs.MitoROS_Module (All Cell Types)",
    title_size        = 24,
    axis_title_size   = 20,
    axis_text_size    = 16,
    legend_title_size = 16,
    legend_text_size  = 14,
    point_size        = 3,
    stroke            = 0.5,
    alpha             = 0.9,
    legend_position   = "bottom",
    equal_aspect      = TRUE,
    
    
    
    # optional fixed limits; leave NULL for auto
    xlim = NULL,
    ylim = NULL
) {
  # --- deps ---
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install 'tidyr'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  dplyr <- asNamespace("dplyr"); tidyr <- asNamespace("tidyr"); ggplot2 <- asNamespace("ggplot2")
  sym <- function(x) rlang::sym(x)
  
  # --- sanity checks ---
  needed <- c(cell_col, condition_col, sig_col, score_col, sample_col)
  miss <- setdiff(needed, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  sigs_present <- unique(df[[sig_col]])
  if (!all(c(x_sig, y_sig) %in% sigs_present)) {
    stop("Requested signatures not found: ",
         paste(setdiff(c(x_sig, y_sig), sigs_present), collapse = ", "))
  }
  
  # --- wide table: one row per sample/cell/condition (x = x_sig score, y = y_sig score) ---
  plot_df <- df |>
    dplyr::filter(.data[[sig_col]] %in% c(x_sig, y_sig)) |>
    dplyr::select(
      !!sym(sample_col), !!sym(cell_col), !!sym(condition_col),
      dplyr::any_of(dataset_col), !!sym(sig_col), !!sym(score_col)
    ) |>
    tidyr::pivot_wider(names_from = !!sym(sig_col), values_from = !!sym(score_col)) |>
    dplyr::rename(x = !!sym(x_sig), y = !!sym(y_sig)) |>
    dplyr::filter(!is.na(x) & !is.na(y)) |>
    dplyr::mutate(
      .cond = factor(.data[[condition_col]], levels = names(condition_palette)),
      .cell = factor(.data[[cell_col]],      levels = names(shape_map))
    )
  
  if (nrow(plot_df) == 0) stop("No complete cases to plot after reshaping.")
  
  # --- FIXED symmetric limits 0..1 on both axes ---
  if (is.null(xlim)) xlim <- c(0, 1)
  if (is.null(ylim)) ylim <- c(0, 1)
  
  # --- resolve axis labels (default -> overrides via vector) ---
  x_lab <- paste0(x_sig, " score")
  y_lab <- paste0(y_sig, " score")
  if (!is.null(axis_labels)) {
    if (!is.null(names(axis_labels))) {
      if (!is.null(axis_labels[["x"]])) x_lab <- axis_labels[["x"]]
      if (!is.null(axis_labels[["y"]])) y_lab <- axis_labels[["y"]]
    } else {
      if (length(axis_labels) >= 1) x_lab <- axis_labels[[1]]
      if (length(axis_labels) >= 2) y_lab <- axis_labels[[2]]
    }
  }
  
  # --- plot ---
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = x, y = y, fill = .cond, shape = .cell)
  ) +
    ggplot2::geom_point(color = "black", size = point_size, stroke = stroke, alpha = alpha) +
    ggplot2::geom_vline(xintercept = origin, colour = "grey40", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = origin, colour = "grey40", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = condition_palette, name = "Condition") +
    ggplot2::scale_shape_manual(values = shape_map, name = "Cell Type") +
    {
      if (equal_aspect)
        ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = TRUE)
      else
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE)
    } +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.title       = ggplot2::element_text(size = axis_title_size),
      axis.text        = ggplot2::element_text(size = axis_text_size),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = legend_position,
      legend.title     = ggplot2::element_text(size = legend_title_size),
      legend.text      = ggplot2::element_text(size = legend_text_size),
      plot.margin      = ggplot2::margin(10, 20, 10, 10)
    ) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(override.aes = list(shape = 21, colour = "black",
                                                        stroke = stroke, size = point_size + 2)),
      shape = ggplot2::guide_legend(override.aes = list(fill = "grey85", colour = "black",
                                                        stroke = stroke, size = point_size + 2))
    )
  
  return(p)
}




#───────────────────────────────────────────────────────────────────────────────
#  shared DEGs heatmap
#───────────────────────────────────────────────────────────────────────────────


get_msigdb_genes_from_gene_set_name <- function(gs_name,
                             category = NULL,
                             species  = "Homo sapiens",
                             msig     = NULL) {
  # ---- sanity checks on inputs ----
  if (!is.character(gs_name) || length(gs_name) != 1L || !nzchar(gs_name)) {
    stop("`gs_name` must be a non-empty character scalar.")
  }
  if (!is.null(category) && (!is.character(category) || length(category) != 1L)) {
    stop("`category` must be NULL or a character scalar (e.g., 'C2', 'H', 'C5').")
  }
  if (!is.character(species) || length(species) != 1L || !nzchar(species)) {
    stop("`species` must be a non-empty character scalar (e.g., 'Homo sapiens').")
  }
  
  # ---- helper to validate a provided msig data.frame ----
  .validate_msig <- function(x) {
    if (!is.data.frame(x)) stop("`msig` must be a data.frame.")
    req_cols <- c("gs_name", "gene_symbol")
    missing  <- setdiff(req_cols, names(x))
    if (length(missing)) {
      stop("`msig` is missing required column(s): ",
           paste(missing, collapse = ", "),
           ". Expect at least: ", paste(req_cols, collapse = ", "), ".")
    }
    # Coerce factor -> character for key columns
    for (cc in c("gs_name", "gene_symbol", "species_name", "gs_cat")) {
      if (cc %in% names(x) && is.factor(x[[cc]])) x[[cc]] <- as.character(x[[cc]])
    }
    if (!is.character(x$gs_name)) {
      stop("`msig$gs_name` must be character (or coercible from factor).")
    }
    if (!is.character(x$gene_symbol)) {
      stop("`msig$gene_symbol` must be character (or coercible from factor).")
    }
    x
  }
  
  # ---- get or use msig table ----
  if (is.null(msig)) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Package 'msigdbr' is required when `msig` is not provided.")
    }
    msig <- if (is.null(category)) {
      msigdbr::msigdbr(species = species)
    } else {
      msigdbr::msigdbr(species = species, category = category)
    }
  } else {
    msig <- .validate_msig(msig)
    
    # Optional species filter (if column exists)
    if ("species_name" %in% names(msig)) {
      sp_levels <- unique(msig$species_name)
      if (!species %in% sp_levels) {
        warning("`msig` contains species_name(s) {",
                paste(stats::na.omit(sp_levels), collapse = ", "),
                "} but requested `species` is '", species,
                "'. Proceeding with rows matching `species` only (may be empty).")
      }
      msig <- msig[is.na(msig$species_name) | msig$species_name == species, , drop = FALSE]
    }
    
    # Optional category filter (if both provided & column exists)
    if (!is.null(category) && "gs_cat" %in% names(msig)) {
      msig <- msig[msig$gs_cat == category, , drop = FALSE]
    }
  }
  
  # ---- extract genes for the requested gs_name ----
  if (!nrow(msig)) {
    warning("MSigDB table is empty after filtering; returning character(0).")
    return(character(0))
  }
  
  hits <- msig[msig$gs_name == gs_name, , drop = FALSE]
  
  if (nrow(hits) == 0L) {
    # Helpful message for common mismatch
    if (!is.null(category) && !"gs_cat" %in% names(msig)) {
      warning("No rows found for gs_name = '", gs_name,
              "'. (Note: `msig` has no 'gs_cat' column to filter by category.) Returning character(0).")
    } else {
      warning("No rows found for gs_name = '", gs_name, "'. Returning character(0).")
    }
    return(character(0))
  }
  
  # Normalize to HUGO-like uppercase symbols
  syms <- hits$gene_symbol
  syms <- toupper(trimws(as.character(syms)))
  syms <- syms[nzchar(syms) & !is.na(syms)]
  
  unique(syms)
}



get_gene_sets_unique <- function(gene_sets_df, powermatching) {
  # ── sanity checks ───────────────────────────────────────────────────────────
  stopifnot(is.data.frame(gene_sets_df))
  req_cols <- c("gene", "gene_set")
  miss <- setdiff(req_cols, names(gene_sets_df))
  if (length(miss)) stop("`gene_sets_df` missing required columns: ", paste(miss, collapse = ", "))
  if (!is.character(gene_sets_df$gene) || !is.character(gene_sets_df$gene_set)) {
    stop("`gene_sets_df$gene` and `gene_sets_df$gene_set` must be character.")
  }
  if (!is.character(powermatching) || length(powermatching) == 0L) {
    stop("`powermatching` must be a non-empty character vector of set names ordered by priority (high→low).")
  }
  
  # Normalize whitespace
  gene_sets_df$gene     <- trimws(gene_sets_df$gene)
  gene_sets_df$gene_set <- trimws(gene_sets_df$gene_set)
  powermatching <- trimws(unique(powermatching))  # keep first occurrence as the higher rank
  
  # ── cross-checks (both directions) ─────────────────────────────────────────
  all_sets <- sort(unique(gene_sets_df$gene_set))
  missing_in_df <- setdiff(powermatching, all_sets)
  missing_in_pm <- setdiff(all_sets, powermatching)
  if (length(missing_in_df)) {
    warning("Sets in `powermatching` but not in `gene_sets_df`: ",
            paste(missing_in_df, collapse = ", "))
  }
  if (length(missing_in_pm)) {
    warning("Sets in `gene_sets_df` but not in `powermatching` (treated as lowest priority): ",
            paste(missing_in_pm, collapse = ", "))
  }
  
  # ── build rank map (lower number = higher priority) ────────────────────────
  rank_map <- setNames(seq_along(powermatching), powermatching)
  # Any set not in powermatching gets the WORST rank (Inf)
  rank_of <- function(set) ifelse(set %in% names(rank_map), rank_map[[set]], Inf)
  
  # ── unique (gene, set) rows and collect per-gene candidate sets ────────────
  unique_rows  <- unique(gene_sets_df[, c("gene", "gene_set")])
  sets_by_gene <- tapply(unique_rows$gene_set, unique_rows$gene, function(v) sort(unique(v)))
  
  # ── assign each gene to its highest-ranked set (tie → lexicographic) ───────
  assigned <- lapply(names(sets_by_gene), function(g) {
    gs <- sets_by_gene[[g]]
    if (length(gs) == 1L) {
      data.frame(gene = g, gene_set = gs, stringsAsFactors = FALSE)
    } else {
      ranks <- vapply(gs, rank_of, numeric(1))
      best  <- gs[order(ranks, gs)][1]  # smallest rank first; tie → alphabetical
      data.frame(gene = g, gene_set = best, stringsAsFactors = FALSE)
    }
  })
  gene_sets_df_unique <- do.call(rbind, assigned)
  rownames(gene_sets_df_unique) <- NULL
  gene_sets_df_unique
}



# Build a long data frame for ComplexHeatmap from meta tables + a unique gene-set map
# Output columns: cell_type, sample, condition, data_set, gene, median_log2fc, gene_set
get_gene_set_headmap_df <- function(meta_table_list,
                                    gene_sets_df_unique,
                                    tables_to_ignore = NULL) {
  # ── sanity checks ───────────────────────────────────────────────────────────
  if (!is.list(meta_table_list) || is.null(names(meta_table_list)))
    stop("`meta_table_list` must be a *named* list of directories (cell type -> path).")
  if (!is.data.frame(gene_sets_df_unique))
    stop("`gene_sets_df_unique` must be a data.frame with columns: gene, gene_set.")
  req_cols <- c("gene", "gene_set")
  miss <- setdiff(req_cols, names(gene_sets_df_unique))
  if (length(miss))
    stop("`gene_sets_df_unique` is missing required column(s): ", paste(miss, collapse = ", "))
  if (!all(vapply(meta_table_list, dir.exists, logical(1)))) {
    bad <- names(meta_table_list)[!vapply(meta_table_list, dir.exists, logical(1))]
    stop("These directories do not exist: ", paste(bad, collapse = ", "))
  }
  if (!is.null(tables_to_ignore) && !is.character(tables_to_ignore))
    stop("`tables_to_ignore` must be NULL or a character vector of file names/paths to skip.")
  
  # standardize gene mapping (use uppercase keys for robust matching)
  gene_sets_df_unique$gene <- toupper(trimws(as.character(gene_sets_df_unique$gene)))
  gene_sets_df_unique$gene_set <- trimws(as.character(gene_sets_df_unique$gene_set))
  gene_to_set <- setNames(gene_sets_df_unique$gene_set, gene_sets_df_unique$gene)
  
  # ignore list: compare by basename and full path
  ignore_full <- if (is.null(tables_to_ignore)) character(0) else normalizePath(tables_to_ignore, mustWork = FALSE)
  ignore_base <- if (is.null(tables_to_ignore)) character(0) else basename(tables_to_ignore)
  
  out_list <- list()
  n_added  <- 0L
  
  for (cell_type in names(meta_table_list)) {
    dir_path <- meta_table_list[[cell_type]]
    
    files <- list.files(dir_path, pattern = "meta_table\\.qs$", full.names = TRUE)
    if (!length(files)) {
      warning("No meta_table .qs files found under: ", dir_path)
      next
    }
    
    for (fp in files) {
      bn <- basename(fp)
      # skip if in ignore list (by basename or by full path)
      if (bn %in% ignore_base || fp %in% ignore_full) next
      
      mt <- tryCatch(qs::qread(fp), error = function(e) {
        warning("Failed to read: ", fp, " (", conditionMessage(e), ")"); NULL
      })
      if (is.null(mt)) next
      if (!"median_log2FC" %in% colnames(mt)) {
        warning("Skipping (no 'median_log2FC' col): ", fp)
        next
      }
      
      # obtain gene symbols from rownames or 'gene' column
      rn <- rownames(mt)
      if ((is.null(rn) || all(rn == "")) && "gene" %in% colnames(mt)) rn <- as.character(mt$gene)
      if (is.null(rn) || !length(rn)) {
        warning("Skipping (no gene identifiers found): ", fp)
        next
      }
      
      df_small <- data.frame(
        gene = toupper(trimws(rn)),
        median_log2fc = mt$median_log2FC,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      
      # keep only genes that are in the provided unique mapping
      hits <- df_small$gene %in% names(gene_to_set)
      if (!any(hits)) next
      df_small <- df_small[hits, , drop = FALSE]
      df_small$gene_set <- gene_to_set[df_small$gene]
      
      # parse file name to extract condition, data_set, sample
      core <- sub("_meta_table\\.qs$", "", bn)
      parts <- strsplit(core, "_", fixed = TRUE)[[1]]
      condition <- if (length(parts) >= 1) parts[1] else NA_character_
      data_set  <- if (length(parts) >= 2) parts[2] else NA_character_
      sample    <- if (length(parts) >= 3) paste(parts[3:length(parts)], collapse = "_") else NA_character_
      
      out_list[[length(out_list) + 1L]] <- data.frame(
        cell_type     = cell_type,
        sample        = sample,
        condition     = condition,
        data_set      = data_set,
        gene          = df_small$gene,
        median_log2fc = df_small$median_log2fc,
        gene_set      = df_small$gene_set,
        stringsAsFactors = FALSE
      )
      n_added <- n_added + nrow(df_small)
    }
  }
  
  if (!length(out_list)) {
    warning("No matching genes were found across the provided meta tables.")
    return(data.frame(
      cell_type     = character(0),
      sample        = character(0),
      condition     = character(0),
      data_set      = character(0),
      gene          = character(0),
      median_log2fc = numeric(0),
      gene_set      = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  output_df <- do.call(rbind, out_list)
  rownames(output_df) <- NULL
  output_df
}






# Further reorder rows *inside* the two meta-effect groups.
# - "Up in target": by median across CD59 Def. columns (desc)
# - "Down in target":
#     (i) genes with MS median >= 0 AND >50% of MS columns > 0, ordered by MS median (desc)
#     (ii) remainder, ordered by Ab-mediated IDD median (desc)
# (unchanged)
reorder_within_meta_groups <- function(
    mat,
    meta,
    row_order,
    row_split,
    target_condition = "CD59 Def.",
    ms_condition     = "Multiple Sclerosis",
    idd_condition    = "Ab-mediated IDD"
){
  stopifnot(length(row_order) == length(row_split))
  
  cols_target <- meta$col_id[meta$condition == target_condition]
  cols_ms     <- meta$col_id[meta$condition == ms_condition]
  cols_idd    <- meta$col_id[meta$condition == idd_condition]
  
  rmed <- function(m) {
    if (is.null(m) || length(m) == 0L || ncol(m) == 0L) return(rep(NA_real_, nrow(m)))
    matrixStats::rowMedians(m, na.rm = TRUE)
  }
  
  out <- character(0)
  
  ## --- Up in target ----------------------------------------------------------
  up_rows <- row_order[!is.na(row_split) & row_split == "Up in target"]
  if (length(up_rows)) {
    if (length(cols_target)) {
      s <- rmed(mat[up_rows, cols_target, drop = FALSE])
      names(s) <- up_rows
      up_rows  <- names(sort(s, decreasing = TRUE, na.last = TRUE))
    }
    out <- c(out, up_rows)
  }
  
  ## --- Down in target --------------------------------------------------------
  dn_rows <- row_order[!is.na(row_split) & row_split == "Down in target"]
  if (length(dn_rows)) {
    
    # (i) MS median >= 0 AND >50% MS columns > 0 -> order by MS median (desc)
    if (length(cols_ms)) {
      s_ms <- rmed(mat[dn_rows, cols_ms, drop = FALSE]); names(s_ms) <- dn_rows
      p_ms <- rowMeans(mat[dn_rows, cols_ms, drop = FALSE] > 0, na.rm = TRUE)
      p_ms[!is.finite(p_ms)] <- 0
      keep1 <- dn_rows[is.finite(s_ms) & s_ms >= 0 & p_ms > 0.5]
      ord1  <- if (length(keep1)) names(sort(s_ms[keep1], decreasing = TRUE, na.last = TRUE)) else character(0)
    } else {
      ord1 <- character(0)
    }
    
    # (ii) among the remainder, use IDD median (desc) when finite
    remainder <- setdiff(dn_rows, ord1)
    if (length(remainder) && length(cols_idd)) {
      s_idd <- rmed(mat[remainder, cols_idd, drop = FALSE]); names(s_idd) <- remainder
      finite_idx <- remainder[is.finite(s_idd)]
      ord2 <- if (length(finite_idx)) names(sort(s_idd[finite_idx], decreasing = TRUE, na.last = TRUE)) else character(0)
      # (iii) leftovers (NA medians / no columns) -> very bottom
      ord3 <- setdiff(remainder, ord2)
    } else {
      ord2 <- character(0)
      ord3 <- remainder
    }
    
    out <- c(out, ord1, ord2, ord3)
  }
  
  missing <- setdiff(row_order, out)
  c(out, missing)
}





get_gene_level_heatmap <- function(
    df,
    title,
    meta_name_fontsize = 11,
    cell_pal = c(
      "NKs"          = "#6A4C93",
      "Monocytes"    = "#5E60CE",
      "B Cells"      = "#C1121F",
      "CD4+ T Cells" = "#1B998B",
      "CD8+ T Cells" = "#D4A373"
    ),
    condition_order   = c("CD59 Def.", "Multiple Sclerosis", "Ab-mediated IDD"),
    cell_type_order   = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
    min_prop_keep     = 0.85,
    condition_gap_mm  = 4,
    row_gap_mm        = 3,     # controls gap between Up/Down row slices
    target_condition  = "CD59 Def.",
    row_split_by_sign = TRUE,
    top_n             = NULL,
    use_olo           = TRUE,
    show_row_split_labels = FALSE  # <<< NEW: control displaying "Up in target"/"Down in target"
){
  have_seriation <- requireNamespace("seriation", quietly = TRUE)
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("Please install matrixStats.")
  
  if (!("dataset" %in% names(df)) && ("data_set" %in% names(df))) df$dataset <- df$data_set
  required_cols <- c("sample","gene","median_log2fc","dataset",
                     "condition","condition_subtype","age","sex","EDSS","cell_type")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse=", "))
  if (!is.numeric(df$median_log2fc)) stop("'median_log2fc' must be numeric")
  if (!is.numeric(df$age))  stop("'age' must be numeric")
  if (!is.numeric(df$EDSS)) stop("'EDSS' must be numeric")
  if (!is.character(df$gene) && !is.factor(df$gene)) stop("'gene' must be chr/fct")
  
  df <- df |>
    dplyr::mutate(
      sample    = as.character(sample),
      gene      = as.character(gene),
      cell_type = as.character(cell_type),
      col_id    = paste0(sample," | ",cell_type)
    )
  
  filt_meta <- df |>
    dplyr::arrange(sample, cell_type) |>
    dplyr::distinct(col_id, .keep_all = TRUE)
  
  mat_wide <- df |>
    dplyr::select(gene, col_id, median_log2fc) |>
    dplyr::group_by(gene, col_id) |>
    dplyr::summarise(value = mean(median_log2fc, na.rm=TRUE), .groups="drop") |>
    tidyr::pivot_wider(names_from = col_id, values_from = value, values_fill = NA_real_) |>
    dplyr::arrange(gene)
  
  gene_mat <- as.matrix(mat_wide[,-1, drop=FALSE])
  rownames(gene_mat) <- mat_wide$gene
  
  filt_meta$condition <- factor(filt_meta$condition, levels = condition_order)
  filt_meta$cell_type <- factor(filt_meta$cell_type, levels = cell_type_order)
  filt_meta <- filt_meta |> dplyr::arrange(condition, dataset, sample, cell_type)
  common_cols <- filt_meta$col_id
  gene_mat    <- gene_mat[, common_cols, drop=FALSE]
  filt_meta   <- filt_meta[match(common_cols, filt_meta$col_id), , drop=FALSE]
  
  fin_mask  <- is.finite(gene_mat)
  keep_rows <- rowMeans(fin_mask) >= min_prop_keep
  gene_mat  <- gene_mat[keep_rows, , drop = FALSE]
  if (nrow(gene_mat) == 0) stop("After filtering (prop >=", min_prop_keep, ") no genes remain.")
  
  clip_lo <- -2; clip_hi <- 2
  gene_mat <- pmax(pmin(gene_mat, clip_hi), clip_lo)
  
  cond_fac <- droplevels(filt_meta$condition)
  
  cond_levels_present <- levels(cond_fac); if (is.null(cond_levels_present)) cond_levels_present <- unique(cond_fac)
  cond_levels_present <- cond_levels_present[cond_levels_present %in% condition_order]
  cond_levels_present <- factor(cond_levels_present, levels = condition_order)
  cond_means <- lapply(levels(cond_levels_present), function(cc){
    if (is.na(cc)) return(NULL)
    cols <- filt_meta$col_id[cond_fac == cc]
    if (!length(cols)) return(NULL)
    rowMeans(gene_mat[, cols, drop=FALSE], na.rm=TRUE)
  })
  names(cond_means) <- levels(cond_levels_present)
  cond_means <- cond_means[!vapply(cond_means, is.null, logical(1))]
  cond_mean_mat <- if (length(cond_means)) do.call(cbind, cond_means) else NULL
  
  per_ct_delta <- function(mat, meta, target){
    genes <- rownames(mat)
    cts   <- unique(as.character(meta$cell_type))
    deltas  <- matrix(NA_real_, nrow = length(genes), ncol = length(cts),
                      dimnames = list(genes, cts))
    weights <- deltas
    for (ct in cts) {
      cols_ct <- meta$col_id[meta$cell_type == ct]
      tgt_ct  <- meta$col_id[meta$cell_type == ct & meta$condition == target]
      oth_ct  <- setdiff(cols_ct, tgt_ct)
      if (length(tgt_ct) && length(oth_ct)) {
        deltas[, ct]  <- matrixStats::rowMedians(mat[, tgt_ct, drop=FALSE], na.rm=TRUE) -
          matrixStats::rowMedians(mat[, oth_ct, drop=FALSE], na.rm=TRUE)
        weights[, ct] <- length(tgt_ct) + length(oth_ct)
      }
    }
    w <- ifelse(is.finite(weights), weights, 0)
    d <- ifelse(is.finite(deltas),  deltas,  0)
    wsum <- rowSums(w)
    meta_delta <- ifelse(wsum > 0, rowSums(d * w) / wsum, NA_real_)
    sign_same  <- sweep(sign(deltas), 1, sign(meta_delta), FUN = "==")
    avail      <- is.finite(deltas)
    consistency <- rowSums(sign_same & avail, na.rm=TRUE) / pmax(1, rowSums(avail, na.rm=TRUE))
    list(meta_delta = meta_delta, consistency = consistency)
  }
  
  ct_pos_frac <- function(mat, meta, condition){
    genes <- rownames(mat)
    cts   <- unique(as.character(meta$cell_type))
    pos <- numeric(length(genes)); tot <- numeric(length(genes))
    names(pos) <- names(tot) <- genes
    for (ct in cts) {
      cols <- meta$col_id[meta$condition == condition & meta$cell_type == ct]
      if (length(cols)) {
        med_ct <- matrixStats::rowMedians(mat[, cols, drop=FALSE], na.rm=TRUE)
        pos <- pos + as.numeric(med_ct > 0)
        tot <- tot + 1
      }
    }
    ifelse(tot > 0, pos / tot, NA_real_)
  }
  
  pe <- per_ct_delta(gene_mat, filt_meta, target_condition)
  meta_delta  <- pe$meta_delta
  consistency <- pe$consistency
  score       <- abs(meta_delta) * consistency
  
  ord <- order(score, decreasing = TRUE, na.last = NA)
  if (!is.null(top_n)) ord <- ord[seq_len(min(top_n, length(ord)))]
  row_order <- names(meta_delta)[ord]
  
  frac_target_pos <- ct_pos_frac(gene_mat, filt_meta, target_condition)
  row_split <- NULL
  if (row_split_by_sign) {
    up_mask <- (meta_delta[row_order] > 0) & (frac_target_pos[row_order] > 0.5)
    up_mask <- ifelse(is.na(up_mask), FALSE, up_mask)
    row_split <- factor(ifelse(up_mask, "Up in target","Down in target"),
                        levels = c("Up in target","Down in target"))
  }
  
  if (use_olo && have_seriation && !is.null(cond_mean_mat)) {
    smooth_block <- function(idx){
      if (length(idx) <= 2) return(idx)
      sub <- cond_mean_mat[row_order[idx], , drop=FALSE]
      d   <- stats::dist(sub)
      oo  <- seriation::seriate(d, method = "OLO")
      idx[seriation::get_order(oo)]
    }
    if (!is.null(row_split)) {
      grp <- split(seq_along(row_order), row_split)
      row_order <- unlist(lapply(grp, smooth_block), use.names = FALSE) |> (\(i) row_order[i])()
    }
  }
  
  if (!is.null(row_split)) names(row_split) <- row_order
  new_order <- reorder_within_meta_groups(
    mat       = gene_mat,
    meta      = filt_meta,
    row_order = row_order,
    row_split = if (is.null(row_split)) factor(rep("Up in target", length(row_order)),
                                               levels = c("Up in target","Down in target")) else row_split,
    target_condition = target_condition,
    ms_condition     = "Multiple Sclerosis",
    idd_condition    = "Ab-mediated IDD"
  )
  row_split <- if (!is.null(row_split)) row_split[new_order] else NULL
  row_split <- if (!is.null(row_split)) droplevels(row_split) else NULL   # drop empty slices
  row_order <- new_order
  gene_mat  <- gene_mat[row_order, , drop = FALSE]
  
  ## >>> NEW: optionally hide the row-split strip labels while keeping the split/order
  if (!is.null(row_split) && !show_row_split_labels) {
    # keep structure but blank the displayed labels
    row_split <- factor(row_split,
                        levels = levels(row_split),
                        labels = rep("", length(levels(row_split))))
  }
  
  dataset_pal <- c(Esaulova="#E9A5F1", Kang="#FFB4B4", Ramesh="#e76f51",
                   Schafflick="#000000", Mevorach="#ff4d6d")
  subtype_pal <- c(CIS="#bb9457", RRMS="#432818", "CD59 Def."="#fb8500",
                   "MOGAD"="#2C6E49", "NMOSD" ="#B2DF8A")
  condition_pal <- c("Multiple Sclerosis"="#023047",
                     "Ab-mediated IDD"="#8ecae6",
                     "CD59 Def."="#fb8500",
                     "CD55 Def."="#8b8c89")
  sex_pal  <- c(Female="#ff8daa", Male="#32b7ff")
  edss_pal <- circlize::colorRamp2(c(0,6), c("#70e000","#e01e37"))
  age_pal  <- circlize::colorRamp2(c(10,80), c("#48CAE4","#03045e"))
  
  if (is.null(cell_pal)) {
    ct_levels <- levels(filt_meta$cell_type)
    if (is.null(ct_levels)) ct_levels <- unique(as.character(filt_meta$cell_type))
    default_ct_cols <- grDevices::rainbow(length(ct_levels))
    names(default_ct_cols) <- ct_levels
    cell_pal <- default_ct_cols
  }
  
  anno <- ComplexHeatmap::HeatmapAnnotation(
    Condition       = filt_meta$condition,
    `Sub-Condition` = filt_meta$condition_subtype,
    Dataset         = filt_meta$dataset,
    `Cell type`     = filt_meta$cell_type,
    Sex             = filt_meta$sex,
    Age             = filt_meta$age,
    EDSS            = filt_meta$EDSS,
    col = list(
      Condition       = condition_pal,
      `Sub-Condition` = subtype_pal,
      Dataset         = dataset_pal,
      `Cell type`     = cell_pal,
      Sex             = sex_pal,
      Age             = age_pal,
      EDSS            = edss_pal
    ),
    na_col = "grey85",
    annotation_legend_param = list(
      EDSS = list(
        title  = "EDSS",
        at     = seq(0,6,2),
        labels = as.character(seq(0,6,2)),
        legend_gp = grid::gpar(fill = edss_pal(seq(0,6,2)))
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp   = grid::gpar(fontsize = meta_name_fontsize)
  )
  
  brks   <- c(-2,-1,0,1,2)
  labels <- c("≤ -2","-1","0","1","≥ 2")
  col_fun <- circlize::colorRamp2(brks, c("#2166AC","#67A9CF","#FFFFFF","#EF8A62","#B2182B"))
  
  col_split <- cond_fac
  
  ## compute a proper row_gap (0 if only one non-empty slice)
  row_gap_unit <- grid::unit(0, "mm")
  if (!is.null(row_split)) {
    nz <- sum(table(row_split) > 0)
    if (nz >= 2) row_gap_unit <- grid::unit(row_gap_mm, "mm")
  }
  
  ComplexHeatmap::Heatmap(
    gene_mat,
    name              = "Median log2FC",
    top_annotation    = anno,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = FALSE,
    row_split         = row_split,
    row_gap           = row_gap_unit,
    cluster_columns   = FALSE,
    column_order      = common_cols,
    column_split      = col_split,
    column_gap        = grid::unit(condition_gap_mm, "mm"),
    col               = col_fun,
    na_col            = "grey85",
    heatmap_legend_param = list(
      at = brks, labels = labels,
      title = "Median log2FC", title_position = "leftcenter-rot"
    ),
    column_title    = title,
    column_title_gp = grid::gpar(fontface="bold", fontsize=14)
  )
}



get_mulitpule_gene_level_heatmap <- function(
    df,
    title,
    meta_name_fontsize = 11,
    cell_pal = c(
      "NKs"          = "#6A4C93",
      "Monocytes"    = "#5E60CE",
      "B Cells"      = "#C1121F",
      "CD4+ T Cells" = "#1B998B",
      "CD8+ T Cells" = "#D4A373"
    ),
    condition_order   = c("CD59 Def.", "Multiple Sclerosis", "Ab-mediated IDD"),
    cell_type_order   = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
    min_prop_keep     = 0.85,
    condition_gap_mm  = 4,
    row_gap_mm        = 3,     # gap between row-slices
    target_condition  = "CD59 Def.",
    row_split_by_sign = TRUE,  # keep the Up/Down internal ordering logic
    top_n             = NULL,  # applied per gene set (if not NULL)
    use_olo           = TRUE,
    show_row_split_labels = FALSE,  # if TRUE and row_split_by_sign, append "Up/Down" in slice titles
    gene_set_column_name,            # column in df that contains gene-set labels
    order_of_gene_sets = NULL,       # display order (top→bottom); defaults to appearance order
    rename_gene_sets   = NULL        # named chr vector c(orig="New Label") for display
){
  have_seriation <- requireNamespace("seriation", quietly = TRUE)
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("Please install matrixStats.")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Please install ComplexHeatmap.")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Please install circlize.")
  
  # ── column sanity ───────────────────────────────────────────────────────────
  if (!("dataset" %in% names(df)) && ("data_set" %in% names(df))) df$dataset <- df$data_set
  req_base <- c("sample","gene","median_log2fc","dataset",
                "condition","condition_subtype","age","sex","EDSS","cell_type")
  miss <- setdiff(req_base, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))
  if (!is.character(gene_set_column_name) || length(gene_set_column_name) != 1L)
    stop("`gene_set_column_name` must be a single column name (character).")
  if (!(gene_set_column_name %in% names(df)))
    stop("`gene_set_column_name` not found in `df`.")
  
  if (!is.numeric(df$median_log2fc)) stop("'median_log2fc' must be numeric")
  if (!is.numeric(df$age))  stop("'age' must be numeric")
  if (!is.numeric(df$EDSS)) stop("'EDSS' must be numeric")
  if (!is.character(df$gene) && !is.factor(df$gene)) stop("'gene' must be chr/fct")
  
  # flag for single cell type
  single_ct <- (length(unique(df$cell_type)) == 1L)
  
  # ── normalize & ids ────────────────────────────────────────────────────────
  df <- df |>
    dplyr::mutate(
      sample    = as.character(sample),
      gene      = as.character(gene),
      cell_type = as.character(cell_type),
      col_id    = paste0(sample," | ",cell_type),
      .gene_set = as.character(.data[[gene_set_column_name]]),
      .row_id   = paste0(.gene_set, " ⟷ ", gene)  # make rows unique per (gene_set, gene)
    )
  
  # metadata table for columns
  filt_meta <- df |>
    dplyr::arrange(sample, cell_type) |>
    dplyr::distinct(col_id, .keep_all = TRUE)
  
  # wide matrix (rows = (gene_set, gene), cols = sample|cell_type)
  mat_wide <- df |>
    dplyr::select(.row_id, col_id, median_log2fc) |>
    dplyr::group_by(.row_id, col_id) |>
    dplyr::summarise(value = mean(median_log2fc, na.rm=TRUE), .groups="drop") |>
    tidyr::pivot_wider(names_from = col_id, values_from = value, values_fill = NA_real_) |>
    dplyr::arrange(.row_id)
  
  gene_mat <- as.matrix(mat_wide[,-1, drop=FALSE])
  rownames(gene_mat) <- mat_wide$.row_id
  
  # align column order & factors
  filt_meta$condition <- factor(filt_meta$condition, levels = condition_order)
  filt_meta$cell_type <- factor(filt_meta$cell_type, levels = cell_type_order)
  filt_meta <- filt_meta |> dplyr::arrange(condition, dataset, sample, cell_type)
  common_cols <- filt_meta$col_id
  gene_mat    <- gene_mat[, common_cols, drop=FALSE]
  filt_meta   <- filt_meta[match(common_cols, filt_meta$col_id), , drop=FALSE]
  
  # keep sufficiently observed rows
  fin_mask  <- is.finite(gene_mat)
  keep_rows <- rowMeans(fin_mask) >= min_prop_keep
  gene_mat  <- gene_mat[keep_rows, , drop = FALSE]
  if (nrow(gene_mat) == 0) stop("After filtering (prop >=", min_prop_keep, ") no (gene_set, gene) rows remain.")
  
  # clipping (same palette range)
  clip_lo <- -2; clip_hi <- 2
  gene_mat <- pmax(pmin(gene_mat, clip_hi), clip_lo)
  
  # recover per-row original gene_set & gene
  row_to_set <- df |> dplyr::select(.row_id, .gene_set) |> dplyr::distinct()
  row_to_set <- row_to_set[row_to_set$.row_id %in% rownames(gene_mat), , drop=FALSE]
  set_for_row <- setNames(row_to_set$.gene_set, row_to_set$.row_id)
  
  # column condition split & means
  cond_fac <- droplevels(filt_meta$condition)
  cond_levels_present <- levels(cond_fac); if (is.null(cond_levels_present)) cond_levels_present <- unique(cond_fac)
  cond_levels_present <- cond_levels_present[cond_levels_present %in% condition_order]
  cond_levels_present <- factor(cond_levels_present, levels = condition_order)
  cond_means <- lapply(levels(cond_levels_present), function(cc){
    if (is.na(cc)) return(NULL)
    cols <- filt_meta$col_id[cond_fac == cc]
    if (!length(cols)) return(NULL)
    rowMeans(gene_mat[, cols, drop=FALSE], na.rm=TRUE)
  })
  names(cond_means) <- levels(cond_levels_present)
  cond_means <- cond_means[!vapply(cond_means, is.null, logical(1))]
  cond_mean_mat <- if (length(cond_means)) do.call(cbind, cond_means) else NULL
  
  # helpers (match single-set behavior)
  per_ct_delta <- function(mat, meta, target){
    rows <- rownames(mat)
    cts   <- unique(as.character(meta$cell_type))
    deltas  <- matrix(NA_real_, nrow = length(rows), ncol = length(cts),
                      dimnames = list(rows, cts))
    weights <- deltas
    for (ct in cts) {
      cols_ct <- meta$col_id[meta$cell_type == ct]
      tgt_ct  <- meta$col_id[meta$cell_type == ct & meta$condition == target]
      oth_ct  <- setdiff(cols_ct, tgt_ct)
      if (length(tgt_ct) && length(oth_ct)) {
        deltas[, ct]  <- matrixStats::rowMedians(mat[, tgt_ct, drop=FALSE], na.rm=TRUE) -
          matrixStats::rowMedians(mat[, oth_ct, drop=FALSE], na.rm=TRUE)
        weights[, ct] <- length(tgt_ct) + length(oth_ct)
      }
    }
    w <- ifelse(is.finite(weights), weights, 0)
    d <- ifelse(is.finite(deltas),  deltas,  0)
    wsum <- rowSums(w)
    meta_delta <- ifelse(wsum > 0, rowSums(d * w) / wsum, NA_real_)
    sign_same  <- sweep(sign(deltas), 1, sign(meta_delta), FUN = "==")
    avail      <- is.finite(deltas)
    consistency <- rowSums(sign_same & avail, na.rm=TRUE) / pmax(1, rowSums(avail, na.rm=TRUE))
    list(meta_delta = meta_delta, consistency = consistency)
  }
  ct_pos_frac <- function(mat, meta, condition){
    rows <- rownames(mat)
    cts   <- unique(as.character(meta$cell_type))
    pos <- numeric(length(rows)); tot <- numeric(length(rows))
    names(pos) <- names(tot) <- rows
    for (ct in cts) {
      cols <- meta$col_id[meta$condition == condition & meta$cell_type == ct]
      if (length(cols)) {
        med_ct <- matrixStats::rowMedians(mat[, cols, drop=FALSE], na.rm=TRUE)
        pos <- pos + as.numeric(med_ct > 0)
        tot <- tot + 1
      }
    }
    ifelse(tot > 0, pos / tot, NA_real_)
  }
  
  # global meta-effect stats (per (gene_set,gene) row)
  pe <- per_ct_delta(gene_mat, filt_meta, target_condition)
  meta_delta_all  <- pe$meta_delta
  consistency_all <- pe$consistency
  score_all       <- abs(meta_delta_all) * consistency_all
  frac_target_pos_all <- ct_pos_frac(gene_mat, filt_meta, target_condition)
  
  # ── gene-set ordering & relabels ───────────────────────────────────────────
  sets_present_all <- unique(set_for_row[rownames(gene_mat)])
  if (is.null(order_of_gene_sets)) {
    order_of_gene_sets <- unique(df[[gene_set_column_name]])
  }
  order_of_gene_sets <- intersect(order_of_gene_sets, sets_present_all)
  if (!length(order_of_gene_sets))
    stop("None of the specified `order_of_gene_sets` are present after filtering.")
  
  label_map <- function(x){
    if (!is.null(rename_gene_sets) && length(rename_gene_sets) && (x %in% names(rename_gene_sets))) {
      rename_gene_sets[[x]]
    } else x
  }
  display_levels <- vapply(order_of_gene_sets, label_map, FUN.VALUE = character(1))
  
  # ── build per-set row orders (and optional Up/Down slice) ──────────────────
  final_row_order <- character(0)
  final_set_label_for_row <- character(0)  # displayed label (after rename)
  final_sign_for_row      <- character(0)  # "Up in target"/"Down in target"
  
  # small helper for OLO on index vector
  .olo <- function(idx) {
    if (!use_olo || !have_seriation || is.null(cond_mean_mat) || length(idx) <= 2) return(idx)
    sub <- cond_mean_mat[idx, , drop=FALSE]
    d   <- stats::dist(sub)
    oo  <- seriation::seriate(d, method = "OLO")
    seriation::get_order(oo)
  }
  
  for (gs in order_of_gene_sets) {
    rows_in_set <- names(set_for_row)[set_for_row == gs]
    rows_in_set <- intersect(rows_in_set, rownames(gene_mat))
    if (!length(rows_in_set)) next
    
    # score-based order (desc), then optional top_n per set
    ord <- rows_in_set[order(score_all[rows_in_set], decreasing = TRUE, na.last = NA)]
    if (!is.null(top_n)) ord <- utils::head(ord, min(top_n, length(ord)))
    if (!length(ord)) next
    
    # within-set Up/Down classification (same logic)
    up_mask <- (meta_delta_all[ord] > 0) & (frac_target_pos_all[ord] > 0.5)
    up_mask <- ifelse(is.na(up_mask), FALSE, up_mask)
    row_split_sign <- factor(ifelse(up_mask, "Up in target", "Down in target"),
                             levels = c("Up in target","Down in target"))
    names(row_split_sign) <- ord
    
    # OLO smoothing within sign slices (robust to NAs)
    if (use_olo && have_seriation && !is.null(cond_mean_mat)) {
      grp <- split(ord, row_split_sign)
      ord <- unlist(lapply(grp, function(xx){
        if (!length(xx)) return(xx)
        rel <- match(xx, rownames(cond_mean_mat))
        good <- !is.na(rel)
        if (sum(good) <= 2) return(xx)
        ord_idx <- .olo(rel[good])
        c(xx[good][ord_idx], xx[!good])
      }), use.names = FALSE)
    }
    
    # meta reordering rule within sign groups
    ord <- reorder_within_meta_groups(
      mat       = gene_mat,
      meta      = filt_meta,
      row_order = ord,
      row_split = row_split_sign[ord],
      target_condition = target_condition,
      ms_condition     = "Multiple Sclerosis",
      idd_condition    = "Ab-mediated IDD"
    )
    # accumulate
    final_row_order <- c(final_row_order, ord)
    final_set_label_for_row <- c(final_set_label_for_row, setNames(rep(label_map(gs), length(ord)), ord))
    final_sign_for_row      <- c(final_sign_for_row, as.character(row_split_sign[ord]))
  }
  
  # if something left unordered (unlikely), append at bottom
  missing_rows <- setdiff(rownames(gene_mat), final_row_order)
  if (length(missing_rows)) {
    final_row_order <- c(final_row_order, missing_rows)
    final_set_label_for_row <- c(final_set_label_for_row,
                                 setNames(label_map(set_for_row[missing_rows]), missing_rows))
    leftover_sign <- ifelse(meta_delta_all[missing_rows] > 0 & frac_target_pos_all[missing_rows] > 0.5,
                            "Up in target","Down in target")
    final_sign_for_row <- c(final_sign_for_row, setNames(leftover_sign, missing_rows))
  }
  
  gene_mat <- gene_mat[final_row_order, , drop=FALSE]
  
  # ── row split structure ────────────────────────────────────────────────────
  # Outer = gene set (renamed), Inner = (optional) Up/Down
  outer_fac <- factor(final_set_label_for_row[final_row_order], levels = display_levels)
  outer_fac <- droplevels(outer_fac)  # drop empty sets
  if (row_split_by_sign) {
    inner_fac <- factor(final_sign_for_row[final_row_order], levels = c("Up in target","Down in target"))
    inner_fac <- droplevels(inner_fac)
    row_split <- data.frame(`Gene set` = outer_fac, Direction = inner_fac, check.names = FALSE)
  } else {
    row_split <- outer_fac
  }
  
  
  # ── row titles (one per slice), hiding any ", NA" suffix ───────────────────
  if (is.data.frame(row_split)) {
    # include NA as an explicit (but hidden) level so #titles == #slices
    inner_levels <- if (row_split_by_sign) addNA(row_split[[2]]) else row_split[[2]]
    key <- interaction(row_split[[1]], inner_levels, drop = TRUE, lex.order = TRUE, sep = "\t")
    lev <- levels(key)
    
    row_titles <- vapply(lev, function(k){
      parts <- strsplit(k, "\t", fixed = TRUE)[[1]]
      gs_lab <- parts[1]
      # Only append Direction when requested AND not NA
      if (isTRUE(show_row_split_labels) && length(parts) >= 2L && !identical(parts[2], "NA")) {
        paste0(gs_lab, " — ", parts[2])
      } else {
        gs_lab
      }
    }, FUN.VALUE = character(1))
    
    n_slices <- length(lev)
  } else {
    row_titles <- levels(row_split)
    n_slices <- length(row_titles)
  }
  row_gap <- if (n_slices >= 2) grid::unit(row_gap_mm, "mm") else grid::unit(0, "mm")
  
  # ── annotations & palettes (identical to single-set function) ──────────────
  dataset_pal <- c(Esaulova="#E9A5F1", Kang="#FFB4B4", Ramesh="#e76f51",
                   Schafflick="#000000", Mevorach="#ff4d6d")
  subtype_pal <- c(CIS="#bb9457", RRMS="#432818", "CD59 Def."="#fb8500",
                   "MOGAD"="#2C6E49", "NMOSD" ="#B2DF8A")
  condition_pal <- c("Multiple Sclerosis"="#023047",
                     "Ab-mediated IDD"="#8ecae6",
                     "CD59 Def."="#fb8500",
                     "CD55 Def."="#8b8c89")
  sex_pal  <- c(Female="#ff8daa", Male="#32b7ff")
  edss_pal <- circlize::colorRamp2(c(0,6), c("#70e000","#e01e37"))
  age_pal  <- circlize::colorRamp2(c(10,80), c("#48CAE4","#03045e"))
  
  if (is.null(cell_pal)) {
    ct_levels <- levels(filt_meta$cell_type)
    if (is.null(ct_levels)) ct_levels <- unique(as.character(filt_meta$cell_type))
    default_ct_cols <- grDevices::rainbow(length(ct_levels))
    names(default_ct_cols) <- ct_levels
    cell_pal <- default_ct_cols
  }
  
  # build annotation arguments, optionally omitting Cell type
  anno_args <- list(
    Condition       = filt_meta$condition,
    `Sub-Condition` = filt_meta$condition_subtype,
    Dataset         = filt_meta$dataset,
    Sex             = filt_meta$sex,
    Age             = filt_meta$age,
    EDSS            = filt_meta$EDSS,
    col = list(
      Condition       = condition_pal,
      `Sub-Condition` = subtype_pal,
      Dataset         = dataset_pal,
      Sex             = sex_pal,
      Age             = age_pal,
      EDSS            = edss_pal
    ),
    na_col = "grey85",
    annotation_legend_param = list(
      EDSS = list(
        title  = "EDSS",
        at     = seq(0,6,2),
        labels = as.character(seq(0,6,2)),
        legend_gp = grid::gpar(fill = edss_pal(seq(0,6,2)))
      )
    ),
    annotation_name_side = "left",
    annotation_name_gp   = grid::gpar(fontsize = meta_name_fontsize)
  )
  
  if (!single_ct) {
    anno_args$`Cell type` <- filt_meta$cell_type
    anno_args$col$`Cell type` <- cell_pal
  }
  
  anno <- do.call(ComplexHeatmap::HeatmapAnnotation, anno_args)
  
  brks   <- c(-2,-1,0,1,2)
  labels <- c("≤ -2","-1","0","1","≥ 2")
  col_fun <- circlize::colorRamp2(brks, c("#2166AC","#67A9CF","#FFFFFF","#EF8A62","#B2182B"))
  
  col_split <- cond_fac
  
  ComplexHeatmap::Heatmap(
    gene_mat,
    name              = "Median log2FC",
    top_annotation    = anno,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    cluster_rows      = FALSE,
    row_split         = row_split,
    row_gap           = row_gap,                 # single unit; ComplexHeatmap recycles
    cluster_columns   = FALSE,
    column_order      = common_cols,
    column_split      = col_split,
    column_gap        = grid::unit(condition_gap_mm, "mm"),
    col               = col_fun,
    na_col            = "grey85",
    heatmap_legend_param = list(
      at = brks, labels = labels,
      title = "Median log2FC", title_position = "leftcenter-rot"
    ),
    column_title    = title,
    column_title_gp = grid::gpar(fontface="bold", fontsize=14),
    row_title       = row_titles,
    row_title_gp    = grid::gpar(fontface = "bold")
  )
}



#───────────────────────────────────────────────────────────────────────────────
#  AUC boxplots for whole gene sets of the phenotype
#───────────────────────────────────────────────────────────────────────────────
 


get_one_boxplot_whole_gene_sets <- function(
    df,
    cell_col        = "cell_type",
    sig_col         = "sig_name",
    condition_col   = "condition",
    value_col       = "sig_score",
    
    # optional single-signature mode
    sig_name        = NULL,
    
    # explicit ordering controls (optional)
    row_order       = NULL,   # order of signatures (rows)
    col_order       = NULL,   # order of cell types (columns)
    row_decreasing  = FALSE,  # reverse row order
    col_decreasing  = FALSE,  # reverse column order
    
    # rename rows (facet strip labels); named character vector old -> new
    row_renames     = c(
      "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS" =
        "Extrinsic Apoptosis",
      "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS" =
        "Apoptosis Inhibitors",
      "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammatory Response",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFα / NF-κB",
      "WP_NRF2_PATHWAY" = "NRF2 Pathway",
      "mitochondria"    = "ETC CII-CIII"
    ),
    
    condition_order = c("CD59.Def","MS","Ab.mediated.IDD"),
    palette         = c("MS" = "#012A3A",
                        "Ab.mediated.IDD" = "#2C7FB8",
                        "CD59.Def" = "#E66100"),
    
    # legend condition label renaming (old level -> new label)
    condition_renames =  c(
      "CD59.Def"        = "CD59 Def.",
      "MS"              = "Multiple Sclerosis",
      "Ab.mediated.IDD" = "Ab-mediated IDD"
    ),
    
    show_points   = FALSE,
    
    # ---- Plot title & legend title controls ----------------------------------
    title               = "Signature AUC Scores from Shared DEGs Across Conditions", 
    legend_title        = "Condition",
    title_size          = 20,
    cell_type_size      = 17,
    sig_size            = 17,
    y_text_size         = 17,
    legend_title_size   = 17,
    legend_text_size    = 17,
    title_face          = "bold",
    title_hjust         = 0.5,
    
    row_gap       = 2.5,        # gap between rows
    row_gap_unit  = "lines",    # gap unit
    
    # horizontal guide at y = 0.5 in every panel
    add_row_y_line      = TRUE,
    row_y_line          = 0.5,
    row_y_line_color    = "grey60",
    row_y_line_size     = 0.9,
    row_y_line_linetype = "dashed"
) {
  stopifnot(is.data.frame(df))
  req_cols <- c(cell_col, sig_col, condition_col, value_col)
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))
  if (!is.numeric(df[[value_col]]))
    stop(sprintf("`%s` must be numeric.", value_col))
  
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(forcats)
  })
  
  # ----- Optional: focus on a single signature --------------------------------
  if (!is.null(sig_name)) {
    df <- df %>% dplyr::filter(.data[[sig_col]] == sig_name)
    if (nrow(df) == 0) {
      stop("No rows found for sig_name = ", sig_name)
    }
  }
  
  # ----- Resolve desired orders ----------------------------------------------
  desired_row <- row_order
  desired_col <- col_order
  
  if (isTRUE(row_decreasing) && !is.null(desired_row)) desired_row <- rev(desired_row)
  if (isTRUE(col_decreasing) && !is.null(desired_col)) desired_col <- rev(desired_col)
  
  if (!is.null(desired_row)) {
    not_in_rows <- setdiff(desired_row, unique(df[[sig_col]]))
    if (length(not_in_rows))
      warning("These 'row_order' entries were not found in ", sig_col, ": ",
              paste(not_in_rows, collapse = ", "))
  }
  if (!is.null(desired_col)) {
    not_in_cols <- setdiff(desired_col, unique(df[[cell_col]]))
    if (length(not_in_cols))
      warning("These 'col_order' entries were not found in ", cell_col, ": ",
              paste(not_in_cols, collapse = ", "))
  }
  
  # ----- Factorize with requested orders -------------------------------------
  df <- df %>%
    mutate(
      !!sig_col       := if (is.null(desired_row)) .data[[sig_col]] else
        fct_relevel(.data[[sig_col]], desired_row),
      !!cell_col      := if (is.null(desired_col)) .data[[cell_col]] else
        fct_relevel(.data[[cell_col]], desired_col),
      !!condition_col := if (is.null(condition_order)) .data[[condition_col]] else
        fct_relevel(.data[[condition_col]], condition_order)
    )
  
  # ----- Optional row renaming (facet strip labels) ---------------------------
  # Skip in single-signature mode (no row facets)
  if (!is.null(row_renames) && is.null(sig_name)) {
    if (!is.character(row_renames) || is.null(names(row_renames)) || any(names(row_renames) == "")) {
      stop("'row_renames' must be a NAMED character vector like c('old'='New Label', ...).")
    }
    sig_levels <- levels(factor(df[[sig_col]]))
    if (!setequal(sig_levels, names(row_renames))) {
      stop("`row_renames` must have names for ALL current signature names. ",
           "Expected names: ", paste(sig_levels, collapse = ", "))
    }
    renames_in_order <- unname(row_renames[sig_levels])
    df[[sig_col]] <- forcats::fct_relabel(factor(df[[sig_col]], levels = sig_levels),
                                          function(x) renames_in_order[match(x, sig_levels)])
  }
  
  # ----- Build plot -----------------------------------------------------------
  p <- ggplot(
    df,
    aes(
      x      = .data[[condition_col]],
      y      = .data[[value_col]],
      colour = .data[[condition_col]]
    )
  )
  
  if (isTRUE(add_row_y_line)) {
    p <- p +
      geom_hline(
        yintercept = row_y_line,
        colour     = row_y_line_color,
        linewidth  = row_y_line_size,
        linetype   = row_y_line_linetype,
        inherit.aes = FALSE
      )
  }
  
  p <- p +
    geom_boxplot(outlier.shape = NA, width = 0.6, fill = NA, size = 0.8)
  
  if (show_points) {
    p <- p + geom_point(
      position = position_jitter(width = 0.07, height = 0, seed = 1),
      size = 0.6, alpha = 0.6, color = "black"
    )
  }
  
  # Determine legend breaks in use
  legend_breaks <- if (!is.null(condition_order)) condition_order else levels(df[[condition_col]])
  
  # Build legend labels (apply renames if provided)
  legend_labels <- legend_breaks
  if (!is.null(condition_renames)) {
    if (!is.character(condition_renames) || is.null(names(condition_renames)) ||
        any(names(condition_renames) == "")) {
      stop("`condition_renames` must be a NAMED character vector like c('old'='New Label', ...).")
    }
    missing_names <- setdiff(legend_breaks, names(condition_renames))
    if (length(missing_names)) {
      stop("`condition_renames` must provide labels for ALL legend entries. Missing: ",
           paste(missing_names, collapse = ", "))
    }
    legend_labels <- unname(condition_renames[legend_breaks])
  }
  
  # ----- Faceting: only in multi-signature mode -------------------------------
  if (is.null(sig_name)) {
    p <- p +
      facet_grid(
        rows   = vars(.data[[sig_col]]),
        cols   = vars(.data[[cell_col]]),
        scales = "fixed",
        switch = "y",
        drop   = FALSE,
        labeller = labeller(.rows = label_wrap_gen(width = 25))
      )
  }
  
  # color scale (always used)
  p <- p +
    scale_color_manual(values = palette,
                       breaks = legend_breaks,
                       labels = legend_labels,
                       drop   = FALSE,
                       name   = legend_title)
  
  # ---- Y scale: fixed 0–1, 0.2 steps, ensure row_y_line included ------------
  p <- p +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = sort(unique(c(seq(0, 1, by = 0.2), row_y_line))),
      expand = expansion(mult = c(0.02, 0.05))
    )
  
  # ---- Title: in single-signature mode, override to "AUC scores for <label>" -
  if (!is.null(sig_name)) {
    sig_raw <- unique(df[[sig_col]])[1]
    sig_label <- sig_raw
    if (!is.null(row_renames) && sig_raw %in% names(row_renames)) {
      sig_label <- row_renames[[sig_raw]]
    }
    title <- paste("AUC scores for", sig_label)
  }
  
  p <- p +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_line(color = "#efefef"),
      panel.grid.minor   = element_blank(),
      strip.placement    = "outside",
      strip.background   = element_blank(),
      strip.text.y.left  = element_text(size = sig_size, face = "bold", angle = 0, hjust = 1),
      strip.text.x       = element_text(size = cell_type_size, face = "bold"),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = y_text_size),
      panel.spacing.y    = grid::unit(row_gap, row_gap_unit),
      panel.spacing.x    = grid::unit(0.8, "lines"),
      legend.position    = "right",
      legend.box         = "vertical",
      legend.title       = element_text(size = legend_title_size, face = "bold"),
      legend.text        = element_text(size = legend_text_size),
      plot.title         = element_text(size = title_size, face = title_face, hjust = title_hjust)
    )
  
  return(p)
}

get_all_boxplot_whole_gene_sets <- function(
    df,
    sig_order = NULL,   # order of signatures (controls subplot order)
    ncol      = 3,      # number of columns in the grid
    ...
) {
  stopifnot(is.data.frame(df))
  
  
  
  sig_col  <- "sig_name"
  all_sigs <- unique(df[[sig_col]])
  
  # If user didn't specify an order, use all signatures as they appear
  if (is.null(sig_order)) {
    sig_order <- all_sigs
  } else {
    missing <- setdiff(sig_order, all_sigs)
    if (length(missing)) {
      stop("These sig_name values are not present in df: ",
           paste(missing, collapse = ", "))
    }
  }
  
  # Build list of individual plots (one per signature)
  plots <- lapply(sig_order, function(s) {
    get_one_boxplot_whole_gene_sets(
      df       = df,
      sig_name = s,
      ...
    )
  })
  
  # Combine into a grid, collect one shared legend, add A/B/C... tags
  combined <- wrap_plots(
    plotlist = plots,
    ncol     = ncol,
    guides   = "collect"
  ) +
    plot_annotation(tag_levels = "A") &
    theme(
      legend.position = "right",
      plot.tag        = element_text(face = "bold", size = 18)
    )
  
  return(combined)
}










