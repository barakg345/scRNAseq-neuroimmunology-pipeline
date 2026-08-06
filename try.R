
seu_hao_3 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/cite-seq_benchmark/seu_hao_3p.qs")
seu_hao_5 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/cite-seq_benchmark/seu_hao_5p.qs")

metric_column <- "is_CD8_T_cell_try"
ground_truth_label <- "CD8 T"
metric_label <- "CD8+ T cells"
ucell_column <-  "cd_8_signature_by_ucell"

plan(multisession, workers = 8)         # or multicore on Linux console
options(future.globals.maxSize = 100 * 1024^3) # 100 GiB



#get_monocytes for seu_hao_3 and seu_hao_5
{
  seu_hao_3 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/cite-seq_benchmark/seu_hao_3p.qs")
  seu_hao_3 <- get_CD4_cells(seu_hao_3)
  
  seu_hao_5 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/cite-seq_benchmark/seu_hao_5p.qs")
  seu_hao_5 <- get_CD8_cells(seu_hao_5)
  
  
  
  
  
  
  
  
  seu_hao_5$SingleR.labels.blueprint_try <- factor(
    case_when(
      seu_hao_5$SingleR.labels.blueprint == "CD8+ T cells" &
        seu_hao_5$is_CD4_T_cell == "CD4+ T cells" ~ "CD4+ T cells",
      TRUE ~ as.character(seu_hao_5$SingleR.labels.blueprint)
    ),
    levels = bp_levels
  )
  
  seu_hao_5$SingleR.labels.monaco_try <- factor(
    case_when(
      seu_hao_5$SingleR.labels.monaco == "CD8+ T cells" &
        seu_hao_5$is_CD4_T_cell == "CD4+ T cells" ~ "CD4+ T cells",
      TRUE ~ as.character(seu_hao_5$SingleR.labels.monaco)
    ),
    levels = mo_levels
  )
  
  
  positive_cd8_signature <- c("CD8A+","CD8B+","CD8B2+", "TRAC+","LCK+","RUNX3+"                              ,"CD3E+","EOMES+"
                              ,"CD3E+","EOMES+"
                              # ,"TRDC+", "IKZF2+"
                              )
  
  
  positive_cd_4_burden_signature <- c("CD4+",
    "CCR7+", "IL7R")
  
  positive_nk_burden_signature <- c("NCR1", "FCGR3A", "TYROBP")
  
  signatures <- list(
    positive_cd8_signature = positive_cd8_signature,
    positive_cd_4_burden_signature = positive_cd_4_burden_signature,
    positive_nk_burden_signature = positive_nk_burden_signature
  )
  
  
  # This appends a column named  'positive_cd_4_signature_UCell'
  seu_hao_5 <- AddModuleScore_UCell(
    seu_hao_5, features = signatures, name = NULL)
  
  # 3. Compute a composite score and threshold
  seu_hao_5$CD8_to_CD4_score_primary     <- seu_hao_5$positive_cd8_signature - seu_hao_5$positive_cd_4_burden_signature
  seu_hao_5$CD8_to_NK_score_primary     <- seu_hao_5$positive_cd8_signature - seu_hao_5$positive_nk_burden_signature

  
  seu_hao_5$cd_8_signature_by_ucell <- ifelse(
    (seu_hao_5$CD8_to_CD4_score_primary >= 0.3) |
      (seu_hao_5$CD8_to_NK_score_primary >= 0.3),
    "CD8+ T cells",
    "other"
  )
  
  seu_hao_5$cd_8_signature_by_ucell <- ifelse((seu_hao_5$cd_8_signature_by_ucell == "CD8+ T cells" & 
                                                 seu_hao_5$is_CD4_T_cell == "CD4+ T cells"), "other",seu_hao_5$cd_8_signature_by_ucell)
  
  table(seu_hao_5$cd_8_signature_by_ucell, seu_hao_5$celltype.l1)
  
  
  
  seu_hao_5$is_CD8_T_cell_try <- dplyr::case_when(
    (seu_hao_5$cd_8_signature_by_ucell == "CD8+ T cells" &
       seu_hao_5$SingleR.labels.blueprint_try == "CD8+ T cells") |
      (seu_hao_5$cd_8_signature_by_ucell == "CD8+ T cells" &
         seu_hao_5$SingleR.labels.monaco_try    == "CD8+ T cells") |
      (seu_hao_5$SingleR.labels.blueprint_try   == "CD8+ T cells" &
         seu_hao_5$SingleR.labels.monaco_try    == "CD8+ T cells") ~ "CD8+ T cells",
    TRUE ~ seu_hao_5$is_CD8_T_cell
  )
    
  
  seu_hao_5$is_CD8_T_cell <- dplyr::case_when(
    (seu_hao_5$cd_8_signature_by_ucell == "CD8+ T cells" &
       seu_hao_5$SingleR.labels.blueprint == "CD8+ T cells") |
      (seu_hao_5$cd_8_signature_by_ucell == "CD8+ T cells" &
         seu_hao_5$SingleR.labels.monaco    == "CD8+ T cells") |
      (seu_hao_5$SingleR.labels.blueprint   == "CD8+ T cells" &
         seu_hao_5$SingleR.labels.monaco    == "CD8+ T cells") ~ "CD8+ T cells",
    TRUE ~ "other"
  )  
  
  
  seu_hao_5$is_CD8_T_cell_try <- case_when(
    seu_hao_5$is_CD8_T_cell == "CD8+ T cells" & seu_hao_5$is_CD4_T_cell == "CD4+ T cells" ~ "other",
    seu_hao_5$is_CD8_T_cell == "CD8+ T cells" & seu_hao_5$is_nk == "NK cells" ~ "other",    
    TRUE ~ seu_hao_5$is_CD8_T_cell_try) 
  
  seu_hao_5$confusion_matrix <- dplyr::case_when(
    seu_hao_5$celltype.l1 == "CD4 T" ~ "CD4+ T cells" ,
    (seu_hao_5$is_CD8_T_cell == "other" &
       seu_hao_5$celltype.l1 == "CD8 T") ~ "FN" ,
    TRUE ~ "other"
  )
   
  Idents(seu_hao_5) <- "confusion_matrix"
  
  fn_to_cd4_markers <- FindMarkers(seu_hao_5, ident.1 = "FN", ident.2 = "CD4+ T cells")
  
  
  
  
  
  
  
  hao_5_confusion_matrix_and_idenfication_umap<- get_confusion_matrix_and_idenfication_4_umaps(seu_hao_5,metric_column, metric_label, ground_truth_column = "celltype.l1", ground_truth_label, ucell_column = ucell_column)
  
  
  votes <- cbind(
    seu_hao_5$cd_8_signature_by_ucell        == "CD8+ T cells",
    seu_hao_5$SingleR.labels.blueprint_try   == "CD8+ T cells",
    seu_hao_5$SingleR.labels.monaco_try      == "CD8+ T cells"
  )
  seu_hao_5$is_CD8_T_cell_try <- ifelse(rowSums(votes) >= 2,
                                        "CD8+ T cells", "other")
  
  
  
}

#get plots
{
  
  #hao_3
  {
    hao_3_confusion_matrix_and_idenfication_umap<- get_confusion_matrix_and_idenfication_4_umaps(seu_hao_3,metric_column, metric_label, ground_truth_column = "celltype.l1", ground_truth_label, ucell_column = ucell_column)
    
    p_seu_hao_3_ground_truth <- dittoDimPlot(seu_hao_3, "celltype.l1", reduction.use = "wnn.umap", size = 0.4) +
      ggplot2::ggtitle(label = "<b>Ground Truth</b> lables for seu_hao_3" ) +
      theme(plot.title = element_markdown())
    
    
  }
  
  #hao_5
  {
    hao_5_confusion_matrix_and_idenfication_umap<- get_confusion_matrix_and_idenfication_4_umaps(seu_hao_5,metric_column, metric_label, ground_truth_column = "celltype.l1", ground_truth_label, ucell_column = ucell_column)
    
    
    p_seu_hao_5_ground_truth <- dittoDimPlot(seu_hao_5, "celltype.l1", reduction.use = "wnn.umap", size = 0.4) +
      ggplot2::ggtitle(label = "<b>Ground Truth</b> lables for seu_hao_5" ) +
      theme(plot.title = element_markdown())
    
    
    
  }
}

#save as pdfs
{
  #hao 3
  {
    pdf("/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp_hao3.pdf", width = 26, height = 20)
    print(hao_3_confusion_matrix_and_idenfication_umap$p_cell_type_idenfication)
    print(p_seu_hao_3_ground_truth )
    dev.off()
    
  }
  
  
  #hao 5
  
  {
    pdf("/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp_hao5.pdf", width = 26, height = 20)
    print(hao_5_confusion_matrix_and_idenfication_umap$p_cell_type_idenfication)
    print(p_seu_hao_5_ground_truth )
    dev.off()
  }
}



