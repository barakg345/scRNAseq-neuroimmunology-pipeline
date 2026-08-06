#Libs and functions
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
    library(circlize)
    library(furrr) 
    library(msigdbdf)
    library(ComplexUpset)
    library(AnnotationDbi)
    library(SBGNview)
    # library(reactome.db)      # offline Reactome mappings :contentReference[oaicite:0]{index=0}
    library(GO.db)
    library(BiocParallel)
    library(forcats)
    library(scMetabolism)
    library(broom)
    # library(ReactomeContentService4R)
    library(Cairo)
    library(MetBrewer)
    
   
  })
  source("/bigdata/bgutman/Collabs/Dror mevorach/scripts/functions.R")
  
  rm(mapped.ids, pathways.info, SBGNhub.id.mapping.tables)
  gc()
  
}

#Color Platte
{
  group_palette  = c(
    "Proliferation"                     = "#E9A5F1",
    "Metabolism / Cellular Respiration" = "#B2DF8A",
    "Stress Adaptation"                 = "#48CAE4",
    "Replication-Stress"                = "#fb8500",
    "Oxidative-Stress"                  = "#FFB4B4",
    "p53-Senescence"                    = "#023e8a",
    "Pro Inflammatory"                  = "#ff4d6d",
    "Pro Inhibitory"                    = "#ffb703",
    "Migration / Chemotaxis"            = "#2C6E49",
    "Other"                             = "grey"
  )
  
  c("#8ecae6", "#219ebc", "#023047", "#ffb703", "#fb8500")
  
  c(
    "NKs"          = "#8ecae6",
    "Monocytes"    = "#219ebc",
    "B Cells"      = "#023047",
    "CD4+ T Cells" = "#ffb703",
    "CD8+ T Cells" = "#fb8500"
  )
  
}


#Basic pipeline
{
  
  #Create healthy list for each cell type
  {
    #Subset to healthy controls only
    {
      #Hong
      seu_hong <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/pSS/Hong/seu_hong.qs")
      
      seu_hong <- subset(seu_hong, subset = sample  %in% c("HC-1", "HC-2",
                                                           "HC-3",
                                                           "HC-4",
                                                           "HC-5"))
      
      #Kang
      seu_kang <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS_and_other_disorders/Kang/seu_kang.qs")
      
      seu_kang <- subset(seu_kang, subset = sample  %in% c("JSB",
                                                           "KYO"))
      
      
      #Ramesh
      seu_ramesh <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS_and_other_disorders/Ramesh/seu_ramesh_5.qs")
      
      seu_ramesh <- subset(seu_ramesh, subset = sample  %in% c("2",
                                                               "3", "32"))
      
      #Mevorach
      seu_cd59 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/seu_cd59.qs")
      
      seu_cd59 <- subset(seu_cd59, subset = condition  == "Healthy")
      
      seu_cd59$sample <-"CD59-control"
      
      
      seu_cd55 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/seu_cd55.qs")
      
      seu_cd55 <- subset(seu_cd55, subset = condition  == "Healthy")
      
      seu_cd55$sample <-"CD55-control"
      
      
      healthy_list <- list(seu_cd55 = seu_cd55,
                           seu_cd59 = seu_cd59,
                           seu_hong = seu_hong,
                           seu_kang = seu_kang,
                           seu_ramesh = seu_ramesh)
      
      rm(seu_cd55, seu_cd59, seu_hong, seu_kang, seu_ramesh)  
      
      
      
    }
    
    # get one healthy_list for all cell types
    {
      plan(multisession, workers = 8)         # or multicore on Linux console
      options(future.globals.maxSize = 320 * 1024^3) # 100 GiB
      
      healthy_list <- lapply(healthy_list, get_all_cell_types)
      
    }
    
    #I manually removed actT from seu_cd55 and seu_cd59, it's in the garbage code
    qsave(healthy_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list.qs")
    
    healthy_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list.qs")
    
    #get five healthy_list, one for each cell type
    #run get_full_seurat on each object in each list
    #and qsave
    {
      
      subset_by_type <- function(seu, cell_type) {
        subset(seu, subset = final_cell_type == cell_type)
      }
      
      ## define the mapping: cell‑type → new‑list name -----------------------------
      cell_map <- c(
        "NK cells"      = "healthy_list_nk",
        "B cells"       = "healthy_list_b_cells",
        "Monocytes"     = "healthy_list_monocytes",
        "CD4+ T cells"  = "healthy_list_cd4_t_cells",
        "CD8+ T cells"  = "healthy_list_cd8_t_cells"
      )
      
      ## generate five lists in one loop -------------------------------------------
      for (ct in names(cell_map)) {
        assign(
          cell_map[[ct]],
          lapply(healthy_list, subset_by_type, cell_type = ct)
        )
      }
      
      
      ## helper ─ retry get_full_seurat until it works ------------------------------
      process_with_retry <- function(seu, start_nbin = 23, min_nbin = 5) {
        for (nb in seq(start_nbin, min_nbin, by = -1)) {
          result <- tryCatch(
            get_full_seurat(seu, nbin_for_cell_cycle = nb),
            error = function(e) {
              message("  » nbin = ", nb, " failed: ", e$message)
              NULL
            }
          )
          if (!is.null(result)) return(result)
        }
        stop("get_full_seurat failed for all nbin values")
      }
      
      ## the five cell‑type lists ---------------------------------------------------
      celltype_lists <- c(
        "healthy_list_nk",
        "healthy_list_b_cells",
        "healthy_list_monocytes",
        "healthy_list_cd4_t_cells",
        "healthy_list_cd8_t_cells"
      )
      
      ## main loop ------------------------------------------------------------------
      for (lst_name in celltype_lists) {
        message("▶ Processing ", lst_name)
        
        lst_in  <- get(lst_name)                            # original list
        lst_out <- lapply(lst_in, process_with_retry)       # retry logic per object
        
        assign(lst_name, lst_out, envir = .GlobalEnv)       # keep the updated list
        
        qsave(
          lst_out,
          sprintf(
            "/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/%s.qs",
            lst_name
          )
        )
      }
      
    }
    
    #run get_healthy_list on each list
    {
      ## names of the five cell‑type lists you want to rebuild ---------------------
      celltype_lists <- c(
        "healthy_list_nk",
        "healthy_list_b_cells",
        "healthy_list_monocytes",
        "healthy_list_cd4_t_cells",
        "healthy_list_cd8_t_cells"
      )
      
      ## loop: rebuild each list in‑place via get_healthy_list ---------------------
      for (lst_name in celltype_lists) {
        message("▶ Re‑creating ", lst_name)
        
        raw_five    <- get(lst_name)          # the list that still contains the 5 Seurat objects
        new_hlist   <- get_healthy_list(raw_five)
        
        assign(lst_name, new_hlist, envir = .GlobalEnv)   # overwrite with the new healthy_list
      }
      
      
    }
  }
  
  #Create subsets of seu_cd59
  {
    
    ## helper ─ retry get_full_seurat_object with shrinking nbin ------------------
    retry_full_pipeline <- function(seu, start_nbin = 23, min_nbin = 5) {
      for (nb in seq(start_nbin, min_nbin, by = -1)) {
        res <- tryCatch(
          get_full_seurat(seu, nbin_for_cell_cycle = nb),
          error = function(e) {
            message("   » nbin = ", nb, " failed: ", e$message)
            NULL
          }
        )
        if (!is.null(res)) return(res)        # succeed → stop retrying
      }
      stop("get_full_seurat_object failed for all nbin values")
    }
    
    ## mapping: cell‑type → object name + output folder --------------------------
    cell_map <- list(
      "NK cells"      = list(obj = "seu_cd59_nk",          folder = "NK"),
      "Monocytes"     = list(obj = "seu_cd59_monocytes",   folder = "Monocytes"),
      "B cells"       = list(obj = "seu_cd59_b_cells",     folder = "B_cells"),
      "CD4+ T cells"  = list(obj = "seu_cd59_cd4_t_cells", folder = "CD4_T_cells"),
      "CD8+ T cells"  = list(obj = "seu_cd59_cd8_t_cells", folder = "CD8_T_cells")
    )
    
    ## main loop ------------------------------------------------------------------
    for (ctype in names(cell_map)) {
      info <- cell_map[[ctype]]
      message("▶ Processing ", ctype, " …")
      
      # 1. subset
      seu_sub  <- subset(seu_cd59, subset = final_cell_type == ctype)
      
      # 2. run full pipeline with automatic retry on nbin
      seu_proc <- retry_full_pipeline(seu_sub)
      
      # 3. keep in workspace
      assign(info$obj, seu_proc, envir = .GlobalEnv)
      
      # 4. save
      out_path <- sprintf(
        "/bigdata/bgutman/Collabs/Dror mevorach/data/%s/%s.qs",
        info$folder, info$obj
      )
      dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
      qsave(seu_proc, out_path)
    }
    
    
    
  }
  
  #list_for_one_on_one
  {
    ## ---------------------------------------------------------------------------
    ## Build `list_for_one_on_one`
    ##   • five top‑level entries named by cell type
    ##   • each entry contains: healthy_list, patient_list, folder
    ## ---------------------------------------------------------------------------
    ## ---------------------------------------------------------------------------
    ## Build `list_for_one_on_one`  (folder names now compact)
    ## ---------------------------------------------------------------------------
    type_map <- list(
      "NK cells"      = list(healthy = "healthy_list_nk",
                             seu     = "seu_cd59_nk",
                             folder  = "NK"),
      "B cells"       = list(healthy = "healthy_list_b_cells",
                             seu     = "seu_cd59_b_cells",
                             folder  = "B_cells"),
      "Monocytes"     = list(healthy = "healthy_list_monocytes",
                             seu     = "seu_cd59_monocytes",
                             folder  = "Monocytes"),
      "CD4+ T cells"  = list(healthy = "healthy_list_cd4_t_cells",
                             seu     = "seu_cd59_cd4_t_cells",
                             folder  = "CD4_T_cells"),
      "CD8+ T cells"  = list(healthy = "healthy_list_cd8_t_cells",
                             seu     = "seu_cd59_cd8_t_cells",
                             folder  = "CD8_T_cells")
    )
    
    list_for_one_on_one <- setNames(
      lapply(names(type_map), function(ct) {
        info <- type_map[[ct]]
        
        list(
          healthy_list = get(info$healthy, envir = .GlobalEnv),
          patient_list = list(
            seu_mevorach_cd59 = list(
              seu       = get(info$seu, envir = .GlobalEnv),
              condition = "CD59.Def",
              dataset   = "Mevorach",
              sample    = "CD59_Patient"
            )
          ),
          folder = info$folder
        )
      }),
      names(type_map)
    )
    
    ## list_for_one_on_one is now ready and uses folder names:
    ## "NK", "Monocytes", "B_cells", "CD4_T_cells", "CD8_T_cells"
    
    
    
  }
  
  #run one on one
  {
    options(future.globals.maxSize = 320 * 1024 * 1024 * 1024)  # 320 GB
    
    ## ---------------------------------------------------------------------------
    ## loop through the five cell‑type bundles in list_for_one_on_one
    ## ---------------------------------------------------------------------------
    for (ct in names(list_for_one_on_one)) {
      
      bundle       <- list_for_one_on_one[[ct]]
      healthy_list <- bundle$healthy_list
      patient_list <- bundle$patient_list      # always length‑1 but loop keeps it generic
      folder       <- bundle$folder            # "NK", "Monocytes", …
      
      message("▶ Cell type: ", ct, "  (folder = ", folder, ")")
      
      for (current_seu_patient in patient_list) {
        
        seu_patient <- current_seu_patient$seu
        output_list <- list()
        
        for (seu in healthy_list) {
          
          seu_control <- seu$seu
          
          # run your one‑on‑one comparison
          output <- one_on_one(
            seu_patient, seu_control,
            patient_type        = current_seu_patient$condition,
            control_type        = "Healthy",
            patient_orig_ident  = current_seu_patient$dataset
          )
          
          # merge into output_list by control sample (orig.ident)
          orig_ident_value <- as.character(unique(seu_control$orig.ident)[1])
          
          if (orig_ident_value %in% names(output_list)) {
            output_list[[orig_ident_value]] <-
              append(output_list[[orig_ident_value]], output[[orig_ident_value]])
          } else {
            output_list[[orig_ident_value]] <- output[[orig_ident_value]]
          }
        }
        
        ## save --------------------------------------------------------------------
        out_dir <- sprintf(
          "/bigdata/bgutman/Collabs/Dror mevorach/data/%s/output_lists/neuroimmunological-healthy",
          folder
        )
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        
        out_path <- sprintf(
          "%s/%s_%s_%s_outputlist.qs",
          out_dir,
          current_seu_patient$condition,
          current_seu_patient$dataset,
          current_seu_patient$sample
        )
        
        qsave(output_list, out_path)
        message("   ✔ saved ", out_path)
        gc()
      }
    }
    
  }
  
  #run create_meta_table
  {
    outputlist_nk <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/output_lists/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_outputlist.qs")
    outputlist_monocytes <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/output_lists/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_outputlist.qs")
    outputlist_b_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/output_lists/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_outputlist.qs")
    outputlist_cd4_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/output_lists/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_outputlist.qs")
    outputlist_cd8_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/output_lists/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_outputlist.qs")
    
    list_of_outputlists <- list(outputlist_nk = outputlist_nk,
                                outputlist_monocytes = outputlist_monocytes,
                                outputlist_b_cells = outputlist_b_cells,
                                outputlist_cd4_t_cells = outputlist_cd4_t_cells,
                                outputlist_cd8_t_cells = outputlist_cd8_t_cells)
    
    rm(outputlist_nk, outputlist_monocytes, outputlist_b_cells, outputlist_cd4_t_cells, outputlist_cd8_t_cells)
    
    ## --- 1 · run create_meta_table() on each output list -------------------------
    meta_tables_list <- lapply(list_of_outputlists, create_meta_table)
    
    ## --- 2 · rename the list elements so they read ‘meta_table_*’ -----------------
    names(meta_tables_list) <- sub("^outputlist_", "meta_table_", names(meta_tables_list))
    
    qsave(meta_tables_list$meta_table_nk, "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    qsave(meta_tables_list$meta_table_monocytes, "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    qsave(meta_tables_list$meta_table_b_cells, "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    qsave(meta_tables_list$meta_table_cd4_t_cells, "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    qsave(meta_tables_list$meta_table_cd8_t_cells, "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    
  }
  
  #run get_integration
  {
    seu_cd59 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/seu_cd59_only_patient_full_annotation.qs")
    healthy_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list.qs")
    
    options(future.globals.maxSize = 320 * 1024 * 1024 * 1024)  # 320 GB
    
    integrated_list <- get_integration(seu_patient = seu_cd59, healthy_list_for_integration = healthy_list, patient_condition = "CD59.Def")
    
    # qsave(return_lists$NK, "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    # qsave(return_lists$B_cells, "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    # qsave(return_lists$Monocytes, "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    # qsave(return_lists$CD4_T_cells, "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    # qsave(return_lists$CD8_T_cells, "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    # 
    
    
    
  }
  
  #run get_meta_int_comparison_pds_for_all_the_samples for all cell types
  {
    
    
    integrated_path =  "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/integrated_objects/neuroimmunulogical-healthy"
    meta_path = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/"
    get_meta_int_comparison_pds_for_all_the_samples(
      integrated_path = integrated_path,
      meta_path = meta_path,
      files_to_exclude_int = NULL,
      files_to_exclude_meta = NULL,
      cell_type = "NKs",
      save_dir_gsea = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/full_gsea_ora_results/neuroimmunological-healthy/GSEA"
    )
    
    
    
    
    integrated_path =  "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/integrated_objects/neuroimmunulogical-healthy"
    meta_path = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/"
    get_meta_int_comparison_pds_for_all_the_samples(
      integrated_path = integrated_path,
      meta_path = meta_path,
      files_to_exclude_int = NULL,
      files_to_exclude_meta = NULL,
      cell_type = "Monocytes",
      save_dir_gsea = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/full_gsea_ora_results/neuroimmunological-healthy/GSEA"
      
    )
    
    
    integrated_path =  "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/integrated_objects/neuroimmunulogical-healthy"
    meta_path = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/"
    get_meta_int_comparison_pds_for_all_the_samples(
      integrated_path = integrated_path,
      meta_path = meta_path,
      files_to_exclude_int = NULL,
      files_to_exclude_meta = NULL,
      cell_type = "B_cells",
      save_dir_gsea = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA")
    
    
    integrated_path =  "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/integrated_objects/neuroimmunulogical-healthy"
    meta_path = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/"
    get_meta_int_comparison_pds_for_all_the_samples(
      integrated_path = integrated_path,
      meta_path = meta_path,
      files_to_exclude_int = NULL,
      files_to_exclude_meta = NULL,
      cell_type = "CD4_T_cells",
      save_dir_gsea = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA"
      
    )
    
    integrated_path =  "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/integrated_objects/neuroimmunulogical-healthy"
    meta_path = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
    get_meta_int_comparison_pds_for_all_the_samples(
      integrated_path = integrated_path,
      meta_path = meta_path,
      files_to_exclude_int = NULL,
      files_to_exclude_meta = NULL,
      cell_type = "CD8_T_cells",
      save_dir_gsea = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA"
      
    )
    
  }
  
  
}

#Basic pipeline for neuroimmunulogical samples
{
  #removing activated T cells from seu_cd55_patient
  {
    
    seu_cd55_patient <- qread("/bigdata/barakgutman/Collabs/Dror_mevorach/CD59/data/cd55/CD55_data/count/seu_cd55.qs")
    
    seu_cd55_patient <- subset(seu_cd55_patient, subset = condition  == "patient")
    
    
    
    seu_cd55_patient <- subset(
      seu_cd55_patient,
      subset = !(seu_cd55_patient$sample %in% c("actT-pt"))
    )
    
    
    seu_cd55_patient$sample <-"CD55-patient"
    
    
    
    
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
        "sample",     "condition"
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
    
    
    seu_cd55_patient <- get_slim_seu(seu_cd55_patient)
    
    
    
    cells    <- colnames(seu_cd55_patient)        # all cell barcodes
    prefixes <- unique(sub("_.*$", "", cells))    # drop everything after first “_”
    sort(prefixes)    
    
    
    
    seu_cd55_patient <- get_full_seurat(seu_cd55_patient)
    
    
    seu_cd55_patient <- get_all_cell_types(seu_cd55_patient, n_cores = 8)
    seu_cd55_patient$condition <- "CD55.Def"
    seu_cd55_patient$orig.ident <- "Mevorach"
    seu_cd55_patient$sample <- "CD55-patient"
    qsave(seu_cd55_patient, "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial/seu_CD55.Def_Mevorach_CD55-patient.qs")
    
  }
  
  
  # create neuroimmunulogical list
  {
    seu_kang <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS_and_other_disorders/Kang/seu_kang.qs")
    seu_kang <- subset(seu_kang, subset = condition  %in% c("MS", "Ab-mediated IDD"))
    
    seu_esaulova <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS_and_other_disorders/Esaulova/seu_esaulova_ms_mog.qs")
    seu_esaulova$orig.ident <- "Esaulova"
    seu_esaulova$condition <- ifelse(seu_esaulova$Sample %in% c("ms1_blood", "ms2_blood"), "MS", "Ab-mediated IDD")
    seu_esaulova$sample <- seu_esaulova$Sample
    seu_esaulova$Sample <- NULL
    
    seu_ramesh_3 <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS_and_other_disorders/Ramesh/seu_ramesh_3.qs")
    seu_ramesh_3 <- subset(seu_ramesh_3, subset = sample  == "16")
    seu_ramesh_3$orig.ident <- "Ramesh"
    seu_ramesh_3$condition <- "MS"
    
    
    seu_ramesh_5 <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS_and_other_disorders/Ramesh/seu_ramesh_5.qs")
    seu_ramesh_5 <- subset(seu_ramesh_5, subset = condition  %in% c("MS", "CIS", "NMO"))
    seu_ramesh_5$condition <- ifelse(
      seu_ramesh_5$condition == "CIS", "MS",
      ifelse(seu_ramesh_5$condition == "NMO", "Ab-mediated IDD", seu_ramesh_5$condition)
    )
    
    # seu_schafflick <- qread("/bigdata/barakgutman/neuroimmunulogical_project/scRNA-seq/scRNA_PBMCs/data/MS/Schafflick/seu_schafflick.qs")
    # seu_schafflick$orig.ident <- "Schafflick"
    # seu_schafflick <- subset(seu_schafflick, subset = condition  == "MS")
    
    neuroimmunulogical_list <- list(seu_kang = seu_kang,
                                    seu_esaulova = seu_esaulova,
                                    seu_ramesh_3 = seu_ramesh_3,
                                    seu_ramesh_5 = seu_ramesh_5)
    
    rm(seu_kang, seu_esaulova, seu_ramesh_3, seu_ramesh_5)
    
    for (nm in names(neuroimmunulogical_list)) {
      seu <- neuroimmunulogical_list[[nm]]
      
      if ("condition" %in% colnames(seu@meta.data)) {
        seu$condition[seu$condition == "Ab-mediated IDD"] <- "Ab.mediated.IDD"
      }
      
      neuroimmunulogical_list[[nm]] <- seu
    }
    
    
    
    
    
  }
  
  #run get_all_cell_types for each sample in neuroimmunulogical list
  {
    
    options(future.globals.maxSize = 320 * 1024 * 1024 * 1024)  # 320 GB
    
    # assumes: neuroimmunulogical_list, get_slim_seu, get_full_seurat, get_all_cell_types are defined
    # and library(qs) is available for qsave()
    
    base_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial"
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
    
    # helper: try get_full_seurat with decreasing nbin (23 → … → 1)
    run_full_with_retry <- function(seu_obj) {
      bins <- c(seq(23, 3, by = -4), 1)  # 23,19,15,11,7,3,1
      last_err <- NULL
      for (nb in bins) {
        message("  • Trying get_full_seurat(nbin_for_cell_cycle = ", nb, ")")
        res <- try(get_full_seurat(seu_obj, nbin_for_cell_cycle = nb, compute_cell_cycle = FALSE), silent = TRUE)
        if (!inherits(res, "try-error")) {
          if (!is.null(last_err)) message("    ↳ Succeeded after previous error: ", last_err)
          return(res)
        } else {
          err_msg <- as.character(res)
          last_err <- err_msg
          message("    ↳ Failed: ", err_msg)
        }
      }
      return(structure(NULL, last_error = last_err))
    }
    
    for (nm in names(neuroimmunulogical_list)) {
      seu <- neuroimmunulogical_list[[nm]]
      
      # 1) Slim
      seu <- get_slim_seu(seu)
      
      # 2) Sanity: need a 'sample' column
      if (!"sample" %in% colnames(seu@meta.data)) {
        message("Skipping '", nm, "' (no 'sample' column).")
        next
      }
      
      # 3) For each sample
      for (smp in sort(unique(seu$sample))) {
        message("Processing dataset '", nm, "' sample '", smp, "'")
        
        # subset one sample (avoid name collision with base::sample)
        seu_sample <- subset(seu, subset = sample == smp)
        
        # rebuild + annotate with retry for cell-cycle bins
        seu_sample_full <- run_full_with_retry(seu_sample)
        if (is.null(seu_sample_full)) {
          last_err <- attr(seu_sample_full, "last_error")
          message("  ✖ Skipping sample '", smp, "' after all retries failed.",
                  if (!is.null(last_err)) paste0(" Last error: ", last_err) else "")
          next
        }
        
        # downstream cell-type annotation (unchanged)
        seu_sample_full <- get_all_cell_types(seu_sample_full, n_cores = 8)
        
        # extract single values for naming (fallbacks if missing/ambiguous)
        cond  <- unique(seu_sample_full$condition);  cond  <- if (length(cond)  >= 1) cond[1]  else "NA"
        oid   <- unique(seu_sample_full$orig.ident); oid   <- if (length(oid)   >= 1) oid[1]   else "NA"
        smpl  <- unique(seu_sample_full$sample);     smpl  <- if (length(smpl)  >= 1) smpl[1]  else "NA"
        
        # sanitize any spaces or slashes
        sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
        
        saving_address <- file.path(
          base_dir,
          paste0("seu_", sanitize(cond), "_", sanitize(oid), "_", sanitize(smpl), ".qs")
        )
        
        qsave(seu_sample_full, saving_address)
        message("  ✓ Saved: ", saving_address)
      }
    }
  }
  
  #create cell_type_lists from neuroimmunulogical list
  {
    
    
    base_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial"
    
    # --- 1) Initialize empty lists ---
    nk_list           <- list()
    monocytes_list    <- list()
    b_cells_list      <- list()
    cd4_t_cells_list  <- list()
    cd8_t_cells_list  <- list()
    
    # --- 2) Gather .qs files ---
    qs_files <- list.files(base_dir, pattern = "\\.qs$", full.names = TRUE)
    if (length(qs_files) == 0) {
      stop(sprintf("No .qs files found in: %s", base_dir))
    }
    
    # --- 3) Process each Seurat object ---
    for (f in qs_files) {
      base_name <- file_path_sans_ext(basename(f))
      
      # Load the Seurat object
      seu <- tryCatch(qread(f), error = function(e) {
        message(sprintf("Skipping '%s': failed to qread(). Error: %s", base_name, e$message))
        return(NULL)
      })
      if (is.null(seu)) next
      
      # Sanity check
      if (!"final_cell_type" %in% colnames(seu@meta.data)) {
        message(sprintf("Skipping '%s': no 'final_cell_type' column in meta.data.", base_name))
        next
      }
      
      # Subsets (note: CD8 fixed; your snippet had CD4 twice)
      seu_nk          <- tryCatch(subset(seu, subset = final_cell_type == "NK cells"),        error = function(e) NULL)
      seu_monocytes   <- tryCatch(subset(seu, subset = final_cell_type == "Monocytes"),       error = function(e) NULL)
      seu_b_cells     <- tryCatch(subset(seu, subset = final_cell_type == "B cells"),         error = function(e) NULL)
      seu_cd4_t_cells <- tryCatch(subset(seu, subset = final_cell_type == "CD4+ T cells"),    error = function(e) NULL)
      seu_cd8_t_cells <- tryCatch(subset(seu, subset = final_cell_type == "CD8+ T cells"),    error = function(e) NULL)
      
      # Helper to get cell counts safely
      n_cells <- function(x) { if (is.null(x)) 0 else length(Cells(x)) }
      
      # Append if non-empty
      if (n_cells(seu_nk)          > 0) nk_list[[base_name]]          <- seu_nk
      if (n_cells(seu_monocytes)   > 0) monocytes_list[[base_name]]   <- seu_monocytes
      if (n_cells(seu_b_cells)     > 0) b_cells_list[[base_name]]     <- seu_b_cells
      if (n_cells(seu_cd4_t_cells) > 0) cd4_t_cells_list[[base_name]] <- seu_cd4_t_cells
      if (n_cells(seu_cd8_t_cells) > 0) cd8_t_cells_list[[base_name]] <- seu_cd8_t_cells
      
      # Progress message
      message(sprintf(
        "Finished: %s  (NK=%d, Monocytes=%d, B=%d, CD4=%d, CD8=%d)",
        base_name,
        n_cells(seu_nk),
        n_cells(seu_monocytes),
        n_cells(seu_b_cells),
        n_cells(seu_cd4_t_cells),
        n_cells(seu_cd8_t_cells)
      ))
    }
    
    # --- 4) Optional: quick summary ---
    message(sprintf("Summary: NK=%d | Monocytes=%d | B=%d | CD4=%d | CD8=%d",
                    length(nk_list),
                    length(monocytes_list),
                    length(b_cells_list),
                    length(cd4_t_cells_list),
                    length(cd8_t_cells_list)))
    
    
    # qsave(nk_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/nk_neuroimmunulogocial_list.qs")
    # qsave(monocytes_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/monocytes_neuroimmunulogocial_list.qs")
    # qsave(b_cells_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/b_cells_neuroimmunulogocial_list.qs")
    # qsave(cd4_t_cells_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/cd4_t_cells_list_neuroimmunulogocial_list.qs")
    # qsave(cd8_t_cells_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/cd8_t_cells_list_neuroimmunulogocial_list.qs")
    # 
    
  }
  
  #run one_one_one on each cell type
  {
    
    # # nk_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/nk_neuroimmunulogocial_list.qs")
    # # keep_nks <- c(
    # #   "seu_Ab.mediated.IDD_Esaulova_mog_blood",
    # #   "seu_Ab.mediated.IDD_Kang_CYG",
    # #   "seu_CD55.Def_Mevorach_CD55-patient",
    # #   "seu_MS_Ramesh_10",
    # #   "seu_MS_Ramesh_16",
    # #   "seu_MS_Schafflick_MS49131"
    # # ) 
    # # 
    # # nk_list <- nk_list[intersect(names(nk_list), keep_nks)]
    # 
    # 
    # nk_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/nk_neuroimmunulogocial_list.qs")
    # 
    # # keep_nks <- c(
    # #   "seu_Ab.mediated.IDD_Esaulova_mog_blood",
    # #   "seu_Ab.mediated.IDD_Kang_CYG",
    # #   "seu_CD55.Def_Mevorach_CD55-patient",
    # #   "seu_MS_Ramesh_10",
    # #   "seu_MS_Ramesh_16",
    # #   "seu_MS_Schafflick_MS49131"
    # # )
    # # 
    # # # Remove all objects listed in keep_nks
    # # nk_list <- nk_list[setdiff(names(nk_list), keep_nks)]
    # 
    # file_basenames <- basename(files_nk)
    # 
    # # Remove "_outputlist.qs" suffix
    # file_basenames <- sub("_outputlist\\.qs$", "", file_basenames)
    # 
    # # Prepend "seu_" to match nk_list names
    # done_samples <- paste0("seu_", file_basenames)
    # 
    # # Remove them from nk_list
    # nk_list <- nk_list[setdiff(names(nk_list), done_samples)]
    
    
    # monocytes_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/monocytes_neuroimmunulogocial_list.qs")
    # b_cells_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/b_cells_neuroimmunulogocial_list.qs")
    # 
    
    cd4_t_cells_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/cd4_t_cells_list_neuroimmunulogocial_list.qs")
    cd8_t_cells_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/other_conditions/neuroimmunulogocial_lists_per_cell_type/cd8_t_cells_list_neuroimmunulogocial_list.qs")
    
    
    # healthy_list_nk <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list_nk.qs")
    # healthy_list_monocytes <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list_monocytes.qs")
    # healthy_list_b_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list_b_cells.qs")
    #
    
    healthy_list_cd4_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list_cd4_t_cells.qs")
    healthy_list_cd8_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/healthy controls/healthy_list_cd8_t_cells.qs")
    
    
    
    options(future.globals.maxSize = 400 * 1024 * 1024 * 1024)  # 320 GB
    
    # --- Helpers (unchanged logic, minimal edits) ---
    get_meta_scalar <- function(seu, fields, default = NA_character_) {
      for (fld in fields) {
        if (!is.null(seu@meta.data[[fld]])) {
          vals <- unique(as.character(seu@meta.data[[fld]]))
          vals <- vals[!is.na(vals) & nzchar(vals)]
          if (length(vals) >= 1) return(vals[1])
        }
      }
      default
    }
    
    parse_name_fallback <- function(x) {
      nm <- sub("^seu_", "", x)
      parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
      condition <- if (length(parts) >= 1) parts[1] else nm
      dataset   <- if (length(parts) >= 2) parts[2] else "unknownDataset"
      sample    <- if (length(parts) >= 3) paste(parts[3:length(parts)], collapse = "_") else "unknownSample"
      list(condition = condition, dataset = dataset, sample = sample)
    }
    
    safe_control_key <- function(seu_control, fallback_name) {
      key <- get_meta_scalar(seu_control, c("orig.ident", "orig_ident", "dataset"))
      if (is.na(key) || !nzchar(key)) key <- fallback_name
      as.character(key)
    }
    
    patient_labels <- function(seu_patient, patient_name) {
      condition <- get_meta_scalar(seu_patient, c("condition", "patient_type", "disease", "group"))
      dataset   <- get_meta_scalar(seu_patient, c("dataset", "orig.ident", "orig_ident", "source"))
      sample    <- get_meta_scalar(seu_patient, c("sample", "sample_id", "sample.name", "sample_name"))
      if (is.na(condition) || is.na(dataset) || is.na(sample) ||
          !nzchar(condition) || !nzchar(dataset) || !nzchar(sample)) {
        parsed <- parse_name_fallback(patient_name)
        if (is.na(condition) || !nzchar(condition)) condition <- parsed$condition
        if (is.na(dataset)   || !nzchar(dataset))   dataset   <- parsed$dataset
        if (is.na(sample)    || !nzchar(sample))    sample    <- parsed$sample
      }
      list(condition = condition, dataset = dataset, sample = sample)
    }
    
    # --- Runner for a cell type (tiny tweak: don't stop() on empty lists) ---
    run_one_on_one_for_celltype <- function(patient_list, healthy_list, folder = "NK") {
      if (length(patient_list) == 0) {
        message("▶ Cell type: ", folder, " — patient_list is empty (skipping).")
        return(invisible(NULL))
      }
      if (length(healthy_list) == 0) {
        message("▶ Cell type: ", folder, " — healthy_list is empty (skipping).")
        return(invisible(NULL))
      }
      
      message("▶ Cell type: ", folder)
      
      for (patient_name in names(patient_list)) {
        tryCatch({
          seu_patient <- patient_list[[patient_name]]
          if (!inherits(seu_patient, "Seurat")) {
            message("Skipping '", patient_name, "': not a Seurat object.")
            next
          }
          
          ids <- patient_labels(seu_patient, patient_name)
          message(sprintf("  • Patient: %s (condition=%s, dataset=%s, sample=%s)",
                          patient_name, ids$condition, ids$dataset, ids$sample))
          
          output_list <- list()
          
          for (control_name in names(healthy_list)) {
            tryCatch({
              seu_control <- healthy_list[[control_name]]
              if (!inherits(seu_control, "Seurat")) {
                message("    - Skipping control '", control_name, "': not a Seurat object.")
                next
              }
              
              control_key <- safe_control_key(seu_control, control_name)
              
              output <- tryCatch(
                one_on_one(
                  seu_patient, seu_control,
                  patient_type       = ids$condition,
                  control_type       = "Healthy",
                  patient_orig_ident = ids$dataset
                ),
                error = function(e) {
                  message("    - one_on_one failed for control '", control_name, "': ", e$message)
                  NULL
                }
              )
              if (is.null(output)) next
              
              if (control_key %in% names(output_list)) {
                output_list[[control_key]] <- append(output_list[[control_key]], output[[control_key]])
              } else {
                if (!is.null(output[[control_key]])) {
                  output_list[[control_key]] <- output[[control_key]]
                } else if (length(output) >= 1) {
                  output_list[[control_key]] <- output[[1]]
                } else {
                  message("    - Warning: unexpected output structure for control '", control_name, "'.")
                }
              }
            }, error = function(e) {
              message("    - Comparison crashed for control '", control_name, "': ", e$message, " (continuing)")
            })
          } # controls
          
          out_dir <- sprintf(
            "/bigdata/bgutman/Collabs/Dror mevorach/data/%s/output_lists/neuroimmunological-healthy",
            folder
          )
          dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
          
          out_path <- sprintf(
            "%s/%s_%s_%s_outputlist.qs",
            out_dir, ids$condition, ids$dataset, ids$sample
          )
          
          tryCatch({
            qsave(output_list, out_path)
            message("     ✔ saved ", out_path)
          }, error = function(e) {
            message("     ✖ failed to save ", out_path, ": ", e$message, " (continuing)")
          })
          
          gc()
        }, error = function(e) {
          message("  ✖ Patient '", patient_name, "' crashed: ", e$message, " (continuing to next patient)")
        })
      }
    }
    
    # --- Execute for all three lists ---
    # # NK
    # run_one_on_one_for_celltype(
    #   patient_list = nk_list,
    #   healthy_list =  healthy_list_nk,
    #   folder       = "NK"
    # )
    # 
    # gc()
    # 
    # # Monocytes
    # run_one_on_one_for_celltype(
    #   patient_list = monocytes_list,
    #   healthy_list = healthy_list_monocytes,
    #   folder       = "Monocytes"
    # )
    # 
    # gc()
    # 
    # 
    # # B cells
    # run_one_on_one_for_celltype(
    #   patient_list = b_cells_list,
    #   healthy_list = healthy_list_b_cells,
    #   folder       = "B_cells"
    # )
    # 
    # gc()
    
    # CD4_T_CELLs
    run_one_on_one_for_celltype(
      patient_list = cd4_t_cells_list,
      healthy_list = healthy_list_cd4_t_cells,
      folder       = "CD4_T_cells"
    )
    
    gc()
    
    
    #CD8_T_cells
    run_one_on_one_for_celltype(
      patient_list = cd8_t_cells_list,
      healthy_list = healthy_list_cd8_t_cells,
      folder       = "CD8_T_cells"
    )
    
    gc()
    
    
    
    
  }
  
  
  #make sure we have all 12 original_tables in each output_list
  {
    
    
    # Folders per cell type
    cell_paths <- list(
      "CD4 + T cells"       = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/output_lists/neuroimmunological-healthy/",
      "CD8 + T cells" = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/output_lists/neuroimmunological-healthy/"
    )
    
    # Expected datasets → samples (12 total)
    expected <- list(
      "Mevorach" = c("CD55-control", "CD59-control"),
      "Hong"     = paste0("HC-", 1:5),
      "Kang"     = c("JSB", "KYO"),
      "Ramesh"   = c("2", "3", "32")
    )
    
    is_valid_table <- function(x) {
      if (is.null(x)) return(FALSE)
      if (length(x) == 1L && is.na(x)) return(FALSE)
      if (is.data.frame(x)) return(nrow(x) > 0)
      FALSE
    }
    
    check_file <- function(file, cell_type, expected_map) {
      obj <- tryCatch(qread(file), error = function(e) NULL)
      
      missing <- character(0)
      for (ds in names(expected_map)) {
        for (smp in expected_map[[ds]]) {
          ot <- tryCatch({
            if (is.null(obj) || is.null(obj[[ds]]) || is.null(obj[[ds]][[smp]])) {
              NA
            } else {
              obj[[ds]][[smp]][["original_table"]]
            }
          }, error = function(e) NA)
          
          if (!is_valid_table(ot)) {
            missing <- c(missing, paste(ds, smp, sep = "/"))
          }
        }
      }
      
      if (length(missing) > 0) {
        data.frame(
          cell_type        = cell_type,
          object           = sub("\\.qs$", "", basename(file)),
          num_of_missing   = length(missing),
          names_of_missing = paste(missing, collapse = ", "),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }  
    # ── Run the checks ─────────────────────────────────────────────────────────────
    res_list <- list()
    
    for (ct in names(cell_paths)) {
      files <- list.files(cell_paths[[ct]], pattern = "\\.qs$", full.names = TRUE)
      if (length(files) == 0L) next
      ct_rows <- lapply(files, check_file, cell_type = ct, expected_map = expected)
      ct_rows <- Filter(Negate(is.null), ct_rows)
      if (length(ct_rows) > 0L) res_list[[ct]] <- do.call(rbind, ct_rows)
    }
    
    missing_df <- if (length(res_list) > 0L) {
      do.call(rbind, res_list)
    } else {
      data.frame(
        cell_type = character(),
        object = character(),
        num_of_missing = integer(),
        names_of_missing = character(),
        stringsAsFactors = FALSE
      )
    }
    
    # Optional: sort rows
    if (nrow(missing_df)) {
      missing_df <- missing_df[order(missing_df$cell_type, missing_df$object), ]
    }
    
    # View the result
    missing_df
    
    
  }
  
  # run create_meta_table
  {
    library(qs)
    
    cell_paths <- list(
      "NKs"            = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/output_lists/neuroimmunological-healthy/",
      "Monocytes"      = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/output_lists/neuroimmunological-healthy/",
      "B Cells"        = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/output_lists/neuroimmunological-healthy/",
      "CD4 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/output_lists/neuroimmunological-healthy/",
      "CD8 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/output_lists/neuroimmunological-healthy/"
    )
    
    for (ct in names(cell_paths)) {
      in_dir  <- cell_paths[[ct]]
      out_dir <- sub("output_lists", "meta_tables", in_dir)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      files <- list.files(in_dir, pattern = "\\.qs$", full.names = TRUE)
      files <- files[!grepl("CD59.Def_Mevorach_CD59_Patient_outputlist\\.qs$", files)]
      
      for (f in files) {
        obj_name <- sub("\\.qs$", "", basename(f))
        
        message("▶ Starting: ", obj_name, " (", ct, ")")
        
        assign(obj_name, qread(f))
        meta_table <- create_meta_table(get(obj_name))
        
        out_file <- file.path(out_dir, paste0(obj_name, "_meta_table.qs"))
        qsave(meta_table, out_file)
        
        message("✓ Finished: ", obj_name, " → saved to ", out_file)
      }
    }
    
    
  }
  
 
}


#GSEA dotplot
{

  #This is an old code, new version is at /bigdata/bgutman/Collabs/Dror mevorach/scripts/old_and_irelevant/NK_anlysis_one_on_one_manually_doing-stuff_and_old_functions.Rmd
  
  NK_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
  monocytes_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
  B_cells_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
  CD4_T_cells_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
  CD8_T_cells_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
  
  gsea_list <- list(
    "NK" = NK_gsea,
    "Monocytes" = monocytes_gsea,
    "B Cells" = B_cells_gsea, 
    "CD4+ T cells" = CD4_T_cells_gsea,
    "CD8+ T cells" = CD8_T_cells_gsea
  )
  
  ID_vector <- c(
    "GO:0002250",
    "hsa04062",
    "HALLMARK_APOPTOSIS",
    # "HALLMARK_E2F_TARGETS",
    # "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_HYPOXIA",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    # "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    # "HALLMARK_TGF_BETA_SIGNALING",
    "GO:0042743",
    "GO:0019814",
    # "GO:0042775",
    "GO:0042101"
  )
  
  complex_lolipop_df <- data.frame(
    ID          = character(),
    Description = character(),
    Cell_type   = character(),
    dataset     = character(),
    NES         = numeric(),
    p.adjust    = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  complex_lolipop_df <- map2_dfr(
    gsea_list,                       # the list elements (data frames)
    names(gsea_list),                # their names (cell types)
    ~ .x %>%                         # .x = current data frame
      filter(ID %in% ID_vector) %>%# keep only the rows we want
      mutate(Cell_type = .y) %>%   # .y = current cell-type name
      dplyr::select(ID, Description, Cell_type, dataset, NES, p.adjust)
  )
  
  complex_lolipop_df <- complex_lolipop_df %>% 
    filter(
      !(
        ## 1) Monocytes or NK with GO:0070098
        (Cell_type %in% c("Monocytes", "NK") & ID == "GO:0002250") |
          ## 2) Monocytes or NK with GO:0042101
          (Cell_type %in% c("Monocytes", "NK","B Cells") & ID == "GO:0042101") |
          ## 3) Monocytes, NK, CD4+ T cells, or CD8+ T cells with GO:0019814
          (Cell_type %in% c("Monocytes", "NK", "CD4+ T cells", "CD8+ T cells") & ID == "GO:0019814")
      )
    )
  
  
  Changing_names_vector <- list(HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
                                # HALLMARK_E2F_TARGETS = "E2F Targets",
                                # HALLMARK_G2M_CHECKPOINT = "G2M Checkpoint",
                                # "mitochondrial ATP synthesis coupled electron transport" = 
                                #   "Mitochondrial ATP Synthesis Coupled Electron Transport",
                                HALLMARK_MTORC1_SIGNALING = "mTORC1 Signaling",
                                # HALLMARK_TGF_BETA_SIGNALING = "TGF Beta Signaling",
                                "chemokine-mediated signaling pathway"= "Chemokine Mediated Signaling Pathway")
  row_order <- c("Hydrogen Peroxide Metabolic Process",
                 "TNFA Signaling via NFKB",
                 "Inflammatory Response",
                 "Apoptosis",
                 "Hypoxia",
                 "Chemokine Signaling Pathway",
                 "Adaptive Immune Response",
                 "T Cell Receptor Complex",
                 "Immunoglobulin Complex"
                 )
  
  dotplot <- get_pathway_dotplot(
    df                       = complex_lolipop_df,
    manuall_changing_names   = Changing_names_vector,
    row_order                = row_order,                # or pass your vector
    celltype_order           = c("CD4+ T cells", "CD8+ T cells",
                                 "NK", "B Cells", "Monocytes")
  )
  
  

  # ggsave(
  #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/CD59.Def_pathway_dotplot.jpg",
  #   plot     = dotplot,
  #   width    = 25,
  #   height   = 20,
  #   units    = "in",
  #   dpi      = 600,
  #   device   = "jpeg",
  #   type     = "cairo",         # <– forces Cairo backend
  #   quality  = 100
  # )
  
}


#Individual cell type pathway lolipops
{
  #NKS
  {
    
    cd_59_nk_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
    
    
    
    list_of_highlight_pathways <- c(
      "GO:0042743",
      "R-HSA-500792",
      "GO:0034703",
      "GO:0099023",
      "HALLMARK_APOPTOSIS",
      "R-HSA-9663891",
      "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
      "R-HSA-448424",
      "HALLMARK_IL2_STAT5_SIGNALING",
      "HALLMARK_MTORC1_SIGNALING",
      "R-HSA-166058",
      "GO:0000118",
      "BIOCARTA_ETS_PATHWAY",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "GO:0042773",
      "GO:0001819",
      "hsa04612"
    )
    
    manuall_changing_names <- list(
      HALLMARK_IL2_STAT5_SIGNALING = "IL2 STAT5 Signaling",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
      "GPCR ligand binding" = "GPCR ligand binding",
      "HALLMARK_MTORC1_SIGNALING" = "mTORC1 Signaling",
      "ATP synthesis coupled electron transport" = "ATP Synthesis Coupled Electron Transport")
    
    cell_type ="NKs"
    condition = "CD59 def."
    sample_name = "CD59 patient"
    dataset_name = "Mevorach"
    

    
    try_lollipop <- get_pathway_lollipop(df = cd_59_nk_gsea,list_of_highlight_pathways = list_of_highlight_pathways,manuall_changing_names =manuall_changing_names,
                                         title = "PDE-Based GSEA: Key Pathways Enriched in CD59-Deficient NKs",lollipop_color = "#B62055")
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/NKs/CD59.Def_Mevorach_CD59_Patient/CD59.Def_NK_highlight_lolipop_plot.jpg",
      plot     = try_lollipop$p_pathway_lollipop,
      width    = 26,
      height   = 20,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",         # <– forces Cairo backend
      quality  = 100
    )
    
    
  }
  
  #monocytes
  
  {
    

    full_meta_gsea_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
    
    
    list_of_highlight_pathways <- c(
      "GO:0042775",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_TGF_BETA_SIGNALING",
      "HALLMARK_APOPTOSIS",
      "hsa04610",
      "GO:0045954",
      "hsa04062",
      "HALLMARK_HYPOXIA")
    
    
    # df <- try_pathways$meta_gsea_top
    
    manuall_changing_names <- list(
      HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory Response",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
      "HALLMARK_TGF_BETA_SIGNALING" = "TGF Beta Signaling",
      "mitochondrial ATP synthesis coupled electron transport" = "Mitochondrial ATP Synthesis Coupled Electron Transport")
    
    cell_type ="Monocytes"
    condition = "CD59 def."
    sample_name = "CD59 patient"
    dataset_name = "Mevorach"
    
    
    try_lollipop <- get_pathway_lollipop(full_meta_gsea_df,list_of_highlight_pathways,manuall_changing_names =manuall_changing_names, 
                                         title = "PDE-Based GSEA: Key Pathways Enriched in CD59-Deficient Monocytes",lollipop_color = "#56B4E9")
    

    
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/Monocytes/CD59.Def_Mevorach_CD59_Patient/CD59.Def_monocytes_meta_gsea_highlight_lolipop.jpg",
      plot     = try_lollipop$p_pathway_lollipop,
      width    = 26,
      height   = 20,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",         # <– forces Cairo backend
      quality  = 100
    )
    
  }
  
  #B cells
  {
    B_cells_full_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
    
    list_of_highlight_pathways <- c(
      "GO:0042775",
      "HALLMARK_E2F_TARGETS",
      "HALLMARK_G2M_CHECKPOINT",
      "GO:0019814",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_APOPTOSIS",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "hsa04062")
    
    
    
    manuall_changing_names <- list(
      HALLMARK_E2F_TARGETS = "E2F Targets",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
      "mitochondrial ATP synthesis coupled electron transport" ="Mitochondrial ATP Synthesis Coupled Electron Transport",
      "HALLMARK_G2M_CHECKPOINT" = "G2M Checkpoint")
    
    cell_type ="B Cells"
    condition = "CD59 def."
    sample_name = "CD59 patient"
    dataset_name = "Mevorach"
    
    
    try_lollipop <- get_pathway_lollipop(B_cells_full_gsea,list_of_highlight_pathways,manuall_changing_names =manuall_changing_names, 
                                         title = "PDE-Based GSEA: Key Pathways Enriched in CD59-Deficient B cells",lollipop_color = "#009E73")
    
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/B_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_B_cells_meta_gsea_highlight_lolipop.jpg",
      plot     = try_lollipop$p_pathway_lollipop,
      width    = 26,
      height   = 20,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",         # <– forces Cairo backend
      quality  = 100
    )
    
  }
  
  # CD4_T_cells
  {
    
    CD4_T_full_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
    
    list_of_highlight_pathways <- c(
      "HALLMARK_G2M_CHECKPOINT",
      "HALLMARK_E2F_TARGETS",
      "HALLMARK_HYPOXIA",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "GO:0042775",
      "HALLMARK_MTORC1_SIGNALING",
      "GO:0042101",
      "GO:0002250")
    
    # df <- try_pathways$meta_gsea_top
    
    manuall_changing_names <- list(
      HALLMARK_E2F_TARGETS = "E2F Targets",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
      "HALLMARK_G2M_CHECKPOINT" = "G2M Checkpoint",
      "mitochondrial ATP synthesis coupled electron transport" ="Mitochondrial ATP Synthesis Coupled Electron Transport",
      "HALLMARK_MTORC1_SIGNALING" ="mTORC1 Signaling",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB")
    
    cell_type ="CD4_T_cells"
    condition = "CD59 def."
    sample_name = "CD59 patient"
    dataset_name = "Mevorach"
    
    highlight_df <- CD4_T_full_gsea %>%
      filter(ID %in% list_of_highlight_pathways)
    
    highlight_df$p.adjust <- round(highlight_df$p.adjust, 2)
    
    
    try_lollipop <- get_pathway_lollipop(highlight_df,list_of_highlight_pathways,manuall_changing_names =manuall_changing_names, 
                                         title = "PDE-Based GSEA: Key Pathways Enriched in CD59-Deficient CD4+ T Cells",lollipop_color = "#F0E442")
    
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD4_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_CD4_T_cells_meta_gsea_highlight_lolipop.jpg",
      plot     = try_lollipop$p_pathway_lollipop,
      width    = 26,
      height   = 20,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",         # <– forces Cairo backend
      quality  = 100
    )
  }
  
  #CD8_T_cells
  {
    CD8_T_full_gsea <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/full_gsea_ora_results/neuroimmunological-healthy/GSEA/CD59.Def_Mevorach_CD59_Patient_full_meta_gsea_df.qs")
    
    
    list_of_highlight_pathways <- c(
      "GO:0042743",
      "HALLMARK_G2M_CHECKPOINT",
      "HALLMARK_E2F_TARGETS",
      "HALLMARK_HYPOXIA",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "GO:0042775",
      "HALLMARK_MTORC1_SIGNALING",
      "GO:0042101",
      "GO:0002250")
    
    # df <- try_pathways$meta_gsea_top
    
    manuall_changing_names <- list(
      HALLMARK_E2F_TARGETS = "E2F Targets",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
      "HALLMARK_G2M_CHECKPOINT" = "G2M Checkpoint",
      "mitochondrial ATP synthesis coupled electron transport" ="Mitochondrial ATP Synthesis Coupled Electron Transport",
      "HALLMARK_MTORC1_SIGNALING" ="mTORC1 Signaling",
      HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB")
    
    cell_type ="CD8_T_cells"
    condition = "CD59 def."
    sample_name = "CD59 patient"
    dataset_name = "Mevorach"
    
    
    try_lollipop <- get_pathway_lollipop(CD8_T_full_gsea,list_of_highlight_pathways,manuall_changing_names =manuall_changing_names,
                                         title = "PDE-Based GSEA: Key Pathways Enriched in CD59-Deficient CD8+ T Cells",lollipop_color = "#F37707")
    
    
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD8_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_CD8_T_cells_meta_gsea_highlight_lolipop.jpg",
      plot     = try_lollipop$p_pathway_lollipop,
      width    = 26,
      height   = 20,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",         # <– forces Cairo backend
      quality  = 100
    )
  }
  
}

#dimplot and dotplot
{
  font_sizes <- list(
    title        = 22,   # plot titles
    legend_title = 16,   # legend headers
    legend_text  = 14,   # legend entries
    axis_text    = 14    # axis tick & gene / cell-type labels
  )
  
  #dimplot
  {
    
    options(future.globals.maxSize = 3200 * 1024 * 1024 * 1024)  # 320 GB
    
    seu_cd59 <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/seu_cd59_only_patient_full_annotation.qs")
    
    dimplot <- DimPlot_scCustom(
      seurat_object = seu_cd59,
      group.by      = "final_cell_type",
      figure_plot   = TRUE,
      colors_use    = c(
        "NK cells"      = "#B62055",
        "Monocytes"     = "#56B4E9",
        "B cells"       = "#009E73",
        "CD4+ T cells"  = "#F0E442",
        "CD8+ T cells"  = "#F37707",
        "Other"         = "#999999"
      )
    )
    

    dimplot <- dimplot & labs(title = NULL)
    
    # add your own, centred, *not* bold
    dimplot <- dimplot +                                   # keep your patchwork chain
      plot_annotation(title = "Cell-Type Annotation of CD59-Deficient PBMCs") &
      theme(
        plot.title    = element_text(size = font_sizes$title,
                                     hjust = 0.5,
                                     face  = "bold"),
        legend.title  = element_text(size = font_sizes$legend_title),
        legend.text   = element_text(size = font_sizes$legend_text)
      )
    
  }
  
  #dotplot
  {
    
    DefaultAssay(object = seu_cd59) <- "SCT"
    Idents(seu_cd59) <- seu_cd59$final_cell_type
    
    
    solarExtra = c("5"='#3361A5', "7"='#248AF3', "1"='#14B3FF', "8"='#88CEEF', "9"='#C1D5DC', "4"='#EAD397', "3"='#FDB31A',"2"= '#E42A2A', "6"='#A31D1D')
    
    
    desired_order <- c("B cells", "CD4+ T cells", "CD8+ T cells",
                       "Monocytes", "NK cells", "Other")
    
    seu_cd59$final_cell_type <- factor(seu_cd59$final_cell_type,
                                       levels = desired_order)   # sets row order
    
    ## 2) Make it the active identity (or keep using group.by)
    Idents(seu_cd59) <- "final_cell_type"
    
    ## 3) Draw DotPlot
    dotplot <- DotPlot(
      object   = seu_cd59,
      features = c("MS4A1","CD79A", "CD19",#B
                   "CCR7", "IL7R", "CD40LG",#CD4
                   "CD8B", "CD8A","ITGAE",#CD8
                   "LYZ", "S100A8", "CD14",#mon
                   "KLRD1", "KLRF1","KLRC1"), #NK
      dot.min  = 0.01
    ) +
      scale_color_gradientn(colours = solarExtra) +
      theme(axis.text.x = element_text(angle = 45)) +
      ## OPTIONAL: put the first level at the very top (ggplot puts it at bottom)
      scale_y_discrete(limits = rev(desired_order)) +
      labs(
        x = "Features",
        y = "Cell Type"
      )
    
    dotplot <- dotplot +                                   # continue your existing chain
      plot_annotation(title = "Selected Genes Expression in CD59-Deficient PBMCs") &
      theme(
        plot.title    = element_text(size = font_sizes$title,
                                     hjust = 0.5,
                                     face  = "bold"),
        legend.title  = element_text(size = font_sizes$legend_title),
        legend.text   = element_text(size = font_sizes$legend_text),
        axis.text.x   = element_text(size = font_sizes$axis_text, angle = 45, hjust = 1),
        axis.text.y   = element_text(size = font_sizes$axis_text)
      )
    
  }
  
  #combined plot
  {
    # tag them if you like
    dimplot_tagged <- dimplot + labs(tag = "A")
    dotplot_tagged <- dotplot + labs(tag = "B")
    
    # now *tile* them
    combined <- plot_grid(
      dimplot, dotplot,
      ncol   = 2,
      labels = c("A", "B"),
      rel_widths = c(1,1))
    
    ggsave(
      "/bigdata/bgutman/Collabs/Dror mevorach/plots/cell_annotation_dimplot_and_dotplot.jpg",
      combined,
      width  = 25, height = 13, units = "in", dpi = 600,
      device  = "jpeg", type = "cairo", quality = 100
    )
    
    # ggsave(
    #   "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
    #   combined,
    #   width  = 25, height = 13, units = "in", dpi = 600,
    #   device  = "jpeg", type = "cairo", quality = 100
    # )
  }
}

#CD59.Def_cd8_integrated_vs_meta_scatter_plot
{
  #regular version
  {
    meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs") 
    cd8_integration_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs") 
    de_table <- cd8_integration_list$de_table 
    
    # selected_genes <- merged %>%
    #   dplyr::filter(median_log2FC < -5, avg_log2FC > -5) %>%
    #   dplyr::pull(gene) 
    # highlight_genes <- c("MAP3K8","S100B", "LAT") 
    # highlight_genes <- c("ATF3", "GADD45G", "GSTM1", "RGS16", "TMEM2", "WASHC1","S100B") 
    # highlight_genes <- c("TMEM2", "WASHC1","S100B") 
    # highlight_genes <- c(
    #   "MAP3K8", "S100B", "LAT",
    #   "ATF3", "GADD45G", "GSTM1",
    #   "RGS16", "TMEM2", "WASHC1"
    # ) 
    # highlight_genes <- c("S100B", "GABARAP", "PPBP") 
    
    highlight_genes <- c("MPV17L", "CD8B2", "ID1", "NCF1") 
    de_table$gene <- rownames(de_table) 
    
    # try <- get_linear_regression(
    #   meta_table = nk_meta_table,
    #   de_table = de_table,
    #   highlight_genes = c("S100B", "JUNB", "CXCR4", "KLRC3", "LAT"),
    #   title_suffix = "yo bro"
    # ) 
    
    meta_table$avg_log2FC <- NULL 
    
    # Inner join on gene
    merged <- merge(meta_table, de_table, by = "gene") 
    
    # # # Flag genes that start with "AC[digits]", "AL[digits]", "AP[digits]", "MT", or "RP" 
    # merged$special <- ifelse(
    #   grepl("^(AC[0-9]+|AL[0-9]+|AP[0-9]+|MT|RP)", merged$gene),
    #   "Non-Coding, Mitochondrial and Ribosomal Genes",
    #   "Other Genes"
    # ) 
    
    # 1) Fetch gene -> biotype from Ensembl
    coding_map <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      
      syms <- unique(merged$gene)
      mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters = "hgnc_symbol",
        values = syms,
        mart = mart
      )
      
      # keep one row per symbol
      annot[!duplicated(annot$hgnc_symbol), , drop = FALSE]
    }, error = function(e) {
      warning("biomaRt annotation failed: ", conditionMessage(e))
      NULL
    })
    
    # 2) Join biotype onto merged
    if (!is.null(coding_map)) {
      merged <- suppressWarnings(
        merge(merged, coding_map,
              by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
      )
    } else {
      merged$gene_biotype <- NA_character_
    }
    
    # 3) Define buckets (smart rules)
    is_noncoding <- !is.na(merged$gene_biotype) & merged$gene_biotype != "protein_coding"
    is_mito <- grepl("^MT-", merged$gene)     # mitochondrial genome–encoded
    is_ribo <- grepl("^(RPS|RPL|MRPS|MRPL)", merged$gene) # ribosomal proteins
    
    # 4) Assign 'special'
    merged$special <- ifelse(
      is_noncoding | is_mito | is_ribo,
      "Non-Coding, Mitochondrial and Ribosomal Genes",
      "Other Genes"
    )
    
    selected_genes <- unique(
      merged$gene[
        merged$median_log2FC < -5 &
          merged$avg_log2FC > -5 &
          merged$special == "Other Genes"
      ]
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
      meta_table = meta_table,
      de_table = de_table,
      meta_effect_size_column = "median_log2FC",
      int_effect_size_column = "avg_log2FC"
    )
    
    lm_fit <- lm(avg_log2FC ~ median_log2FC, data = merged)
    
    # label_data <- label_data %>%
    #   mutate(
    #     y_hat = predict(lm_fit, newdata = .),
    #     nudge_y = ifelse(avg_log2FC >= y_hat, 0.3, -0.3)
    #   )
    
    label_data <- label_data %>%
      mutate(
        y_hat   = predict(lm_fit, newdata = .),
        nudge_y = ifelse(avg_log2FC >= y_hat, 0.3, -0.3),
        # push CD8B2 and ID1 labels to at most y = -0.6 (below the x-axis)
        nudge_y = ifelse(gene %in% c("CD8B2","ID1") & (avg_log2FC + nudge_y) > -0.6,
                         -0.6 - avg_log2FC, nudge_y)
      )
    
    p <- ggplot(merged, aes(x = median_log2FC, y = avg_log2FC, fill = special)) +
      geom_point(shape = 21, colour = "black", size = 3, stroke = 0.5, alpha = 0.9) +
      
      ## —— highlight genes, now with per-point nudges ——
      geom_label_repel(
        data = label_data,
        aes(label = gene, nudge_y = nudge_y), # ← the key line
        colour = "white",
        fontface = "bold",
        size = 5,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey40",
        segment.size = 0.4,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      
      ## regression line (unchanged)
      geom_smooth(
        aes(group = 1), # ← guarantees a single fit
        method = "lm", se = FALSE,
        colour = "black", size = 0.8
      ) +
      
      geom_vline(xintercept = 0, colour = "darkgrey", size = 0.8) +
      geom_hline(yintercept = 0, colour = "darkgrey", size = 0.8) +
      
      annotate(
        "text", x = -Inf, y = Inf, hjust = 0, vjust = 1.1,
        label = paste0(
          "Linear Regression Slope: ", round(coef(lm_fit)[2], 2),
          "\nEffect Size Correlation (Pearson): ",
          round(effect_size_correlation$effect_size_pearson_correlation, 2)
        ),
        size = 6, fontface = "bold"
      ) +
      
      scale_fill_manual(
        values = c(
          "Non-Coding, Mitochondrial and Ribosomal Genes" = "#56B4E9",
          "Other Genes" = "#009E73"
        ),
        name = "Gene Group",
        guide = guide_legend(
          override.aes = list(shape = 21, colour = "black",
                              stroke = 0.5, size = 5, linetype = 0)
        )
      ) +
      
      labs(
        title = "PDE Median log2FC vs. Integration log2FC in CD59-Deficient CD8+ T Cells",
        x = "PDE Median log2FC",
        y = "Integration log2FC"
      ) +
      
      theme_minimal(base_size = 20) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(10, 20, 10, 10) # give labels a bit more air
      )
    
    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD8_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_CD8_integrated_vs_meta_scatter_plot.jpg",
    #   plot = p,
    #   width = 16, height = 13, units = "in",
    #   dpi = 600, device = "jpeg", type = "cairo", # <– forces Cairo backend
    #   quality = 100
    # )
    
    
  }
  
  #SD version
  {
    
    meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    cd8_integration_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    de_table <- cd8_integration_list$de_table
    
    
    highlight_genes <- c("MPV17L", "CD8B2", "ID1", "NCF1") 
    
    de_table$gene <- rownames(de_table) 
    # try <- get_linear_regression(meta_table = nk_meta_table, de_table = de_table, highlight_genes = c("S100B", "JUNB", "CXCR4", "KLRC3", "LAT"),
    # title_suffix = "yo bro") 
    
    meta_table$avg_log2FC <- NULL 
    
    # Inner join on gene 
    merged <- merge(meta_table, de_table, by = "gene") 
    
    # # Flag genes that start with "AC[digits]", "AL[digits]", "AP[digits]", "MT", or "RP" 
    # merged$special <- ifelse(
    #   grepl("^(AC[0-9]+|AL[0-9]+|AP[0-9]+|MT|RP)", merged$gene),
    #   "Non-Coding, Mitochondrial and Ribosomal Genes",
    #   "Other Genes"
    # ) 
    
    # 1) Fetch gene -> biotype from Ensembl 
    coding_map <- tryCatch({ 
      if (!requireNamespace("biomaRt", quietly = TRUE)) 
        stop("Package 'biomaRt' not installed.") 
      
      syms <- unique(merged$gene) 
      mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl") 
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters = "hgnc_symbol",
        values = syms,
        mart = mart
      ) 
      
      # keep one row per symbol 
      annot[!duplicated(annot$hgnc_symbol), , drop = FALSE] 
    }, error = function(e) { 
      warning("biomaRt annotation failed: ", conditionMessage(e)) 
      NULL 
    }) 
    
    # 2) Join biotype onto merged 
    if (!is.null(coding_map)) { 
      merged <- suppressWarnings(
        merge(merged, coding_map, by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
      ) 
    } else { 
      merged$gene_biotype <- NA_character_ 
    } 
    
    # 3) Define buckets (smart rules) 
    is_noncoding <- !is.na(merged$gene_biotype) & merged$gene_biotype != "protein_coding" 
    is_mito <- grepl("^MT-", merged$gene) # mitochondrial genome–encoded 
    is_ribo <- grepl("^(RPS|RPL|MRPS|MRPL)", merged$gene) # cytosolic + mitochondrial ribosomal proteins 
    
    # 4) Assign 'special' 
    merged$special <- ifelse(
      is_noncoding | is_mito | is_ribo,
      "Non-Coding, Mitochondrial and Ribosomal Genes",
      "Other Genes"
    ) 
    
    selected_genes <- unique(
      merged$gene[
        merged$median_log2FC < -5 &
          merged$avg_log2FC > -5 &
          merged$special == "Other Genes"
      ]
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
      meta_table = meta_table,
      de_table = de_table,
      meta_effect_size_column= "median_log2FC",
      int_effect_size_column = "avg_log2FC"
    ) 
    
    lm_fit <- lm(avg_log2FC ~ median_log2FC, data = merged) 
    
    # ---- SD bands (parallel to the regression line) ---- 
    # Residual SD (a.k.a. RMSE) of the model 
    resid_sd <- sigma(lm_fit) 
    
    # X grid to draw smooth lines across the plot range 
    xgrid <- data.frame(
      median_log2FC = seq(min(merged$median_log2FC, na.rm = TRUE),
                          max(merged$median_log2FC, na.rm = TRUE),
                          length.out = 200)
    ) 
    
    # Predicted line and ±1σ / ±2σ offsets 
    xgrid$yhat <- predict(lm_fit, newdata = xgrid) 
    band_df <- transform(
      xgrid,
      y_p1 = yhat + resid_sd, # +1 SD
      y_m1 = yhat - resid_sd, # -1 SD
      y_p2 = yhat + 2 * resid_sd, # +2 SD
      y_m2 = yhat - 2 * resid_sd # -2 SD
    ) 
    
    # label_data <- label_data %>%
    #   mutate(
    #     y_hat    = predict(lm_fit, newdata = .),
    #     nudge_y  = ifelse(avg_log2FC >= y_hat, 0.3, -0.3),
    #     nudge_y  = ifelse(gene %in% c("CD8B2","ID1") & (avg_log2FC + nudge_y) > -0.6,
    #                       -0.6 - avg_log2FC, nudge_y),
    #     y_label  = avg_log2FC + nudge_y
    #   )
    
    
    # label_data <- label_data %>%
    #   mutate(
    #     y_hat    = predict(lm_fit, newdata = .),
    #     nudge_y  = ifelse(avg_log2FC >= y_hat, 0.3, -0.3),
    #     nudge_y  = ifelse(gene %in% c("CD8B2","ID1") & (avg_log2FC + nudge_y) > -0.6,
    #                       -0.6 - avg_log2FC, nudge_y),
    #     y_label  = avg_log2FC + nudge_y,
    #     nudge_x  = ifelse(gene %in% c("CD8B2","ID1"), 0.2, 0)   # <— add this
    #   )
    
    label_data <- label_data %>%
      mutate(
        y_hat    = predict(lm_fit, newdata = .),
        nudge_y  = ifelse(avg_log2FC >= y_hat, 0.3, -0.3),
        nudge_y  = ifelse(gene %in% c("CD8B2","ID1") & (avg_log2FC + nudge_y) > -0.6,
                          -0.6 - avg_log2FC, nudge_y),
        # NEW: also push these labels below the -1σ band (with 0.25 margin)
        nudge_y  = ifelse(gene %in% c("CD8B2","ID1"),
                          pmin(nudge_y, (y_hat - resid_sd - 0.25) - avg_log2FC),
                          nudge_y),
        y_label  = avg_log2FC + nudge_y,
        nudge_x  = ifelse(gene %in% c("CD8B2","ID1"), 0.2, 0)
      )
    
    
    
    p <- ggplot(merged, aes(x = median_log2FC, y = avg_log2FC, fill = special)) + 
      geom_point(shape = 21, colour = "black", size = 3, stroke = 0.5, alpha = 0.9) + 
      ## —— highlight genes, now with per-point nudges —— 
      geom_label_repel(
        data = label_data,
        inherit.aes = FALSE,
        aes(x = median_log2FC + nudge_x, y = y_label, label = gene, fill = special),
        colour        = "white",
        fontface      = "bold",
        size          = 5,
        box.padding   = 0.5,
        point.padding = 0.6,     # a touch more clearance from points
        segment.color = "grey40",
        segment.size  = 0.4,
        max.overlaps  = Inf,
        show.legend   = FALSE
      ) +
      
      ## regression line (unchanged) 
      geom_smooth(
        aes(group = 1), # ← guarantees a single fit
        method = "lm",
        se = FALSE,
        colour = "black",
        size = 0.8
      ) + 
      # ---- add SD lines ---- 
    geom_line(data = band_df, aes(x = median_log2FC, y = y_p1),
              linetype = "dashed", size = 0.6, colour = "grey20",
              inherit.aes = FALSE) +
      geom_line(data = band_df, aes(x = median_log2FC, y = y_m1),
                linetype = "dashed", size = 0.6, colour = "grey20",
                inherit.aes = FALSE) +
      geom_line(data = band_df, aes(x = median_log2FC, y = y_p2),
                linetype = "dotted", size = 0.6, colour = "grey20",
                inherit.aes = FALSE) +
      geom_line(data = band_df, aes(x = median_log2FC, y = y_m2),
                linetype = "dotted", size = 0.6, colour = "grey20",
                inherit.aes = FALSE) + 
      geom_vline(xintercept = 0, colour = "darkgrey", size = 0.8) + 
      geom_hline(yintercept = 0, colour = "darkgrey", size = 0.8) + 
      annotate(
        "text", x = -Inf, y = Inf, hjust = 0, vjust = 1.1,
        label = paste0(
          "Linear Regression Slope: ", round(coef(lm_fit)[2], 2),
          "\nEffect Size Correlation (Pearson): ", 
          round(effect_size_correlation$effect_size_pearson_correlation, 2),
          "\nResidual SD (σ): ", round(resid_sd, 2)
        ),
        size = 6, fontface = "bold"
      ) + 
      scale_fill_manual(
        values = c(
          "Non-Coding, Mitochondrial and Ribosomal Genes" = "#56B4E9",
          "Other Genes" = "#009E73"
        ),
        name = "Gene Group",
        guide = guide_legend(
          override.aes = list(shape = 21, colour = "black",
                              stroke = 0.5, size = 5, linetype = 0)
        )
      ) + 
      labs(
        title = "PDE Median log2FC vs. Integration log2FC in CD59-Deficient CD8+ T Cells",
        x = "PDE Median log2FC",
        y = "Integration log2FC"
      ) + 
      theme_minimal(base_size = 20) + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(10, 20, 10, 10) # give labels a bit more air
      ) 
    
    
    
    ## ── Percent of "Non-Coding, Mitochondrial and Ribosomal Genes" by SD bands ──
    ## Paste this AFTER your existing code (it uses `merged`, `lm_fit`, and `resid_sd`)
    
    # 1) Compute residuals and assign SD-distance bands (absolute distance from LM)
    banded_df <- merged %>%
      dplyr::filter(is.finite(avg_log2FC), is.finite(median_log2FC)) %>%
      dplyr::mutate(
        y_hat     = predict(lm_fit, newdata = .),
        resid     = avg_log2FC - y_hat,
        abs_resid = abs(resid),
        band = cut(
          abs_resid,
          breaks  = c(-Inf, resid_sd, 2 * resid_sd, Inf),
          labels  = c("LM_to_SD1", "SD1_to_SD2", "SD2_plus"),
          right   = TRUE
        ),
        is_special = special == "Non-Coding, Mitochondrial and Ribosomal Genes",
        side       = ifelse(resid >= 0, "Above_LM", "Below_LM")
      )
    
    # 2) Summary: within each band, what % are special?
    n_all             <- nrow(banded_df)
    n_special_total   <- sum(banded_df$is_special, na.rm = TRUE)
    
    band_summary <- banded_df %>%
      dplyr::group_by(band) %>%
      dplyr::summarise(
        n_total   = dplyr::n(),
        n_special = sum(is_special, na.rm = TRUE),
        pct_special_within_band   = 100 * n_special / n_total,          # % of genes in this band that are special
        pct_of_all_genes_in_band  = 100 * n_total   / n_all,            # optional: % of ALL genes that lie in this band
        pct_of_all_specials_in_band = 100 * n_special / n_special_total,# optional: % of ALL SPECIAL genes that lie in this band
        .groups = "drop"
      ) %>%
      dplyr::arrange(factor(band, levels = c("LM_to_SD1","SD1_to_SD2","SD2_plus")))
    
    # 3) (Optional) Same breakdown but split by side of the line (above/below)
    band_by_side <- banded_df %>%
      dplyr::group_by(side, band) %>%
      dplyr::summarise(
        n_total   = dplyr::n(),
        n_special = sum(is_special, na.rm = TRUE),
        pct_special_within_band   = 100 * n_special / n_total,
        pct_of_all_genes_in_band  = 100 * n_total   / n_all,
        pct_of_all_specials_in_band = 100 * n_special / n_special_total,
        .groups = "drop"
      ) %>%
      dplyr::arrange(side, factor(band, levels = c("LM_to_SD1","SD1_to_SD2","SD2_plus")))
    
    # 4) Inspect results
    band_summary
    band_by_side
    
    
    # 1) Pull the % we want (one decimal) for each band
    vals <- band_summary |>
      dplyr::transmute(
        band = as.character(band),
        pct  = round(pct_special_within_band, 1)
      )
    
    # 2) Choose a right-side x position inside the plot area
    xr <- range(merged$median_log2FC, na.rm = TRUE)
    x_pos <- xr[2] - 0.03 * diff(xr)  # 3% in from the right border
    
    # Predicted y at that x; then midlines for the three gaps
    yhat_right <- predict(lm_fit, newdata = data.frame(median_log2FC = x_pos))
    label_pos <- tibble::tibble(
      band = c("LM_to_SD1", "SD1_to_SD2", "SD2_plus"),
      x    = x_pos,
      y    = c(
        yhat_right + 0.5 * resid_sd,   # between LM and +1σ
        yhat_right + 1.5 * resid_sd,   # between +1σ and +2σ
        yhat_right + 2.35 * resid_sd   # "≥ +2σ" space (slightly above +2σ)
      )
    ) |>
      dplyr::left_join(vals, by = "band") |>
      dplyr::mutate(label = sprintf("%.1f%%", pct))
    
    # 3) Add the labels to your existing plot `p`
    p <- p +
      geom_text(
        data = label_pos,
        inherit.aes = FALSE,
        aes(x = x, y = y, label = label),
        color = "#56B4E9", fontface = "bold", size = 6, hjust = 0,
        position = position_nudge(x = -4, y = -3.5)  # tweak these
      )
    
    
    p_sd <- p
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD8_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_CD8_integrated_vs_meta_scatter_SD_plot.jpg",
      plot = p,
      width = 16,
      height = 13,
      units = "in",
      dpi = 600,
      device = "jpeg",
      type = "cairo", # <– forces Cairo backend
      quality = 100
    )
    
    
  
    
  }
  
  #horizontal lines version
  {
    meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs") 
    cd8_integration_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs") 
    de_table <- cd8_integration_list$de_table 
    
    highlight_genes <- c("MPV17L", "CD8B2", "ID1", "NCF1") 
    de_table$gene <- rownames(de_table) 
    meta_table$avg_log2FC <- NULL 
    
    # Inner join on gene
    merged <- merge(meta_table, de_table, by = "gene") 
    
    # biotypes
    coding_map <- tryCatch({
      if (!requireNamespace("biomaRt", quietly = TRUE))
        stop("Package 'biomaRt' not installed.")
      syms <- unique(merged$gene)
      mart <- biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
      annot <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "gene_biotype"),
        filters = "hgnc_symbol",
        values = syms,
        mart = mart
      )
      annot[!duplicated(annot$hgnc_symbol), , drop = FALSE]
    }, error = function(e) {
      warning("biomaRt annotation failed: ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(coding_map)) {
      merged <- suppressWarnings(
        merge(merged, coding_map,
              by.x = "gene", by.y = "hgnc_symbol", all.x = TRUE)
      )
    } else {
      merged$gene_biotype <- NA_character_
    }
    
    is_noncoding <- !is.na(merged$gene_biotype) & merged$gene_biotype != "protein_coding"
    is_mito <- grepl("^MT-", merged$gene)
    is_ribo <- grepl("^(RPS|RPL|MRPS|MRPL)", merged$gene)
    
    merged$special <- ifelse(
      is_noncoding | is_mito | is_ribo,
      "Non-Coding, Mitochondrial and Ribosomal Genes",
      "Other Genes"
    )
    
    selected_genes <- unique(
      merged$gene[
        merged$median_log2FC < -5 &
          merged$avg_log2FC > -5 &
          merged$special == "Other Genes"
      ]
    )
    
    lm_model <- lm(avg_log2FC ~ median_log2FC, data = merged)
    linear_regression_slope <- coef(lm_model)["median_log2FC"]
    
    if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
      label_data <- subset(merged, gene %in% highlight_genes)
    } else {
      label_data <- merged[FALSE, ]
    }
    
    effect_size_correlation <- get_effect_size_correlation(
      meta_table = meta_table,
      de_table = de_table,
      meta_effect_size_column = "median_log2FC",
      int_effect_size_column = "avg_log2FC"
    )
    
    lm_fit <- lm(avg_log2FC ~ median_log2FC, data = merged)
    
    label_data <- label_data %>%
      mutate(
        y_hat   = predict(lm_fit, newdata = .),
        nudge_y = ifelse(avg_log2FC >= y_hat, 0.3, -0.3),
        nudge_y = ifelse(gene %in% c("CD8B2","ID1") & (avg_log2FC + nudge_y) > -0.6,
                         -0.6 - avg_log2FC, nudge_y)
      )
    
    p <- ggplot(merged, aes(x = median_log2FC, y = avg_log2FC, fill = special)) +
      geom_point(shape = 21, colour = "black", size = 3, stroke = 0.5, alpha = 0.9) +
      geom_label_repel(
        data = label_data,
        aes(label = gene, nudge_y = nudge_y),
        colour = "white", fontface = "bold", size = 5,
        box.padding = 0.5, point.padding = 0.3,
        segment.color = "grey40", segment.size = 0.4,
        max.overlaps = Inf, show.legend = FALSE
      ) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                  colour = "black", size = 0.8) +
      # Axes and the two dashed reference lines at x=4, y=4
      geom_vline(xintercept = 0, colour = "darkgrey", size = 0.8) +
      geom_hline(yintercept = 0, colour = "darkgrey", size = 0.8) +
      geom_vline(xintercept = 4, colour = "black", size = 0.8, linetype = "dashed") +
      geom_hline(yintercept = 4, colour = "black", size = 0.8, linetype = "dashed") +
      # Small “x=4” / “y=4” labels (keep or delete these two lines if you want only lines)
      { 
        xr <- range(merged$median_log2FC, na.rm = TRUE)
        yr <- range(merged$avg_log2FC,    na.rm = TRUE)
        dx <- diff(xr); dy <- diff(yr)
        list(
          coord_cartesian(clip = "off"),
          theme(plot.margin = margin(20, 28, 20, 20)),
          annotate("text", x = 4 + 0.02*dx, y = yr[2] - 0.03*dy,
                   label = "x = 4", size = 6, fontface = "bold", hjust = 0, vjust = 1),
          annotate("text", x = xr[2] - 0.03*dx, y = 4 + 0.02*dy,
                   label = "y = 4", size = 6, fontface = "bold", hjust = 1, vjust = 0)
        )
      } +
      annotate(
        "text", x = -Inf, y = Inf, hjust = 0, vjust = 1.1,
        label = paste0(
          "Linear Regression Slope: ", round(coef(lm_fit)[2], 2),
          "\nEffect Size Correlation (Pearson): ",
          round(effect_size_correlation$effect_size_pearson_correlation, 2)
        ),
        size = 6, fontface = "bold"
      ) +
      scale_fill_manual(
        values = c(
          "Non-Coding, Mitochondrial and Ribosomal Genes" = "#56B4E9",
          "Other Genes" = "#009E73"
        ),
        name = "Gene Group",
        guide = guide_legend(
          override.aes = list(shape = 21, colour = "black",
                              stroke = 0.5, size = 5, linetype = 0)
        )
      ) +
      labs(
        title = "PDE Median log2FC vs. Integration log2FC in CD59-Deficient CD8+ T Cells",
        x = "PDE Median log2FC",
        y = "Integration log2FC"
      ) +
      theme_minimal(base_size = 20) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(10, 20, 10, 10)
      )
    
    p_thres <- p
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD8_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_CD8_integrated_vs_meta_scatter_lines_plot.jpg",
      plot = p,
      width = 16, height = 13, units = "in",
      dpi = 600, device = "jpeg", type = "cairo",
      quality = 100
    )
  }
  
  
  #CD59.Def_CD8_integrated_vs_meta_scatter_combined.jpg
  {
    
    # we've created p_sd and p_thres in the  "#SD version" and the "#horizontal lines version" sections before
    
    # Combine plots: p_sd on the left, p_thres on the right
    # p_scatter <- (p_sd + p_thres) +
    #   plot_annotation(tag_levels = 'A') &          # automatic A, B, C… labels
    #   theme(plot.tag = element_text(face = "bold", size = 18))  # style of A/B labels
    # 
    
    p_scatter <- (p_sd / p_thres) +
      plot_annotation(tag_levels = 'A') &
      theme(plot.tag = element_text(face = "bold", size = 18))
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD8_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59.Def_CD8_integrated_vs_meta_scatter_sd_and_threshold.jpg",
      plot = p_scatter,
      width = 16, height = 26, units = "in",
      dpi = 600, device = "jpeg", type = "cairo",
      quality = 100
    )
    
  }
  


}


#Volcano, scatter and forset mega plot
{
  #volcano
  {
    # meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # nk_integration_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    # de_table <- nk_integration_list$de_table
    # 
    
    meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    cd8_integration_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/integrated_objects/neuroimmunulogical-healthy/CD59.Def_Mevorach_CD59_Patient_integrated.qs")
    de_table <- cd8_integration_list$de_table
    
    meta_table$avg_log2FC <- NULL
    de_table$gene <- rownames(de_table)
    
    # highlight_genes <- c("S100B", "MAP3K8", "LAT")
    # highlight_genes <- c("ATF3", "GADD45G", "GSTM1", "RGS16", "TMEM2", "WASHC1","S100B")
    # highlight_genes <- c("TMEM2", "WASHC1","S100B")
    highlight_genes <- c("MPV17L",
                         "CD8B2", "ID1", "NCF1")
    
    volcano_plot <- get_volcano(de_results = de_table,key_genes = highlight_genes, plottitle = "Integrated Volcano: CD59-Def. CD8+ T Cells", group_a = "CD59 Def.", group_b = "Healthy")
    
    
    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
    #   plot     = volcano_plot,
    #   width    = 16,
    #   height   = 13,
    #   units    = "in",
    #   dpi      = 600,
    #   device   = "jpeg",
    #   type     = "cairo",         # <– forces Cairo backend
    #   quality  = 100
    # )

    
  }
  
  #scatter
  {
    
    meta_table
    x_variable               = "dsi"
    y_variable               = "median_log2FC"
    upper_threshold          = 0.75
    below_threshold          = -0.75
    upper_y_line             = 1
    below_y_line             = -1
    full_x_variable_name     = "Differential Stability Index (DSI)"
    full_y_variable_name     = "Median log2 Fold Change"
    display_house_keeping_labels = TRUE
    # costume_de_genes         = c("S100B", "MAP3K8", "LAT")
    # costume_de_genes <- c("ATF3", "GADD45G", "GSTM1", "RGS16", "TMEM2", "WASHC1","S100B")
    # costume_de_genes <- c("TMEM2", "WASHC1","S100B")
    costume_de_genes <- c("MPV17L",
                         "CD8B2", "ID1", "NCF1")
    
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
        title = "PDE Plot: DSI vs. Median Log2FC of CD59-Def. CD8+ T Cells",
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
    
    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
    #   plot     = p_scatter_avg_log2fc,
    #   width    = 16,
    #   height   = 13,
    #   units    = "in",
    #   dpi      = 600,
    #   device   = "jpeg",
    #   type     = "cairo",         # <– forces Cairo backend
    #   quality  = 100
    # )
    
  }
  
  #forest
  {
    
    # highlight_genes <- c("S100B", "MAP3K8", "LAT")
    # highlight_genes <- c("S100B", "TMEM2", "WASHC1")
    highlight_genes <- c("MPV17L",
                         # "CD8B2",
                         "ID1", "NCF1")
    
    
    output_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/output_lists/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_outputlist.qs")
    
    
    get_3_forest_plots <- function(output_list, three_genes_vector) {
      # --------------------------- sanity checks --------------------------- #
      if (!is.vector(three_genes_vector) || length(three_genes_vector) != 3) {
        stop(sprintf(
          "three_genes_vector must contain exactly three genes – you supplied %d.",
          length(three_genes_vector)
        ))
      }
      
      # --------------------------- build plots ----------------------------- #
      tags <- paste0("C", 1:3)
      plot_list <- setNames(vector("list", 3), tags)
      
      for (i in seq_along(three_genes_vector)) {
        gene         <- three_genes_vector[i]
        forest_table <- get_table_for_forest(output_list, gene)
        forest_plot  <- create_forest_plot(forest_table, gene)
        
        plot_list[[ tags[i] ]] <- forest_plot
      }
      
      # --------------------------- tag panels ------------------------------ #
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
        plot_list, names(plot_list)
      )
      
      # --------------------------- assemble grid --------------------------- #
      # Stacked vertically: C1 / C2 / C3
      combined_plot <- plot_list[["C1"]] /
        plot_list[["C2"]] /
        plot_list[["C3"]]
      
      combined_plot
    }
    
    three_forset_plot <- get_3_forest_plots(output_list, c("MPV17L", "ID1", "NCF1"))
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
      plot     = three_forset_plot,
      width    = 20,
      height   = 15,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",         # <– forces Cairo backend
      quality  = 100
    )
  }
  
  
  #───────────────────────────────────────────────────────────────────────────────
  p_scatter_labeled <- p_scatter_avg_log2fc +
    annotate("text",
             x       = -Inf,  y =  Inf,        # top-left corner
             label   = "A",
             hjust   = -0.6,  vjust =  1.2,
             fontface = "bold", size = 6)
  
  volcano_labeled   <- volcano_plot +
    annotate("text",
             x       = -Inf,  y =  Inf,
             label   = "B",
             hjust   = -0.6,  vjust =  1.2,
             fontface = "bold", size = 6)
  
  #───────────────────────────────────────────────────────────────────────────────
  # 2. build layout  -------------------------------------------------------------
  #───────────────────────────────────────────────────────────────────────────────
  right_column <- p_scatter_labeled /
    volcano_labeled +
    plot_layout(heights = c(1, 1))      # stack 50 % / 50 %
  
  final_plot <- right_column|                  # left half
    three_forset_plot +                       # right half
    plot_layout(widths = c(1, 1))        # 50 % / 50 %
  
  
  # ggsave(
  #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/NKs/CD59.Def_Mevorach_CD59_Patient/CD59_def_Nk_volcano_scatter_and_forset_mega_plot.jpg",
  #   plot     = final_plot,
  #   width    = 25,
  #   height   = 20,
  #   units    = "in",
  #   dpi      = 600,
  #   device   = "jpeg",
  #   type     = "cairo",         # <– forces Cairo backend
  #   quality  = 100
  # )
  
  
  ggsave(
    filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/CD8_T_cells/CD59.Def_Mevorach_CD59_Patient/CD59_def_CD8_volcano_scatter_and_forset_mega_plot.jpg",
    plot     = final_plot,
    width    = 25,
    height   = 20,
    units    = "in",
    dpi      = 600,
    device   = "jpeg",
    type     = "cairo",         # <– forces Cairo backend
    quality  = 100
  )
  
}

# forset plot
{
  highlight_genes <- c("S100B", "TMEM2", "WASHC1")
  
  
  output_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/output_lists/neuroimmunological-healthy/cd59_output_list.qs")
  
  create_forest_plot <- function(forest_table, gene,
                                 panel_ratio = 1) {   # keep if you use it
    forest_table <- forest_table[order(forest_table$avg_log2FC), ]
    forest_table$sample <- factor(forest_table$sample, levels = forest_table$sample)
    
    ggplot(forest_table, aes(y = avg_log2FC, x = sample)) +
      geom_point(aes(size = pct_2, color = dataset), shape = 16, alpha = 0.8) +
      geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
      
      # ­­scales ---------------------------------------------------------------
    scale_size_continuous(range = c(3, 8)) +
      scale_y_continuous(breaks = scales::pretty_breaks(10)) +
      scale_color_manual(values = c("#ffc857", "#e9724c", "#c5283d",
                                    "#481d24", "#976391", "#255f85",
                                    "#2a9d8f", "#283618", "#606c38")) +
      
      # legends ---------------------------------------------------------------
    guides(
      size  = guide_legend(
        title = "pct.2",
        order = 1),                       # ← FIRST (on top)
      color = guide_legend(
        title = "Dataset",
        order = 2,                        # ← SECOND
        override.aes = list(size = 7, shape = 16))
    ) +
      
      # labels ----------------------------------------------------------------
    ylab("Log2 Fold Change") +
      ggtitle(sprintf("%d Comparisons of <b>%s</b> Expression (Pct.1 = %.3f)",
                      nrow(forest_table), gene, unique(forest_table$pct_1)[1])) +
      
      # theme -----------------------------------------------------------------
    theme_minimal() +
      theme(
        aspect.ratio   = panel_ratio,
        plot.title     = element_markdown(size = 24, hjust = 0.5),
        axis.title.y   = element_text(size = 18),
        axis.title.x   = element_text(size = 18),
        axis.text.x    = element_text(size = 18, angle = 45, hjust = 1),
        axis.text.y    = element_text(size = 18),
        legend.title   = element_text(size = 20, face = "bold"),
        legend.text    = element_text(size = 20),
        legend.key.size = grid::unit(1.1, "cm")
      )
  }
  
  
  get_3_forest_plots <- function(output_list,
                                 three_genes_vector,
                                 panel_ratio = 1,
                                 tag_size    = 10   # ← change this to control C-label size
  ) {
    if (!is.vector(three_genes_vector) || length(three_genes_vector) != 3)
      stop("three_genes_vector must contain exactly three genes.")
    
    tags <- paste0("C", 1:3)
    plot_list <- vector("list", 3)
    names(plot_list) <- tags
    
    for (i in seq_along(three_genes_vector)) {
      gene         <- three_genes_vector[i]
      forest_table <- get_table_for_forest(output_list, gene)
      forest_plot  <- create_forest_plot(forest_table, gene,
                                         panel_ratio = panel_ratio)
      
      # add corner tag (C1, C2, C3 …) with bigger text
      plot_list[[ tags[i] ]] <- forest_plot +
        ggplot2::annotate(
          "text",
          x = -Inf, y = Inf,
          label    = tags[i],
          hjust    = -0.5, vjust = 1.5,
          size     = tag_size,    # ← bigger
          fontface = "bold"
        )
    }
    
    # stack vertically and KEEP every legend
    (plot_list[["C1"]] /
        plot_list[["C2"]] /
        plot_list[["C3"]]) +
      patchwork::plot_layout(guides = "keep")
  }
  
  
  three_forest_plot <- get_3_forest_plots(
    output_list,
    c("S100B", "TMEM2", "WASHC1"),
    panel_ratio = 1  # tweak if you want them even taller/narrower
  )
  
  overall_width  <- 14   # inches – adjust for your slide
  overall_height <- 24   # 3 × panels → 3 × 8" = 24"
  
  ggsave(
    filename = "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
    plot     = three_forest_plot,
    width    = overall_width,
    height   = overall_height,
    units    = "in",
    dpi      = 600,
    device   = "jpeg",
    type     = "cairo",
    quality  = 100
  )
  
}


  
#Upset Plot DE genes
{
  #create positive_gene_upset_df and negative_gene_upset_df
  {
    # 
    # NK_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # monocytes_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # B_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # CD4_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # CD8_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # 
    # 
    # meta_table_list <- list(NK_meta_table = NK_meta_table,
    #                         monocytes_meta_table = monocytes_meta_table,
    #                         B_cells_meta_table = B_cells_meta_table,
    #                         CD4_T_cells_meta_table = CD4_T_cells_meta_table,
    #                         CD8_T_cells_meta_table = CD8_T_cells_meta_table)
    # 
    # 
    # rm(NK_meta_table, monocytes_meta_table, B_cells_meta_table, CD4_T_cells_meta_table, CD8_T_cells_meta_table)
    # 
    # 
    # ## ---- 1.  “Positive” genes  --------------------------------------------------
    # meta_table_list_positive <- lapply(
    #   meta_table_list,
    #   function(df) df %>%                         # if dplyr is loaded, %>% works; otherwise remove it
    #     filter(median_log2FC >  1.25, q_val_fisher <= 0.05)
    # )
    # 
    # ## ---- 2.  “Negative” genes  --------------------------------------------------
    # meta_table_list_negative <- lapply(
    #   meta_table_list,
    #   function(df) df %>%
    #     filter(median_log2FC < -1.25, q_val_fisher <= 0.05)
    # )
    # 
    # ## (Optional) make sure the new lists carry the same element names
    # names(meta_table_list_positive) <- names(meta_table_list)
    # names(meta_table_list_negative) <- names(meta_table_list)
    # 
    # 
    # positive_gene_upset_df <- build_presence_df(
    #   meta_table_list_positive,   # the list you built earlier
    #   id_col    = "gene",          # auto-detect would also work
    #   filter_nocoding = TRUE)
    # 
    # 
    # 
    # 
    # negative_gene_upset_df <- build_presence_df(
    #   meta_table_list_negative,
    #   id_col = "gene",
    #   filter_nocoding = TRUE)
    # 
    
    
  }
  
  #add pathway group annotations using REACTOME
  {

    
    # param <- MulticoreParam(workers = 6)   # or SnowParam on Windows
    # positive_gene_upset_df <- add_pathway_group_reactome_and_goslim (positive_gene_upset_df,  BPPARAM = param)
    # negative_gene_upset_df <- add_pathway_group_reactome_and_goslim (negative_gene_upset_df,  BPPARAM = param)
    # 
    # positive_gene_upset_df$group_reactome_copy <- positive_gene_upset_df$group_reactome
    # negative_gene_upset_df$group_reactome_copy <- negative_gene_upset_df$group_reactome
    # 
    # qsave(positive_gene_upset_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
    # qsave(negative_gene_upset_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")
    # 
    
    positive_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
    negative_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")
    
    # 
    # positive_gene_upset_df$group_reactome <- positive_gene_upset_df$group_reactome_copy
    # negative_gene_upset_df$group_reactome <- negative_gene_upset_df$group_reactome_copy
    # 
    # 
    # positive_gene_upset_df <- positive_gene_upset_df %>%
    #   mutate(
    #     group_reactome = case_when(
    #       group_reactome %in% c("Metabolism", "Immune System") ~ group_reactome,
    #       group_reactome %in% c("Metabolism of proteins", "Metabolism of RNA") ~
    #         "RNA/Proteins Metabolism",
    #       TRUE ~ "Other/No Pathway"
    #     )
    #   )
    # 
    # negative_gene_upset_df <- negative_gene_upset_df %>%
    #   mutate(
    #     group_reactome = case_when(
    #       group_reactome %in% c("Metabolism", "Immune System") ~ group_reactome,
    #       group_reactome %in% c("Metabolism of proteins", "Metabolism of RNA") ~ "RNA/Proteins Metabolism",
    #       TRUE ~ "Other/No Pathway"
    #     )
    #   )
    

  
  }
  
  
  #Old and Irrelevant parts
  {
    #OLD: Grouping the results
    {
      
      # 
      # positive_gene_upset_df <- positive_gene_upset_df %>%
      #   mutate(group_go = case_when(
      #     group_go %in% c(
      #       "hydrolase activity",
      #       "catalytic activity",
      #       "sulfur compound metabolic process",
      #       "transporter activity",
      #       "vesicle-mediated transport",
      #       "tRNA metabolic process",
      #       "protein catabolic process",
      #       "amino acid metabolic process",
      #       "mRNA metabolic process",
      #       "generation of precursor metabolites and energy",
      #       "ribosome biogenesis",
      #       "lipid binding",
      #       "lipid metabolic process",
      #       "cellular modified amino acid metabolic process",
      #       "mitochondrial gene expression",
      #       "nucleobase-containing small molecule metabolic process",
      #       "carbohydrate metabolic process",
      #       "vitamin metabolic process",
      #       "cytoplasmic translation",
      #       "protein folding",
      #       "carbohydrate derivative metabolic process",
      #       "snRNA metabolic process",
      #       "protein glycosylation")  ~ "Metabolism",
      #     group_go %in% c("Immune Processes", "cell differentiation", "cytoskeleton organization", "immune system process",
      #                     "defense response to other organism", "cell motility", "inflammatory response","wound healing")  ~ "Immune Processes",
      #     group_go %in% c(
      #       "programmed cell death",
      #       "DNA repair",
      #       "autophagy",
      #       "mitochondrion organization",
      #       "DNA binding",
      #       "DNA-templated transcription",
      #       "Transcription regulator activity",
      #       "Chromatin organization",
      #       "mitotic cell cycle",
      #       "DNA replication",
      #       "telomere organization",
      #       "nucleocytoplasmic transport",
      #       "transcription regulator activity",
      #       "regulation of DNA-templated transcription",
      #       "regulatory ncRNA-mediated gene silencing",
      #       "chromosome segregation ",
      #       "chromatin organization"
      #     ) ~ "Stress and Genome Regulation",
      #     group_go %in% c("nervous system process", "reproductive process",
      #                     "structural molecule activity","intracellular protein transport",
      #                     "anatomical structure development", "signaling", "RNA binding",
      #                     "molecular function regulator activity",
      #                     "cytoskeletal protein binding",
      #                     "GTPase activity",
      #                     "membrane organization",
      #                     "molecular transducer activity",
      #                     "muscle system process",
      #                     "translation regulator activity",
      #                     "cell adhesion",
      #                     "transmembrane transport",
      #                     "circulatory system process",
      #                     "cytokinesis",
      #                     "molecular adaptor activity",
      #                     "protein localization to plasma membrane")  ~ "Other / No Pathway",
      #     
      #     # TRUE ~ negative_gene_upset_df$group_go       # keep any pre-existing label
      #     TRUE ~ "Other / No Pathway"      # label as "other"
      #     
      #   ))
      # 
      # 
      # 
      # 
      # negative_gene_upset_df <- negative_gene_upset_df %>%
      #   mutate(group_go = case_when(
      #     group_go %in% c(
      #       "hydrolase activity",
      #       "catalytic activity",
      #       "sulfur compound metabolic process",
      #       "transporter activity",
      #       "vesicle-mediated transport",
      #       "tRNA metabolic process",
      #       "protein catabolic process",
      #       "amino acid metabolic process",
      #       "mRNA metabolic process",
      #       "generation of precursor metabolites and energy",
      #       "ribosome biogenesis",
      #       "lipid binding",
      #       "lipid metabolic process",
      #       "cellular modified amino acid metabolic process",
      #       "mitochondrial gene expression",
      #       "nucleobase-containing small molecule metabolic process",
      #       "carbohydrate metabolic process",
      #       "vitamin metabolic process",
      #       "cytoplasmic translation",
      #       "protein folding",
      #       "carbohydrate derivative metabolic process",
      #       "snRNA metabolic process",
      #       "protein glycosylation")  ~ "Metabolism",
      #     group_go %in% c("Immune Processes", "cell differentiation", "cytoskeleton organization", "immune system process",
      #                     "defense response to other organism", "cell motility", "inflammatory response","wound healing")  ~ "Immune Processes",
      #     group_go %in% c(
      #       "programmed cell death",
      #       "DNA repair",
      #       "autophagy",
      #       "mitochondrion organization",
      #       "DNA binding",
      #       "DNA-templated transcription",
      #       "Transcription regulator activity",
      #       "Chromatin organization",
      #       "mitotic cell cycle",
      #       "DNA replication",
      #       "telomere organization",
      #       "nucleocytoplasmic transport",
      #       "transcription regulator activity",
      #       "regulation of DNA-templated transcription",
      #       "regulatory ncRNA-mediated gene silencing",
      #       "chromosome segregation ",
      #       "chromatin organization"
      #     ) ~ "Stress and Genome Regulation",
      #     group_go %in% c("nervous system process", "reproductive process",
      #                     "structural molecule activity","intracellular protein transport",
      #                     "anatomical structure development", "signaling", "RNA binding",
      #                     "molecular function regulator activity",
      #                     "cytoskeletal protein binding",
      #                     "GTPase activity",
      #                     "membrane organization",
      #                     "molecular transducer activity",
      #                     "muscle system process",
      #                     "translation regulator activity",
      #                     "cell adhesion",
      #                     "transmembrane transport",
      #                     "circulatory system process",
      #                     "cytokinesis",
      #                     "molecular adaptor activity",
      #                     "protein localization to plasma membrane")  ~ "Other / No Pathway",
      #     
      #     # TRUE ~ negative_gene_upset_df$group_go       # keep any pre-existing label
      #     TRUE ~ "Other / No Pathway"      # label as "other"
      #     
      #   ))
      
    }
    
    
    
    #OLD: experimenting with always_true_df
    {
      
      #ORA (for annotating always_true_dfs later)
      {
        
        #positive ORA
        {
          # 
          # cols <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
          # 
          # 
          # positive_ora <- get_ora(meta_table_list, log2FC_threshold_for_gene = 1.25, positive = TRUE ,q_value_threshold_for_gene = 0.05, cols, q_value_threshold_for_pathway = 0.1)
          # 
          # positive_ora_df <- get_full_gsea_and_ora_dfs(ora_return_list  = positive_ora,
          #                                              save_dir_ora     = "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/",
          #                                              condition = "CD59.Def", dataset_name = "Mevorach", sample_name = "CD59_Patient")
          # 
          # excluded_ids <- c(
          #   "R-HSA-192823",
          #   "R-HSA-975956",
          #   "R-HSA-975957",
          #   "R-HSA-9735869",
          #   "R-HSA-9754678"
          # )
          # 
          # positive_ora_df <- positive_ora_df$meta_ora_top %>%
          #   filter(! ID %in% excluded_ids)
          
        }
        
        
        #negative ORA
        {                                                                                                                                                                                                                         
          # cols <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
          # 
          # 
          # negative_ora <- get_ora(meta_table_list, log2FC_threshold_for_gene = 1.25, positive = FALSE ,q_value_threshold_for_gene = 0.05, cols)
          # 
          # negative_ora_df <- get_full_gsea_and_ora_dfs(ora_return_list  = negative_ora)
          # 
          # negative_ora_df <- negative_ora_df$meta_ora_top
          # 
          # excluded_ids <- c(
          #   "R-HSA-168255",
          #   "R-HSA-192823",
          #   "hsa05134",
          #   "hsa05171",
          #   "GO:0032496",
          #   "GO:0002237",
          #   "GO:0071222",
          #   "GO:0071219",
          #   "GO:0071216",
          #   "GO:0045444",
          #   "GO:0030217",
          #   "GO:0045580",
          #   "GO:0043367",
          #   "GO:0035710",
          #   "GO:0046632",
          #   "GO:1904019",
          #   "GO:0043370",
          #   "GO:0042093",
          #   "GO:0046637",
          #   "GO:0002294",
          #   "GO:0002287",
          #   "GO:0002293",
          #   "GO:2000514",
          #   "GO:0002292",
          #   "GO:2000319",
          #   "GO:2000316",
          #   "GO:0072539",
          #   "GO:0045622",
          #   "GO:0072538",
          #   "GO:0050673",
          #   "GO:0050678",
          #   "R-HSA-975956",
          #   "GO:0019915")
          # 
          # 
          # negative_ora_df <- negative_ora_df %>%
          #   filter(! ID %in% excluded_ids)
          
        }
        
      }
      
      # positive_gene_upset_df <- positive_gene_upset_df %>%
      #   mutate(group_reactome = case_when(
      #     group_reactome %in% c("Metabolism of proteins", "Metabolism of RNA")  ~ "Metabolism",
      #     group_reactome %in% c("Disease", "Immune System")  ~ "Immunopathology",
      #     group_reactome %in% c("DNA Repair", "Cell Cycle","Organelle biogenesis and maintenance",
      #                           "Autophagy", "Programmed Cell Death", "Chromatin organization")  ~ "Stress Response & Quality Control",
      #     TRUE ~ "Other / No Pathway"       # keep any pre-existing label
      #   ))
      # 
      
      
      
      
      # always_true_df positive
      
      {
        
        # columns_to_remove <- c(
        #   "NKs",
        #   "Monocytes",
        #   "B Cells",
        #   "CD4+ T Cells",
        #   "CD8+ T Cells"
        # )
        # 
        # always_true_positive <- positive_gene_upset_df %>%
        #   filter(if_all(c(NKs, Monocytes, `B Cells`, `CD4+ T Cells`, `CD8+ T Cells`), ~ .)) %>%
        #   dplyr:: select(-all_of(columns_to_remove))
        # 
        # cols <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
        # 
        # always_true_in_atleast_4 <- positive_gene_upset_df %>%
        #   filter(rowSums(across(all_of(cols), ~ .x)) >3 ) %>%   # exactly 4 TRUEs
        #   dplyr::select(-all_of(columns_to_remove))
        # always_true_in_atleast_4 %>%
        #   count(group_reactome, sort = TRUE)
        # 
        # sort(table(always_true_in_atleast_4$group_reactome), decreasing = TRUE)
        # 
        # 
        # always_negative_in_atleast_4 <- negative_gene_upset_df %>%
        #   filter(rowSums(across(all_of(cols), ~ .x)) >3 ) %>%   # exactly 4 TRUEs
        #   dplyr::select(-all_of(columns_to_remove))
        # 
        # 
        # sort(table(always_negative_in_atleast_4$group_reactome), decreasing = TRUE)
        # 
        # 
        # 
        # 
        # # 1. Bind the needed columns from all meta tables into one long table
        # meta_long <- meta_table_list %>%
        #   imap_dfr(function(tbl, mt_name) {
        #     # ensure gene column exists
        #     if (!"gene" %in% colnames(tbl)) {
        #       tbl <- tibble::rownames_to_column(tbl, "gene")
        #     }
        #     tbl %>%
        #       dplyr::select(gene,
        #                     median_log2FC,
        #                     q_val_fisher,
        #                     dsi,
        #                     comparisons_num) %>%
        #       mutate(meta_source = mt_name)
        #   })
        # 
        # # 2. For each gene, compute minima across the 5 meta tables
        # meta_mins <- meta_long %>%
        #   group_by(gene) %>%
        #   summarize(
        #     min_median_log2FC = min(median_log2FC, na.rm = TRUE),
        #     min_q_val_fisher  = min(q_val_fisher,  na.rm = TRUE),
        #     min_dsi           = min(dsi,           na.rm = TRUE),
        #     min_comparisons_num = min(comparisons_num, na.rm = TRUE),
        #     .groups = "drop"
        #   )
        # 
        # # 3. Join back onto always_true_df
        # always_true_positive <- always_true_positive %>%
        #   left_join(meta_mins, by = "gene")
        # 
        # 
        # always_true_positive <- annotate_with_positive_ORA(
        #   always_true_positive,
        #   positive_ora_df
        # )
        # 
        # 
        # always_true_positive_important_genes <- always_true_positive %>%
        #   filter(positive_ORA == TRUE)
        # 
        # sort(table(always_true_positive$group_go), decreasing = TRUE)
        
      }
      
      # always_true_df negative
      {
        
        #  columns_to_remove <- c(
        #    "NKs",
        #    "Monocytes",
        #    "B Cells",
        #    "CD4+ T Cells",
        #    "CD8+ T Cells"
        #  )
        #  
        #  always_true_negative <- negative_gene_upset_df %>% 
        #    filter(if_all(c(NKs, Monocytes, `B Cells`, `CD4+ T Cells`, `CD8+ T Cells`), ~ .)) %>% 
        #    dplyr:: select(-all_of(columns_to_remove))
        #  
        #  
        #  
        #  # 1. Bind the needed columns from all meta tables into one long table
        #  meta_long <- meta_table_list %>%
        #    imap_dfr(function(tbl, mt_name) {
        #      # ensure gene column exists
        #      if (!"gene" %in% colnames(tbl)) {
        #        tbl <- tibble::rownames_to_column(tbl, "gene")
        #      }
        #      tbl %>%
        #        dplyr::select(gene,
        #                      median_log2FC,
        #                      q_val_fisher,
        #                      dsi,
        #                      comparisons_num) %>%
        #        mutate(meta_source = mt_name)
        #    })
        #  
        #  # 2. For each gene, compute minima across the 5 meta tables
        #  meta_mins <- meta_long %>%
        #    group_by(gene) %>%
        #    summarize(
        #      min_median_log2FC = min(median_log2FC, na.rm = TRUE),
        #      min_q_val_fisher  = min(q_val_fisher,  na.rm = TRUE),
        #      min_dsi           = min(dsi,           na.rm = TRUE),
        #      min_comparisons_num = min(comparisons_num, na.rm = TRUE),
        #      .groups = "drop"
        #    )
        #  
        #  # 3. Join back onto always_true_df
        #  always_true_negative <- always_true_negative %>%
        #    left_join(meta_mins, by = "gene")
        #  
        #  four_of_five <- negative_gene_upset_df %>%
        #    dplyr::filter({
        #      n_true <- rowSums(dplyr::across(c(NKs, Monocytes, `B Cells`, `CD4+ T Cells`, `CD8+ T Cells`)))
        #      n_true >= 4 & n_true < 5
        #    }) %>%
        #    dplyr::select(-all_of(columns_to_remove))
        #  
        #  
        # sort(table(always_true_negative$group_go), decreasing = TRUE)
        # sort(table(four_of_five$group_go), decreasing = TRUE)
        
        
      }
      
      
    }
    
  }
  
  #upset plots themselves
  {
    #positive_gene_upset
    {
      positive_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
      
      positive_gene_upset <- get_upset_plot(
        df = positive_gene_upset_df,
        plot_title = "Overlap of Positive DE Protein Coding Genes between Cell Types",
        sets           = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"),
        group_column    = "group_reactome",        
        group_palette  =  c(
          "Metabolism"                        = "#219ebc",
          "RNA/Proteins Metabolism"       = "#ffb703",
          "Immune System"                     = "#fb8500",
          "Other/No Pathway"                = "#023047"),
        # legend_order   = c("Metabolism / Cellular Respiration","Oxidative-Stress", "Proliferation",
        #                    "Replication-Stress", "p53-Senescence", "Migration / Chemotaxis", "Other / No Pathway"),
        mode       = "inclusive_intersection",
        legend_title = "Reactome Based Groups",
        big_title  = 32,
        big_axis   = 30,
        big_txt    = 22,
        big_legend = 30
      )
      
      # ggsave(
      #   filename =       "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/upsets/positive_gene_upset_plot_reactome.jpg",
      #   plot     = positive_gene_upset,
      #   width    = 34,
      #   height   = 24,
      #   units    = "in",
      #   dpi      = 600
      # )


    }
    
    #negative_gene_upset
    {
      negative_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")
      
      negative_gene_upset <- get_upset_plot(
        df =  negative_gene_upset_df,
        plot_title = "Overlap of Negative DE Protein Coding Genes between Cell Types",
        sets           = c("B Cells","NKs", "Monocytes",  "CD4+ T Cells", "CD8+ T Cells"),
        group_column    = "group_reactome",        
        group_palette  =  c(
          "Metabolism"                        = "#219ebc",
          "RNA/Proteins Metabolism"       = "#ffb703",
          "Immune System"                     = "#fb8500",
          "Other/No Pathway"                = "#023047"),
        # legend_order   = c("Pro Inflammatory", "Stress Adaptation", "Migration / Chemotaxis",
        #                    "Pro Inhibitory", "Proliferation", "Other / No Pathway"),
        mode       = "inclusive_intersection",
        legend_title = "Reactome Based Groups",
        big_title  = 32,
        big_axis   = 30,
        big_txt    = 22,
        big_legend = 30
      )
      

      
      # ggsave(
      #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/upsets/negative_gene_upset_plot_reactome.jpg",
      #   plot     = negative_gene_upset,
      #   width    = 34,
      #   height   = 24,
      #   units    = "in",
      #   dpi      = 600
      # )
    }
    
    #combined upsetplot
    {
      # Combine plots: p_sd on the left, p_thres on the right
      p_whole_upset <- (positive_gene_upset / negative_gene_upset) +
        plot_annotation(tag_levels = 'A') &          # automatic A, B, C… labels
        theme(plot.tag = element_text(face = "bold", size = 45))  # style of A/B labels
      
      
      ggsave(
        filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/upsets/whole_gene_upset_plot_reactome.jpg",
        plot     = p_whole_upset,
        width    = 34,     # or any value >50
        height   = 60,
        units    = "in",
        dpi      = 400,
        limitsize = FALSE   # ⬅️ this bypasses the 50-inch limit
      )
      
 
    }
    
  }
  
}

#Significance Barplot of Overlap genes 
{
  
  #Calculate foldchange of independence
  {
    NK_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    monocytes_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    B_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    CD4_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    CD8_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    
    
    meta_table_list <- list(NK_meta_table = NK_meta_table,
                            monocytes_meta_table = monocytes_meta_table,
                            B_cells_meta_table = B_cells_meta_table,
                            CD4_T_cells_meta_table = CD4_T_cells_meta_table,
                            CD8_T_cells_meta_table = CD8_T_cells_meta_table)
    
    
    rm(NK_meta_table, monocytes_meta_table, B_cells_meta_table, CD4_T_cells_meta_table, CD8_T_cells_meta_table)
    
    overlap_genes_significance <- get_overlap_genes_significance(obj_list = meta_table_list)
    
    
    
  }
  
  #barplot_list 
  {
    barplot_list <- overlap_genes_significance$EF_list %>%
      map(~ tibble(
        name    = names(.x),
        log2ef  = log2(map_dbl(.x, "ef")),
        q_val   = map_dbl(.x, "q_val")
      ))
    

  }
  
  #call get_overlap_genes_significance_barplot
  {
    p_final <-get_overlap_genes_significance_barplot(barplot_list = barplot_list)
    
    
    # ggsave(
    #   filename =       "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/upsets/overlap_genes_significance_barplot.jpg",
    #   plot     = p_final,
    #   width    = 8,
    #   height   = 6,
    #   units    = "in",
    #   dpi      = 600
    # )
  }
  
 
  
}

#Composition plots
{
  #create positive_gene_upset_df and negative_gene_upset_df
  {
    # 
    # NK_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # monocytes_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # B_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # CD4_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # CD8_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    # 
    # 
    # meta_table_list <- list(NK_meta_table = NK_meta_table,
    #                         monocytes_meta_table = monocytes_meta_table,
    #                         B_cells_meta_table = B_cells_meta_table,
    #                         CD4_T_cells_meta_table = CD4_T_cells_meta_table,
    #                         CD8_T_cells_meta_table = CD8_T_cells_meta_table)
    # 
    # 
    # rm(NK_meta_table, monocytes_meta_table, B_cells_meta_table, CD4_T_cells_meta_table, CD8_T_cells_meta_table)
    # 
    # 
    # ## ---- 1.  “Positive” genes  --------------------------------------------------
    # meta_table_list_positive <- lapply(
    #   meta_table_list,
    #   function(df) df %>%                         # if dplyr is loaded, %>% works; otherwise remove it
    #     filter(median_log2FC >  1.25, q_val_fisher <= 0.05)
    # )
    # 
    # ## ---- 2.  “Negative” genes  --------------------------------------------------
    # meta_table_list_negative <- lapply(
    #   meta_table_list,
    #   function(df) df %>%
    #     filter(median_log2FC < -1.25, q_val_fisher <= 0.05)
    # )
    # 
    # ## (Optional) make sure the new lists carry the same element names
    # names(meta_table_list_positive) <- names(meta_table_list)
    # names(meta_table_list_negative) <- names(meta_table_list)
    # 
    # 
    # positive_gene_upset_df <- build_presence_df(
    #   meta_table_list_positive,   # the list you built earlier
    #   id_col    = "gene",          # auto-detect would also work
    #   filter_nocoding = TRUE)
    # 
    # 
    # 
    # 
    # negative_gene_upset_df <- build_presence_df(
    #   meta_table_list_negative,
    #   id_col = "gene",
    #   filter_nocoding = TRUE)
    
    
    
  }

  #add pathway group annotations using GOslim
  {
    positive_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
    negative_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")
    # 
    # param <- MulticoreParam(workers = 6)   # or SnowParam on Windows
    # positive_gene_upset_df <- add_pathway_group_reactome_and_goslim (positive_gene_upset_df,  BPPARAM = param)
    # negative_gene_upset_df <- add_pathway_group_reactome_and_goslim (negative_gene_upset_df,  BPPARAM = param)
    # 
    # positive_gene_upset_df$group_reactome_copy <- positive_gene_upset_df$group_reactome
    # negative_gene_upset_df$group_reactome_copy <- negative_gene_upset_df$group_reactome
    # 
    # qsave(positive_gene_upset_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
    # qsave(negative_gene_upset_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")
    # 
# 
#     positive_gene_upset_df$group_reactome <- positive_gene_upset_df$group_reactome_copy
#     negative_gene_upset_df$group_reactome <- negative_gene_upset_df$group_reactome_copy
#     
#     
#     positive_gene_upset_df <- positive_gene_upset_df %>%
#       mutate(
#         group_reactome = case_when(
#           group_reactome %in% c("Metabolism", "Immune System") ~ group_reactome,
#           group_reactome %in% c("Metabolism of proteins", "Metabolism of RNA") ~
#             "RNA/Proteins Metabolism",
#           TRUE ~ "Other/No Pathway"
#         )
#       )
#     
#     negative_gene_upset_df <- negative_gene_upset_df %>%
#       mutate(
#         group_reactome = case_when(
#           group_reactome %in% c("Metabolism", "Immune System") ~ group_reactome,
#           group_reactome %in% c("Metabolism of proteins", "Metabolism of RNA") ~ "RNA/Proteins Metabolism",
#           TRUE ~ "Other/No Pathway"
#         )
#       )
    
    
  }
  
  #get_intersection_composition_df
  {
    positive_intersection_composition_df <- get_intersection_composition_df(df = positive_gene_upset_df,
                                                                            group_col   = "group_reactome",
                                                                            groups_map  =    c("Immune System"             = "Immune_percent",
                                                                                               "Metabolism"                   = "Metabolism_percent",
                                                                                               "RNA/Proteins Metabolism"      = "RNA_and_Proteins_Metabolism_percent",
                                                                                               "Other/No Pathway"           = "Other_percent"))

    negative_intersection_composition_df <- get_intersection_composition_df(df = negative_gene_upset_df,
                                                                            group_col   = "group_reactome",
                                                                            groups_map  =    c("Immune System"             = "Immune_percent",
                                                                                               "Metabolism"                   = "Metabolism_percent",
                                                                                               "RNA/Proteins Metabolism"      = "RNA_and_Proteins_Metabolism_percent",
                                                                                               "Other/No Pathway"           = "Other_percent"))
                                                                                               
  }    
  
  # calling get_intersection_composition_stack_barplot
  {
    
    group_name_map <- c(
      Immune_percent               = "Immune System",
      Metabolism_percent           = "Metabolism",
      RNA_and_Proteins_Metabolism_percent  = "RNA/Proteins Metabolism",
      Other_percent                = "Other/No Pathway")
    
    group_palette <-  c(
      "Metabolism"                        = "#219ebc",
      "Immune System"                     = "#fb8500",
      "RNA/Proteins Metabolism"       = "#ffb703",
      "Other/No Pathway"                = "#023047")
    
    composition_plot_of_positive_overlap_genes <- get_intersection_composition_stack_barplot(
      df              = positive_intersection_composition_df,
      cell_cols       = c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells"),
      percent_cols    = c("Immune_percent","Metabolism_percent","RNA_and_Proteins_Metabolism_percent","Other_percent"),
      group_name_map  = group_name_map,
      group_palette   = group_palette,
      legend_title    = "Gene Label Group",
      group_order    = c("Other/No Pathway", "Immune System","RNA/Proteins Metabolism", "Metabolism"),
      legend_order   = c("Metabolism", "RNA/Proteins Metabolism","Immune System","Other/No Pathway"),
      plot_title = "Composition of Shared Positive Genes Across Immune Cell Types",
      y_line =10
    )                                                                                            
    
    # ggsave(
    #   filename =       "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/composition_stacked_bar_plots/positive_overlap_genes_composition_plot_all_labels.jpg",
    #   plot     = composition_plot_of_positive_overlap_genes$plot,
    #   width    = 17,
    #   height   = 12,
    #   units    = "in",
    #   dpi      = 600
    # )
    
    composition_plot_of_negative_overlap_genes <- get_intersection_composition_stack_barplot(
      df              = negative_intersection_composition_df,
      cell_cols       = c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells"),
      percent_cols    = c("Immune_percent","Metabolism_percent","RNA_and_Proteins_Metabolism_percent","Other_percent"),
      group_name_map  = group_name_map,
      group_palette   = group_palette,
      legend_title    = "Gene Label Group",
      group_order    = c("Other/No Pathway", "Immune System","Metabolism", "RNA/Proteins Metabolism"),
      legend_order   = c("Metabolism", "RNA/Proteins Metabolism","Immune System","Other/No Pathway"),
      plot_title = "Composition of Shared Negative Genes Across Immune Cell Types",
      y_line =13)                                                                                            
    
    # ggsave(
    #   filename =       "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/composition_stacked_bar_plots/negative_overlap_genes_composition_plot_all_labels.jpg",
    #   plot     = composition_plot_of_negative_overlap_genes$plot,
    #   width    = 17,
    #   height   = 12,
    #   units    = "in",
    #   dpi      = 600
    # )
    
    
  }

  #overlaps one on top of the other
  {
    # Combine plots: p_sd on the left, p_thres on the right
    # p_whole_composition <-
    #   (composition_plot_of_positive_overlap_genes$plot |
    #      composition_plot_of_negative_overlap_genes$plot) +
    #   plot_annotation(tag_levels = 'A') &
    #   theme(plot.tag = element_text(face = "bold", size = 18))
    
    p_whole_composition <-
      (composition_plot_of_positive_overlap_genes$plot /
         composition_plot_of_negative_overlap_genes$plot) +
      plot_annotation(tag_levels = 'A') &
      theme(plot.tag = element_text(face = "bold", size = 18))
    
    ggsave(
      filename =       "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/composition_stacked_bar_plots/whole_overlap_genes_composition_plot_all_labels.jpg",
      plot     = p_whole_composition,
      width    = 17,     # or any value >50
      height   = 24,
      units    = "in",
      dpi      = 600
    )
    
  }


}

#ORA
{

  #ORA itself
  {
    NK_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    monocytes_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    B_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    CD4_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    CD8_T_cells_meta_table <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    
    
    meta_table_list <- list(NK_meta_table = NK_meta_table,
                            monocytes_meta_table = monocytes_meta_table,
                            B_cells_meta_table = B_cells_meta_table,
                            CD4_T_cells_meta_table = CD4_T_cells_meta_table,
                            CD8_T_cells_meta_table = CD8_T_cells_meta_table)
    
    
    rm(NK_meta_table, monocytes_meta_table, B_cells_meta_table, CD4_T_cells_meta_table, CD8_T_cells_meta_table)
    
    #positive ORA
    {
      
      cols <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
      
      
      # positive_ora <- get_ora(meta_table_list, log2FC_threshold_for_gene = 1.25, positive = TRUE ,q_value_threshold_for_gene = 0.05, cols, q_value_threshold_for_pathway = 0.1)
      # 
      # positive_ora_df <- get_full_gsea_and_ora_dfs(ora_return_list  = positive_ora,
      #                                              save_dir_ora     = "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/",
      #                                              condition = "CD59.Def", dataset_name = "Mevorach", sample_name = "CD59_Patient")

      # qsave(positive_ora_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/positive_ora_df.qs")
      positive_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/positive_ora_df.qs")
      
      positive_ora_df <- positive_ora_df$meta_ora_top
      excluded_ids <- c(
        "R-HSA-192823",
        "R-HSA-975956",
        "R-HSA-975957",
        "R-HSA-9735869",
        "R-HSA-9754678"
      )
      
      positive_ora_df <- positive_ora_df %>%
        filter(! ID %in% excluded_ids) %>% 
        filter(Count >3)
      
    }
    
    
    #negative ORA
    {                                                                                                                                                                                                                         
      cols <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
      
      
      # negative_ora <- get_ora(meta_table_list, log2FC_threshold_for_gene = 1.25, positive = FALSE ,q_value_threshold_for_gene = 0.05, cols)
      # 
      # negative_ora_df <- get_full_gsea_and_ora_dfs(ora_return_list  = negative_ora)
      
      # qsave(negative_ora_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/negative_ora_df.qs")
      negative_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/negative_ora_df.qs")
      
      negative_ora_df <- negative_ora_df$meta_ora_top
      
      excluded_ids <- c(
        "R-HSA-168255",
        "R-HSA-192823",
        "hsa05134",
        "hsa05171",
        "GO:0032496",
        "GO:0002237",
        "GO:0071222",
        "GO:0071219",
        "GO:0071216",
        "GO:0045444",
        "GO:0030217",
        "GO:0045580",
        "GO:0043367",
        "GO:0035710",
        "GO:0046632",
        "GO:1904019",
        "GO:0043370",
        "GO:0042093",
        "GO:0046637",
        "GO:0002294",
        "GO:0002287",
        "GO:0002293",
        "GO:2000514",
        "GO:0002292",
        "GO:2000319",
        "GO:2000316",
        "GO:0072539",
        "GO:0045622",
        "GO:0072538",
        "GO:0050673",
        "GO:0050678",
        "R-HSA-975956",
        "GO:0019915")
      
      
      negative_ora_df <- negative_ora_df %>%
        filter(! ID %in% excluded_ids) %>% 
        filter(Count >3)
      
    }
    
  }
  
  #lolipops
  {
      #positive ORA lolipop
      {
        positive_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/positive_ora_df.qs")
        positive_ora_df <- positive_ora_df$meta_ora_top 
        # list_of_highlight_pathways <- c("R-HSA-71291", "R-HSA-376176", "R-HSA-2408557", "R-HSA-2408522", "R-HSA-9633012", "R-HSA-927802", "R-HSA-156827")
        list_of_highlight_pathways <- c("R-HSA-71291", "R-HSA-2408557", "R-HSA-9633012")
        
        # try <- positive_ora_df %>%
        #   filter(ID %in% c("R-HSA-71291", "R-HSA-376176", "R-HSA-2408557", "R-HSA-2408522", "R-HSA-9633012", "R-HSA-927802", "R-HSA-156827")) %>% 
        #   pull(Description)
        
        
        
        manuall_changing_names <- list(
          "Selenocysteine synthesis" = "Selenocysteine Synthesis",
          # "Selenoamino acid metabolism" = "Selenoamino Acid Metabolism",
          "Response of EIF2AK4 (GCN2) to amino acid deficiency" ="GCN2 Response to Amino Acid Deficiency"
          # "L13a-mediated translational silencing of Ceruloplasmin expression" = "L13a-Mediated CP Translational Silencing",
          # "Nonsense-Mediated Decay (NMD)" = "Nonsense-Mediated Decay (NMD)",
          # "Signaling by ROBO receptors" = "Signaling by ROBO Receptors"
          )
        

        
    positive_ORA_lolipop <- get_pathway_lollipop(
    df =positive_ora_df ,
    list_of_highlight_pathways = list_of_highlight_pathways,
    manuall_changing_names = manuall_changing_names,
    title = "ORA of Positive DEGs Shared Across All Five Cell Types",
    lollipop_color = "#219ebc")
        


        # ggsave(
        #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/ORA/positive_ORA_lolipop.jpg",
        #   plot     = positive_ORA_lolipop$p_pathway_lollipop,
        #   width    = 26,
        #   height   = 20,
        #   units    = "in",
        #   dpi      = 600,
        #   device   = "jpeg",
        #   type     = "cairo",
        #   quality  = 100
        # )
        
      }
      
      #negative ORA lolipop
      {
        negative_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/negative_ora_df.qs")
        negative_ora_df <-  negative_ora_df$meta_ora_top

        
        # list_of_highlight_pathways <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
        #                                 "GO:0002757", "R-HSA-9711097", "R-HSA-9633012",
        #                                 "R-HSA-2408557", "R-HSA-2408522", "GO:0097191")
        
        list_of_highlight_pathways <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                                        "GO:0002757", "R-HSA-9633012",
                                        "R-HSA-2408557", "GO:0097191")
        
        

        # try <- negative_ora_df %>%
        #   filter(ID %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
        #                    "GO:0002757", "R-HSA-9711097", "R-HSA-9633012",
        #                    "R-HSA-2408557", "R-HSA-2408522", "GO:0097191")) %>% 
        #   pull(Description)
        
        
        
        manuall_changing_names <- list(
          "Selenocysteine synthesis" = "Selenocysteine Synthesis",
          # "Selenoamino acid metabolism" = "Selenoamino Acid Metabolism",
          "Response of EIF2AK4 (GCN2) to amino acid deficiency" ="GCN2 Response to Amino Acid Deficiency",
          HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA Signaling via NFKB",
          HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory Response",
          "immune response-activating signaling pathway" = "Immune Response-Activating Pathway",
          "extrinsic apoptotic signaling pathway" = "Extrinsic Apoptotic Signaling")
        

        
        positive_ORA_lolipop <- get_pathway_lollipop(
          df =negative_ora_df ,
          list_of_highlight_pathways = list_of_highlight_pathways,
          manuall_changing_names = manuall_changing_names,
          title = "ORA of Negative DEGs Shared Across All Five Cell Types",
          lollipop_color = "#fb8500")
        
        


        ggsave(
          filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/ORA/negative_ORA_lolipop.jpg",
          plot     = positive_ORA_lolipop$p_pathway_lollipop,
          width    = 26,
          height   = 20,
          units    = "in",
          dpi      = 600,
          device   = "jpeg",
          type     = "cairo",
          quality  = 100
        )

      }
    }
    
  #scatter
  {
      positive_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/positive_ora_df.qs")
      positive_ora_df <- positive_ora_df$meta_ora_top
      positive_ora_df <- positive_ora_df %>%
        filter(Count >3)
      negative_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/negative_ora_df.qs")
      negative_ora_df <-  negative_ora_df$meta_ora_top
      negative_ora_df <- negative_ora_df %>%
        filter(Count >3)
      
    p  <- get_scatter_ORA(
        positive_ora_df, negative_ora_df,
        title      = "ORA of Shared DEGs in CD59-Deficient PBMCs",
        big_title  = 30,
        big_axis   = 26,
        big_txt    = 20,
        big_legend = 20
      )
      ggsave(
        filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/ORA/scatter_ORA_plot.jpg",
        plot     = p,
        width    = 20,
        height   = 17,
        units    = "in",
        dpi      = 100,
        device   = "jpeg",
        type     = "cairo",
        quality  = 600
      )

      
    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
    #   plot     = p,
    #   width    = 20,
    #   height   = 17,
    #   units    = "in",
    #   dpi      = 100,
    #   device   = "jpeg",
    #   type     = "cairo",
    #   quality  = 600
    # )

      
      
      
      
    }

}



#pathview plots
{
  data("sbgn.xmls")          # SBGN-ML files for Reactome, SMPDB, PANTHER, etc.
  data("pathways.info")      # metadata for every pathway in the collection
  data("pathways.stats")     # counts per database (optional, for curiosity)
  
  hits <- SBGNview::findPathways("ISR")
  
  #relevant pathwatts
  {
    # KEGG: TNF signaling (hsa04668)
    # Reactome: TNFR1-induced NF-κB signaling (R-HSA-5357956)
    # Reactome: TNFR2 non-canonical NF-κB (R-HSA-5668541)
    # Reactome: Detoxification of Reactive Oxygen Species (R-HSA-3299685)	
  }
  
  #Get the median of the genes you need
  {
    
    meta_table_nk          <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_monocytes   <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_b_cells     <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_cd4_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_cd8_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    
    
    meta_table_list <- list(meta_table_nk = meta_table_nk,
                            meta_table_monocytes = meta_table_monocytes,
                            meta_table_b_cells = meta_table_b_cells,
                            meta_table_cd4_t_cells = meta_table_cd4_t_cells,
                            meta_table_cd8_t_cells = meta_table_cd8_t_cells)
    
    rm(meta_table_nk, meta_table_monocytes, meta_table_b_cells, meta_table_cd4_t_cells, meta_table_cd8_t_cells)
    

    gene_data_tnf_kegg <- get_SBGNview_df(pathway_id = "R-HSA-5357956", meta_table_list = meta_table_list)
    
    # # find a Reactome SBGN for TNF/NF-κB
    # rx <- findPathways("TNF", keyword.type = "pathway.name")
    # rx <- subset(rx, sub.database == "reactome")
    # pid <- rx$pathway.id["hsa04668"]   
    
    # SBGNview.obj <- SBGNview(
    #   gene.data         = gene_data_tnf_kegg,
    #   gene.id.type      = "SYMBOL",         
    #   input.sbgn        = "R-HSA-5357956", 
    #   sbgn.gene.id.type = "pathwayCommons",  # glyph ID type in the bundled SBGN files
    #   output.file       = "tnf_meta_scores",
    #   output.formats    = c("svg","png")
    # )
    
  
    
  }
  
  
  # pathview tnf_kegg
  {
    #get data
    {
      meta_table_nk          <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
      meta_table_monocytes   <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
      meta_table_b_cells     <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
      meta_table_cd4_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
      meta_table_cd8_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
      
      
      meta_table_list <- list(meta_table_nk = meta_table_nk,
                              meta_table_monocytes = meta_table_monocytes,
                              meta_table_b_cells = meta_table_b_cells,
                              meta_table_cd4_t_cells = meta_table_cd4_t_cells,
                              meta_table_cd8_t_cells = meta_table_cd8_t_cells)
      
      rm(meta_table_nk, meta_table_monocytes, meta_table_b_cells, meta_table_cd4_t_cells, meta_table_cd8_t_cells)
      
      
      gene_data_tnf_kegg <- get_SBGNview_df(pathway_id = "hsa04668", meta_table_list = meta_table_list)
      
      gene.vec <- setNames(as.numeric(gene_data_tnf_kegg[, "meta_score"]),
                           rownames(gene_data_tnf_kegg))
      
    }
    
    #get pathview
    {
      gene_disc5 <- cut(
        gene.vec,
        breaks = c(-Inf, -0.5, -0.25, 0.25, 0.5, Inf),
        labels = c(-2L, -1L, 0L, 1L, 2L)
      )
      
      # Convert factor to numeric and keep gene names
      gene_disc5 <- setNames(as.numeric(as.character(gene_disc5)), names(gene.vec))
      
      # Run Pathview with 5 discrete bins
      pv <- pathview::pathview(
        gene.data   = gene_disc5,
        pathway.id  = "hsa04668",
        species     = "hsa",
        gene.idtype = "SYMBOL",
        kegg.native = TRUE,
        node.sum    = "mean",      # <- mean across isoforms in each node
        out.suffix  = "tnf_blocked_pm5_discrete",
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
    
    # get legend
    {
      brks <- c(min(gene.vec, na.rm = TRUE), -0.5, -0.25, 0.25, 0.5, max(gene.vec, na.rm = TRUE))
      cols <- c("red", "#FF4646", "gray", "#B2DF8A", "green")
      labs <- format(round(brks, 2), trim = TRUE)

      draw_legend(type = "png", file = "legend_discrete_5bins.png")
      
      
      
    }
    
 

  }
  
  
  # calling get_pathview_full on many pathways
  {

    meta_table_nk          <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_monocytes   <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_b_cells     <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_cd4_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_cd8_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    
    
    meta_table_list <- list(meta_table_nk = meta_table_nk,
                            meta_table_monocytes = meta_table_monocytes,
                            meta_table_b_cells = meta_table_b_cells,
                            meta_table_cd4_t_cells = meta_table_cd4_t_cells,
                            meta_table_cd8_t_cells = meta_table_cd8_t_cells)
    
    rm(meta_table_nk, meta_table_monocytes, meta_table_b_cells, meta_table_cd4_t_cells, meta_table_cd8_t_cells)
    
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa00480",
                      out.suffix = "glut_met",
                                  breaks = c(-Inf, -0.5, -0.25, 0.25, 0.5, Inf),
                                  labels = c(-2L, -1L, 0L, 1L, 2L))
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa00190",
                      out.suffix = "oxphos")
    
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa04670",
                      out.suffix = "trans_migration")
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa05208",
                      out.suffix = "ROS")  
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa04218",
                      out.suffix = "senescence") 
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa04010",
                      out.suffix = "mapk") 
    
    
    get_pathview_full(meta_table_list = meta_table_list,
                      pathway_id = "hsa04610",
                      out.suffix = "complement") 
    
    
    
  }
  
  # calling get_pathview_per_cell_type 
  {
    get_pathview_per_cell_type(meta_table_list = meta_table_list,
                               pathway_id = "hsa04010",
                               out.suffix = "mapk")
  }
  

  
    
}


#CD59 only AUC heatmaps
{
  #getting meta_table_list
  {
    meta_table_nk          <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_monocytes   <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_b_cells     <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_cd4_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    meta_table_cd8_t_cells <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/CD59.Def_Mevorach_CD59_Patient_meta_table.qs")
    
    
    meta_table_list <- list(meta_table_nk = meta_table_nk,
                            meta_table_monocytes = meta_table_monocytes,
                            meta_table_b_cells = meta_table_b_cells,
                            meta_table_cd4_t_cells = meta_table_cd4_t_cells,
                            meta_table_cd8_t_cells = meta_table_cd8_t_cells)
    
    rm(meta_table_nk, meta_table_monocytes, meta_table_b_cells, meta_table_cd4_t_cells, meta_table_cd8_t_cells)
    
  }
  
  #mitochondrial complexes auc heatmap
  {
    vec_path_ids <- c(
      "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_I",
      "REACTOME_COMPLEX_I_BIOGENESIS",
      
      "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
      "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
      
      "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
      "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
      # "REACTOME_COMPLEX_III_ASSEMBLY",
      # "HP_DECREASED_ACTIVITY_OF_MITOCHONDRIAL_COMPLEX_III",
      
      "WP_MITOCHONDRIAL_COMPLEX_IV_ASSEMBLY",
      "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_IV",
      
      "GOBP_MITOCHONDRIAL_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX_ASSEMBLY",
      'GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS')
    
    resp_df <- get_meta_table_auc_score_for_many_cell_types_and_pathways(meta_table_list = meta_table_list,
                                                                         vec_path_ids = vec_path_ids)
    
    

    ## ── 1. Prepare matrix for heatmap ─────────────────────────────
    # drop the num_genes column
    mat <- as.matrix(resp_df[, c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells")])
    
    # rows are pathways
    rownames(mat) <- rownames(resp_df)
    
    ## ── 2. Rename pathways if you want shorter labels ─────────────
    # Example mapping (replace with your own edits)
    pathway_names <- c(
      "Kegg Medicus Reference Electron Transfer in Complex i" = "KEGG-MR: ETC in Complex I",
      "Reactome Complex i Biogenesis" = "Reactome: Complex I Biogenesis",
      "Wp Mitochondrial Complex Ii Assembly" = "WP: Complex II Assembly",
      "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "GO-BP:ETC in Complex II",
      "Wp Mitochondrial Complex Iii Assembly" = "WP: Complex III Assembly",
      "Kegg Medicus Reference Electron Transfer in Complex Iii" = "KEGG-MR: ETC in Complex III",
      "Wp Mitochondrial Complex Iv Assembly" = "WP: Complex IV Assembly",
      "Kegg Medicus Reference Electron Transfer in Complex Iv" = "KEGG-MR: ETC in Complex IV",
      "Gobp Mitochondrial Proton Transporting Atp Synthase Complex Assembly" = "GO-BP: Complex V Assembly",
      "Gobp Proton Motive Force Driven Atp Synthesis" = "GO-BP: PMF-Driven ATP Synthesis"
      # "Wp Oxidative Damage Response" = "WP: Oxidative Damage Response",
    )
    
    rownames(mat) <- pathway_names[rownames(mat)]
    
    ## ── 3. Control order of rows and columns ─────────────────────
    # Row order (define as character vector in desired sequence)
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
      # "WP: Oxidative Damage Response",
    )
    
    # Column order
    col_order <- c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells")
    
    ## ── 4. Define color scale (blue → grey → red, with grey at 0.5) ─────────────
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "grey", "red"))
    
    ## ── 5. Column annotation (legend for cell types) ─────────────
    # make the annotation a factor to control legend order
    cell_types <- factor(colnames(mat), levels = col_order)
    
    cell_colors <- c(
      "NKs"          = "#8ecae6",
      "Monocytes"    = "#219ebc",
      "B Cells"      = "#023047",
      "CD4+ T Cells" = "#ffb703",
      "CD8+ T Cells" = "#fb8500"
    )
    
    ha_col <- HeatmapAnnotation(
      CellType = cell_types,
      col = list(CellType = cell_colors),
      annotation_legend_param = list(title = "Cell Types", at = col_order),
      annotation_name_gp = gpar(fontsize = 0)  # hide the "CellType" label next to columns
    )
    
    rn_gp <- gpar(fontsize = 11)  # tweak as you like
    
    ht <- Heatmap(
      mat,
      name = "AUC",
      col = col_fun,
      column_title = "Mitochondrial Complexes AUC Heatmap Across Cell Types",
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_order = row_order,
      column_order = col_order,
      top_annotation = ha_col,
      show_column_names = FALSE,
      row_names_side = "left",
      column_names_side = "top",
      column_names_rot = 45,
      row_names_gp = rn_gp,
      row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), gp = rn_gp),
      heatmap_legend_param = list(
        title  = "AUC score",
        at     = seq(0, 1, by = 0.1),
        labels = sprintf("%.2f", seq(0, 1, by = 0.1))
      ),
      cell_fun = function(j, i, x, y, width, height, fill) {
        val <- mat[i, j]
        lab <- sprintf("%.2f", trunc(val * 100) / 100)  # truncate to 2 decimals, no rounding
        grid.text(lab, x, y, gp = gpar(col = "black", fontsize = 9, fontface = "bold"))
      }
    )
    
    
    
    jpeg(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress/mitochondrial_complexes_auc_heatmap.jpg",
      width    = 11, height = 8, units = "in", res = 600, quality = 100
    )

    # add extra left padding (top, right, bottom, left)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         padding = unit(c(5, 5, 5, 25), "mm"))

    dev.off()
    



    
    

    
    
  }
  
  #tca auc heatmap
  
  {
    tca_sets <- c(
      "WP_THIAMINE_METABOLIC_PATHWAYS",
      "GOBP_PYRUVATE_METABOLIC_PROCESS",
      "GOMF_PYRUVATE_DEHYDROGENASE_ACETYL_TRANSFERRING_PHOSPHATASE_ACTIVITY",
      "REACTOME_MALATE_ASPARTATE_SHUTTLE",
      "GOBP_GLYCEROL_3_PHOSPHATE_METABOLIC_PROCESS",
      "GOBP_GLUTAMATE_METABOLIC_PROCESS",
      "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
      "KEGG_PROPANOATE_METABOLISM",
      "GOBP_SUCCINYL_COA_CATABOLIC_PROCESS",
      "GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE",
      "GOBP_UNSATURATED_FATTY_ACID_METABOLIC_PROCESS",
      # "GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_OXIDASE",
      "GOBP_KETONE_METABOLIC_PROCESS",
      "GOMF_ACETYL_COA_C_ACYLTRANSFERASE_ACTIVITY",
      "REACTOME_CARNITINE_SHUTTLE",
      "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION"
    )
    
    

    tca_df <- get_meta_table_auc_score_for_many_cell_types_and_pathways(meta_table_list = meta_table_list,
                                                                        vec_path_ids = tca_sets)
    
    
    
    
    
    ## ── 1. Prepare matrix for heatmap ─────────────────────────────
    # drop the num_genes column
    mat <- as.matrix(tca_df[, c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells")])
    
    # rows are pathways
    rownames(mat) <- rownames(tca_df)
    
    ## ── 2. Rename pathways if you want shorter labels ─────────────
    # Example mapping (replace with your own edits)
    pathway_names <- c(
      "Wp Thiamine Metabolic Pathways" = "WP: Thiamine Metabolism",
      "Gobp Pyruvate Metabolic Process" = "GO-BP: Pyruvate Metabolism",
      "Gomf Pyruvate Dehydrogenase Acetyl Transferring Phosphatase Activity" = "GO-MF: PDH Phosphatase Activity",
      "Reactome Malate Aspartate Shuttle" = "Reactome: Malate-Aspartate Shuttle",
      "Gobp Glycerol 3 Phosphate Metabolic Process" = "GO-BP: Glycerol-3P Metabolism",
      "Gobp Glutamate Metabolic Process" = "GO-BP: Glutamate Metabolism",
      "Kegg Valine Leucine and Isoleucine Degradation" = "KEGG: BCAA (Val/Leu/Ile) Degradation",
      "Kegg Propanoate Metabolism" = "KEGG: Propanoate Metabolism",
      "Gobp Succinyl Coa Catabolic Process" = "GO-BP: Succinyl-CoA Catabolism",
      "Gobp Fatty Acid Beta Oxidation using Acyl Coa Dehydrogenase" = "GO-BP: FA β-Ox (ACAD)",
      "Gobp Unsaturated Fatty Acid Metabolic Process" = "GO-BP: Unsaturated FA Metabolism",
      # "Gobp Fatty Acid Beta Oxidation using Acyl Coa Oxidase" = "GO-BP: FA β-Ox (ACOX)",
      "Gobp Ketone Metabolic Process" = "GO-BP: Ketone Metabolism",
      "Gomf Acetyl Coa c Acyltransferase Activity" = "GO-MF: Thiolase Activity",
      "Reactome Carnitine Shuttle" = "Reactome: Carnitine Shuttle",
      "Reactome Mitochondrial Fatty Acid Beta Oxidation" = "Reactome:FA β-Ox (Mito)"
      
    )
    
    rownames(mat) <- pathway_names[rownames(mat)]
    
    ## ── 3. Control order of rows and columns ─────────────────────
    # Row order (define as character vector in desired sequence)
    row_order <- c(
      "WP: Thiamine Metabolism",
      "GO-BP: Pyruvate Metabolism",
      "GO-MF: PDH Phosphatase Activity",
      "GO-BP: Glutamate Metabolism",
      "GO-BP: Ketone Metabolism",
      "GO-BP: Unsaturated FA Metabolism",
      "Reactome: Malate-Aspartate Shuttle",
      "GO-BP: Glycerol-3P Metabolism",
      "KEGG: BCAA (Val/Leu/Ile) Degradation",
      "KEGG: Propanoate Metabolism",
      "GO-BP: Succinyl-CoA Catabolism",
      "GO-BP: FA β-Ox (ACAD)",
      "GO-MF: Thiolase Activity",
      "Reactome:FA β-Ox (Mito)",
      "Reactome: Carnitine Shuttle"
    )
    
    # Column order
    col_order <- c("NKs","Monocytes","B Cells","CD4+ T Cells","CD8+ T Cells")
    
    ## ── 4. Define color scale (blue → grey → red, with grey at 0.5) ─────────────
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "grey", "red"))
    
    ## ── 5. Column annotation (legend for cell types) ─────────────
    # make the annotation a factor to control legend order
    cell_types <- factor(colnames(mat), levels = col_order)
    
    cell_colors <- c(
      "NKs"          = "#8ecae6",
      "Monocytes"    = "#219ebc",
      "B Cells"      = "#023047",
      "CD4+ T Cells" = "#ffb703",
      "CD8+ T Cells" = "#fb8500"
    )
    
    ha_col <- HeatmapAnnotation(
      CellType = cell_types,
      col = list(CellType = cell_colors),
      annotation_legend_param = list(title = "Cell Types", at = col_order),
      annotation_name_gp = gpar(fontsize = 0)  # hide the "CellType" label next to columns
    )
    
    rn_gp <- gpar(fontsize = 11)  # tweak as you like
    
    ht <- Heatmap(
      mat,
      name = "AUC",
      col = col_fun,
      column_title = "TCA Pathways AUC Heatmap Across Cell Types",
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),  # <-- make title bold
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_order = row_order,
      column_order = col_order,
      top_annotation = ha_col,
      row_names_side = "left",
      column_names_side = "top",
      column_names_rot = 45,
      show_column_names = FALSE,
      row_names_gp = rn_gp,
      row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), gp = rn_gp),
      heatmap_legend_param = list(
        title  = "AUC score",
        at     = seq(0.1, 1, by = 0.1),
        labels = sprintf("%.1f", seq(0.1, 1, by = 0.1))
      )
    )
    
    
    jpeg(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress/tca_auc_heatmap.jpg",
      width    = 11, height = 8, units = "in", res = 600, quality = 100
    )
    
    # add extra left padding (top, right, bottom, left)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         padding = unit(c(5, 5, 5, 25), "mm"))
    
    dev.off()
    
    
    
  }
  
  msig_try <- find_msig_sets_with_genes(c("GLRX5","BOLA1"),
                                        require_all = TRUE) 
  
  # DUSP3, PTPN7, PTPRR, PTPN5
  # "MAPK14","TP53", "DUSP3"

#manually checking the genes of ETC in complexes 4 and 5
  {
    
    
    # ---- inputs ----
    pathways <- c(
      IV = "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_IV",
      V = "GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS"
    )
    
    # Optional: nice column names for your meta tables
    cell_name_map <- c(
      meta_table_nk          = "NK",
      meta_table_monocytes   = "Monocytes",
      meta_table_b_cells     = "B_cells",
      meta_table_cd4_t_cells = "CD4_T_cells",
      meta_table_cd8_t_cells = "CD8_T_cells"
    )
    
    # ---- pull gene sets from msigdbr ----
    msig <- msigdbr(species = "Homo sapiens") %>%
      dplyr::select(gs_name, gene_symbol)
    
    get_genes_for <- function(gs_name) {
      msig %>%
        filter(gs_name == !!gs_name) %>%
        pull(gene_symbol) %>%
        unique() %>%
        sort()
    }
    
    # ---- helper: build a gene x cell-type table of median_log2FC ----
    build_matrix_for_set <- function(genes, meta_table_list, cell_name_map = NULL) {
      # Start with all pathway genes to ensure rows exist even if all values are NA
      out <- tibble(gene = genes)
      
      # For each cell table: join its median_log2FC on 'gene'
      for (nm in names(meta_table_list)) {
        colname <- if (!is.null(cell_name_map) && nm %in% names(cell_name_map)) {
          cell_name_map[[nm]]
        } else nm
        
        add_col <- meta_table_list[[nm]] %>%
          dplyr::select(gene, median_log2FC) %>%
          distinct() %>%
          rename(!!colname := median_log2FC)
        
        out <- out %>% left_join(add_col, by = "gene")
      }
      
      # Put genes as rownames, return a plain data.frame
      out %>%
        arrange(match(gene, genes)) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("gene")
    }
    
    # ---- build the two tables ----
    genes_gocc <- get_genes_for(pathways[["IV"]])
    genes_gomf <- get_genes_for(pathways[["V"]])
    
    complex_v_list <- list(
      IV = build_matrix_for_set(genes_gocc, meta_table_list, cell_name_map),
      V = build_matrix_for_set(genes_gomf, meta_table_list, cell_name_map)
    )
    

    
    
  }


  
  # carbs:
  #   Wp Thiamine Metabolic Pathways - 0.4
  # "Gobp Pyruvate Metabolic Process" - 0.4
  # "Gomf Pyruvate Dehydrogenase Acetyl Transferring Phosphatase Activity" - 0.4
  # 
  # redox import:
  #   Reactome Malate Aspartate Shuttle - 0.4,
  # Gobp Glycerol 3 Phosphate Metabolic Process - 0.3,
  # 
  # anaplerosis shifts:
  # GOBP_GLUTAMATE_METABOLIC_PROCESS - 0.25
  # Kegg Valine Leucine and Isoleucine Degradation - 0.5
  # Kegg Propanoate Metabolism - 0.5
  # Gobp Succinyl Coa Catabolic Process - 0.6
  # 
  # mitochondrial fats
  # Gobp Fatty Acid Beta Oxidation using Acyl Coa Dehydrogenase - 0.55
  # Gobp Unsaturated Fatty Acid Metabolic Process - 0.2
  # Gobp Fatty Acid Beta Oxidation using Acyl Coa Oxidase - 0.3
  # 
  # Acetyl-CoA finishing step is primed, but ketones aren’t the fallback:
  # Gobp Ketone Metabolic Process-  0.3,
  # Gomf Acetyl Coa c Acyltransferase Activity -0.6


  
}

#AUC scores on neuroimuunulogical samples
{
  #AUC mitchondrial complexes
  {
    #create neuroimmunulogical_meta_table_lists for AUC scores
    {
      cell_paths <- list(
        "NKs"            = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/",
        "Monocytes"      = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/",
        "B Cells"        = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/",
        "CD4 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/",
        "CD8 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
      )
      
      # Map cell types to element names inside meta_table_list
      element_name_map <- c(
        "NKs"            = "meta_table_nk",
        "Monocytes"      = "meta_table_monocytes",
        "B Cells"        = "meta_table_b_cells",
        "CD4 + T cells"  = "meta_table_cd4_t_cells",
        "CD8 + T cells"  = "meta_table_cd8_t_cells"
      )
      
      # Helper: safe qread
      safe_qread <- function(path) {
        tryCatch(qread(path), error = function(e) NULL)
      }
      
      # ── Gather filenames per cell type ─────────────────────────────────────────────
      files_by_ct <- lapply(cell_paths, function(dir) {
        if (!dir.exists(dir)) return(character(0))
        f <- list.files(dir, pattern = "\\.qs$", full.names = TRUE)
        setNames(f, basename(f))  # names are basenames; values are full paths
      })
      
      # Union of sample filenames across all folders
      all_samples <- sort(unique(unlist(lapply(files_by_ct, names))))
      
      # ── Build the list-of-lists ───────────────────────────────────────────────────
      neuroimmunulogical_meta_table_lists <- list()
      
      for (sample_file in all_samples) {
        # Store under name without "_meta_table.qs"
        sample_key <- sub("_meta_table$", "", sub("\\.qs$", "", sample_file))
        
        meta_table_list <- list()
        for (ct in names(cell_paths)) {
          paths_vec <- files_by_ct[[ct]]
          if (length(paths_vec) && sample_file %in% names(paths_vec)) {
            df <- safe_qread(paths_vec[[sample_file]])
            if (!is.null(df)) {
              el_name <- element_name_map[[ct]]
              meta_table_list[[el_name]] <- df
            }
          }
        }
        
        if (length(meta_table_list) > 0) {
          neuroimmunulogical_meta_table_lists[[sample_key]] <- meta_table_list
        }
      }
      
      # Optional sanity check:
      # str(neuroimmunulogical_meta_table_lists, max.level = 2)
      
      # Optional save:
      # qsave(neuroimmunulogical_meta_table_lists,
      #       "/bigdata/bgutman/Collabs/Dror mevorach/data/neuroimmunulogical_meta_table_lists.qs")
    }

    #resp AUC scores for each neuroimmunulogical condition
    {
      # ## ── Vector of pathway IDs ─────────────────────────────────────────────────────
      # vec_path_ids <- c(
      #   "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_I",
      #   "REACTOME_COMPLEX_I_BIOGENESIS",
      #   
      #   "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
      #   "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
      #   
      #   "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
      #   "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
      #   # "REACTOME_COMPLEX_III_ASSEMBLY",
      #   # "HP_DECREASED_ACTIVITY_OF_MITOCHONDRIAL_COMPLEX_III",
      #   
      #   "WP_MITOCHONDRIAL_COMPLEX_IV_ASSEMBLY",
      #   "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_IV",
      #   
      #   "GOBP_MITOCHONDRIAL_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX_ASSEMBLY",
      #   'GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS'
      # )
      # 
      # ## ── Iterate and compute AUC per sample ───────────────────────────────────────
      # resp_auc_list <- list()
      # 
      # for (nm in names(neuroimmunulogical_meta_table_lists)) {
      #   that_meta_table_list <- neuroimmunulogical_meta_table_lists[[nm]]
      #   res <- tryCatch(
      #     get_meta_table_auc_score_for_many_cell_types_and_pathways(
      #       meta_table_list = that_meta_table_list,
      #       vec_path_ids    = vec_path_ids
      #     ),
      #     error = function(e) {
      #       message("AUC calc failed for ", nm, ": ", conditionMessage(e))
      #       NULL
      #     }
      #   )
      #   if (!is.null(res)) resp_auc_list[[nm]] <- res
      # }
      # 
      # qsave(resp_auc_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/oxidative_stress/neuroimmunulogical_resp_auc_list.qs")
      # 
      # 
      
    }
  
    #get resp_auc_cell_type_list
    {
      resp_auc_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/oxidative_stress_auc_scores/neuroimmunulogical_resp_auc_list.qs")
      
      ## Build a list of 5 data frames (one per cell type) from resp_auc_list
      ## Each row = pathway; columns = pathway, AUC (for that cell type), condition, dataset, sample
      
      cell_types <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
      
      resp_auc_cell_type_list <- setNames(vector("list", length(cell_types)), cell_types)
      
      for (ct in cell_types) {
        # Collect per-sample data frames for this cell type
        pieces <- lapply(names(resp_auc_list), function(nm) {
          mat <- resp_auc_list[[nm]]
          if (is.null(mat) || !is.data.frame(mat) || !(ct %in% colnames(mat))) return(NULL)
          
          # Parse condition / dataset / sample from the name: cond_dataset_sample(with underscores)
          parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
          condition <- if (length(parts) >= 1) parts[1] else NA_character_
          dataset   <- if (length(parts) >= 2) parts[2] else NA_character_
          sample    <- if (length(parts) >= 3) paste(parts[3:length(parts)], collapse = "_") else NA_character_
          
          # Build tidy block for this sample & cell type
          df <- data.frame(
            pathway   = rownames(mat),
            AUC       = as.numeric(mat[[ct]]),
            condition = condition,
            dataset   = dataset,
            sample    = sample,
            stringsAsFactors = FALSE
          )
          
          # (Optional) keep the found_genes column—uncomment if you want it
          # if ("found_genes" %in% colnames(mat)) {
          #   df$found_genes <- as.character(mat[["found_genes"]])
          # }
          
          df
        })
        
        pieces <- Filter(Negate(is.null), pieces)
        resp_auc_cell_type_list[[ct]] <- if (length(pieces)) do.call(rbind, pieces) else
          data.frame(pathway=character(), AUC=numeric(), condition=character(),
                     dataset=character(), sample=character(), stringsAsFactors=FALSE)
        
        rownames(resp_auc_cell_type_list[[ct]]) <- NULL
      }
      
      # Quick glance:
      # lapply(resp_auc_cell_type_list, head)
      
      #turn "CD55-patient" into "CD55_Patient"
      {
        # Fix sample naming inconsistency across all cell type dfs
        for (ct in names(resp_auc_cell_type_list)) {
          df <- resp_auc_cell_type_list[[ct]]
          
          # Replace only the exact match "CD55-patient"
          df$sample[df$sample == "CD55-patient"] <- "CD55_Patient"
          
          resp_auc_cell_type_list[[ct]] <- df
        }
        
        # Quick check
        lapply(resp_auc_cell_type_list, function(df) unique(df$sample))
      }
      
      
      
      # add meta data from information_df_neuroimmunulogical_samples
      {
        information_df_neuroimmunulogical_samples <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/information_df_neuroimmunulogical_samples.qs")
        
        cols_to_add <- c("condition", "condition_subtype", "age", "sex", "EDSS")
        
        # Loop through each cell type in the list and merge
        resp_auc_cell_type_list <- lapply(
          resp_auc_cell_type_list,
          function(df) {
            if ("condition" %in% names(df)) df$condition <- NULL  # drop old condition
            merge(
              df,
              information_df_neuroimmunulogical_samples[, c("sample", cols_to_add)],
              by = "sample",
              all.x = TRUE
            )
          }
        )
      }
      
      
      # Turn "CDMS" into "RRMS" under condition_subtype
      {
        for (ct in names(resp_auc_cell_type_list)) {
          df <- resp_auc_cell_type_list[[ct]]
          
          # Replace only the exact match "CDMS"
          df$condition_subtype[df$condition_subtype == "CDMS"] <- "RRMS"
          
          resp_auc_cell_type_list[[ct]] <- df
        }
        
        # Quick check
        lapply(resp_auc_cell_type_list, function(df) unique(df$condition_subtype))
      }
      
    }
    
    #AUC sanity checks
    {
      #Number of NAs, per disease per cell type
      {
        
        ## ── Config ───────────────────────────────────────────────────────────────────
        cell_types <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
        
        # Extract disease from list element name (prefix before first "_")
        get_disease <- function(nm) sub("_.*$", "", nm)
        
        # All pathways across all samples
        all_paths <- sort(unique(unlist(lapply(resp_auc_list, function(x) rownames(x)))))
        
        # Disease per sample and order by first appearance
        diseases_vec <- vapply(names(resp_auc_list), get_disease, FUN.VALUE = character(1))
        disease_levels <- unique(diseases_vec)
        
        # Safe extractor: value for (pathway, cell_type) or NA if missing
        get_cell <- function(mat, pathway, ct) {
          if (is.null(mat) || (!is.matrix(mat) && !is.data.frame(mat))) return(NA_real_)
          if (!(ct %in% colnames(mat))) return(NA_real_)
          if (!(pathway %in% rownames(mat))) return(NA_real_)
          as.numeric(mat[pathway, ct])
        }
        
        ## ── Compute NA counts ────────────────────────────────────────────────────────
        na_stats_long <- do.call(rbind, lapply(disease_levels, function(d) {
          sample_ids <- names(resp_auc_list)[diseases_vec == d]
          mats <- resp_auc_list[sample_ids]
          
          do.call(rbind, lapply(cell_types, function(ct) {
            na_counts <- vapply(all_paths, function(pw) {
              vals <- vapply(mats, get_cell, FUN.VALUE = numeric(1), pathway = pw, ct)
              sum(is.na(vals))
            }, FUN.VALUE = integer(1))
            
            data.frame(
              disease  = d,
              cell_type = ct,
              pathway  = names(na_counts),
              num_NA   = as.integer(na_counts),
              row.names = NULL,
              stringsAsFactors = FALSE
            )
          }))
        }))
        
        ## ── Optional: wide format (one row per disease+pathway, columns per cell type) ─
        na_stats_wide <- reshape(
          na_stats_long,
          idvar   = c("disease", "pathway"),
          timevar = "cell_type",
          direction = "wide"
        )
        
        # Tidy column names: num_NA.NKs -> NKs, num_NA.Monocytes -> Monocytes, num_NA.B Cells -> B_Cells
        names(na_stats_wide) <- sub("^num_NA\\.", "", names(na_stats_wide))
        names(na_stats_wide) <- sub("B Cells$", "B_Cells", names(na_stats_wide))
      }
      
      #pathways median, per cell type in each conditon
      {
        # ── Config ───────────────────────────────────────────────────────────────────
        cell_types <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")  # source column names
        out_cols   <- c("NKs", "Monocytes", "B_cells", "CD4+ T Cells", "CD8+ T Cells")  # desired output column names
        
        # Extract disease from list element name (prefix before first "_")
        get_disease <- function(nm) sub("_.*$", "", nm)
        
        # All pathways across all samples (union)
        all_paths <- sort(unique(unlist(lapply(resp_auc_list, rownames))))
        
        # Disease per sample and order by first appearance
        diseases_vec   <- vapply(names(resp_auc_list), get_disease, FUN.VALUE = character(1))
        disease_levels <- unique(diseases_vec)
        
        # Safe extractor: value for (pathway, cell_type) or NA if missing
        get_cell <- function(mat, pathway, ct) {
          if (is.null(mat) || (!is.matrix(mat) && !is.data.frame(mat))) return(NA_real_)
          if (!(ct %in% colnames(mat))) return(NA_real_)
          if (!(pathway %in% rownames(mat))) return(NA_real_)
          as.numeric(mat[pathway, ct])
        }
        
        # Desired ID order (as given)
        vec_path_ids <- c(
          "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_I",
          "REACTOME_COMPLEX_I_BIOGENESIS",
          "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
          "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
          "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
          "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
          "WP_MITOCHONDRIAL_COMPLEX_IV_ASSEMBLY",
          "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_IV",
          "GOBP_MITOCHONDRIAL_PROTON_TRANSPORTING_ATP_SYNTHASE_COMPLEX_ASSEMBLY",
          "GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS"
        )
        
        # Convert IDs -> display names (matches your rowname style incl. lower-case roman numerals)
        get_pathway_name_by_id <- function(path_id) {
          nm <- gsub("_", " ", toupper(path_id), fixed = TRUE)
          nm <- tools::toTitleCase(tolower(nm))
          nm <- sub("Complex I\\b",   "Complex i",   nm)
          nm <- sub("Complex Ii\\b",  "Complex ii",  nm)
          nm <- sub("Complex Iii\\b", "Complex iii", nm)
          nm <- sub("Complex Iv\\b",  "Complex iv",  nm)
          nm
        }
        desired_row_order <- vapply(vec_path_ids, get_pathway_name_by_id, character(1))
        
        reorder_rows <- function(df, desired) {
          present <- desired[desired %in% rownames(df)]
          extras  <- setdiff(rownames(df), present)
          df[c(present, extras), , drop = FALSE]
        }
        
        # ── Compute MEDIAN (ignore NAs) per disease × cell type × pathway ────────────
        median_auc_by_condition <- setNames(vector("list", length(disease_levels)), disease_levels)
        
        for (d in disease_levels) {
          sample_ids <- names(resp_auc_list)[diseases_vec == d]
          mats <- resp_auc_list[sample_ids]
          
          # preallocate matrix [pathways × out_cols]
          M <- matrix(NA_real_, nrow = length(all_paths), ncol = length(out_cols),
                      dimnames = list(all_paths, out_cols))
          
          for (pw in all_paths) {
            for (j in seq_along(cell_types)) {
              ct <- cell_types[j]
              vals <- vapply(mats, get_cell, FUN.VALUE = numeric(1), pathway = pw, ct = ct)
              M[pw, out_cols[j]] <- if (all(is.na(vals))) NA_real_ else median(vals, na.rm = TRUE)
            }
          }
          
          df <- as.data.frame(M, check.names = FALSE, stringsAsFactors = FALSE)
          df <- reorder_rows(df, desired_row_order)   # ← reorder HERE
          median_auc_by_condition[[d]] <- df
        }
        
        # Quick peek
        # lapply(median_auc_by_condition, head)
        
        
      }
      
    }
    
    resp_auc_cell_type_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/oxidative_stress_auc_scores/resp_auc_cell_type_list.qs")
    resp_auc_cell_type_list <- lapply(resp_auc_cell_type_list, function(df) {
      df %>%
        filter(sample != "CD55_Patient") %>%      # remove CD55_Patient rows
        mutate(age = ifelse(sample == "CD59_Patient", 13, age))  # set age=13 for CD59_Patient
    })
    
    
    #heatmaps
    {
      pathway_names <- c(
        "Kegg Medicus Reference Electron Transfer in Complex i"       = "KEGG-MR: ETC in Complex I",
        "Reactome Complex i Biogenesis"                               = "Reactome: Complex I Biogenesis",
        "Wp Mitochondrial Complex Ii Assembly"                        = "WP: Complex II Assembly",
        "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "GO-BP:ETC in Complex II",
        "Wp Mitochondrial Complex Iii Assembly"                       = "WP: Complex III Assembly",
        "Kegg Medicus Reference Electron Transfer in Complex Iii"     = "KEGG-MR: ETC in Complex III",
        "Wp Mitochondrial Complex Iv Assembly"                        = "WP: Complex IV Assembly",
        "Kegg Medicus Reference Electron Transfer in Complex Iv"      = "KEGG-MR: ETC in Complex IV",
        "Gobp Mitochondrial Proton Transporting Atp Synthase Complex Assembly" = "GO-BP: Complex V Assembly",
        "Gobp Proton Motive Force Driven Atp Synthesis"               = "GO-BP: PMF-Driven ATP Synthesis"
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
      # sanitize cell-type names for filenames
      sanitize <- function(x) {
        x <- gsub("[^A-Za-z0-9]+", "_", x)
        x <- gsub("^_|_$", "", x)
        tolower(x)
      }
      
      out_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      for (ct in names(resp_auc_cell_type_list)) {
        df <- resp_auc_cell_type_list[[ct]]
        
        # build title and output path
        ttl <- paste0("AUC of Mitochondrial Respiratory Pathways in ", ct)
        outfile <- file.path(out_dir, paste0("auc_heatmap_", sanitize(ct), ".jpg"))
        
        # make heatmap object
        ht <- get_auc_heatmap(df, title = ttl, pathway_names = pathway_names,row_order = row_order)
        
        # open device, draw, close
        CairoJPEG(
          filename = outfile,
          width    = 15,
          height   = 10,
          units    = "in",
          dpi      = 600,
          quality  = 100
        )
        draw(ht, annotation_legend_list = list(legend_NA("NA (missing)")))
        dev.off()
        
        message("Saved: ", outfile)
      }
      
    }
    
    #auc_heatmap_4_cell_types.jpg
    {

      
      # order and panel letters
      ct_order  <- c("NKs", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
      panel_tag <- c("A", "B", "C", "D")
      
      ht_list <- NULL
      for (i in seq_along(ct_order)) {
        ct <- ct_order[i]
        df <- resp_auc_cell_type_list[[ct]]
        ttl <- paste0(panel_tag[i], ") AUC of Mitochondrial Respiratory Pathways in ", ct)
        
        ht_i <- get_auc_heatmap(
          df,
          title         = ttl,
          pathway_names = pathway_names,
          row_order     = row_order
        )
        
        ht_list <- c(ht_list, ht_i)
      }
      
      # Save the combined figure
      outfile <- file.path(out_dir, "auc_heatmap_4_cell_types.jpg")
      CairoJPEG(filename = outfile, width = 12, height = 20, units = "in", dpi = 600, quality = 100)
      draw(
        ht_list,
        merge_legend = TRUE,                   # keeps legends unified
        heatmap_legend_side = "right",
        annotation_legend_side = "right"
      )      dev.off()
      message("Saved: ", outfile)
    }
    
    
    #boxplots per cell type
    {
      
      # where to save
      # out_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress"
      out_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress/"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      # helper to make safe filenames
      sanitize <- function(x) {
        x <- gsub("[^A-Za-z0-9]+", "_", x)
        x <- gsub("^_|_$", "", x)
        tolower(x)
      }
      
      for (ct in names(resp_auc_cell_type_list)) {
        df_ct <- resp_auc_cell_type_list[[ct]]
        
        # title per cell type
        ttl <- paste0("AUC of Mitochondrial Respiratory Pathways in ", ct)
        
        # make plot
        p <- get_auc_boxplot(df = df_ct, title = ttl)
        
        # file per cell type
        outfile <- file.path(out_dir, paste0("auc_boxplot_", sanitize(ct), ".jpg"))
        
        # save
        ggsave(
          filename = outfile,
          plot     = p,
          width    = 30,
          height   = 20,
          units    = "in",
          dpi      = 600,
          device   = "jpeg",
          type     = "cairo",
          quality  = 100
        )
        
        message("Saved: ", outfile)
      }
      
      
      
      
    }
    
    # Boxplots for pathway (complexes II & III)
    {
      # Pathway mappings: raw_name -> short display name
      pathway_names <- c(
        "Wp Mitochondrial Complex Ii Assembly"                          = "WP: Complex II Assembly",
        "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "GO-BP:ETC in Complex II",
        "Wp Mitochondrial Complex Iii Assembly"                         = "WP: Complex III Assembly",
        "Kegg Medicus Reference Electron Transfer in Complex Iii"       = "KEGG-MR: ETC in Complex III"
      )
      
      # File names for saving (parallel order to pathway_names)
      file_names <- c(
        "Wp Mitochondrial Complex Ii Assembly"                          = "auc_box_plot_complex_ii_assembly",
        "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "auc_box_plot_complex_ii_etc",
        "Wp Mitochondrial Complex Iii Assembly"                         = "auc_box_plot_complex_iii_assembly",
        "Kegg Medicus Reference Electron Transfer in Complex Iii"       = "auc_box_plot_complex_iii_etc"
      )
      
      # Output directory
      out_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress/"
      
      failed <- character(0)
      
      for (raw_name in names(pathway_names)) {
        short_name <- pathway_names[raw_name]   # for plot title
        outfile    <- file.path(out_dir, paste0(file_names[raw_name], ".jpg"))  # for file name
        
        p <- try(
          get_auc_boxplot_for_pathway(
            auc_list     = resp_auc_cell_type_list,
            pathway_name = raw_name,
            title        = paste0(short_name, " across Cell Types")
          ),
          silent = TRUE
        )
        
        if (inherits(p, "try-error")) {
          failed <- c(failed, raw_name)
          next
        }
        
        ggsave(
          filename = outfile,
          plot     = p,
          width    = 30,
          height   = 20,
          units    = "in",
          dpi      = 600,
          device   = "jpeg",
          type     = "cairo",
          quality  = 100
        )
      }
      
      if (length(failed)) {
        message("Skipped (no rows / error): ", paste(failed, collapse = ", "))
      } else {
        message("Saved all boxplots to ", out_dir)
      }
    }
    
    
    #auc_4_box_plots
    {
      # Boxplots for pathway (complexes II & III)
      {
        # Pathway mappings: raw_name -> short display name
        pathway_names <- c(
          "Wp Mitochondrial Complex Ii Assembly"                          = "WP: Complex II Assembly",
          "Gobp Mitochondrial Electron Transport Succinate to Ubiquinone" = "GO-BP:ETC in Complex II",
          "Wp Mitochondrial Complex Iii Assembly"                         = "WP: Complex III Assembly",
          "Kegg Medicus Reference Electron Transfer in Complex Iii"       = "KEGG-MR: ETC in Complex III"
        )
        
        out_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress/"
        
        # Build the four plots (keep titles)
        p_list <- lapply(names(pathway_names), function(raw_name) {
          short_name <- pathway_names[raw_name]
          get_auc_boxplot_for_pathway(
            auc_list     = resp_auc_cell_type_list,
            pathway_name = raw_name,
            title        = paste0(short_name, " across Cell Types")
          )
        })
        
        # Name them for clarity
        p1 <- p_list[[1]]
        p2 <- p_list[[2]]
        p3 <- p_list[[3]]
        p4 <- p_list[[4]]
        
        # Arrange 2x2, collect a single legend, tag A–D at upper-right of each subplot
        auc_4_box_plots <-
          ((p1 + p2) / (p3 + p4)) +
          plot_layout(guides = "collect") +
          plot_annotation(tag_levels = "A") &
          theme(
            legend.position = "right",
            plot.tag.position = c(0.02, 0.98),       # upper-left corner
            plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 26)
          )
        
        # Save one combined figure
        ggsave(
          filename = file.path(out_dir, "auc_4_box_plots.jpg"),
          plot     = auc_4_box_plots,
          width    = 30,
          height   = 20,
          units    = "in",
          dpi      = 600,
          device   = "jpeg",
          type     = "cairo",
          quality  = 100
        )
        
        message("Saved combined plot to ", file.path(out_dir, "auc_4_box_plots.jpg"))
      }
      
    }
    
  }
  

  
  #TCA   
  {
    #create neuroimmunulogical_meta_table_lists for AUC scores
    {
      cell_paths <- list(
        "NKs"            = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/",
        "Monocytes"      = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/",
        "B Cells"        = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/",
        "CD4 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/",
        "CD8 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
      )
      
      # Map cell types to element names inside meta_table_list
      element_name_map <- c(
        "NKs"            = "meta_table_nk",
        "Monocytes"      = "meta_table_monocytes",
        "B Cells"        = "meta_table_b_cells",
        "CD4 + T cells"  = "meta_table_cd4_t_cells",
        "CD8 + T cells"  = "meta_table_cd8_t_cells"
      )
      
      # Helper: safe qread
      safe_qread <- function(path) {
        tryCatch(qread(path), error = function(e) NULL)
      }
      
      # ── Gather filenames per cell type ─────────────────────────────────────────────
      files_by_ct <- lapply(cell_paths, function(dir) {
        if (!dir.exists(dir)) return(character(0))
        f <- list.files(dir, pattern = "\\.qs$", full.names = TRUE)
        setNames(f, basename(f))  # names are basenames; values are full paths
      })
      
      # Union of sample filenames across all folders
      all_samples <- sort(unique(unlist(lapply(files_by_ct, names))))
      
      # ── Build the list-of-lists ───────────────────────────────────────────────────
      neuroimmunulogical_meta_table_lists <- list()
      
      for (sample_file in all_samples) {
        # Store under name without "_meta_table.qs"
        sample_key <- sub("_meta_table$", "", sub("\\.qs$", "", sample_file))
        
        meta_table_list <- list()
        for (ct in names(cell_paths)) {
          paths_vec <- files_by_ct[[ct]]
          if (length(paths_vec) && sample_file %in% names(paths_vec)) {
            df <- safe_qread(paths_vec[[sample_file]])
            if (!is.null(df)) {
              el_name <- element_name_map[[ct]]
              meta_table_list[[el_name]] <- df
            }
          }
        }
        
        if (length(meta_table_list) > 0) {
          neuroimmunulogical_meta_table_lists[[sample_key]] <- meta_table_list
        }
      }
      
      # Optional sanity check:
      # str(neuroimmunulogical_meta_table_lists, max.level = 2)
      
      # Optional save:
      # qsave(neuroimmunulogical_meta_table_lists,
      #       "/bigdata/bgutman/Collabs/Dror mevorach/data/neuroimmunulogical_meta_table_lists.qs")
    }
    
    
    #TCA AUC scores for each neuroimmunulogical condition
    {
      # ## ── Vector of pathway IDs ─────────────────────────────────────────────────────
      # vec_path_ids <-  c(
      #   "WP_THIAMINE_METABOLIC_PATHWAYS",
      #   "GOBP_PYRUVATE_METABOLIC_PROCESS",
      #   "GOMF_PYRUVATE_DEHYDROGENASE_ACETYL_TRANSFERRING_PHOSPHATASE_ACTIVITY",
      #   "REACTOME_MALATE_ASPARTATE_SHUTTLE",
      #   "GOBP_GLYCEROL_3_PHOSPHATE_METABOLIC_PROCESS",
      #   "GOBP_GLUTAMATE_METABOLIC_PROCESS",
      #   "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
      #   "KEGG_PROPANOATE_METABOLISM",
      #   "GOBP_SUCCINYL_COA_CATABOLIC_PROCESS",
      #   "GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE",
      #   "GOBP_UNSATURATED_FATTY_ACID_METABOLIC_PROCESS",
      #   "GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_OXIDASE",
      #   "GOBP_KETONE_METABOLIC_PROCESS",
      #   "GOMF_ACETYL_COA_C_ACYLTRANSFERASE_ACTIVITY"
      # )
      # 
      # ## ── Iterate and compute AUC per sample ───────────────────────────────────────
      # tca_auc_list <- list()
      # 
      # for (nm in names(neuroimmunulogical_meta_table_lists)) {
      #   that_meta_table_list <- neuroimmunulogical_meta_table_lists[[nm]]
      #   res <- tryCatch(
      #     get_meta_table_auc_score_for_many_cell_types_and_pathways(
      #       meta_table_list = that_meta_table_list,
      #       vec_path_ids    = vec_path_ids
      #     ),
      #     error = function(e) {
      #       message("AUC calc failed for ", nm, ": ", conditionMessage(e))
      #       NULL
      #     }
      #   )
      #   if (!is.null(res)) tca_auc_list[[nm]] <- res
      # }
      # 
      # qsave(tca_auc_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/oxidative_stress/neuroimmunulogical_tca_auc_list.qs")
      # 
      # 
      
    }
    
    
    #get tca_auc_cell_type_list
    {
      tca_auc_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/oxidative_stress/neuroimmunulogical_tca_auc_list.qs")
      
      ## Build a list of 5 data frames (one per cell type) from resp_auc_list
      ## Each row = pathway; columns = pathway, AUC (for that cell type), condition, dataset, sample
      
      cell_types <- c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells")
      
      resp_auc_cell_type_list <- setNames(vector("list", length(cell_types)), cell_types)
      
      for (ct in cell_types) {
        # Collect per-sample data frames for this cell type
        pieces <- lapply(names(tca_auc_list), function(nm) {
          mat <- tca_auc_list[[nm]]
          if (is.null(mat) || !is.data.frame(mat) || !(ct %in% colnames(mat))) return(NULL)
          
          # Parse condition / dataset / sample from the name: cond_dataset_sample(with underscores)
          parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
          condition <- if (length(parts) >= 1) parts[1] else NA_character_
          dataset   <- if (length(parts) >= 2) parts[2] else NA_character_
          sample    <- if (length(parts) >= 3) paste(parts[3:length(parts)], collapse = "_") else NA_character_
          
          # Build tidy block for this sample & cell type
          df <- data.frame(
            pathway   = rownames(mat),
            AUC       = as.numeric(mat[[ct]]),
            condition = condition,
            dataset   = dataset,
            sample    = sample,
            stringsAsFactors = FALSE
          )
          
          # (Optional) keep the found_genes column—uncomment if you want it
          # if ("found_genes" %in% colnames(mat)) {
          #   df$found_genes <- as.character(mat[["found_genes"]])
          # }
          
          df
        })
        
        pieces <- Filter(Negate(is.null), pieces)
        resp_auc_cell_type_list[[ct]] <- if (length(pieces)) do.call(rbind, pieces) else
          data.frame(pathway=character(), AUC=numeric(), condition=character(),
                     dataset=character(), sample=character(), stringsAsFactors=FALSE)
        
        rownames(resp_auc_cell_type_list[[ct]]) <- NULL
      }
      
      # Quick glance:
      # lapply(resp_auc_cell_type_list, head)
      
      #turn "CD55-patient" into "CD55_Patient"
      {
        # Fix sample naming inconsistency across all cell type dfs
        for (ct in names(resp_auc_cell_type_list)) {
          df <- resp_auc_cell_type_list[[ct]]
          
          # Replace only the exact match "CD55-patient"
          df$sample[df$sample == "CD55-patient"] <- "CD55_Patient"
          
          resp_auc_cell_type_list[[ct]] <- df
        }
        
        # Quick check
        lapply(resp_auc_cell_type_list, function(df) unique(df$sample))
      }
      
      
      
      # add meta data from information_df_neuroimmunulogical_samples
      {
        information_df_neuroimmunulogical_samples <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/information_df_neuroimmunulogical_samples.qs")
        
        cols_to_add <- c("condition", "condition_subtype", "age", "sex", "EDSS")
        
        # Loop through each cell type in the list and merge
        resp_auc_cell_type_list <- lapply(
          resp_auc_cell_type_list,
          function(df) {
            if ("condition" %in% names(df)) df$condition <- NULL  # drop old condition
            merge(
              df,
              information_df_neuroimmunulogical_samples[, c("sample", cols_to_add)],
              by = "sample",
              all.x = TRUE
            )
          }
        )
      }
      
      
      # Turn "CDMS" into "RRMS" under condition_subtype
      {
        for (ct in names(resp_auc_cell_type_list)) {
          df <- resp_auc_cell_type_list[[ct]]
          
          # Replace only the exact match "CDMS"
          df$condition_subtype[df$condition_subtype == "CDMS"] <- "RRMS"
          
          resp_auc_cell_type_list[[ct]] <- df
        }
        
        # Quick check
        lapply(resp_auc_cell_type_list, function(df) unique(df$condition_subtype))
      }
      
    }
    
    tca_auc_cell_type_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/oxidative_stress/tca_auc_cell_type_list.qs")
    
    #heatmaps
    {
      pathway_names <- c(
        "Wp Thiamine Metabolic Pathways" = "WP: Thiamine Metabolism",
        "Gobp Pyruvate Metabolic Process" = "GO-BP: Pyruvate Metabolism",
        "Gomf Pyruvate Dehydrogenase Acetyl Transferring Phosphatase Activity" = "GO-MF: PDH Phosphatase Activity",
        "Reactome Malate Aspartate Shuttle" = "Reactome: Malate-Aspartate Shuttle",
        "Gobp Glycerol 3 Phosphate Metabolic Process" = "GO-BP: Glycerol-3P Metabolism",
        "Gobp Glutamate Metabolic Process" = "GO-BP: Glutamate Metabolism",
        "Kegg Valine Leucine and Isoleucine Degradation" = "KEGG: BCAA (Val/Leu/Ile) Degradation",
        "Kegg Propanoate Metabolism" = "KEGG: Propanoate Metabolism",
        "Gobp Succinyl Coa Catabolic Process" = "GO-BP: Succinyl-CoA Catabolism",
        "Gobp Fatty Acid Beta Oxidation using Acyl Coa Dehydrogenase" = "GO-BP: FA β-Ox (ACAD)",
        "Gobp Unsaturated Fatty Acid Metabolic Process" = "GO-BP: Unsaturated FA Metabolism",
        # "Gobp Fatty Acid Beta Oxidation using Acyl Coa Oxidase" = "GO-BP: FA β-Ox (ACOX)",
        "Gobp Ketone Metabolic Process" = "GO-BP: Ketone Metabolism",
        "Gomf Acetyl Coa c Acyltransferase Activity" = "GO-MF: Thiolase Activity"
        # "Wp Oxidative Damage Response" = "WP: Oxidative Damage Response",
      )
      

      row_order <- c(
        "WP: Thiamine Metabolism",
        "GO-BP: Pyruvate Metabolism",
        "GO-MF: PDH Phosphatase Activity",
        "GO-BP: Glutamate Metabolism",
        "GO-BP: Ketone Metabolism",
        "GO-BP: Unsaturated FA Metabolism",
        "Reactome: Malate-Aspartate Shuttle",
        "GO-BP: Glycerol-3P Metabolism",
        "KEGG: BCAA (Val/Leu/Ile) Degradation",
        "KEGG: Propanoate Metabolism",
        "GO-BP: Succinyl-CoA Catabolism",
        "GO-BP: FA β-Ox (ACAD)",
        "GO-MF: Thiolase Activity"
        # "GO-BP: FA β-Ox (ACOX)",
      )
      # sanitize cell-type names for filenames
      sanitize <- function(x) {
        x <- gsub("[^A-Za-z0-9]+", "_", x)
        x <- gsub("^_|_$", "", x)
        tolower(x)
      }
      
      out_dir <- "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/oxidative_stress"
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      for (ct in names(tca_auc_cell_type_list)) {
        df <- tca_auc_cell_type_list[[ct]]
        
        # build title and output path
        ttl <- paste0("AUC of TCA Pathways in ", ct)
        outfile <- file.path(out_dir, paste0("tca_heatmap_", sanitize(ct), ".jpg"))
        
        # make heatmap object
        ht <- get_auc_heatmap(df, title = ttl, pathway_names = pathway_names,row_order = row_order)
        
        # open device, draw, close
        CairoJPEG(
          filename = outfile,
          width    = 15,
          height   = 10,
          units    = "in",
          dpi      = 600,
          quality  = 100
        )
        draw(ht, annotation_legend_list = list(legend_NA("NA (missing)")))
        dev.off()
        
        message("Saved: ", outfile)
      }
      
    }

  }


}



# Shared DEGs signatures AUC Scores
{
  
  # creating shared_degs_signature_list
  {
    #Shared nDEGs signatures
    {
      #getting sig_gene_hugo_long: 
      {


        #subsetting to only shred DEGs
        {
          # positive_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
          # negative_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")
          # 
          # 
          # # Subset genes with TRUE in all cell types
          # positive_all_true <- positive_gene_upset_df %>%
          #   filter(if_all(c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"), ~ .x == TRUE))
          # 
          # negative_all_true <- negative_gene_upset_df %>%
          #   filter(if_all(c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"), ~ .x == TRUE))



        }


        #ORA pathways for the signature
        {

# 
#           # --- NF-κB/Immune-Activation -----------------------------------------------
#           # HALLMARK_TNFA_SIGNALING_VIA_NFKB
#           # immune response-activating signaling pathway
#           # immune response-regulating signaling pathway
#           # regulation of inflammatory response
#           # canonical NF-kappaB signal transduction
#           # HALLMARK_INFLAMMATORY_RESPONSE
#           # immune response-activating cell surface receptor signaling pathway
#           # immune response-regulating cell surface receptor signaling pathway
#           # regulation of canonical NF-kappaB signal transduction
#           # leukocyte activation involved in immune response
#           # cell activation involved in immune response
#           # innate immune response-activating signaling pathway
#           # HALLMARK_TGF_BETA_SIGNALING
#           # negative regulation of NF-kappaB transcription factor activity
#           # chronic inflammatory response
#           sig_nfkB_immune_activation <- c(
#             "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#             # "GO:0002757",
#             # "GO:0002764",
#             # "GO:0050727",
#             # "GO:0007249",
#             "HALLMARK_INFLAMMATORY_RESPONSE"
#             # "GO:0002429",
#             # "GO:0002768",
#             # "GO:0043122",
#             # "GO:0002366",
#             # "GO:0002263",
#             # "GO:0002758",
#             # "HALLMARK_TGF_BETA_SIGNALING",
#             # "GO:0032088",
#             # "GO:0002544"
#           )
# 
#           # --- Metabolic Adaptation ---------------------------------------------------
#           # Cellular response to starvation
#           # Response of EIF2AK4 (GCN2) to amino acid deficiency
#           # sig_metabolic_adaptation <- c(
#           #   "R-HSA-9711097",
#           #   "R-HSA-9633012"
#           # )
# 
#           # --- Seleno-Antioxidant Response -------------------------------------------
#           # Selenocysteine synthesis
#           # Selenoamino acid metabolism
#           sig_seleno_antioxidant <- c(
#             "R-HSA-2408557",
#             "R-HSA-2408522"
#           )
# 
#           # --- Extrinsic Apoptosis Regulation ----------------------------------------
#           # extrinsic apoptotic signaling pathway
#           # extrinsic apoptotic signaling pathway via death domain receptors
#           sig_extrinsic_apoptosis <- c(
#             # "GO:0097191",
#             "GO:0008625"
#           )
# 
# 
# 
#           # (Optional) Put them together if helpful:
#           signatures <- list(
#             NFkB_Immune_Activation = sig_nfkB_immune_activation,
#             # Metabolic_Adaptation    = sig_metabolic_adaptation,
#             Seleno_Antioxidant      = sig_seleno_antioxidant,
#             Extrinsic_Apoptosis     = sig_extrinsic_apoptosis
#           )


        }


        #creating sig_gene_long
        {
          # positive_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/positive_ora_df.qs")
          # positive_ora_df <- positive_ora_df$meta_ora_top
          # negative_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/negative_ora_df.qs")
          # negative_ora_df <- negative_ora_df$meta_ora_top
          # 
          # 
          # # --- Example usage ------------------------------------------------------------
          # res <- extract_signature_genes(negative_ora_df, signatures)
          # 
          # # List of unique gene tokens per signature:
          # genes_by_signature <- res$genes_by_signature
          # 
          # # Tidy long DF (two columns: signature, gene):
          # sig_gene_long <- res$tidy
          # 
          # # Peek:
          # lapply(genes_by_signature, function(v) head(v, 10))
          # head(sig_gene_long)

        }

        #creating sig_gene_hugo_long
        {

# 
#           # --- Example usage ------------------------------------------------------------
#           # Make sure org.Hs.eg.db is available for your helper:
#           if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#             stop("Please install org.Hs.eg.db: BiocManager::install('org.Hs.eg.db')")
#           }
#           if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
#             stop("Please install AnnotationDbi: BiocManager::install('AnnotationDbi')")
#           }
# 
#           hugo_res <- to_hugo_by_signature(negative_ora_df, signatures)
# 
#           # Per-signature unique HUGO symbols (duplicates allowed BETWEEN signatures)
#           genes_by_signature_hugo <- hugo_res$genes_by_signature_hugo
# 
#           # Tidy long DF: signature | gene
#           sig_gene_hugo_long <- hugo_res$tidy_hugo
# 
#           # Quick peek
#           lapply(genes_by_signature_hugo, head)
#           head(sig_gene_hugo_long)



        }

        #adding metabolisim genes based on general reactome categories
        {

          #
          # # --- libraries ---------------------------------------------------------------
          # # (Assumes dplyr, tibble, stringr are installed)
          # # library(dplyr)
          # # library(tibble)
          # # library(stringr)
          #
          # # --- helpers -----------------------------------------------------------------
          # get_sig_genes <- function(df, sig) {
          #   df %>%
          #     dplyr::filter(signature == sig) %>%
          #     dplyr::pull(gene) %>%
          #     unique() %>%
          #     trimws()
          # }
          #
          # add_genes_to_signature <- function(sig_df, sig_name, new_gene_vec) {
          #   new_gene_vec <- unique(trimws(new_gene_vec))
          #   if (length(new_gene_vec) == 0) return(sig_df)
          #
          #   added_block <- tibble::tibble(signature = sig_name, gene = new_gene_vec)
          #   sig_df %>%
          #     dplyr::bind_rows(added_block) %>%
          #     dplyr::distinct(signature, gene, .keep_all = TRUE) # no dups within signature
          # }
          #
          #
          # metab_genes <- negative_all_true %>%
          #   dplyr::filter(group_reactome_copy %in% c("Metabolism", "Metabolism of proteins"
          #   )) %>%
          #   dplyr::pull(gene) %>%
          #   unique() %>%
          #   trimws() %>%
          #   stats::na.omit() %>%
          #   as.character()
          #
          # existing_metab <- get_sig_genes(sig_gene_hugo_long, "Metabolic_Adaptation")
          # to_add_metab <- setdiff(metab_genes, existing_metab)
          #
          # cat("Metabolism/RNA-Proteins → Metabolic_Adaptation: adding", length(to_add_metab), "new genes\n")
          # if (length(to_add_metab)) print(to_add_metab)
          #
          # sig_gene_hugo_long <- add_genes_to_signature(
          #   sig_gene_hugo_long, "Metabolic_Adaptation", to_add_metab
          # )
          #
          # # --- 3) sanity: uniqueness within each signature -----------------------------
          # sig_gene_hugo_long <- sig_gene_hugo_long %>%
          #   dplyr::mutate(signature = trimws(signature), gene = trimws(gene)) %>%
          #   dplyr::distinct(signature, gene, .keep_all = TRUE)
          #
          # # --- Optional quick summary ---------------------------------------------------
          # summary_tbl <- sig_gene_hugo_long %>%
          #   dplyr::count(signature, name = "n_genes") %>%
          #   dplyr::arrange(dplyr::desc(n_genes)) %>%
          #   tibble::as_tibble()
          #
          # # tibble::print(summary_tbl, n = Inf)
        }



# 
#         my_genes <- toupper(trimws(unique(negative_all_true$gene)))
# 
# 
# 
#         # ====================================================================== #
#         # 1) C2 ----------------------------------------------------------------- #
#         overlaps_C2 <- overlaps_for_category("C2", my_genes)
# 
#         map_C2 <- tibble(
#           gs_name = c(
#             "GRAESSMANN_APOPTOSIS_BY_DOXORUBICIN_UP",
#             "CONCANNON_APOPTOSIS_BY_EPOXOMICIN_UP",
#             "KRIEG_HYPOXIA_NOT_VIA_KDM3A"
#             # "REACTOME_RNA_POLYMERASE_II_TRANSCRIPTION"
#           ),
#           signature = c(
#             "Intrinsic_Apoptosis",
#             "Intrinsic_Apoptosis",
#             "Hypoxia"
#             # "Transcriptional_Machinery"
#           )
#         )
# 
#         sig_gene_hugo_long <- add_mappings(overlaps_C2, map_C2, sig_gene_hugo_long)
# 
#         # 2) C5 ----------------------------------------------------------------- #
#         overlaps_C5 <- overlaps_for_category("C5", my_genes)
# 
#         map_C5 <- tibble(
#           gs_name = c(
#             "GOBP_INTRACELLULAR_TRANSPORT",
#             "GOBP_PROTEIN_TRANSPORT",
#             "GOBP_INTRACELLULAR_PROTEIN_TRANSPORT",
#             "GOBP_VESICLE_MEDIATED_TRANSPORT",
#             "GOCC_MITOCHONDRION"
#           ),
#           signature = c(rep("Endomembrane_Trafficking", 4), "Mito_Function_Maintenance")
#         )
# 
#         sig_gene_hugo_long <- add_mappings(overlaps_C5, map_C5, sig_gene_hugo_long)
# 
#         # 3) H (Hallmark) ------------------------------------------------------- #
#         overlaps_H <- overlaps_for_category("H", my_genes)
# 
#         map_H <- tibble(
#           gs_name   = "HALLMARK_HYPOXIA",
#           signature = "Hypoxia"
#         )
# 
#         sig_gene_hugo_long <- add_mappings(overlaps_H, map_H, sig_gene_hugo_long)
# 
#         # 4) Antioxidant modules (multi-collection keys) ------------------------ #
#         set_keys_antiox <- c(
#           "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS","KEGG_GLUTATHIONE_METABOLISM",
#           "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS",
#           "GOBP_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS","GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
#           "GOBP_GLUTATHIONE_CATABOLIC_PROCESS","GOBP_GLUTATHIONE_METABOLIC_PROCESS",
#           "KEGG_GLUTATHIONE_METABOLISM","KEGG_MEDICUS_REFERENCE_GLUTATHIONE_BIOSYNTHESIS",
#           "REACTOME_GLUTATHIONE_CONJUGATION","REACTOME_GLUTATHIONE_SYNTHESIS_AND_RECYCLING",
#           "WP_GLUTATHIONE_METABOLISM","MORF_ATRX","BIOCARTA_ARENRF2_PATHWAY",
#           "IBRAHIM_NRF2_DOWN","IBRAHIM_NRF2_UP","KEGG_MEDICUS_REFERENCE_KEAP1_NRF2_SIGNALING_PATHWAY",
#           "WP_NRF2_PATHWAY","WP_NRF2ARE_REGULATION"
#         )
# 
#         overlaps_antiox <- compute_overlaps_for_sets(my_genes, set_keys_antiox)
# 
#         map_antiox <- tibble(
#           gs_name   = c("WP_NRF2_PATHWAY","IBRAHIM_NRF2_DOWN","WP_GLUTATHIONE_METABOLISM"),
#           signature = "Seleno_Antioxidant"
#         )
# 
#         sig_gene_hugo_long <- add_mappings(overlaps_antiox, map_antiox, sig_gene_hugo_long) |>
#           mutate(signature = dplyr::case_when(
#             signature == "Seleno_Antioxidant"        ~ "Antioxidant_response",
#             # signature == "Transcriptional_Machinery" ~ "Transcriptional_Control",
#             TRUE ~ signature
#           ))
# 
#         # 5) Transcriptional_Control (multi-collection keys) -------------------- #
#         # overlaps_tx <- compute_overlaps_for_sets(
#         #   my_genes,
#         #   c("GOMF_SEQUENCE_SPECIFIC_DNA_BINDING","GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY")
#         # )
#         #
#         # map_tx <- tibble(
#         #   gs_name   = c("GOMF_SEQUENCE_SPECIFIC_DNA_BINDING","GOMF_TRANSCRIPTION_REGULATOR_ACTIVITY"),
#         #   signature = "Transcriptional_Control"
#         # )
#         #
#         # sig_gene_hugo_long <- add_mappings(overlaps_tx, map_tx, sig_gene_hugo_long)
# 
#         # (optional) quick counts
#         table(sig_gene_hugo_long$signature)
# 
# 
#         # Ensure character columns
#         sig_gene_hugo_long$signature <- as.character(sig_gene_hugo_long$signature)
#         sig_gene_hugo_long$gene      <- as.character(sig_gene_hugo_long$gene)
# # # 
# # #         ## 1) Merge Metabolic_Adaptation into Transcriptional_Control
# # #         sig_gene_hugo_long$signature[
# # #           sig_gene_hugo_long$signature == "Metabolic_Adaptation"
# # #         ] <- "Transcriptional_Control"
# # # 
# # #         ## 2) (Optional) Rename the merged bucket to include ISR/translation
# # #         sig_gene_hugo_long$signature[
# # #           sig_gene_hugo_long$signature == "Transcriptional_Control"
# # #         ] <- "Transcription/Translation_Control"
# # # 
# # #         ## 3) (Optional) Tighten the mito label
# # #         sig_gene_hugo_long$signature[
# # #           sig_gene_hugo_long$signature == "Mito_Function_Maintenance"
# # #         ] <- "Mitochondrial_Homeostasis"
# # 
# #         ## 4) Ensure uniqueness within signature after merges/renames
# #         sig_gene_hugo_long <- unique(sig_gene_hugo_long)
# # 
# #         ## 5) Quick check
# #         table(sig_gene_hugo_long$signature)
# 

      }
      
      #removing some signatures, observing gene overlap between the rest:
      {
        # 
        # sig_gene_hugo_long <- sig_gene_hugo_long[
        #   !sig_gene_hugo_long$signature %in% c("Antioxidant_response", "Intrinsic_Apoptosis", "Hypoxia"),
        # ]
        # 
        # 
        # # ensure uniqueness & clean casing
        # sg <- sig_gene_hugo_long |>
        #   mutate(signature = as.character(signature),
        #          gene = toupper(trimws(gene))) |>
        #   distinct(signature, gene)
        # 
        # # list of genes per signature
        # sig_list <- sg |>
        #   group_by(signature) |>
        #   summarise(genes = list(sort(unique(gene))), n = length(unique(gene)), .groups = "drop")
        # 
        # sig_names <- sig_list$signature
        # gene_sets <- sig_list$genes
        # sizes     <- sig_list$n
        # names(gene_sets) <- sig_names
        # names(sizes)     <- sig_names
        # 
        # # ---------- (A) Long/tidy pairwise overlap table ----------
        # pairwise_tbl <- {
        #   idx <- combn(sig_names, 2, simplify = FALSE)
        #   tibble(
        #     sig1 = map_chr(idx, 1),
        #     sig2 = map_chr(idx, 2)
        #   ) |>
        #     mutate(
        #       n_sig1   = sizes[sig1],
        #       n_sig2   = sizes[sig2],
        #       shared   = map2(sig1, sig2, ~ intersect(gene_sets[[.x]], gene_sets[[.y]])),
        #       n_shared = map_int(shared, length),
        #       union    = map2(sig1, sig2, ~ union(gene_sets[[.x]], gene_sets[[.y]])),
        #       n_union  = map_int(union, length),
        #       jaccard  = if_else(n_union > 0, n_shared / n_union, NA_real_),
        #       overlap_coeff = if_else(pmin(n_sig1, n_sig2) > 0, n_shared / pmin(n_sig1, n_sig2), NA_real_),
        #       genes_shared = map_chr(shared, ~ paste(.x, collapse = ", "))
        #     ) |>
        #     arrange(desc(n_shared), desc(jaccard), sig1, sig2) |>
        #     dplyr::select(sig1, sig2, n_sig1, n_sig2, n_shared, jaccard, overlap_coeff, genes_shared)
        # }
        # 
        # # Example: write to TSV (without the large gene list, if you prefer)
        # # write.table(select(pairwise_tbl, -genes_shared), "signature_overlaps.tsv",
        # #             sep = "\t", row.names = FALSE, quote = FALSE)
        # 
        # # ---------- (B) Symmetric overlap matrix (counts) ----------
        # overlap_mat <- {
        #   m <- outer(sig_names, sig_names, Vectorize(function(a,b) length(intersect(gene_sets[[a]], gene_sets[[b]]))))
        #   dimnames(m) <- list(sig_names, sig_names)
        #   as.data.frame(m)  # nice to view; or keep as matrix 'm'
        # }
        # 
        # # Optional: Jaccard matrix instead of counts
        # jaccard_mat <- {
        #   U <- outer(sig_names, sig_names, Vectorize(function(a,b) length(union(gene_sets[[a]], gene_sets[[b]]))))
        #   J <- ifelse(U > 0, overlap_mat / U, NA_real_)
        #   as.data.frame(J)
        # }
        # 
        # # Quick peek
        # # head(pairwise_tbl, 10)
        # # overlap_mat
        # # jaccard_mat
        # 

      }
      
      #powermatching, creating sig_gene_hugo_long_sig
      {
        
        # sig_gene_hugo_long <- power_match_sig_list(sig_gene_hugo_long)
        # 
        # sig_gene_hugo_long_negative <- sig_gene_hugo_long

      }
      
    }
    
    
    #Shared pDEGs signature:
    {
      #subsetting to only shred DEGs
      {
        positive_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/positive_gene_upset_df.qs")
        negative_gene_upset_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/upset/negative_gene_upset_df.qs")


        # Subset genes with TRUE in all cell types
        positive_all_true <- positive_gene_upset_df %>%
          filter(if_all(c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"), ~ .x == TRUE))

        negative_all_true <- negative_gene_upset_df %>%
          filter(if_all(c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"), ~ .x == TRUE))



      }

    
    positive_ora_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/ORA/positive_ora_df.qs")
    positive_ora_df <- positive_ora_df$meta_ora_top
    
    
    my_genes <- toupper(trimws(unique(positive_all_true$gene)))
    
    sig_gene_hugo_long <- data.frame(signature = character(), gene = character(), stringsAsFactors = FALSE)
    
    overlaps_all <- overlaps_for_category(genes = my_genes)
    
    
    
    map_all <- tibble(
      gs_name = c(
        # "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
        # "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
        # "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
        # "REACTOME_METABOLISM_OF_RNA",
        # "REACTOME_SELENOAMINO_ACID_METABOLISM",
        # "REACTOME_METABOLISM_OF_LIPIDS",
        # "REACTOME_GLYCOSAMINOGLYCAN_METABOLISM",
        # "REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM",
        "GRAESSMANN_APOPTOSIS_BY_DOXORUBICIN_DN",
        "HOLLMANN_APOPTOSIS_VIA_CD40_DN"
        # "PID_RETINOIC_ACID_PATHWAY",
    
        # "GOBP_LIPID_METABOLIC_PROCESS",
        # "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
        # "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
        # "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
        # "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
        # "GOBP_HYDROGEN_PEROXIDE_METABOLIC_PROCESS"

    
      ),
      signature = c(
        # "Amino-Starvation ISR",
        # "Amino-Starvation ISR",
        # "Amino-Starvation ISR",
        # "Amino-Starvation ISR",
        # "MitoROS_Module",
        # "MitoROS_Module",
        # "MitoROS_Module",
        # "MitoROS_Module",
        "Apoptosis_Resistance",
        "Apoptosis_Resistance"
        # "MitoROS_Module",
        # 
        # 
        # "MitoROS_Module",
        # "MitoROS_Module",
        # "MitoROS_Module",
        # "MitoROS_Module",
        # "MitoROS_Module",
        # "MitoROS_Module"

      )
    )
    
    
    sig_gene_hugo_long <- add_mappings(overlaps_all, map_all, sig_gene_hugo_long)
    
    
    # sig_gene_hugo_long <- sig_gene_hugo_long[
    #   !(sig_gene_hugo_long$signature %in% c("Apoptosis_Resistance")),
    #   , drop = FALSE
    # ]
    
    
    
    # try <- find_msig_sets_with_genes(genes = c("PSEN1", "RARB", "RBP1", "RBP4", "RARA"),
    #                                  require_all = FALSE)
    
    
    
    
    
    
    sig_gene_hugo_long_positive <- power_match_sig_list(sig_gene_hugo_long)


      
    }
    
    # creating and saving the sig_list
    {

      # sig_list <- list(
      #   positive_sigs = sig_gene_hugo_long_positive,
      #   negative_sigs = sig_gene_hugo_long_negative
      # )
      # 
      # names(sig_list$positive_sigs)[names(sig_list$positive_sigs) == "Amino-Starvation ISR"] <- "Amino_Starvation_ISR"
      # 
      # qsave(sig_list, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_genes_signatures/shared_degs_signature_list.qs")


    }
    
  }

  # signature_composition_barplot
  {
    sig_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_genes_signatures/shared_degs_signature_list.qs")
    
    signature_composition_barplot <- get_deg_signature_composition_plot(sig_list = sig_list)
    

    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/shared_degs_signatures/signature_composition_barplot.jpg",
    #   plot     = signature_composition_barplot,
    #   width    = 13,
    #   height   = 10,
    #   units    = "in",
    #   dpi      = 600,
    #   device   = "jpeg",
    #   type     = "cairo",         # <– forces Cairo backend
    #   quality  = 100
    # )
  }
  
  # creating sig_score_long_df
  {
    
    sig_list <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_genes_signatures/shared_degs_signature_list.qs")
    
    
    
    #create neuroimmunulogical_meta_table_lists for AUC scores
    {
      cell_paths <- list(
        "NKs"            = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/",
        "Monocytes"      = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/",
        "B Cells"        = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/",
        "CD4 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/",
        "CD8 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
      )
      
      # Map cell types to element names inside meta_table_list
      element_name_map <- c(
        "NKs"            = "meta_table_nk",
        "Monocytes"      = "meta_table_monocytes",
        "B Cells"        = "meta_table_b_cells",
        "CD4 + T cells"  = "meta_table_cd4_t_cells",
        "CD8 + T cells"  = "meta_table_cd8_t_cells"
      )
      
      # Helper: safe qread
      safe_qread <- function(path) {
        tryCatch(qread(path), error = function(e) NULL)
      }
      
      # ── Gather filenames per cell type ─────────────────────────────────────────────
      files_by_ct <- lapply(cell_paths, function(dir) {
        if (!dir.exists(dir)) return(character(0))
        f <- list.files(dir, pattern = "\\.qs$", full.names = TRUE)
        setNames(f, basename(f))  # names are basenames; values are full paths
      })
      
      # Union of sample filenames across all folders
      all_samples <- sort(unique(unlist(lapply(files_by_ct, names))))
      
      # ── Build the list-of-lists ───────────────────────────────────────────────────
      neuroimmunulogical_meta_table_lists <- list()
      
      for (sample_file in all_samples) {
        # Store under name without "_meta_table.qs"
        sample_key <- sub("_meta_table$", "", sub("\\.qs$", "", sample_file))
        
        meta_table_list <- list()
        for (ct in names(cell_paths)) {
          paths_vec <- files_by_ct[[ct]]
          if (length(paths_vec) && sample_file %in% names(paths_vec)) {
            df <- safe_qread(paths_vec[[sample_file]])
            if (!is.null(df)) {
              el_name <- element_name_map[[ct]]
              meta_table_list[[el_name]] <- df
            }
          }
        }
        
        if (length(meta_table_list) > 0) {
          neuroimmunulogical_meta_table_lists[[sample_key]] <- meta_table_list
        }
      }
      
      #remove CD55
      {
        neuroimmunulogical_meta_table_lists <- neuroimmunulogical_meta_table_lists[
          names(neuroimmunulogical_meta_table_lists) != "CD55.Def_Mevorach_CD55-patient"
        ]
      }
    }
    
    
    
    #create sig_score_results
    {
      
      # # ---- Setup -------------------------------------------------------------------
      # save_path <- "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_degs_signature/sig_score_results.qs"
      # 
      # # Ensure qs is available
      # if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required.")
      # 
      # # Helper for timestamped logs
      # ts_msg <- function(...) {
      #   msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""))
      #   cat(msg, "\n")
      #   flush.console()
      # }
      # 
      # # ---- Initialize nested results container -------------------------------------
      # sig_score_results <- setNames(vector("list", length(sig_list)), names(sig_list))
      # 
      # # Pre-allocate second level to keep structure stable from the start
      # for (sig_type in names(sig_list)) {
      #   sig_score_results[[sig_type]] <- setNames(
      #     vector("list", length(neuroimmunulogical_meta_table_lists)),
      #     names(neuroimmunulogical_meta_table_lists)
      #   )
      # }
      # 
      # # Initial save of empty structure for safety
      # qs::qsave(sig_score_results, save_path)
      # ts_msg("Initialized empty results and saved to: ", save_path)
      # 
      # # ---- Main loops --------------------------------------------------------------
      # for (sig_type in names(sig_list)) {
      #   vec_paths <- sig_list[[sig_type]]
      #   ts_msg("Starting signature set: '", sig_type, "' (", length(vec_paths), " pathways)")
      # 
      #   for (mt_name in names(neuroimmunulogical_meta_table_lists)) {
      #     ts_msg("  Processing meta_table: '", mt_name, "' for sig_type: '", sig_type, "'")
      # 
      #     mt <- neuroimmunulogical_meta_table_lists[[mt_name]]
      # 
      #     # Compute with safety net
      #     res <- tryCatch(
      #       {
      #         get_meta_table_auc_score_for_many_cell_types_and_pathways(
      #           meta_table_list = mt,
      #           vec_path_ids    = vec_paths,
      #           filter_nocoding = TRUE,
      #           is_msigdbr      = FALSE
      #         )
      #       },
      #       error = function(e) {
      #         ts_msg("  ERROR in '", sig_type, "' / '", mt_name, "': ", conditionMessage(e))
      #         # Store the error so you can inspect later, but do not crash the loop
      #         list(error = TRUE, message = conditionMessage(e))
      #       }
      #     )
      # 
      #     # Write result immediately into the nested structure
      #     sig_score_results[[sig_type]][[mt_name]] <- res
      # 
      #     # Save after each update for crash-safety
      #     qs::qsave(sig_score_results, save_path)
      #     ts_msg("  Saved progress after '", sig_type, "' / '", mt_name, "'")
      #   }
      # 
      #   ts_msg("Completed signature set: '", sig_type, "'")
      # }
      # 
      # # ---- Final save --------------------------------------------------------------
      # qs::qsave(sig_score_results, save_path)
      # ts_msg("All done. Final results saved to: ", save_path)
      
    }
    
    
    #create sig_score_long_df
    {
      # sig_score_results <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_degs_signature/sig_score_results.qs")
      # 
      # # Flatten sig_score_results -> long data frame
      # sig_score_long_df <- imap_dfr(sig_score_results, function(sample_list, sig_type_name) {
      #   imap_dfr(sample_list, function(mat, full_name) {
      #     df <- as.data.frame(mat, stringsAsFactors = FALSE)
      #     df$sig_name <- rownames(df)
      # 
      #     has_fg   <- "found_genes" %in% names(df)
      #     cell_cols <- setdiff(names(df), c("sig_name", if (has_fg) "found_genes"))
      # 
      #     out <- df %>%
      #       pivot_longer(
      #         cols = all_of(cell_cols),
      #         names_to = "cell_type",
      #         values_to = "sig_score"
      #       ) %>%
      #       mutate(
      #         condition = sub("_.*$", "", full_name),
      #         dataset   = sub("^[^_]*_(.*?)_.*$", "\\1", full_name),
      #         sample    = sub("^[^_]*_[^_]*_(.*)$", "\\1", full_name),
      #         sig_type  = ifelse(sig_type_name == "positive_sigs", "positive", "negative")
      #       )
      # 
      #     if (!has_fg) out$found_genes <- NA_character_
      # 
      #     out %>%
      #       dplyr::select(cell_type, sample, condition, dataset, sig_name, sig_type, sig_score, found_genes)
      #   })
      # })
      # 
      # qsave(sig_score_long_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_degs_signature/sig_score_long_df.qs")
    }
    
    
    #create sig_condition_medians (for inspection and observation)
    {
      # sig_condition_medians <- sig_score_long_df %>%
      #   group_by(cell_type, sig_name, sig_type, condition) %>%
      #   summarise(median = median(sig_score, na.rm = TRUE), .groups = "drop") %>%
      #   mutate(condition = recode(condition,
      #                             "CD59.Def"        = "CD59.Def_median",
      #                             "MS"              = "MS_median",
      #                             "Ab.mediated.IDD" = "Ab.mediated.IDD_median")) %>%
      #   pivot_wider(
      #     names_from  = condition,
      #     values_from = median
      #   ) %>%
      #   # ---- ordering --------------------------------------------------------------
      # mutate(
      #   sig_type  = factor(sig_type, levels = c("positive", "negative")),
      #   cell_type = factor(cell_type,
      #                      levels = c("NKs", "Monocytes", "B Cells", "CD4+ T Cells", "CD8+ T Cells"))
      # ) %>%
      #   arrange(sig_type, sig_name, cell_type) %>%  # positive first; then by signature; then by cell type order
      #   dplyr::select(cell_type, sig_name, sig_type,
      #                 CD59.Def_median, MS_median, Ab.mediated.IDD_median)
      # 
      
    }
    
  }
  
  # creating shared_degs_signatures_boxplots
  {
    sig_score_long_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_degs_signature/sig_score_long_df.qs")
    
    shared_degs_signatures_boxplots <- plot_shared_degs_signatures_boxplots(df = sig_score_long_df)
    

    
    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/shared_degs_signatures/shared_degs_signatures_boxplots_across_conditions.jpg",
    #   plot     = shared_degs_signatures_boxplots,
    #   width    = 25,
    #   height   = 20,
    #   units    = "in",
    #   dpi      = 600,
    #   device   = "jpeg",
    #   type     = "cairo",         # <– forces Cairo backend
    #   quality  = 100
    # )
  }
  
  # creating scatter plots
  {
    
    
    sig_score_long_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/shared_degs_signature/sig_score_long_df.qs")
    
    
    
    
    p <- auc_scores_scatter_across_condtions(df = sig_score_long_df,    x_sig = "NFkB_Immune_Activation",
                                     y_sig = "Apoptosis_Resistance",
                                     axis_labels = c("NF-κB Pathway (AUC)", "Anti-apoptotic Signaling (AUC)"), 
                                     title  = "NF-κB Pathway vs. Anti-apoptotic Signaling AUC Scores Across Immune Cell Types")
    
# 
#     ggsave(
#       filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/shared_degs_signatures/NFkB_Apoptosis_Resistance_scatter.jpg",
#       plot     = p,
#       width    = 25,
#       height   = 20,
#       units    = "in",
#       dpi      = 600,
#       device   = "jpeg",
#       type     = "cairo",         # <– forces Cairo backend
#       quality  = 100
#     )
    
    
   
    
    p <- auc_scores_scatter_across_condtions(df = sig_score_long_df,
                                     x_sig = "Extrinsic_Apoptosis", y_sig = "Apoptosis_Resistance",
                                     axis_labels = c(  "Extrinsic Apoptosis (AUC)", "Anti-apoptotic Signaling (AUC)"), 
                                     title  =     "Extrinsic Apoptosis Signaling vs. Anti-apoptotic Signaling AUC Scores Across Immune Cell Types")
    
    
    # ggsave(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/shared_degs_signatures/Extrinsic_Apoptosis_Apoptosis_Resistance_scatter.jpg",
    #   plot     = p,
    #   width    = 25,
    #   height   = 20,
    #   units    = "in",
    #   dpi      = 600,
    #   device   = "jpeg",
    #   type     = "cairo",         # <– forces Cairo backend
    #   quality  = 100
    # )
    # 
    
    

    
    
    
  }


}

# heatmaps of whole gene sets
{
  #gene_sets_df: all genes from all gene sets (with duplicates)
  {
    msig <- msigdbr(species = "Homo sapiens")
    
    
    gene_sets <- c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      
      "WP_NRF2_PATHWAY",
      # "WP_GLUTATHIONE_METABOLISM",
      # "REACTOME_SELENOAMINO_ACID_METABOLISM",  # R-HSA-2408522
      
      # "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
      "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
      "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
      # "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
      
      # "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
      # "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      
      "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
      "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
      "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
      "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY"
    )
    
    
    # Build (gene, gene_set) table from your gene_sets using the preloaded `msig`
    # Assumes: `gene_sets` (character vector) and `msig` (msigdbr table) already exist,
    # and your function get_msigdb_genes_from_gene_set_name() is defined.
    {
      stopifnot(is.character(gene_sets), length(gene_sets) > 0)
      
      rows_list <- vector("list", length(gene_sets))
      not_found <- character(0)
      
      for (i in seq_along(gene_sets)) {
        gs <- gene_sets[i]
        message("Fetching genes for: ", gs)
        
        genes <- tryCatch(
          get_msigdb_genes_from_gene_set_name(gs_name = gs, msig = msig),
          error = function(e) {
            warning("Error fetching '", gs, "': ", conditionMessage(e))
            character(0)
          }
        )
        
        if (length(genes) == 0L) {
          not_found <- c(not_found, gs)
          rows_list[[i]] <- NULL
        } else {
          rows_list[[i]] <- data.frame(
            gene     = genes,
            gene_set = gs,
            stringsAsFactors = FALSE
          )
        }
      }
      
      gene_sets_df <- do.call(rbind, rows_list)
      if (is.null(gene_sets_df)) {
        gene_sets_df <- data.frame(gene = character(0),
                                   gene_set = character(0),
                                   stringsAsFactors = FALSE)
      } else {
        rownames(gene_sets_df) <- NULL
      }
      
      # Optional: quick peek / summary
      head(gene_sets_df)
      if (length(not_found)) {
        message("No genes found for ", length(not_found), " set(s): ",
                paste(unique(not_found), collapse = ", "))
      }
      
    }
    
  }
  
  #gene_sets_df_unique: gene_sets_df with power matching (and hence no duplications)
  {
    #power matching:
    
    
    # "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
    # "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
    # "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
    # "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
    # 
    # "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    # "HALLMARK_INFLAMMATORY_RESPONSE",
    # 
    # "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
    # 
    # "WP_NRF2_PATHWAY",
    # "WP_GLUTATHIONE_METABOLISM",
    # 
    # 
    #   "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
    # 
    # "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY"
    
    
    
    powermatching <- c(
      "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
      "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
      "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
      "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      
      "WP_NRF2_PATHWAY",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
      "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS")
    
    
    
    gene_sets_df_unique <- get_gene_sets_unique(gene_sets_df, powermatching)
    
    
  }
  

  
  #Create heatmap_df
  {
    
      meta_table_nk          <- "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/"
  meta_table_monocytes   <- "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/"
  meta_table_b_cells     <- "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/"
  meta_table_cd4_t_cells <- "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/"
  meta_table_cd8_t_cells <- "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
  
  
  meta_table_list <- list('NKs' = meta_table_nk,
                          'Monocytes' = meta_table_monocytes,
                          'B Cells' = meta_table_b_cells,
                          'CD4+ T Cells' = meta_table_cd4_t_cells,
                          'CD8+ T Cells' = meta_table_cd8_t_cells)
    
    heatmap_df <- get_gene_set_headmap_df(meta_table_list = meta_table_list,
                                          gene_sets_df_unique = gene_sets_df_unique,
                                          tables_to_ignore = c("CD55.Def_Mevorach_CD55-patient_meta_table.qs"))
    
    information_to_add_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/information_df_neuroimmunulogical_samples.qs")
    
    
    # ensure character
    heatmap_df$sample <- as.character(heatmap_df$sample)
    information_to_add_df$sample <- as.character(information_to_add_df$sample)
    
    # keep one row per sample; drop dataset & number_of_cells
    info_clean <- information_to_add_df %>%
      dplyr::select(-dataset, -number_of_cells) %>%
      distinct(sample, .keep_all = TRUE)
    
    # left-join; replace condition with info's condition; drop the helper column
    heatmap_df <- heatmap_df %>%
      left_join(info_clean, by = "sample", suffix = c("", ".info")) %>%
      mutate(condition = coalesce(condition.info, condition)) %>%
      dplyr::select(-condition.info)
    
    heatmap_df$gene_set_merged <- ifelse(
      heatmap_df$gene_set %in% c(
        "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
        "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
        "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY"
      ),
      "mitochondria",
      heatmap_df$gene_set
    )
    
  }
  

 #עשיתי שטיק ותכלס עיצבתי את זה בצייר
  #כי הייתה לי בעיה עם הכותרת של האפופטוזיס, רציתי שזה יתפוס לשני הבלוקים
  
  # Heatmap itself and save
  {
    ht <- get_mulitpule_gene_level_heatmap(
      df = heatmap_df,
      title = "Gene-Level Profiles of Full Gene Sets Across Conditions",
      gene_set_column_name = "gene_set_merged",
      order_of_gene_sets = c(
        "mitochondria",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "WP_NRF2_PATHWAY",
        "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
        "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS"
      ),
      rename_gene_sets = c(
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "Extrinsic Apoptosis & Inhibitors",
        "WP_NRF2_PATHWAY" = "NRF2 Pathway",
        "HALLMARK_INFLAMMATORY_RESPONSE" = "Inflammatory Response",
        "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" = "Cellular Starvation Response",
        "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS" =
          " ",
        "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS" =
          "",
        "mitochondria" = "ETC CII–CIII"
      ),
      show_row_split_labels = FALSE  # set TRUE to append “— Up/Down in target”
    )
    
    CairoJPEG(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/scripts/tryp.jpg",
      width    = 15,
      height   = 20,
      units    = "in",
      dpi      = 200,
      quality  = 100
    )
    draw(ht, annotation_legend_list = list(legend_NA("NA (missing)")))
    dev.off()
    
    
    # CairoJPEG(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/full_gene_sets_comparisons/gene_level_comparison_of_full_gene_sets_across_conditions_heatmap.jpg",
    #   width    = 15,
    #   height   = 20,
    #   units    = "in",
    #   dpi      = 600,
    #   quality  = 100
    # )
    # draw(ht, annotation_legend_list = list(legend_NA("NA (missing)")))
    # dev.off()
  }
  
}


# heatmaps of effector_vs_exhuasted_cd8_t_cells whole gene sets
{
  #gene_sets_df: all genes from all gene sets (with duplicates)
  {
    
    msig <- msigdbr(species = "Homo sapiens")
    
    gene_sets <- c(
      "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
      "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"
    )
    
    
    # Build (gene, gene_set) table from your gene_sets using the preloaded `msig`
    # Assumes: `gene_sets` (character vector) and `msig` (msigdbr table) already exist,
    # and your function get_msigdb_genes_from_gene_set_name() is defined.
    {
      stopifnot(is.character(gene_sets), length(gene_sets) > 0)
      
      rows_list <- vector("list", length(gene_sets))
      not_found <- character(0)
      
      for (i in seq_along(gene_sets)) {
        gs <- gene_sets[i]
        message("Fetching genes for: ", gs)
        
        genes <- tryCatch(
          get_msigdb_genes_from_gene_set_name(gs_name = gs, msig = msig),
          error = function(e) {
            warning("Error fetching '", gs, "': ", conditionMessage(e))
            character(0)
          }
        )
        
        if (length(genes) == 0L) {
          not_found <- c(not_found, gs)
          rows_list[[i]] <- NULL
        } else {
          rows_list[[i]] <- data.frame(
            gene     = genes,
            gene_set = gs,
            stringsAsFactors = FALSE
          )
        }
      }
      
      gene_sets_df <- do.call(rbind, rows_list)
      if (is.null(gene_sets_df)) {
        gene_sets_df <- data.frame(gene = character(0),
                                   gene_set = character(0),
                                   stringsAsFactors = FALSE)
      } else {
        rownames(gene_sets_df) <- NULL
      }
      
      # Optional: quick peek / summary
      head(gene_sets_df)
      if (length(not_found)) {
        message("No genes found for ", length(not_found), " set(s): ",
                paste(unique(not_found), collapse = ", "))
      }
      
    }
    
  }
  
  #gene_sets_df_unique: gene_sets_df with power matching (and hence no duplications)
  {
    
    
    gene_sets_df_unique <- get_gene_sets_unique(gene_sets_df, gene_sets)
    
    
  }
  
  
  
  #Create heatmap_df
  {
    
    meta_table_nk          <- "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/"
    meta_table_monocytes   <- "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/"
    meta_table_b_cells     <- "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/"
    meta_table_cd4_t_cells <- "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/"
    meta_table_cd8_t_cells <- "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
    
    
    meta_table_list <- list('NKs' = meta_table_nk,
                            'Monocytes' = meta_table_monocytes,
                            'B Cells' = meta_table_b_cells,
                            'CD4+ T Cells' = meta_table_cd4_t_cells,
                            'CD8+ T Cells' = meta_table_cd8_t_cells)
    
    heatmap_df <- get_gene_set_headmap_df(meta_table_list = meta_table_list,
                                          gene_sets_df_unique = gene_sets_df_unique,
                                          tables_to_ignore = c("CD55.Def_Mevorach_CD55-patient_meta_table.qs"))
    
    information_to_add_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/information_df_neuroimmunulogical_samples.qs")
    
    
    # ensure character
    heatmap_df$sample <- as.character(heatmap_df$sample)
    information_to_add_df$sample <- as.character(information_to_add_df$sample)
    
    # keep one row per sample; drop dataset & number_of_cells
    info_clean <- information_to_add_df %>%
      dplyr::select(-dataset, -number_of_cells) %>%
      distinct(sample, .keep_all = TRUE)
    
    # left-join; replace condition with info's condition; drop the helper column
    heatmap_df <- heatmap_df %>%
      left_join(info_clean, by = "sample", suffix = c("", ".info")) %>%
      mutate(condition = coalesce(condition.info, condition)) %>%
      dplyr::select(-condition.info)
    
    
  }
  
  #Subset to CD8+ T cells only
  
  heatmap_df_only_cd8 <- heatmap_df %>%
    subset (cell_type == "CD8+ T Cells")
  
  
  # Heatmap itself and save
  {
    ht <- get_mulitpule_gene_level_heatmap(
      df = heatmap_df_only_cd8,
      title = "Effector vs. Exhausted Markers in CD8+ T Cells Across Conditions",
      cell_pal = c(
        "CD8+ T Cells" = "#D4A373"),
      gene_set_column_name = "gene_set",
      order_of_gene_sets = c(
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"
      ),
      rename_gene_sets = c(
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP" = "Effector Markers",
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN" = "Exhaustion Markers"
      ),
      show_row_split_labels = FALSE  # set TRUE to append “— Up/Down in target”
    )
    
    # CairoJPEG(
    #   filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/full_gene_sets_comparisons/effector_vs_exhuasted_cd8_t_cells.jpg",
    #   width    = 15,
    #   height   = 10,
    #   units    = "in",
    #   dpi      = 200,
    #   quality  = 100
    # )
    # draw(ht, annotation_legend_list = list(legend_NA("NA (missing)")))
    # dev.off()
    
    
    
    }
}

#AUC scores of whole gene sets
{
  
  #get sig_list
  {
    #gene_sets_df: all genes from all gene sets (with duplicates)
    {
      msig <- msigdbr(species = "Homo sapiens")
      
      
      gene_sets <- c(
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        
        "WP_NRF2_PATHWAY",
        # "WP_GLUTATHIONE_METABOLISM",
        # "REACTOME_SELENOAMINO_ACID_METABOLISM",  # R-HSA-2408522
        
        # "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
        "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
        "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
        # "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
        
        # "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
        # "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
        
        "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
        "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
        "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
        "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY"
      )
      
      
      # Build (gene, gene_set) table from your gene_sets using the preloaded `msig`
      # Assumes: `gene_sets` (character vector) and `msig` (msigdbr table) already exist,
      # and your function get_msigdb_genes_from_gene_set_name() is defined.
      {
        stopifnot(is.character(gene_sets), length(gene_sets) > 0)
        
        rows_list <- vector("list", length(gene_sets))
        not_found <- character(0)
        
        for (i in seq_along(gene_sets)) {
          gs <- gene_sets[i]
          message("Fetching genes for: ", gs)
          
          genes <- tryCatch(
            get_msigdb_genes_from_gene_set_name(gs_name = gs, msig = msig),
            error = function(e) {
              warning("Error fetching '", gs, "': ", conditionMessage(e))
              character(0)
            }
          )
          
          if (length(genes) == 0L) {
            not_found <- c(not_found, gs)
            rows_list[[i]] <- NULL
          } else {
            rows_list[[i]] <- data.frame(
              gene     = genes,
              gene_set = gs,
              stringsAsFactors = FALSE
            )
          }
        }
        
        gene_sets_df <- do.call(rbind, rows_list)
        if (is.null(gene_sets_df)) {
          gene_sets_df <- data.frame(gene = character(0),
                                     gene_set = character(0),
                                     stringsAsFactors = FALSE)
        } else {
          rownames(gene_sets_df) <- NULL
        }
        
        # Optional: quick peek / summary
        head(gene_sets_df)
        if (length(not_found)) {
          message("No genes found for ", length(not_found), " set(s): ",
                  paste(unique(not_found), collapse = ", "))
        }
        
      }
      
    }
    
    #gene_sets_df_unique: gene_sets_df with power matching (and hence no duplications)
    {
      #power matching:
      
      
      # "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
      # "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
      # "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
      # "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
      # 
      # "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      # "HALLMARK_INFLAMMATORY_RESPONSE",
      # 
      # "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
      # 
      # "WP_NRF2_PATHWAY",
      # "WP_GLUTATHIONE_METABOLISM",
      # 
      # 
      #   "REACTOME_CASPASE_ACTIVATION_VIA_EXTRINSIC_APOPTOTIC_SIGNALLING_PATHWAY",
      # 
      # "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY"
      
      
      
      powermatching <- c(
        "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY",
        "GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_SUCCINATE_TO_UBIQUINONE",
        "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
        "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        
        "WP_NRF2_PATHWAY",
        "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
        "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
        "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS")
      
      
      
      gene_sets_df_unique <- get_gene_sets_unique(gene_sets_df, powermatching)
      
      
    }
    
    # Finally, sig_list itself
    {
      gene_sets_df_unique$gene_set <- ifelse(
        gene_sets_df_unique$gene_set %in% c(
          "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY",
          "KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_III",
          "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY"
        ),
        "mitochondria",
        gene_sets_df_unique$gene_set
      )
      
      sig_list <- gene_sets_df_unique %>%
        filter(!is.na(gene), !is.na(gene_set)) %>%
        mutate(across(c(gene, gene_set), as.character)) %>%
        group_by(gene_set) %>%
        summarise(genes = list(sort(unique(gene))), .groups = "drop") %>%
        { rlang::set_names(.$genes, .$gene_set) }
    }
    
  }
  
  # get sig_score_long_df
  {
    #create neuroimmunulogical_meta_table_lists for AUC scores
    {
      cell_paths <- list(
        "NKs"            = "/bigdata/bgutman/Collabs/Dror mevorach/data/NK/meta_tables/neuroimmunological-healthy/",
        "Monocytes"      = "/bigdata/bgutman/Collabs/Dror mevorach/data/Monocytes/meta_tables/neuroimmunological-healthy/",
        "B Cells"        = "/bigdata/bgutman/Collabs/Dror mevorach/data/B_cells/meta_tables/neuroimmunological-healthy/",
        "CD4 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD4_T_cells/meta_tables/neuroimmunological-healthy/",
        "CD8 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
      )
      
      # Map cell types to element names inside meta_table_list
      element_name_map <- c(
        "NKs"            = "meta_table_nk",
        "Monocytes"      = "meta_table_monocytes",
        "B Cells"        = "meta_table_b_cells",
        "CD4 + T cells"  = "meta_table_cd4_t_cells",
        "CD8 + T cells"  = "meta_table_cd8_t_cells"
      )
      
      # Helper: safe qread
      safe_qread <- function(path) {
        tryCatch(qread(path), error = function(e) NULL)
      }
      
      # ── Gather filenames per cell type ─────────────────────────────────────────────
      files_by_ct <- lapply(cell_paths, function(dir) {
        if (!dir.exists(dir)) return(character(0))
        f <- list.files(dir, pattern = "\\.qs$", full.names = TRUE)
        setNames(f, basename(f))  # names are basenames; values are full paths
      })
      
      # Union of sample filenames across all folders
      all_samples <- sort(unique(unlist(lapply(files_by_ct, names))))
      
      # ── Build the list-of-lists ───────────────────────────────────────────────────
      neuroimmunulogical_meta_table_lists <- list()
      
      for (sample_file in all_samples) {
        # Store under name without "_meta_table.qs"
        sample_key <- sub("_meta_table$", "", sub("\\.qs$", "", sample_file))
        
        meta_table_list <- list()
        for (ct in names(cell_paths)) {
          paths_vec <- files_by_ct[[ct]]
          if (length(paths_vec) && sample_file %in% names(paths_vec)) {
            df <- safe_qread(paths_vec[[sample_file]])
            if (!is.null(df)) {
              el_name <- element_name_map[[ct]]
              meta_table_list[[el_name]] <- df
            }
          }
        }
        
        if (length(meta_table_list) > 0) {
          neuroimmunulogical_meta_table_lists[[sample_key]] <- meta_table_list
        }
      }
      
      #remove CD55
      {
        neuroimmunulogical_meta_table_lists <- neuroimmunulogical_meta_table_lists[
          names(neuroimmunulogical_meta_table_lists) != "CD55.Def_Mevorach_CD55-patient"
        ]
      }
    }
    
    #create sig_score_results
    {
      
      # # ---- Setup -------------------------------------------------------------------
      # save_path <- "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/phenotype_whole_gene_sets_signature/sig_score_results.qs"
      # 
      # # Ensure qs is available
      # if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required.")
      # 
      # # Helper for timestamped logs
      # ts_msg <- function(...) {
      #   msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""))
      #   cat(msg, "\n")
      #   flush.console()
      # }
      # 
      # # ---- Initialize nested results container -------------------------------------
      # sig_score_results <- setNames(vector("list", length(sig_list)), names(sig_list))
      # 
      # # Pre-allocate second level to keep structure stable from the start
      # for (sig_type in names(sig_list)) {
      #   sig_score_results[[sig_type]] <- setNames(
      #     vector("list", length(neuroimmunulogical_meta_table_lists)),
      #     names(neuroimmunulogical_meta_table_lists)
      #   )
      # }
      # 
      # # Initial save of empty structure for safety
      # qs::qsave(sig_score_results, save_path)
      # ts_msg("Initialized empty results and saved to: ", save_path)
      # 
      # # ---- Main loops --------------------------------------------------------------
      # for (sig_name in names(sig_list)) {
      #   # Make sure we always pass a *named list of pathways* to the scorer
      #   entry <- sig_list[[sig_name]]
      #   pathways <- if (is.list(entry) && !is.character(entry)) {
      #     # old structure: already a list of pathways
      #     entry
      #   } else if (is.character(entry)) {
      #     # new structure: a single gene-set vector -> wrap & name it
      #     setNames(list(entry), sig_name)
      #   } else {
      #     stop(sprintf("Unrecognized sig_list entry for '%s'", sig_name))
      #   }
      #   
      #   ts_msg("Starting signature: '", sig_name, "' (", length(pathways), " pathway(s))")
      #   
      #   for (mt_name in names(neuroimmunulogical_meta_table_lists)) {
      #     ts_msg("  Processing meta_table: '", mt_name, "' for signature: '", sig_name, "'")
      #     
      #     mt <- neuroimmunulogical_meta_table_lists[[mt_name]]
      #     
      #     res <- tryCatch(
      #       {
      #         get_meta_table_auc_score_for_many_cell_types_and_pathways(
      #           meta_table_list = mt,
      #           vec_path_ids    = pathways,   # <-- always a named list now
      #           filter_nocoding = TRUE,
      #           is_msigdbr      = FALSE
      #         )
      #       },
      #       error = function(e) {
      #         ts_msg("  ERROR in '", sig_name, "' / '", mt_name, "': ", conditionMessage(e))
      #         list(error = TRUE, message = conditionMessage(e))
      #       }
      #     )
      #     
      #     sig_score_results[[sig_name]][[mt_name]] <- res
      #     qs::qsave(sig_score_results, save_path)
      #     ts_msg("  Saved progress after '", sig_name, "' / '", mt_name, "'")
      #   }
      #   
      #   ts_msg("Completed signature: '", sig_name, "'")
      # }
      # # ---- Final save --------------------------------------------------------------
      # qs::qsave(sig_score_results, save_path)
      # ts_msg("All done. Final results saved to: ", save_path)
      
    }
    
    # sig_score_results <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/phenotype_whole_gene_sets_signature/sig_score_results.qs")
    
    #create sig_score_long_df
    {
      # sig_score_results <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/phenotype_whole_gene_sets_signature/sig_score_results.qs")
      # 
      # # Flatten sig_score_results -> long data frame
      # sig_score_long_df <- imap_dfr(sig_score_results, function(sample_list, sig_type_name) {
      #   imap_dfr(sample_list, function(mat, full_name) {
      #     
      #     # keep original colnames like "CD4+ T Cells"
      #     df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
      #     df$sig_name <- rownames(df)
      #     
      #     has_fg    <- "found_genes" %in% names(df)
      #     cell_cols <- setdiff(names(df), c("sig_name", if (has_fg) "found_genes"))
      #     
      #     out <- df %>%
      #       pivot_longer(
      #         cols = all_of(cell_cols),
      #         names_to  = "cell_type",
      #         values_to = "sig_score"
      #       ) %>%
      #       mutate(
      #         condition = sub("_.*$", "", full_name),
      #         dataset   = sub("^[^_]*_(.*?)_.*$", "\\1", full_name),
      #         sample    = sub("^[^_]*_[^_]*_(.*)$", "\\1", full_name)
      #         # old 'sig_type' (positive/negative) is no longer applicable
      #       )
      #     
      #     if (!has_fg) out$found_genes <- NA_character_
      #     
      #     out %>%
      #       dplyr::select(cell_type, sample, condition, dataset, sig_name, sig_score, found_genes)
      #   })
      # })
      # 
      # # Remove NAs(as there are no B cells under the KHW sample)
      # sig_score_long_df <- sig_score_long_df %>%
      #   {
      #     removed <- dplyr::filter(., is.na(sig_score))
      #     message("Total rows removed: ", nrow(removed))
      #     dplyr::count(removed, sample, name = "n_removed_per_sample") %>% print()
      #     dplyr::filter(., !is.na(sig_score))
      #   }
      # qsave(sig_score_long_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/phenotype_whole_gene_sets_signature/sig_score_long_df.qs")
    }
    
    
    
  }
  
  sig_score_long_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/phenotype_whole_gene_sets_signature/sig_score_long_df.qs")
  
  #Create the 6 boxplots 
  {
    six_boxplots <- get_all_boxplot_whole_gene_sets(
      df       = sig_score_long_df,
      sig_order = c(
        "mitochondria",
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        "HALLMARK_INFLAMMATORY_RESPONSE",
        "WP_NRF2_PATHWAY",
        "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS",
        "GOBP_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS"
      )
    )
    
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/full_gene_sets_comparisons/auc_boxplots_of_full_gene_sets_across_conditions.jpg",
      plot     = six_boxplots,
      width    = 30,
      height   = 15,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",
      quality  = 600
    )

    
  }
  
}

#AUC scores of whole gene sets
{
  
  #get sig_list
  {
    #gene_sets_df: all genes from all gene sets (with duplicates)
    {
      
      msig <- msigdbr(species = "Homo sapiens")
      
      gene_sets <- c(
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"
      )
      
      
      # Build (gene, gene_set) table from your gene_sets using the preloaded `msig`
      # Assumes: `gene_sets` (character vector) and `msig` (msigdbr table) already exist,
      # and your function get_msigdb_genes_from_gene_set_name() is defined.
      {
        stopifnot(is.character(gene_sets), length(gene_sets) > 0)
        
        rows_list <- vector("list", length(gene_sets))
        not_found <- character(0)
        
        for (i in seq_along(gene_sets)) {
          gs <- gene_sets[i]
          message("Fetching genes for: ", gs)
          
          genes <- tryCatch(
            get_msigdb_genes_from_gene_set_name(gs_name = gs, msig = msig),
            error = function(e) {
              warning("Error fetching '", gs, "': ", conditionMessage(e))
              character(0)
            }
          )
          
          if (length(genes) == 0L) {
            not_found <- c(not_found, gs)
            rows_list[[i]] <- NULL
          } else {
            rows_list[[i]] <- data.frame(
              gene     = genes,
              gene_set = gs,
              stringsAsFactors = FALSE
            )
          }
        }
        
        gene_sets_df <- do.call(rbind, rows_list)
        if (is.null(gene_sets_df)) {
          gene_sets_df <- data.frame(gene = character(0),
                                     gene_set = character(0),
                                     stringsAsFactors = FALSE)
        } else {
          rownames(gene_sets_df) <- NULL
        }
        
        # Optional: quick peek / summary
        head(gene_sets_df)
        if (length(not_found)) {
          message("No genes found for ", length(not_found), " set(s): ",
                  paste(unique(not_found), collapse = ", "))
        }
        
      }
      
    }
    
    #gene_sets_df_unique: gene_sets_df with power matching (and hence no duplications)
    {
      
      
      gene_sets_df_unique <- get_gene_sets_unique(gene_sets_df, gene_sets)
      
      
    }
    
    # Finally, sig_list itself
    {
      
      
      sig_list <- gene_sets_df_unique %>%
        filter(!is.na(gene), !is.na(gene_set)) %>%
        mutate(across(c(gene, gene_set), as.character)) %>%
        group_by(gene_set) %>%
        summarise(genes = list(sort(unique(gene))), .groups = "drop") %>%
        { rlang::set_names(.$genes, .$gene_set) }
    }
    
  }
  
  # get sig_score_long_df
  {
    #create neuroimmunulogical_meta_table_lists for AUC scores - only cd8+ t cells this time
    {
      cell_paths <- list(
        "CD8 + T cells"  = "/bigdata/bgutman/Collabs/Dror mevorach/data/CD8_T_cells/meta_tables/neuroimmunological-healthy/"
      )
      
      # Map cell types to element names inside meta_table_list
      element_name_map <- c(
        "CD8 + T cells"  = "meta_table_cd8_t_cells"
      )
      
      # Helper: safe qread
      safe_qread <- function(path) {
        tryCatch(qread(path), error = function(e) NULL)
      }
      
      # ── Gather filenames per cell type ─────────────────────────────────────────────
      files_by_ct <- lapply(cell_paths, function(dir) {
        if (!dir.exists(dir)) return(character(0))
        f <- list.files(dir, pattern = "\\.qs$", full.names = TRUE)
        setNames(f, basename(f))  # names are basenames; values are full paths
      })
      
      # Union of sample filenames across all folders
      all_samples <- sort(unique(unlist(lapply(files_by_ct, names))))
      
      # ── Build the list-of-lists ───────────────────────────────────────────────────
      neuroimmunulogical_meta_table_lists <- list()
      
      for (sample_file in all_samples) {
        # Store under name without "_meta_table.qs"
        sample_key <- sub("_meta_table$", "", sub("\\.qs$", "", sample_file))
        
        meta_table_list <- list()
        for (ct in names(cell_paths)) {
          paths_vec <- files_by_ct[[ct]]
          if (length(paths_vec) && sample_file %in% names(paths_vec)) {
            df <- safe_qread(paths_vec[[sample_file]])
            if (!is.null(df)) {
              el_name <- element_name_map[[ct]]
              meta_table_list[[el_name]] <- df
            }
          }
        }
        
        if (length(meta_table_list) > 0) {
          neuroimmunulogical_meta_table_lists[[sample_key]] <- meta_table_list
        }
      }
      
      #remove CD55
      {
        neuroimmunulogical_meta_table_lists <- neuroimmunulogical_meta_table_lists[
          names(neuroimmunulogical_meta_table_lists) != "CD55.Def_Mevorach_CD55-patient"
        ]
      }
    }
    
    #create sig_score_results
    {
      
      # ---- Setup -------------------------------------------------------------------
      save_path <- "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/exhuastion_cd8_t_cells_signature/sig_score_results.qs"
      
      # Ensure qs is available
      if (!requireNamespace("qs", quietly = TRUE)) stop("Package 'qs' is required.")
      
      # Helper for timestamped logs
      ts_msg <- function(...) {
        msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""))
        cat(msg, "\n")
        flush.console()
      }
      
      # ---- Initialize nested results container -------------------------------------
      sig_score_results <- setNames(vector("list", length(sig_list)), names(sig_list))
      
      # Pre-allocate second level to keep structure stable from the start
      for (sig_type in names(sig_list)) {
        sig_score_results[[sig_type]] <- setNames(
          vector("list", length(neuroimmunulogical_meta_table_lists)),
          names(neuroimmunulogical_meta_table_lists)
        )
      }
      
      # Initial save of empty structure for safety
      qs::qsave(sig_score_results, save_path)
      ts_msg("Initialized empty results and saved to: ", save_path)
      
      # ---- Main loops --------------------------------------------------------------
      for (sig_name in names(sig_list)) {
        # Make sure we always pass a *named list of pathways* to the scorer
        entry <- sig_list[[sig_name]]
        pathways <- if (is.list(entry) && !is.character(entry)) {
          # old structure: already a list of pathways
          entry
        } else if (is.character(entry)) {
          # new structure: a single gene-set vector -> wrap & name it
          setNames(list(entry), sig_name)
        } else {
          stop(sprintf("Unrecognized sig_list entry for '%s'", sig_name))
        }
        
        ts_msg("Starting signature: '", sig_name, "' (", length(pathways), " pathway(s))")
        
        for (mt_name in names(neuroimmunulogical_meta_table_lists)) {
          ts_msg("  Processing meta_table: '", mt_name, "' for signature: '", sig_name, "'")
          
          mt <- neuroimmunulogical_meta_table_lists[[mt_name]]
          
          res <- tryCatch(
            {
              get_meta_table_auc_score_for_many_cell_types_and_pathways(
                meta_table_list = mt,
                vec_path_ids    = pathways,   # <-- always a named list now
                filter_nocoding = TRUE,
                is_msigdbr      = FALSE
              )
            },
            error = function(e) {
              ts_msg("  ERROR in '", sig_name, "' / '", mt_name, "': ", conditionMessage(e))
              list(error = TRUE, message = conditionMessage(e))
            }
          )
          
          sig_score_results[[sig_name]][[mt_name]] <- res
          qs::qsave(sig_score_results, save_path)
          ts_msg("  Saved progress after '", sig_name, "' / '", mt_name, "'")
        }
        
        ts_msg("Completed signature: '", sig_name, "'")
      }
      # ---- Final save --------------------------------------------------------------
      qs::qsave(sig_score_results, save_path)
      ts_msg("All done. Final results saved to: ", save_path)
      
    }
    
    # sig_score_results <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/exhuastion_cd8_t_cells_signature/sig_score_results.qs")
    
    #create sig_score_long_df
    {
      # sig_score_results <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/exhuastion_cd8_t_cells_signature/sig_score_results.qs")
      # 
      # # Flatten sig_score_results -> long data frame
      # sig_score_long_df <- imap_dfr(sig_score_results, function(sample_list, sig_type_name) {
      #   imap_dfr(sample_list, function(mat, full_name) {
      # 
      #     # keep original colnames like "CD4+ T Cells"
      #     df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
      #     df$sig_name <- rownames(df)
      # 
      #     has_fg    <- "found_genes" %in% names(df)
      #     cell_cols <- setdiff(names(df), c("sig_name", if (has_fg) "found_genes"))
      # 
      #     out <- df %>%
      #       pivot_longer(
      #         cols = all_of(cell_cols),
      #         names_to  = "cell_type",
      #         values_to = "sig_score"
      #       ) %>%
      #       mutate(
      #         condition = sub("_.*$", "", full_name),
      #         dataset   = sub("^[^_]*_(.*?)_.*$", "\\1", full_name),
      #         sample    = sub("^[^_]*_[^_]*_(.*)$", "\\1", full_name)
      #         # old 'sig_type' (positive/negative) is no longer applicable
      #       )
      # 
      #     if (!has_fg) out$found_genes <- NA_character_
      # 
      #     out %>%
      #       dplyr::select(cell_type, sample, condition, dataset, sig_name, sig_score, found_genes)
      #   })
      # })
      # 
      # # Remove NAs(as there are no B cells under the KHW sample)
      # sig_score_long_df <- sig_score_long_df %>%
      #   {
      #     removed <- dplyr::filter(., is.na(sig_score))
      #     message("Total rows removed: ", nrow(removed))
      #     dplyr::count(removed, sample, name = "n_removed_per_sample") %>% print()
      #     dplyr::filter(., !is.na(sig_score))
      #   }
      # qsave(sig_score_long_df, "/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/exhuastion_cd8_t_cells_signature/sig_score_long_df.qs")
    }
    
    
    
  }
  
  sig_score_long_df <- qread("/bigdata/bgutman/Collabs/Dror mevorach/data/multi_cell_analysis/exhuastion_cd8_t_cells_signature/sig_score_long_df.qs")
  
  #Create the 2 boxplots 
  {
    two_boxplots <- get_all_boxplot_whole_gene_sets(
      df       = sig_score_long_df,
      ncol      = 2,      # number of columns in the grid
      sig_order = c(
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"),
      row_renames     = c(
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP" =
          "Effectot Markers",
        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN" =
          "Exhaustion Markers"
      )
      
    )
    
    
    ggsave(
      filename = "/bigdata/bgutman/Collabs/Dror mevorach/plots/multi_cell_analysis/full_gene_sets_comparisons/auc_boxplots_of_effector_vs_exhuasted_cd8_t_cells.jpg",
      plot     = two_boxplots,
      width    = 20,
      height   = 8,
      units    = "in",
      dpi      = 600,
      device   = "jpeg",
      type     = "cairo",
      quality  = 600
    )
    
    
  }
  
}


















