build_cellchat_merged <- function(
    seurat_obj,
    subtype_levels = c("LME_1","LME_2","LME_3"),
    group_col      = "fine_cell_type",
    subtype_col    = "subtype",
    assay          = "RNA",
    species        = c("human","mouse"),
    min_cells      = 10,
    workers        = 1
) {
  species <- match.arg(species)
  db_use  <- if (species == "human") CellChatDB.human else CellChatDB.mouse
  ppi_use <- if (species == "human") PPI.human       else PPI.mouse
  
  seurat_obj[[subtype_col]] <- factor(seurat_obj[[subtype_col]][, 1], levels = subtype_levels)
  cell_type_order <- levels(seurat_obj[[group_col]][, 1])
  
  plan("multisession", workers = workers)
  
  obj_list <- vector("list", length(subtype_levels))
  names(obj_list) <- subtype_levels
  
  for (sb in subtype_levels) {
    seu_sb <- subset(seurat_obj, subset = .data[[subtype_col]] == sb)
    if (ncol(seu_sb) == 0) {
      obj_list[[sb]] <- NULL
      next
    }
    
    cc <- createCellChat(object = seu_sb, group.by = group_col, assay = assay)
    cc@DB <- db_use
    
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- smoothData(cc, adj = ppi_use)
    cc <- computeCommunProb(cc, population.size = TRUE, raw.use = FALSE)
    cc <- filterCommunication(cc, min.cells = min_cells)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    
    # harmonize idents order
    cc@idents <- factor(cc@idents, levels = cell_type_order)
    
    # centrality for downstream comparison plots
    cc <- netAnalysis_computeCentrality(cc)
    
    obj_list[[sb]] <- cc
  }
  
  # drop NULLs (in case a subtype had no cells)
  obj_list <- obj_list[!vapply(obj_list, is.null, logical(1))]
  if (length(obj_list) == 0) stop("No CellChat objects were created.")
  
  merged <- mergeCellChat(obj_list, add.names = names(obj_list))
  list(merged = merged, per_subtype = obj_list)
}
