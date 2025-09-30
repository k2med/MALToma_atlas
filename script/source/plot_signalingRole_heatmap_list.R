plot_signalingRole_heatmap_list <- function(object.list, signaling = NULL, pattern = c("outgoing", 
                                                                                 "incoming", "all"),
                                            slot.name = "netP", color.use = NULL, 
                                            color.heatmap = "BuGn", color.bar = "#3b5684",
                                            title = NULL, width = 10, height = 8, 
                                            font.size = 8, font.size.title = 10, cluster.rows = FALSE, 
                                            cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  heatmap.list <- list()
  all.mat.values <- c()
  all.row.strength <- c()
  all.col.strength <- c()
  
  for (object in object.list) {
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (j in 1:length(centr)) {
      outgoing[, j] <- centr[[j]]$outdeg
      incoming[, j] <- centr[[j]]$indeg
    }
    mat <- switch(pattern,
                  outgoing = t(outgoing),
                  incoming = t(incoming),
                  all = t(outgoing + incoming))
    if (!is.null(signaling)) {
      mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
      mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
      idx <- match(rownames(mat1), signaling)
      mat[idx[!is.na(idx)], ] <- mat1
    }
    mat[mat == 0] <- NA
    mat <- -1 / log(mat)
    all.mat.values <- c(all.mat.values, mat[!is.na(mat)])
    all.row.strength <- c(all.row.strength, rowSums(mat, na.rm = TRUE))
    all.col.strength <- c(all.col.strength, colSums(mat, na.rm = TRUE))
  }
  
  global.breaks <- round(range(all.mat.values, na.rm = TRUE), 1)
  global.row.range <- round(range(all.row.strength, na.rm = TRUE), 1)
  global.col.range <- round(range(all.col.strength, na.rm = TRUE), 1)
  
  for (i in seq_along(object.list)) {
    object <- object.list[[i]]
    object.name <- names(object.list)[i]
    
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (j in 1:length(centr)) {
      outgoing[, j] <- centr[[j]]$outdeg
      incoming[, j] <- centr[[j]]$indeg
    }
    
    mat <- switch(pattern,
                  outgoing = t(outgoing),
                  incoming = t(incoming),
                  all = t(outgoing + incoming))
    
    this.title <- if (is.null(title)) {
      paste0(pattern, " signaling - ", object.name)
    } else {
      paste0(pattern, " signaling - ", title, " - ", object.name)
    }
    
    if (!is.null(signaling)) {
      mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
      mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
      idx <- match(rownames(mat1), signaling)
      mat[idx[!is.na(idx)], ] <- mat1
      dimnames(mat) <- list(signaling, colnames(mat1))
    }
    
    mat.ori <- mat
    mat[mat == 0] <- NA
    mat <- -1/log(mat)
    
    this.color.use <- if (is.null(color.use)) scPalette(length(colnames(mat))) else color.use
    color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(100)
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    names(this.color.use) <- colnames(mat)
    
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = this.color.use), which = "column", 
                                        show_legend = FALSE, show_annotation_name = FALSE, 
                                        simple_anno_size = grid::unit(0.2, "cm"))
    
    ha1 <- rowAnnotation(
      Strength = anno_barplot(
        rowSums(mat, na.rm = TRUE),
        border = FALSE,
        bar_width = 0.7,
        lwd = 0.25,
        gp = gpar(col = "white", fill = color.bar),
        ylim = global.row.range,
        axis_param = list(gp = gpar(fontsize = 5,
                                    lwd = 0.33))
      ),
      width = unit(0.6, "cm"), 
      height = unit(height, 'cm'),
      show_annotation_name = FALSE
    )
    
    ha2 <- HeatmapAnnotation(
      Strength = anno_barplot(
        colSums(mat, na.rm = TRUE),
        border = FALSE,
        bar_width = 0.7,
        lwd = 0.25,
        gp = gpar(col = "white", fill = color.bar),
        ylim = global.col.range,
        axis_param = list(gp = gpar(fontsize = 5,
                                    lwd = 0.33))
      ),
      width = unit(width, "cm"),
      height = unit(0.6, 'cm'),
      show_annotation_name = FALSE
    )
    
    ht <- Heatmap(mat,
                  col = color.heatmap.use,
                  na_col = "white",
                  name = "Relative strength",
                  top_annotation = ha2,
                  right_annotation = ha1,
                  cluster_rows = cluster.rows,
                  cluster_columns = cluster.cols,
                  row_names_side = "left",
                  row_names_rot = 0,
                  row_names_gp = gpar(fontsize = font.size),
                  column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"),
                  height = unit(height, "cm"),
                  column_title = this.title,
                  column_title_gp = gpar(fontsize = font.size.title),
                  column_names_rot = 90,
                  heatmap_legend_param = list(
                    title_gp = gpar(fontsize = font.size, fontface = "plain"),
                    title_position = "leftcenter-rot",
                    border = NA,
                    at = global.breaks,
                    legend_height = unit(15, "mm"),
                    labels_gp = gpar(fontsize = font.size),
                    grid_width = unit(2, "mm")
                  )
    )
    heatmap.list[[i]] <- ht
  }
  
  combined.ht <- Reduce(`+`, heatmap.list)
  return(combined.ht)
}
