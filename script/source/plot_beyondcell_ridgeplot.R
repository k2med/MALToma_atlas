bcRidge <- function (bc, signatures, idents = NULL) 
{
  if (class(bc) != "beyondcell") 
    stop("bc must be a beyondcell object.")
  if (!is.character(signatures)) {
    stop("Signatures must be a character vector.")
  }
  if (length(signatures) == 1 & signatures[1] == "all") {
    signatures <- rownames(bc@normalized)
    in.signatures <- rep(TRUE, times = nrow(bc@normalized))
  }
  else {
    in.signatures <- !is.null(signatures) & signatures %in% 
      rownames(bc@normalized)
    if (all(!in.signatures)) {
      stop("None of the specified signatures were found.")
    }
    else if (any(!in.signatures)) {
      warning(paste0("These signatures were not found in bc: ", 
                     paste0(signatures[!in.signatures], collapse = ", "), 
                     "."))
    }
  }
  if (!is.null(idents)) {
    if (length(idents) != 1) {
      stop("Idents must be a single metadata column.")
    }
    if (idents %in% colnames(bc@meta.data)) {
      meta <- bc@meta.data[colnames(bc@normalized), 
                           idents, drop = TRUE]
    }
    else {
      stop("Idents not found.")
    }
  }
  else {
    meta <- rep("", times = ncol(bc@normalized))
  }
  lvls <- levels(as.factor(meta))
  sub.bc <- bc@normalized[signatures[in.signatures], , drop = FALSE]
  limits <- c(min(as.vector(sub.bc), na.rm = TRUE), max(as.vector(sub.bc), 
                                                        na.rm = TRUE))
  info <- FindDrugs(bc, x = signatures[in.signatures])
  p <- lapply(signatures[in.signatures], function(x) {
    sub.df <- na.omit(data.frame(bcscore = sub.bc[x, ], 
                                 condition = as.factor(meta), row.names = colnames(sub.bc)))
    if (x %in% info$IDs) {
      drug.and.MoA <- info[which(info$IDs == x), c("preferred.and.sigs", 
                                                   "MoAs")]
      drug.and.MoA[2] <- ifelse(test = drug.and.MoA[2] == 
                                  "NA", yes = "", no = drug.and.MoA[2])
    }
    else {
      drug.and.MoA <- c(x, "")
    }
    ridge <- ggplot(sub.df, aes(x = bcscore, y = condition, fill = condition)) + 
      geom_density_ridges(alpha = 0.5, scale = 1, size = 0.116,
                          quantile_lines = TRUE,
                          quantile_fun = function(bcscore, ...) mean(bcscore)) +
      scale_y_discrete(limits = rev) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(x = 'Sensitivity',
           y = NULL,
           title = drug.and.MoA[1],
           subtitle = drug.and.MoA[2]) +
      theme(text = element_text(size = 5),
            panel.grid.major = element_line(colour = "grey90", size = 0.116),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks =  element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 5),
            plot.title = element_text(colour = "black", size = 5),
            plot.subtitle = element_text(colour = "black", size = 5),
            legend.title = element_text(colour = "black", size = 5))
    return(ridge)
  })
  if (length(p) == 1) 
    p <- p[[1]]
  return(p)
}
