# Functions for plotting GSEA adapted from Zhang et al. 2025 (DOI: 10.1002/mdr2.70000), with local modifications
gseaNb <- function (object = NULL, subPlot = 3, lineSize = 0.8, geneSetID = NULL, 
          rmSegment = FALSE, termWidth = 40, segCol = "red", addGene = NULL, 
          geneCol = NULL, arrowAngle = 20, arrowLength = 0.2, arrowEnd = "first", 
          arrowType = "closed", curveCol = c("#76BA99", "#EB4747", 
                                             "#996699"), htCol = c("#08519C", "#A50F15"), rankCol = c("#08519C", 
                                                                                                      "white", "#A50F15"), rankSeq = 5000, htHeight = 0.3, 
          force = 20, max.overlaps = 50, geneSize = 4, newGsea = FALSE, 
          addPoint = TRUE, newCurveCol = c("#336699", "white", "#993399"), 
          newHtCol = c("#336699", "white", "#993399"), rmHt = FALSE, 
          addPval = FALSE, pvalX = 0.9, pvalY = 0.9, pvalSize = 4, 
          pCol = "grey30", pHjust = 1, rmPrefix = TRUE, nesDigit = 2, 
          pDigit = 2, markTopgene = FALSE, topGeneN = 5, kegg = FALSE, 
          legend.position = "right", add.geneExpHt = FALSE, exp = NULL, 
          scale.exp = TRUE, sample.order = NULL, exp.col = c("blue", 
                                                             "white", "red"), ht.legend = TRUE, ght.relHight = 0.4, 
          ght.geneText.size = 6, ght.facet = FALSE, ght.facet.scale = "free", 
          termID.order = NULL, rank.gene = NULL, rank.gene.nudgey = 2) 
{
  gsdata <- purrr::map_df(geneSetID, function(setid) {
    gsInfo(object, geneSetID = setid) %>% dplyr::mutate(id = setid)
  })
  if (kegg == FALSE) {
    gsdata1 <- purrr::map_df(unique(gsdata$Description), 
                             function(setid) {
                               tmp <- gsdata %>% dplyr::filter(Description == 
                                                                 setid) %>% dplyr::mutate(gene_name = names(object@geneList)) %>% 
                                 dplyr::filter(position == 1)
                             })
  }
  else {
    gene2Symbol <- object@gene2Symbol %>% data.frame()
    gsdata1 <- purrr::map_df(unique(gsdata$Description), 
                             function(setid) {
                               tmp <- gsdata %>% dplyr::filter(Description == 
                                                                 setid) %>% dplyr::mutate(gene_name = gene2Symbol$.) %>% 
                                 dplyr::filter(position == 1)
                             })
  }
  data_ga <- data.frame(object) %>% dplyr::filter(ID %in% 
                                                    geneSetID)
  data_ga <- data_ga[unique(gsdata$id), ]
  niceTit <- purrr::map_chr(unique(gsdata$Description), function(x) {
    tit <- unlist(strsplit(x, split = "_"))
    if (length(tit) == 1) {
      niceTit <- paste(stringr::str_to_title(tit[1:length(tit)]), 
                       collapse = " ") %>% stringr::str_wrap(., width = termWidth)
    }
    else {
      if (rmPrefix == TRUE) {
        niceTit <- paste(stringr::str_to_title(tit[2:length(tit)]), 
                         collapse = " ") %>% stringr::str_wrap(., width = termWidth)
      }
      else {
        niceTit <- paste(stringr::str_to_title(tit[1:length(tit)]), 
                         collapse = " ") %>% stringr::str_wrap(., width = termWidth)
      }
    }
  })
  if (length(geneSetID) != 1) {
    ledend.t <- niceTit
    niceTit <- ""
  }
  if (length(geneSetID) == 1) {
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~runningScore), 
                               size = lineSize)
    line.col <- ggplot2::scale_color_gradient(low = curveCol[1], 
                                              high = curveCol[2])
    legend.position = "none"
  }
  else {
    mulcol <- curveCol
    names(mulcol) <- unique(gsdata$Description)
    line <- ggplot2::geom_line(ggplot2::aes_(color = ~Description), 
                               size = lineSize)
    line.col <- ggplot2::scale_color_manual(values = mulcol, 
                                            labels = ledend.t, name = "Term Name")
    legend.position = legend.position
  }
  pcurve <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, 
                                                  y = ~runningScore)) + line + line.col + ggplot2::geom_hline(yintercept = 0, 
                                                                                                              size = lineSize, color = "black", lty = "dashed") + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                               0)) + ggplot2::theme(legend.position = legend.position, 
                                                                                                    legend.box.background = ggplot2::element_blank(), plot.title = ggplot2::element_text(hjust = 0.5), 
                                                                                                    panel.grid = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
                                                                                                    axis.text.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                                                                    axis.title.x = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent"), 
                                                                                                    plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, 
                                                                                                                                  l = 0.2, unit = "cm")) + ggplot2::ylab("Running Enrichment Score") + 
    ggplot2::ggtitle(niceTit)
  midpoint <- sum(range(gsdata$runningScore))/2
  pnew <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore, 
                                                color = ~runningScore)) + ggplot2::geom_hline(yintercept = 0, 
                                                                                              size = lineSize, color = "black", lty = "dashed") + 
    ggplot2::geom_line(size = lineSize) + ggplot2::geom_segment(data = gsdata1, 
                                                                ggplot2::aes_(xend = ~x, yend = 0)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::scale_color_gradient2(low = newCurveCol[1], 
                                   mid = newCurveCol[2], high = newCurveCol[3], midpoint = midpoint) + 
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::theme(legend.position = "none", 
                                                                   plot.title = ggplot2::element_text(hjust = 0.5), axis.ticks.x = ggplot2::element_blank(), 
                                                                   axis.text.x = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                                   axis.title.x = ggplot2::element_blank(), legend.background = ggplot2::element_rect(fill = "transparent"), 
                                                                   plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, 
                                                                                                 l = 0.2, unit = "cm")) + ggplot2::ylab("Running Enrichment Score") + 
    ggplot2::ggtitle(niceTit)
  if (addPoint == TRUE) {
    panother <- pnew + ggplot2::geom_point()
  }
  else {
    panother <- pnew
  }
  if (newGsea == FALSE) {
    pcurveRes <- pcurve
  }
  else {
    pcurveRes <- panother
  }
  if (is.null(addGene)) {
    plabel <- pcurveRes
  }
  else {
    if (markTopgene == TRUE) {
      geneLabel <- gsdata1 %>% dplyr::arrange(x) %>% dplyr::slice_head(n = topGeneN)
    }
    else {
      geneLabel <- gsdata1 %>% dplyr::filter(gene_name %in% 
                                               addGene)
    }
    if (nrow(geneLabel) == 0) {
      message("Your gene is not in this pathway! Please choose again!")
    }
    else {
      if (rmSegment == TRUE) {
        if (is.null(geneCol)) {
          plabel <- pcurveRes + ggrepel::geom_text_repel(data = geneLabel, 
                                                         ggplot2::aes_(label = ~gene_name), force = force, 
                                                         max.overlaps = max.overlaps, size = geneSize, 
                                                         fontface = "italic", arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                     length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                     ends = arrowEnd, type = arrowType))
        }
        else {
          plabel <- pcurveRes + ggrepel::geom_text_repel(data = geneLabel, 
                                                         ggplot2::aes_(label = ~gene_name), force = force, 
                                                         max.overlaps = max.overlaps, size = geneSize, 
                                                         fontface = "italic", color = geneCol, arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                                                      length = ggplot2::unit(arrowLength, "cm"), 
                                                                                                                      ends = arrowEnd, type = arrowType))
        }
      }
      else {
        if (is.null(geneCol)) {
          plabel <- pcurveRes + ggplot2::geom_segment(data = geneLabel, 
                                                      ggplot2::aes_(xend = ~x, yend = 0), color = segCol) + 
            ggrepel::geom_text_repel(data = geneLabel, 
                                     ggplot2::aes_(label = ~gene_name), force = force, 
                                     max.overlaps = max.overlaps, size = geneSize, 
                                     fontface = "italic", arrow = ggplot2::arrow(angle = arrowAngle, 
                                                                                 length = ggplot2::unit(arrowLength, 
                                                                                                        "cm"), ends = arrowEnd, type = arrowType))
        }
        else {
          plabel <- pcurveRes + ggplot2::geom_segment(data = geneLabel, 
                                                      ggplot2::aes_(xend = ~x, yend = 0), color = segCol) + 
            ggrepel::geom_text_repel(data = geneLabel, 
                                     ggplot2::aes_(label = ~gene_name), force = force, 
                                     max.overlaps = max.overlaps, size = geneSize, 
                                     fontface = "italic", color = geneCol, 
                                     arrow = ggplot2::arrow(angle = arrowAngle, 
                                                            length = ggplot2::unit(arrowLength, 
                                                                                   "cm"), ends = arrowEnd, type = arrowType))
        }
      }
    }
  }
  if (addPval == TRUE) {
    pLabel <- paste0("NES: ", round(data_ga$NES, digits = nesDigit), 
                     "\n", "P value: ", ifelse(data_ga$pvalue < 0.001, 
                                              "< 0.001", round(data_ga$pvalue, digits = pDigit)), 
                     "\n", "Adjusted P value: ", ifelse(data_ga$p.adjust < 
                                                        0.001, "< 0.001", round(data_ga$p.adjust, digits = pDigit)), 
                     "\n", sep = " ")
    px <- pvalX * nrow(gsdata[which(gsdata$id == geneSetID[1]), 
    ])
    py <- pvalY * sum(abs(range(gsdata$runningScore))) + 
      min(gsdata$runningScore)
    if (length(geneSetID) == 1) {
      pLabelOut <- plabel + ggplot2::annotate(geom = "text", 
                                              x = px, y = py, label = pLabel, size = pvalSize, 
                                              color = pCol, fontface = "italic", hjust = pHjust)
    }
    else {
      mytable <- tibble::tibble(x = px, y = py, table = list(tibble::tibble(NES = round(data_ga$NES, 
                                                                                        digits = nesDigit), `P value` = ifelse(data_ga$pvalue < 
                                                                                                                              0.001, "< 0.001", round(data_ga$pvalue, digits = pDigit)), 
                                                                            `Adjusted P value` = ifelse(data_ga$p.adjust < 
                                                                                                        0.001, "< 0.001", round(data_ga$p.adjust, 
                                                                                                                                digits = pDigit)))))
      pLabelOut <- plabel + ggpp::geom_table(data = mytable, 
                                             ggplot2::aes(px, py, label = table))
    }
  }
  else {
    pLabelOut <- plabel
  }
  if (length(geneSetID) == 1) {
    line.col <- ggplot2::scale_color_manual(values = "black")
  }
  if (add.geneExpHt == TRUE) {
    pseg.b = 0
  }
  else {
    pseg.b = 0.2
  }
  pseg <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore, 
                                                color = ~Description)) + ggplot2::geom_segment(data = gsdata1, 
                                                                                               ggplot2::aes_(x = ~x, xend = ~x, y = 0, yend = 1), show.legend = F) + 
    line.col + ggplot2::scale_x_continuous(expand = c(0, 
                                                      0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                                                       axis.text = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), 
                                                       panel.grid = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(), 
                                                       strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_blank(), 
                                                       panel.spacing = ggplot2::unit(0.1, "cm"), plot.margin = ggplot2::margin(t = 0, 
                                                                                                                               r = 0.2, b = pseg.b, l = 0.2, unit = "cm")) + ggplot2::xlab("Rank in Ordered Dataset") + 
    ggplot2::facet_wrap(~Description, ncol = 1)
  if (subPlot > 2) {
    pseg <- pseg + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  else {
    pseg <- pseg
  }
  d <- purrr::map_df(unique(gsdata$Description), function(setid) {
    tmp <- gsdata %>% dplyr::filter(Description == setid)
    v <- seq(1, sum(tmp$position), length.out = 9)
    inv <- findInterval(rev(cumsum(tmp$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }
    color <- (grDevices::colorRampPalette(c(htCol[1], "white", 
                                            htCol[2])))(10)
    ymin <- 0
    yy <- htHeight
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = color[unique(inv)], Description = setid)
  })
  pseg_ht <- pseg + ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin, 
                                                     xmax = ~xmax, ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), 
                                       data = d, alpha = 0.8, inherit.aes = FALSE)
  pseg_ht1 <- pseg_ht + ggplot2::xlab("") + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
                                                           plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0, 
                                                                                         l = 0.2, unit = "cm"))
  if (add.geneExpHt == TRUE) {
    prank.b = 0
  }
  else {
    prank.b = 0.2
  }
  prank <- ggplot2::ggplot(gsdata[which(gsdata$Description == 
                                          unique(gsdata$Description)[1]), ], ggplot2::aes_(x = ~x, 
                                                                                           y = ~geneList)) + ggplot2::geom_col(ggplot2::aes_(fill = ~geneList), 
                                                                                                                               width = 1, color = NA, show.legend = F) + ggplot2::scale_fill_gradient2(low = rankCol[1], 
                                                                                                                                                                                                       mid = rankCol[2], high = rankCol[3], midpoint = 0) + 
    ggplot2::geom_hline(yintercept = 0, size = 0.8, color = "black", 
                        lty = "dashed") + ggplot2::scale_x_continuous(breaks = seq(0, 
                                                                                   nrow(gsdata), rankSeq)) + ggplot2::theme_bw(base_size = 14) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   plot.margin = ggplot2::margin(t = -0.1, r = 0.2, 
                                                 b = prank.b, l = 0.2, unit = "cm")) + ggplot2::coord_cartesian(expand = 0) + 
    ggplot2::ylab("Ranked List") + ggplot2::xlab("Rank in Ordered Dataset")
  if (add.geneExpHt == TRUE) {
    prank <- prank + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  else {
    prank <- prank
  }
  if (kegg == FALSE) {
    rank.g <- data.frame(logfc = object@geneList, gene_name = names(object@geneList)) %>% 
      dplyr::mutate(x = 1:length(object@geneList))
  }
  else {
    rank.g <- data.frame(logfc = object@geneList, gene_name = object@gene2Symbol) %>% 
      dplyr::mutate(x = 1:length(object@geneList))
  }
  if (!is.null(rank.gene)) {
    target.rank.g <- rank.g %>% dplyr::filter(gene_name %in% 
                                                rank.gene) %>% dplyr::mutate(vjust = ifelse(logfc > 
                                                                                              0, "bottom", "top"), nudge_y = ifelse(logfc > 0, 
                                                                                                                                    -rank.gene.nudgey, rank.gene.nudgey))
    prank <- prank + ggrepel::geom_text_repel(data = target.rank.g, 
                                              ggplot2::aes(x = as.numeric(x), y = 0, label = gene_name, 
                                                           vjust = vjust, nudge_y = nudge_y), max.overlaps = 200, 
                                              direction = "x", angle = 90, fontface = "italic", 
                                              size = geneSize)
  }
  d <- purrr::map_df(unique(d$Description), function(x) {
    tmp <- d %>% dplyr::filter(Description == x)
    htcolor <- (grDevices::colorRampPalette(newHtCol))(nrow(tmp))
    tmp <- tmp %>% dplyr::mutate(htcol = htcolor)
  })
  ht <- ggplot2::ggplot(gsdata, ggplot2::aes_(x = ~x, y = ~runningScore)) + 
    ggplot2::geom_rect(ggplot2::aes_(xmin = ~xmin, xmax = ~xmax, 
                                     ymin = ~ymin, ymax = ~ymax, fill = ~I(htcol)), data = d, 
                       alpha = 0.8, inherit.aes = FALSE) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme_bw(base_size = 14) + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                                                       axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), 
                                                       axis.title = ggplot2::element_blank(), plot.margin = ggplot2::margin(t = 0, 
                                                                                                                            r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
  if (add.geneExpHt == TRUE) {
    target.g <- purrr::map_df(data_ga$ID, function(x) {
      tmp <- data_ga %>% dplyr::filter(ID == x)
      coregene <- unique(unlist(strsplit(tmp$core_enrichment, 
                                         split = "\\/")))
      output <- data.frame(gene_name = coregene, ID = x, 
                           Description = tmp$Description) %>% dplyr::distinct(., 
                                                                              gene_name, .keep_all = TRUE)
    })
    gpos <- if (kegg == TRUE) {
      match(target.g$gene_name, object@gene2Symbol)
    }
    else {
      match(target.g$gene_name, names(object@geneList))
    }
    ginfo <- target.g %>% dplyr::mutate(gpos = gpos) %>% 
      dplyr::arrange(gpos)
    if (scale.exp == TRUE) {
      gexp <- t(scale(t(exp[, 2:ncol(exp)]), scale = TRUE, 
                      center = TRUE)) %>% data.frame()
      gexp$gene_name <- exp[, 1]
    }
    else {
      gexp <- exp
      colnames(gexp)[1] <- "gene_name"
    }
    exp.long <- gexp %>% dplyr::filter(gene_name %in% unique(ginfo$gene_name)) %>% 
      dplyr::left_join(., ginfo[, 1:2], by = "gene_name") %>% 
      reshape2::melt(., id.vars = c("gene_name", "ID"))
    exp.long$gene_name <- factor(exp.long$gene_name, levels = unique(ginfo$gene_name))
    if (!is.null(sample.order)) {
      exp.long$variable <- factor(exp.long$variable, levels = sample.order)
    }
    if (!is.null(termID.order)) {
      exp.long$ID <- factor(exp.long$ID, levels = termID.order)
    }
    ght <- ggplot2::ggplot(exp.long) + ggplot2::geom_tile(ggplot2::aes(x = gene_name, 
                                                                       y = variable, fill = value), color = NA, show.legend = ht.legend) + 
      ggplot2::theme_bw(base_size = 14) + ggplot2::coord_cartesian(expand = 0) + 
      ggplot2::scale_fill_gradient2(low = exp.col[1], 
                                    mid = exp.col[2], high = exp.col[3], midpoint = 0, 
                                    name = "Z-Score") + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                                                                           vjust = 0.5, hjust = 1, size = ght.geneText.size), 
                                                                       axis.text = ggplot2::element_text(color = "black"), 
                                                                       axis.ticks.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), 
                                                                       plot.margin = ggplot2::margin(t = -0.1, r = 0.2, 
                                                                                                     b = 0.2, l = 0.2, unit = "cm")) + ggplot2::scale_y_discrete(position = "right") + 
      ggplot2::xlab("") + ggplot2::ylab("")
    if (ght.facet == TRUE) {
      fght <- ght + ggplot2::facet_wrap(~ID, ncol = 1, 
                                        scales = ght.facet.scale, strip.position = "left") + 
        ggplot2::theme(strip.background = ggplot2::element_rect(color = NA, 
                                                                fill = "grey90"), strip.placement = "outside")
    }
    else {
      fght <- ght
    }
  }
  if (newGsea == FALSE) {
    if (subPlot == 1) {
      pres <- pLabelOut
    }
    else if (subPlot == 2) {
      if (rmHt == FALSE) {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg_ht), ncol = 1, heights = c(0.8, 0.2))
      }
      else {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg), ncol = 1, heights = c(0.8, 0.2))
      }
    }
    else if (subPlot == 3) {
      if (rmHt == FALSE) {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg_ht1, prank), ncol = 1, heights = c(0.5, 
                                                                                       0.2, 0.3))
      }
      else {
        pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                               pseg, prank), ncol = 1, heights = c(0.5, 0.2, 
                                                                                   0.3))
      }
    }
    else {
      message("Please give 1/2/3 parameters!")
    }
  }
  else {
    if (rmHt == FALSE) {
      pres <- aplot::plot_list(gglist = list(pLabelOut, 
                                             ht), ncol = 1, heights = c(0.9, 0.1))
    }
    else {
      pres <- pLabelOut
    }
  }
  if (add.geneExpHt == TRUE) {
    pfinal <- aplot::plot_list(gglist = list(pres + ggplot2::xlab(""), 
                                             fght), ncol = 1, heights = c(1 - ght.relHight, ght.relHight))
  }
  else {
    pfinal <- pres
  }
  return(pfinal)
}