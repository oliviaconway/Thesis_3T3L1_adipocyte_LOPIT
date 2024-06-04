plotTranslocations <- function (params, type = "alluvial", all = FALSE, fcol, col, 
                                labels = TRUE, labels.par = "adj", cex = 1, spacer = 4, table = FALSE, 
                                ...) 
{
  
  stopifnot(inherits(params, "bandleParams") | inherits(params, 
                                                        "list"))
  type <- match.arg(type, c("alluvial", "chord"))
  labels.par <- match.arg(labels.par, c("adj", "repel"))
  if (inherits(params, "bandleParams")) {
    res1 <- summaries(params)[[1]]@posteriorEstimates$bandle.allocation
    res2 <- summaries(params)[[2]]@posteriorEstimates$bandle.allocation
    cl1 <- colnames(summaries(params)[[1]]@bandle.joint)
    cl2 <- colnames(summaries(params)[[2]]@bandle.joint)
  }
  if (inherits(params, "list")) {
    stopifnot(unlist(lapply(params, function(z) inherits(z, 
                                                         "MSnSet"))))
    params <- commonFeatureNames(params)
    params <- list(params[[1]], params[[2]])
    if (missing(fcol)) 
      stop(paste("Missing fcol, please specify feature columns"))
    if (length(fcol) == 1) {
      fcol <- rep(fcol, 2)
      message(paste0(c("------------------------------------------------", 
                       "\nIf length(fcol) == 1 it is assumed that the", 
                       "\nsame fcol is to be used for both datasets", 
                       "\nsetting fcol = c(", fcol[1], ",", fcol[2], 
                       ")", "\n----------------------------------------------")))
    }
    for (i in seq(fcol)) {
      if (!is.null(fcol[i]) && !fcol[i] %in% fvarLabels(params[[i]])) 
        stop("No fcol found in MSnSet, please specify a valid fcol ", 
             immediate. = TRUE)
    }
    res1 <- fData(params[[1]])[, fcol[1]]
    res2 <- fData(params[[2]])[, fcol[2]]
    cl1 <- names(table(res1))
    cl2 <- names(table(res2))
  }
  fct.lev <- union(cl1, cl2)
  res1_lev <- factor(res1, fct.lev)
  res2_lev <- factor(res2, fct.lev)
  dat <- data.frame(x = res1_lev, y = res2_lev)
  dat$z <- 1
  datdf <- dat %>% group_by(x, y, .drop = FALSE) %>% dplyr:::summarise(count = sum(z), 
                                                                       .groups = "keep")
  if (!all) {
    torm <- which(datdf$x == datdf$y)
    datdf <- datdf[-torm, ]
  }
  df <- as.data.frame(datdf)
  if (missing(col)) {
    setStockcol(NULL)
    grid.col <- setNames(getStockcol()[seq(fct.lev)], 
                                   fct.lev)
    if (length(fct.lev) > length(getStockcol())) 
      grid.col <- setNames(rainbow(length(fct.lev)), 
                                     fct.lev)
  }
  else {
    if (length(fct.lev) > length(col)) 
      stop(message("Not enough colours specified for subcellular classes"))
    grid.col <- col
    if (is.null(names(col))) {
      names(grid.col) <- fct.lev
      warning(message("Please specify a named character of colours\n to ensure matching of correct colours to subcellular classes "))
    }
  }
  if (type == "chord") {
    circos.clear()
    par(mar = c(1, 1, 1, 1) * spacer, cex = 1, xpd = NA)
    circos.par(gap.degree = 4)
    chordDiagram(df, annotationTrack = "grid", preAllocateTracks = 1, 
                 grid.col = grid.col, directional = 1, direction.type = c("diffHeight", 
                                                                          "arrows"), link.arr.type = "big.arrow", ...)
    if (labels) {
      circos.trackPlotRegion(track.index = 1, panel.fun = function(x, 
                                                                   y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(CELL_META$xcenter, ylim[1] + cm_h(2), 
                    sector.name, facing = "clockwise", niceFacing = TRUE, 
                    adj = c(0, 0.5), cex = cex, col = grid.col[sector.name], 
                    font = 2)
        circos.axis(h = "bottom", labels.cex = 0.6, sector.index = sector.name)
      }, bg.border = NA)
    }
    else {
      circos.trackPlotRegion(track.index = 1, panel.fun = function(x, 
                                                                   y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.axis(h = "top", labels.cex = 0.6, major.tick.length = 1, 
                    sector.index = sector.name, track.index = 2)
      }, bg.border = NA)
    }
    circos.clear()
  }
  if (type == "alluvial") {
    torm <- which(df$count == 0)
    df <- df[-torm, ]
    names(df) <- c("Condition1", "Condition2", "value")
    levs1 <- levels(df$Condition1)
    levs2 <- levels(df$Condition2)
    res1 <- unique(df$Condition1)
    res2 <- unique(df$Condition2)
    cond1_col <- grid.col[levs1[levs1 %in% res1]]
    cond2_col <- grid.col[levs2[levs2 %in% res2]]
    columncol <- c(cond1_col, cond2_col)
    stratcol <- c(rev(cond1_col), rev(cond2_col))
    df_expanded <- df[rep(row.names(df), df$value), ]
    df_expanded <- df_expanded %>% mutate(id = row_number()) %>% 
      pivot_longer(-c(value, id), names_to = "Condition", 
                   values_to = "label")
    if (table == TRUE) {
      return(df = df)
      stop("Returning translocation data.frame")
    }
    q <- ggplot(df_expanded, aes(x = Condition, stratum = label, 
                                 alluvium = id, fill = label)) + geom_flow(width = 0) + 
      scale_fill_manual(values = columncol) + scale_color_manual(values = stratcol) + 
      geom_stratum(width = 1/8, color = "white") + scale_x_discrete(expand = c(0.25, 
                                                                               0.25)) + scale_y_continuous(breaks = NULL) + theme_minimal() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(size = 12), panel.grid.major.y = element_blank(), 
            panel.grid.major.x = element_blank(), axis.title.x = element_blank()) + 
      theme(legend.position = "none") + ylab(NULL)
    if (labels == "TRUE") {
      if (labels.par == "adj") {
        q <- q + geom_text(aes(label = after_stat(stratum), 
                               hjust = ifelse(Condition == "Condition1", 1, 
                                              0), x = as.numeric(factor(Condition)) + 0.075 * 
                                 ifelse(Condition == "Condition1", -1, 1), 
                               color = after_stat(stratum)), stat = "stratum", 
                           fontface = "bold", size = 4)
      }
      if (labels.par == "repel") {
        q <- q + ggrepel::geom_text_repel(aes(label = ifelse(after_stat(x) == 
                                                               1, as.character(after_stat(stratum)), "")), 
                                          stat = "stratum", size = 4, direction = "y", 
                                          nudge_x = -0.6) + ggrepel::geom_text_repel(aes(label = ifelse(after_stat(x) == 
                                                                                                          2, as.character(after_stat(stratum)), "")), 
                                                                                     stat = "stratum", size = 4, direction = "y", 
                                                                                     nudge_x = 0.6)
      }
      labels.par <- match.arg(labels.par, c("repel", "adj"))
    }
    q
  }
}
