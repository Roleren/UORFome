# TODO: Option to remove 3' UTR tail, something is strange

library(ggplot2)
library(gridExtra)

#' Get coverage window plot of reads
#' 
#' Spanning a region like a transcripts, plot how the reads distribute.
#' 
#' If you return this function without assigning it and output is NULL,
#' it will automaticly plot the figure in your session. If output is assigned,
#' no plot will be shown in session. 
#' @param coverage a data.table, output of scaledWindowCoverage
#' @param output character vector (NULL), if set, saves the plot as pdf
#' to path given.
#' @param scoring character vector (zscore), either of zScore, 
#' transcriptNormalized, sum, mean
#' @param colors character vector colors to use in plot
#' @return a ggplot object of the coverage plot, NULL if output is set
windowCoveragePlot <- function(coverage, output = NULL, scoring = "zscore",
                               colors = c('skyblue4', 'orange')) {

  coverage$feature  <- factor(coverage$feature, levels = unique(coverage$feature), 
                              labels = unique(coverage$feature))
  coverage$fraction <- factor(coverage$fraction, levels = unique(coverage$fraction), 
                              labels = unique(coverage$fraction))
  
  # Keep only genes with coverage > 1 on SSU
  # ifelse("SSU" %in% unique(coverage$fraction), hit <- "SSU", hit <- unique(coverage$fraction))
  
  # valid_genes <- unique(coverage[feature == "",]$genes)
  # unique_valid_genes <- coverage[, .(unique = unique(genes)), by = list(fraction, feature)]
  # unique_valid_genes <- unique_valid_genes[, .(unique %in% unique), by = fraction]
  # coverage_sum_subset <- coverage[genes %in% unique_valid_genes,]
  
  coverage[, `:=` (gene_sum = sum(count)), by = list(genes, fraction)]
  coverage_sum_subset <- coverage
  # unique_valid_genes <- unique(coverage[fraction == hit &
  #                                       gene_sum > 0 ,]$genes)
  
  
  if (scoring == "zscore") {
    # z score over transcript per fraction
    coverage_sum_subset[, `:=` (windowSD = sd(count), windowMean = mean(count)), by = list(genes, fraction)] 
    coverage_sum_subset[, zscore := (count-windowMean)/windowSD]
    # create mean and sum scores per position, per feature
    coverage_score <- coverage_sum_subset[, .(score = mean(zscore, na.rm = TRUE)),
                                           by = list(fraction, position, feature)]
  } else if (scoring == "transcriptNormalized") {
    coverage_score <- coverage_sum_subset[, .(score = sum(count / gene_sum, na.rm = TRUE)),
                                          by = list(fraction, position, feature)]
    
  } else if (scoring == "mean") {
    coverage_score <- coverage_sum_subset[, .(score = mean(count, na.rm = TRUE)),
                                          by = list(fraction, position, feature)]
  } else if (scoring == "sum") {
    coverage_score <- coverage_sum_subset[, .(score = sum(count, na.rm = TRUE)),
                                          by = list(fraction, position, feature)]
  } else stop(paste("Invalid scoring: ", scoring))
  
  coverage_score[, `:=` (fraction_min=min(score)), by = list(fraction)]
  
  plot <- ggplot(data=as.data.frame(coverage_score), aes(x=position, ymax=score, ymin=fraction_min, y=score, colour = as.factor(fraction))) +
    geom_ribbon(stat="identity", position = "identity", aes(fill= as.factor(fraction), alpha=0.5)) +
    geom_line() +
    theme_bw() +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    ggtitle(label = "TCP-seq", subtitle = paste0("Genes n=", length(unique(coverage_sum_subset$genes))) ) +
    xlab("Scaled position in transcript") + ylab(paste0(scoring, " over transcript")) +
    theme(legend.position="none") +
    facet_grid(fraction ~ feature, scales = "free")

  if (!is.null(output)) {
    if(is.character(output) && dir.exists(dirname(outName))) {
      if (tools::file_ext(output) != "pdf") output <- paste0(output, ".pdf")
      ggsave(output, plot = plot, width=200, height=150, units="mm",
             dpi=100, limitsize = FALSE)
    } else {
      stop("output does not name a valid directory")
    }
    return(NULL)
  } else return(plot)
}


