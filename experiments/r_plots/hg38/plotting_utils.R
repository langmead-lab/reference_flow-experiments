# Borrowed from: https://rpubs.com/sjackman/grid_arrange_shared_legend
# Thanks to Shaun Jackman
grid_arrange_shared_legend <- function(show_legend, num_rows, legend_plot_id, ...) {
  plots <- list(...)
  if (show_legend) {
    g <- ggplotGrob(plots[[legend_plot_id]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    plots_arg <- lapply(plots, function(x)  x + theme(legend.position="none"))
    plots_arg$nrow <- num_rows
    grid.arrange(
      do.call(arrangeGrob, plots_arg),
      legend,
      ncol = 1,
      heights = unit.c(unit(1, "npc") - lheight, lheight))
  }
  else {
    g <- ggplotGrob(plots[[legend_plot_id]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    grid.arrange(
      do.call(arrangeGrob, lapply(plots, function(x)
        x + theme(legend.position="none"))),
      ncol = 1,
      heights = unit.c(unit(1, "npc")))
  }
}

# Filter allelic bias dfs by number of overlapping reads
filter_df <- function (df, cov_thrsd = 16){
  df <- df[df$NUM_READS >= cov_thrsd,]
  df$REFBIAS <- df$REF_COUNT / (df$REF_COUNT + df$ALT_COUNT)
  df$AVG_MAPQ <- df$SUM_MAPQ / df$NUM_READS
  # df <- df[df$AVG_MAPQ > 10,]
  return (df)
}

# Get bias stats for an allelic bias data frame
get_bias_stats <- function(df, diff=0.2){
  return (c(
    nrow(df[df$REFBIAS >= 1-diff,]), 
    nrow(df[df$REFBIAS <= diff,]), 
    nrow(df[df$REFBIAS <= diff,]) + nrow(df[df$REFBIAS >= 1-diff,]), 
    sum(df$NUM_READS) / nrow(df),
    sum(df$REF_COUNT) / sum(df$ALT_COUNT)
  ))
}

# Build histograms
build_hist <- function(df, hist_step){
  half_step <- hist_step / 2
  histogram <- hist(df$REFBIAS, breaks=seq(-half_step, 1+half_step, by=hist_step))$counts
  histogram <- c(histogram, histogram[length(histogram)])
  return (histogram)
}

# Plot contour-only histograms
plot_contour_hist <- function (hist_targets, hist_names, prefix, ylimits, y_tran=F){
  hist_bins <- rep(adjusted_bins, length(hist_names))
  df_hist_onepass <- data.frame(xx = hist_targets, yy = rep(hist_names, each=length(hist_targets)/length(hist_names)))
  
  if (y_tran == T){
    p_hist <- ggplot() + 
      scale_x_continuous(trans = squish_trans(0.2, 0.8, 10)) +
      # scale_y_continuous(limits=c(0,3000)) +
      scale_y_continuous(trans = 'log10', limits=ylimits) +
      theme_bw() + 
      xlab("Allelic balance") + ylab("Count") + labs(color='Method') + 
      geom_step(mapping=aes(x=hist_bins, y=df_hist_onepass$xx, colour=df_hist_onepass$yy), size=1) +
      theme(legend.position = "bottom", aspect.ratio=NULL)
  }
  else{
    p_hist <- ggplot() + 
      # scale_x_continuous(trans = squish_trans(0.2, 0.8, 10)) +
      # scale_y_continuous(limits=c(0,3000)) +
      scale_y_continuous(trans = 'log10', limits=ylimits) +
      theme_bw() + 
      xlab("Allelic balance") + ylab("Count") + labs(color='Method') + 
      geom_step(mapping=aes(x=hist_bins, y=df_hist_onepass$xx, colour=df_hist_onepass$yy), size=1) +
      theme(legend.position = "bottom", aspect.ratio=NULL)
  }
  # if (save){
  #   ggsave(paste(prefix, ".jpg", sep=""), p_hist, width = 12, height = 8)
  #   ggsave(paste(prefix, ".pdf", sep=""), p_hist, width = 12, height = 8)
  # }
  return (p_hist)
}

library(scales)
squish_trans <- function(from, to, factor) {
  
  trans <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    
    return(x)
  }
  
  inv <- function(x) {
    
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    
    return(x)
  }
  
  # return the transformation
  return(trans_new("squished", trans, inv))
}