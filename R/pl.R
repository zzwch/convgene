#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Use dotplot/bubbleplot to visualize average/median and percentage expression of features (eg. genes) in groups
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.
#' @param features features to be visualized
#' @param group_by a vector including cell annotations
#' @param groups only these selected groups to be visualized
#' @param method using `base::mean` or `stats::median` or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If `TRUE` log_base will be ignored.
#' @param log_base `base::exp(1)` or `2` or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set `0.99` to exclude top 1% cells. (default: 1)
#' @param cap_value the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)
#' @param swap_axes swap x-axis and y-axis or not. (default: False)
#' @param plot_gene_force Force plot gene, which is missing in data.
#' @param gradientn_colors the color vector for visualization.
#' @param title title to display in the top.
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
pl_dotplot <- function(
  data,
  features,
  group_by,
  groups = NULL,
  method = base::mean,
  by_log = F,
  log_base = base::exp(1),
  exclude_zero = F,
  outlier_cutoff = 1,
  cap_value = NULL,
  swap_axes = F,
  plot_gene_force = F,
  gradientn_colors = c("grey90","grey", "orange", "red","brown"),
  title = NULL){

  # swap axes
  if(!swap_axes){
    features <- rev(features)
  }

  # do summarise
  ret_tmp <- SummariseDataByGroup(data, features,
                                  group_by, groups,
                                  method, by_log, log_base,
                                  exclude_zero, outlier_cutoff,
                                  cap_value, gene_force = plot_gene_force)
  meanData <- ret_tmp$mean; percData <- ret_tmp$percentage; rm(ret_tmp)

  # prepare data to plot
  dotData <- cbind(reshape2::melt(meanData), reshape2::melt(percData))[,-c(4,5)]
  colnames(dotData) <- c("Group", "Features", "Average", "Proportion")
  if(swap_axes) dotData$Group <- factor(dotData$Group, levels =  rev(levels(dotData$Group)))

  # set cap_value
  if(!IsNULLorNA(cap_value)){
    dotData$Average[dotData$Average > cap_value] <- cap_value
  }else{
    cap_value <- NA
  }

  # set x, y axis
  xx <- ifelse(swap_axes, "Features", "Group")
  yy <- ifelse(swap_axes, "Group", "Features")

  # do plot
  p <- ggplot(dotData) +
    geom_point(mapping = aes_(as.name(xx), as.name(yy),
                              color = as.name("Average"),
                              size = as.name("Proportion"))) +
    scale_color_gradientn(colours = gradientn_colors,
                          limits = c(0, cap_value),
                          breaks = if(IsNULLorNA(cap_value)) waiver() else seq(0,cap_value, 1))+
    #scale_color_gradientn(colours = c("grey90","grey",brewer.pal(n = 9, name = "YlOrRd")[-c(1,2,3,5,7,8)]))+
    scale_size_continuous(range = c(0,5), limits = c(0,1))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(linetype = 1, color = "black", fill = NA)) +
    labs(x = NULL, y = NULL, title = title)
  # return plot
  return(p)
}


#' Voinlin and Boxplot in one plot
#'
#' @param meta_data data.frame
#' @param x a string in colnames of meta_data, to be mapped as x-axis
#' @param y a string in colnames of meta_data, to be mapped as y-axis
#' @param fill_colors fill colors of violins
#' @param box_color border color of boxplot
#' @param box_fill_color fill color of boxplot
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
pl_vioboxplot <- function(
  meta_data,
  x,
  y,
  fill_colors = NULL,
  box_color = "black",
  box_fill_color = "white"){

  p <- ggplot(meta_data, mapping = aes_string(x,y)) +
    geom_violin(mapping = aes_string(fill = x),
                show.legend = F,
                color= 'white',
                scale = "width") +
    geom_boxplot(color = box_color,
                 show.legend = F,
                 fill = box_fill_color,
                 width = 0.25,
                 outlier.shape = NA)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(!IsNULLorNA(fill_colors)){
    p <- p + scale_fill_manual(values = fill_colors)
  }

  return(p)
}


#' stacked Violin Plots in one plot
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.
#' @param features features to be visualized
#' @param group_by a vector including cell annotations
#' @param groups only these selected groups to be visualized
#' @param ncol number of columns to plot
#' @param facet_scales scales of facet_wrap
#' @param fill_colors violion fill colors
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
pl_stackedViolinPlot <- function(
  data,
  features,
  group_by,
  groups = NULL,
  fill_colors = NULL,
  ncol = 1,
  facet_scales = "free_y"){

  # data <- as.matrix(data)
  ret_tmp <- SubsetDataAndGroup(data, group_by, groups)
  data <- ret_tmp$data; group_by <- ret_tmp$group_by; rm(ret_tmp)

  ggData <- reshape2::melt(
    cbind(as.data.frame(t(as.matrix(data[features,,drop = F]))),
          Group = group_by),
    id.vars = "Group",
    measure.vars = features,
    variable.name = 'Gene',
    value.name = 'Expression')

  p <- ggplot(ggData) +
    geom_violin(mapping = aes_string(x = "Group", y = 'Expression', fill = "Group"),
                scale = 'width') +
    facet_wrap(~Gene, ncol=ncol, strip.position = "right", scales = facet_scales) +
    #coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(linetype = 1, color = "black", fill = NA),
          legend.position = 'none')
  if(!IsNULLorNA(fill_colors)){
    p <- p + scale_fill_manual(values = fill_colors)
  }
  return(p)
}


#' Heatmap of average expression of select features
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.
#' @param features features to be visualized
#' @param group_by a vector including cell annotations
#' @param groups only these selected groups to be visualized
#' @param method using `base::mean` or `stats::median` or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If `TRUE` log_base will be ignored.
#' @param log_base `base::exp(1)` or `2` or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set `0.99` to exclude top 1% cells. (default: 1)
#' @param heat_color heatmap colors
#' @param scale_feature whether do scale by feature
#' @param cap_max if scale_feature, set cap for max
#' @param cap_min if scale_feature, set cap for min
#' @param cluster_rows row clustering
#' @param cluster_cols column clustering
#' @param show_rownames show rownames
#' @param ...
#'
#' @return a pheatmap object
#' @export
#'
#' @examples
pl_averageHeatmap <- function(
  data, features, group_by, groups = NULL,
  method = base::mean,
  by_log = F, log_base = base::exp(1), exclude_zero = F, outlier_cutoff = 1,
  heat_color = c("darkblue", "royalblue","grey100","tomato", "brown"),
  scale_feature = T, cap_max = 2, cap_min = -2,
  cluster_rows = F, cluster_cols = F,
  show_rownames = F, ...){

  # do summarise
  ret_tmp <- SummariseDataByGroup(data, features,
                                     group_by, groups,
                                     method, by_log, log_base,
                                     exclude_zero, outlier_cutoff,
                                     cap_value = NULL, gene_force = NULL)
  meanData <- t(ret_tmp$mean); rm(ret_tmp)


  if(scale_feature){
    meanData <- t(scale(t(meanData)))
    meanData[meanData > cap_max] <- cap_max
    meanData[meanData < cap_min] <- cap_min
  }

  pheatmap::pheatmap(meanData,
           color = colorRampPalette(heat_color,
                                    interpolate = "spline", alpha = T)(99),
           cluster_cols = cluster_cols, cluster_rows = cluster_rows,
           show_rownames = show_rownames, ...)
}



#' Heatmap of a 2D-table
#'
#' @param tab a table object
#' @param palette palette name in RColorBrewer or a vector of gradient colors
#' @param n length of gradient colors
#' @param min less than `min` will be setted to `NA`
#' @param title plot title
#' @param xlab plot x lab
#' @param ylab plot y lab
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
pl_tableHeatmap <- function(tab, palette = "YlGn", n = 5, min = 0, title = NULL, xlab = NULL, ylab = NULL){
  # ggData
  ggTab <- reshape2::melt(tab)
  ggTab$value[ggTab$value == 0] <- NA
  ggTab[,1] <- factor(ggTab[,1],
                      levels = rev(levels(ggTab[,1])))
  # ggTab[,2] <- factor(ggTab[,2],
  #                     levels = levels(ggTab[,2]))

  # background colors
  colours <- CheckBrewerPal(palette, n = n)
  # axis
  xx <- colnames(ggTab)[2]
  yy <- colnames(ggTab)[1]

  p <- ggplot(ggTab, mapping = aes_(x = as.name(xx), y = as.name(yy))) +
    geom_tile(mapping = aes(fill = log(value+1)
                            #, height = sqrt((value+1)/max(value)), width = sqrt((value+1)/max(value))
    ),
    color = "grey",
    show.legend = F) +
    #geom_tile(fill = NA, color = "grey") +
    labs(title = title, x = xlab, y = ylab) +
    geom_text(mapping = aes(label = value),
              data = subset(ggTab, !is.na(ggTab$value)),
              color = ifelse(ggTab$value[!is.na(ggTab$value)] >= min, "black", "grey")) +
    scale_fill_gradientn(colours = colours, na.value = "white") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          panel.border = element_rect(linetype = 1, color = "black", fill = NA))

  return(p)
}

