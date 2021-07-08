#' @include utils.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' plot heatmap using pheatmap
#'
#' @param data feature by cell data matrix
#' @param meta_data cell annotation data frame
#' @param annot_col select colnames of meta_data
#' @param feature_data feature annotation data frame
#' @param annot_row select colnames of feature_data
#' @param annot_colors a list of annotation colors named with annotation categories
#' @param features features to display
#' @param color heatmap colors
#' @param permutation_by permutate and order the cells
#' @param seed permutation seed
#' @param smooth_n smooth the expression over smooth_n cells
#' @param do_scale scale data
#' @param cap_max capping the maximum of data
#' @param cap_min capping the minimum of data
#' @param ... other params passed to pheatmap
#' @param smooth_lowess smooth using lowess regression
#' @param lowess_span lowess span
#' @param rownames_to_show select rownames to show to the right
#'
#' @return a pheatmap
#' @import pheatmap
#' @export
#'
#' @examples
#'
pl_heatmap <- function(
  data, features = NULL, color = gradient_colors$f_BuWtRd(99, interpolate = "linear"),
  meta_data = NA, annot_col = NA,
  permutation_by = NA, seed = 666,
  feature_data = NA, annot_row = NA, smooth_n = NA,
  do_scale = F,
  smooth_lowess = F, lowess_span = 2/3,
  cap_max = NA, cap_min = NA,
  annot_colors = NA,
  rownames_to_show = NA,
  ...){

  if(IsNULLorNA(features)){
    features <- rownames(data)
  }
  if(!IsNULLorNA(meta_data)) {
    if(!IsNULLorNA(permutation_by)) {
      set.seed(seed)
      meta_data <- meta_data[sample(1:nrow(meta_data)),]
      meta_data <- meta_data[do.call(what = order, meta_data[,permutation_by,drop = F]),]
    }
    if(!IsNULLorNA(annot_col)){
      meta_data <- meta_data[, annot_col,drop=F]
    }
  }else{
    meta_data <- NA
  }

  if(!IsNULLorNA(feature_data)){
    if(!IsNULLorNA(annot_row)){
      feature_data <- feature_data[, annot_row, drop= F]
    }
  }else{
    feature_data <- NA
  }

  phData <- data[features, rownames(meta_data)]
  if(do_scale) phData <- Matrix::t(scale(Matrix::t(phData)))
  if(!IsNULLorNA(smooth_n)) phData <- rowSmooth(phData, smooth_n)
  if(smooth_lowess) {
    tmp <- t(apply(phData, 1, function(x) stats::lowess(x, f = lowess_span)$y))
    dimnames(tmp) <- dimnames(phData)
    phData <- tmp
  }

  phData <- Seurat::MinMax(phData, min = cap_min, max = cap_max)

  phLabels_row <- rownames(phData)
  if(IsNULLorNA(rownames_to_show)){
    show_rowname <- F
  }else{
    phLabels_row[! phLabels_row %in% rownames_to_show] <- ""
    show_rowname <- T
  }

  pheatmap(mat = phData,
           color = color,
           annotation_col = meta_data,
           annotation_row = feature_data,
           annotation_colors = annot_colors,
           show_rownames = show_rowname,
           labels_row = phLabels_row,
           ...)
}

#' DimPlot and Slingshot curves
#'
#' @param object Seurat3 object
#' @param ggplot a ggplot object as background
#' @param color_by slingshot curve color_by
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#'
pl_slingshot <- function(object, ggplot, color_by = "ss_lineage", pt.size = 1){
  dr <- grep(pattern = "ss[0-9]+", Reductions(object), value = T, perl = T, ignore.case = T)
  p <- ggplot + ggnewscale::new_scale_color()
  for(i in 1:length(dr)){
    ss_dim <- cbind(na.omit(as.data.frame(Seurat::Embeddings(object, reduction = dr[i]))),
                    ss_lineage = stringr::str_replace(dr[i] , pattern = "ss", "Lineage"))
    ndim <- ncol(ss_dim)-1
    colnames(ss_dim)[1:ndim] <- paste0("ss_", 1:ndim)

    dim1 <- colnames(ss_dim)[1]
    dim2 <- colnames(ss_dim)[2]

    p <- p + geom_line(mapping = aes_(as.name(dim1), as.name(dim2),
                                      color = as.name(color_by)),
              data = cbind(ss_dim, object@meta.data[rownames(ss_dim),]),
              size = 2)
  }

  p <- p + guides(color = guide_legend(override.aes = list(size = 3*pt.size)))
  return(p)
}

#' plot x-y scatters
#'
#' @param meta_data a data.frame including x, y, group_by
#' @param x colname of meta_data
#' @param y colname of meta_data
#' @param group_by identical to point_group_by
#' @param groups identical to groups
#' @param show_points show geom_point
#' @param show_ellipse show stat_ellipse
#' @param ellipse_level set ellipse area level
#' @param ellipse_alpha ellipse aes, a value of 0 to 1
#' @param ellipse_color ellipse aes, a color
#' @param show_center show mean point
#' @param center_method function used to calc mean center
#' @param center_size center point aes, a value
#' @param center_shape center point aes, a value
#' @param center_alpha center point aes, a value
#' @param colors modify the point colors
#' @param point_group_by color the points by which vector
#' @param point_groups subset data to show only those points belonging to groups of group_by
#' @param point_shape_by point aes, a value or colname of meta_data
#' @param point_size_by point aes, a value or colname of meta_data
#' @param ellipse_group_by calc ellipse by which vector
#' @param ellipse_groups subset data to calc ellipse for only those points belonging to groups of group_by
#' @param center_group_by calc center by which vector
#' @param center_groups subset data to calc center for only those points belonging to groups of group_by
#' @param ... params passed to center_method
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#'
pl_scatterplot <- function(
  meta_data, x, y,

  show_points = T, point_shape_by = NULL, point_size_by = NULL,
  group_by = NULL,
  groups = NULL,

  point_group_by = group_by,
  point_groups = groups,

  show_ellipse = F, ellipse_level = 0.5, ellipse_alpha = 0.8, ellipse_color = "white",
  ellipse_group_by = point_group_by,
  ellipse_groups = point_groups,

  show_center = F, center_method = SummariseExpr,
  center_size = 3, center_shape = 1, center_alpha = 1,
  center_group_by = ellipse_group_by,
  center_groups = ellipse_groups,

  colors = scanpy_colors$default_64, ...){

  # POINT
  point_data <- meta_data
  ## point shape
  if(IsNULLorNA(point_shape_by) || (!point_shape_by[[1]] %in% colnames(point_data))){
    point_data <- SetDefaultAes(point_data, point_shape_by, default = "shape", name = "shape")
    point_shape_by <- "shape"
  }else{
    point_shape_by <- point_shape_by[[1]]
  }
  ## point size
  if(IsNULLorNA(point_size_by) || (!point_size_by[[1]] %in% colnames(point_data))){
    point_data <- SetDefaultAes(point_data, point_size_by, default = "size", name = "size")
    point_size_by <- "size"
  }else{
    point_size_by <- point_size_by[[1]]
  }

  ## point group
  if(IsNULLorNA(point_group_by) || (!point_group_by[[1]] %in% colnames(point_data))){
    point_data <- SetDefaultAes(point_data, point_group_by, default = "group", name = "group")
    point_group_by <- "group"
  }else{
    point_group_by <- point_group_by[[1]]
    if(!IsNULLorNA(point_groups)){
      point_data <- subset(point_data, point_data[,point_group_by] %in% point_groups)
    }
  }



  # ELLIPSE
  ellipse_data <- meta_data
  ## ellipse group
  if(IsNULLorNA(ellipse_group_by) || (!ellipse_group_by[[1]] %in% colnames(ellipse_data))){
    ellipse_data <- SetDefaultAes(ellipse_data, ellipse_group_by, default = "group", name = "group")
    ellipse_group_by <- "group"
  }else{
    ellipse_group_by <- ellipse_group_by[[1]]
    if(!IsNULLorNA(ellipse_groups)){
      ellipse_data <- subset(ellipse_data, ellipse_data[,ellipse_group_by] %in% ellipse_groups)
    }
  }

  # CENTER
  center_data <- meta_data
  ## center group
  if(IsNULLorNA(center_group_by) || (!center_group_by[[1]] %in% colnames(center_data))){
    center_data <- SetDefaultAes(center_data, center_group_by, default = "group", name = "group")
    center_group_by <- "group"
  }else{
    center_group_by <- center_group_by[[1]]
    if(!IsNULLorNA(center_groups)){
      center_data <- subset(center_data, center_data[,center_group_by] %in% center_groups)
    }
  }

  p <- ggplot()
  # POINT
  if(show_points){
    p <- p + geom_point(mapping = aes_(x = as.name(x),
                                      y = as.name(y),
                                      color = as.name(point_group_by),
                                      shape = as.name(point_shape_by),
                                      size = as.name(point_size_by)),
                        data = point_data
                        )
  }
  # ELLIPSE
  if(show_ellipse){
    p <- p + stat_ellipse(mapping = aes_(x = as.name(x),
                                        y = as.name(y),
                                        fill = as.name(ellipse_group_by),
                                        group = as.name(ellipse_group_by)),
                          data = ellipse_data,
                          color = ellipse_color,
                          size = 1,
                          geom = "polygon",
                          level = ellipse_level,
                          alpha = ellipse_alpha)
  }

  if(show_center){
    meanData <- as.data.frame(apply(center_data[,c(x,y)], 2, function(x) {
      tapply(x , list(center_data[[center_group_by]]), function(y){
        center_method(y, by_log = T, ...)
      })
    }))

    meanData[,center_group_by] <- rownames(meanData)

    p <- p + geom_point(data = meanData,
                        mapping = aes_(x = as.name(x),
                                      y = as.name(y),
                                      color = as.name(center_group_by)),
               alpha = center_alpha,
               shape = center_shape,
               size = center_size)
  }

  if(!IsNULLorNA(colors)){
    p <- p +
      scale_color_manual(values = colors, na.value = "grey70") +
      scale_fill_manual(values = colors, na.value = "grey70")
  }
  return(p)
}

#' Use dotplot/bubbleplot to visualize average/median and percentage expression of features (eg. genes) in groups
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.
#' @param features features to be visualized
#' @param group_by a vector including cell annotations
#' @param groups only these selected groups to be visualized
#' @param method using base::mean or stats::median or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If TRUE log_base will be ignored.
#' @param log_base base::exp(1) or 2 or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set 0.99 to exclude top 1 percent cells. (default: 1)
#' @param cap_value the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)
#' @param swap_axes swap x-axis and y-axis or not. (default: False)
#' @param plot_gene_force Force plot gene, which is missing in data.
#' @param gradientn_colors the color vector for visualization.
#' @param title title to display in the top.
#' @param duplicate_features whether display duplicate features or not.
#' @param do_scale scale the mean value across groups. Useful to demonstrate the differences among groups.
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#'
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
  do_scale = F,
  plot_gene_force = F,
  duplicate_features = T,
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
                                  cap_value, do_scale,
                                  gene_force = plot_gene_force)
  meanData <- ret_tmp$mean; percData <- ret_tmp$percentage; rm(ret_tmp)

  colnames <- c("Group", "Features", "Average", "Proportion")

  # prepare data to plot
  if(duplicate_features) {
    features_ <- colnames(meanData); colnames(meanData) <- NULL; colnames(percData) <- NULL
    dotData <- cbind(reshape2::melt(meanData), reshape2::melt(percData))[,-c(4,5)]
    colnames(dotData) <- colnames
    dotData$Features <- factor(dotData$Features)
  }else{
    dotData <- cbind(reshape2::melt(meanData), reshape2::melt(percData))[,-c(4,5)]
    colnames(dotData) <- colnames
  }
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

  if(duplicate_features){
    p <- if(swap_axes) p + scale_x_discrete(labels = features_) else p + scale_y_discrete(labels = features_)
  }
  # return plot
  return(p)
}


#' Voinlin and Boxplot in one plot
#'
#' @param meta_data data.frame
#' @param x a string in colnames of meta_data, to be mapped as x-axis
#' @param y a string in colnames of meta_data, to be mapped as y-axis
#' @param fill_by a string in colnames of meta_data
#' @param fill_colors fill colors of violins
#' @param box_color border color of boxplot
#' @param box_fill_color fill color of boxplot
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#'
pl_vioboxplot <- function(
  meta_data,
  x,
  y,
  fill_by = x,
  fill_colors = NULL,
  box_color = "black",
  box_fill_color = "white"){

  if(class(meta_data) == "Seurat") meta_data <- FetchData(meta_data, vars = c(x, y, fill_by))
  p <- ggplot(meta_data, mapping = aes_string(x,y)) +
    geom_violin(mapping = aes_string(fill = fill_by),
                show.legend = F,
                color= NA,
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
#'
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
#' @param method using base::mean or stats::median or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If TRUE log_base will be ignored.
#' @param log_base base::exp(1) or 2 or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set 0.99 to exclude top 1 percent cells. (default: 1)
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
#'
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
#' @param min less than min will be setted to NA
#' @param title plot title
#' @param xlab plot x lab
#' @param ylab plot y lab
#'
#' @return a ggplot object
#' @import ggplot2
#' @export
#'
#' @examples
#'
pl_tableHeatmap <- function(tab, palette = "YlGn", n = 5, min = 0, title = NULL, xlab = NULL, ylab = NULL){
  # ggData
  ggTab <- reshape2::melt(tab)
  ggTab$value[ggTab$value == 0] <- NA

  if(!is.factor(ggTab[,1])) ggTab[,1] <- factor(ggTab[,1])
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


#' Dimensional reduction plot
#'
#' Enhanced DimPlot of Seurat v3
#'
#' @param object Seurat object
#' @param group_by Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
#' You can also plot features with gradient cols.
#' @param ncol Number of columns for display when faceting plots
#' @param subset_by enhanced option, subset cells by the colname.
#' @param subset_groups enhanced option, subset cells of these categories.
#' @param subset_cells A vector of cell names. if not NULL, subset these cells, otherwise use subset_by and subset_groups
#' @param reduction Which dimensionality reduction to use. e.g. umap, tsne, pca
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param pt.size Adjust point size for plotting
#' @param cols Vector of colors, each color corresponds to an identity class.
#' @param slot FetchData of features from which slot.
#' @param shape_by one column name to control point shape
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'
pl_dimplot <- function(object, group_by, shape_by = NULL, ncol = 2, slot = "data",
                      subset_by = NULL, subset_groups = NULL, subset_cells = NULL,
                      reduction = "umap", dims = c(1,2), pt.size = 1,
                      cols = scanpy_colors$default_64){

  if(is.null(subset_cells)){
    if(!is.null(subset_groups) && !is.null(subset_by)){
      subset_cells <- rownames(object@meta.data[object@meta.data[,subset_by] %in% subset_groups,])
    }
  }

  embeddings <- Embeddings(object, reduction = reduction)[,dims]
  dims <- colnames(embeddings)

  ggData <- reshape2::melt(cbind(Seurat::FetchData(object, c(group_by, shape_by), slot = slot), embeddings),
                           id.vars = c(dims, shape_by),
                           measure.vars = group_by,
                           variable.name = "group_by",
                           value.name = "value",
                           factorsAsStrings = F)
  ggData2 <- reshape2::melt(cbind(Seurat::FetchData(object, c(group_by, shape_by), slot = slot), embeddings)[subset_cells,],
                            id.vars = c(dims, shape_by),
                            measure.vars = group_by,
                            variable.name = "group_by",
                            value.name = "value",
                            factorsAsStrings = F)

  if(!is.null(subset_cells)){
    p <- ggplot() +
      geom_point(mapping = aes_string(dims[1], dims[2], shape = shape_by),
                 color = "grey80", data = ggData, size = pt.size) +
      geom_point(mapping = aes_string(dims[1], dims[2], shape = shape_by,
                                      color = "value",group = "group_by"),
                 data = ggData2, size = pt.size)
  }else{
    p <- ggplot() +
      geom_point(mapping = aes_string(dims[1], dims[2],
                                      color = "value",
                                      group = "group_by",
                                      shape = shape_by),
                 data = ggData, size = pt.size)
  }

  if(length(group_by) > 1){
    p <- p + facet_wrap(~group_by, ncol = ncol)
  }

  if(!is.null(cols)) {
    if(group_by %in% colnames(object@meta.data)){
      p <- p + scale_color_manual(values = cols)
    }else{
      p <- p + scale_color_gradientn(colors = cols)
    }
  }

  p <- p + guides(color = guide_legend(override.aes = list(size = 3*pt.size))) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_blank())

  return(p)
}


#' rastr version of partial FeaturePlot
#'
#' @param object seurat object
#' @param features genes
#' @param coord_fixed_ratio coord_fixed
#' @param ncol ncol of cowplot::plot_grid
#' @param align align of cowplot::plot_grid
#' @param cols colors of gradientn
#' @param pt.size point size
#' @param order sort orders of cells
#' @param reduction dimension reduction
#' @param dims dimensions
#'
#' @return ggplot
#' @export
#'
#' @examples
#'
pl_rastrFeaturePlot <- function(object, features, coord_fixed_ratio = 1, ncol = NULL, align = c("hv"),
                                cols = c("grey", "yellow", "red"), pt.size = 1, order = T, rasterise = T,
                                reduction = "umap", dims = c(1,2)){
  dimension <- Embeddings(object, reduction = reduction)[,dims]
  xy <- colnames(dimension)
  ggData <- cbind(dimension, FetchData(object, vars = features)) %>%
    as.data.frame()
  features_ <- gsub("-", "____", features)
  colnames(ggData) <- c(xy, features_)

  cowplot::plot_grid(
    plotlist = lapply(seq_along(features), function(i){
      x <- features_[i]
      ggplot(data = if(order) ggData %>% arrange(.data[[x]]) else ggData ,
             mapping = aes_string(x = xy[1], y = xy[2], color = x)) +
        {if(rasterise) ggrastr::rasterise(geom_point(size = pt.size), dpi = 300) else geom_point(size = pt.size)} +
        scale_color_gradientn(colours = cols) +
        coord_fixed(ratio = coord_fixed_ratio) +
        labs(title = features[i]) +
        guides(color = guide_colorbar(title = NULL)) +
        theme_classic() +
        theme(axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold.italic"),
              panel.border = element_rect(fill = NA, color = "black"))
    }),ncol = ncol, align = align)
}


#' Plot the expression of one gene using bar/vline, or so called barplot
#'
#' @param x numeric vector, such as expression of one gene in all cells.
#' @param cluster character vector of cell annotation, such as cell clusters.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'
pl_segmentPlot <- function(x, cluster){
  ggData <- data.frame(yend = x,
                       Cluster = cluster) %>%
    arrange(Cluster, yend) %>%
    mutate(x = seq_along(Cluster), y = 0)

  ggplot() +
    geom_segment(mapping = aes(x = x, y = y, xend = x, yend = yend, color = Cluster),
                 data = ggData ) +
    geom_vline(mapping = aes(xintercept = c(1, unname(cumsum(table(ggData$Cluster))))),
               color = "black", linetype = "dashed") +
    scale_x_continuous(expand = expansion()) +
    scale_y_continuous(expand = expansion(c(0, 0.1))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
}


#' plot multiple DimPlot for each batch
#'
#' This may be useful for displaying multiple datasets
#'
#' @param object Seurat object
#' @param reduction umap/tsne/pca
#' @param split_by split dimplot by this metadata
#' @param groups_label plot using DimPlot with label = TRUE
#' @param groups plot using pl_dimplot
#' @param coord_fixed_ratio only valid for groups
#' @param pt.size point size
#' @param colors point color
#'
#' @return multiple plots
#' @export
#'
#' @examples
#'
pl_splitDimplot <- function(object, reduction = "umap",
                            split_by = "Dataset",
                   groups_label = c("seurat_clusters"),
                   groups = c("Batch", "Phase"),
                   coord_fixed_ratio = 1, pt.size = 1,
                   colors = scanpy_colors$default_64){

  split_by <- object[[split_by]][[split_by]]
  for(i in sort(unique(split_by))){
    subset_cells <-  colnames(object)[split_by == i]

    for(j in groups_label){
      print(
        DimPlot(object, group.by = j, label = T,
                #cols = cluster_color,
                cols = colors,
                cells = subset_cells,
                pt.size = pt.size, reduction = reduction) +
          NoLegend() + coord_fixed()
      )
    }

    for(k in groups){
      print(
        pl_dimplot(object, group_by = k, subset_cells = subset_cells, reduction = reduction,
                   cols = colors,  pt.size = pt.size) +
          coord_fixed(coord_fixed_ratio)
      )
    }
  }
}


#' visualize genes patterns heatmap using kmeans clustering
#'
#' @param mat gene x cell expression matrix. cell order is kept for heatmap.
#' @param genes_select genes to be display
#' @param cellAnnot cell annotation data frame
#' @param geneAnnot gene annotation, used for annotation bar in heatmap.
#' @param cluster_colname which column to be used for group
#' @param clusters_display subset clusters to be displayed
#' @param clusters_stat subset clusters to be used for filtering genes by statistical test
#' @param filter_gene_by_stat whether to do filter gene
#' @param stat_pvalue threshold for filtering genes
#' @param stat_padj threshold for filtering genes. both pvalue and padjust will be fullfilment.
#' @param do_gene_kmeans use kmeans to determine gene order
#' @param kmeans_k k parameter in kmeans
#' @param kmeans_seed random seeds for kmeans
#' @param gene_scale do gene-wise scale for display
#' @param min cap the scaled value
#' @param max cap the scaled value
#' @param display_gaps do \code{gaps_row} and \code{gaps_col} of pheatmap
#' @param ... parameters used in pheatmap
#' @param return_data return mat and kmeans results
#'
#' @return pheatmap plot
#' @export
#'
#' @examples
#'
pl_genePatternHeatmap <- function(
  mat, genes_select, cellAnnot, geneAnnot,
  cluster_colname = "identity",
  clusters_display,
  clusters_stat = NULL,
  filter_gene_by_stat = F, stat_pvalue = 0.05, stat_padj = 1,
  do_gene_kmeans = T, kmeans_k = 10, kmeans_seed = 1,
  gene_scale = T, min = -2, max = 2,
  display_gaps= c("row", "col"),
  return_data = F,
  ...
){
  # cell order
  cellAnnot <- cellAnnot[colnames(mat), ] # keep column order of mat as is which will be displayed in heatmap

  # gene
  genes_notin <- paste0(setdiff(genes_select, rownames(mat)), collapse = "," )
  genes_select <- intersect(genes_select, rownames(mat))
  if(genes_notin != "") message(stringr::str_glue("NOT-IN-DATA: {genes_notin}"))

  # filter by stat
  if(is.null(clusters_stat)) clusters_stat <- clusters_display
  cellAnnot_stat <- cellAnnot[cellAnnot[[cluster_colname]] %in% clusters_stat, ]
  mat_stat <- mat[genes_select, rownames(cellAnnot_stat), drop = F]
  if(filter_gene_by_stat){
    gene_pvalue <- apply(mat_stat, 1, function(x){
      kruskal.test(x, cellAnnot_stat[[cluster_colname]])$p.value}) %>%
      sort()
    gene_padj <- p.adjust(gene_pvalue)
    genes_select <- names(gene_pvalue)[gene_pvalue <= stat_pvalue & gene_padj <= stat_padj]
  }

  #
  cellAnnot <- cellAnnot[cellAnnot[[cluster_colname]] %in% clusters_display, ]
  geneAnnot <- geneAnnot[genes_select,,drop = F]
  mat <- mat[genes_select, rownames(cellAnnot), drop = F]

  if(gene_scale){
    phMat <- t(MinMax(scale(t(mat)), min, max))
  }

  # kmeans
  if(do_gene_kmeans){
    set.seed(kmeans_seed)
    phk <- kmeans(
      # apply(phMat, 1, function(x)
      #   tapply(x,
      #          INDEX = cellAnnot[[cluster_colname]] %>% droplevels,
      #          FUN = mean)
      # ) %>% t,
      phMat,
      centers = kmeans_k)

    phk_hc <- hclust(((1- (phk$centers %>% t %>% cor))/2) %>% as.dist())

    geneAnnot[names(phk$cluster), "kmean"] <- plyr::mapvalues(
      phk$cluster,
      sort(unique(phk$cluster)),
      phk_hc$order %>% order)
    phMat <- phMat[geneAnnot$kmean %>% order,]
    gaps_genes <- table(geneAnnot$kmean) %>% cumsum()
  }

  gaps_cells <- table(cellAnnot[[cluster_colname]] %>% droplevels()) %>% cumsum()
  # heatmap
  if(return_data){
    pheatmap::pheatmap(
      phMat,
      annotation_col = cellAnnot,
      annotation_row = geneAnnot,
      gaps_row = if ("row" %in% display_gaps) gaps_genes else NULL,
      gaps_col = if ("col" %in% display_gaps) gaps_cells else NULL,
      ...)
    return(list(phMat = phMat,
                cellAnnot = cellAnnot,
                geneAnnot = geneAnnot))
  }else{
    pheatmap::pheatmap(
      phMat,
      annotation_col = cellAnnot,
      annotation_row = geneAnnot,
      gaps_row = if ("row" %in% display_gaps) gaps_genes else NULL,
      gaps_col = if ("col" %in% display_gaps) gaps_cells else NULL,
      ...)
  }
}

#' Visualize gene expression pattern using boxVlnPlot and smooth line
#'
#' @param mat gene x cell matrix
#' @param cell_group a vector of cells' group corresponding to colnames(mat)
#' @param gene_group a vector of genes' group corresponding to rownames(mat)
#' @param return_data if TRUE return summarised data, otherwise return ggplot
#'
#' @return ggplot or data.frame
#' @export
#'
#' @examples
#'
pl_genePatternDiagram <- function(mat, cell_group, gene_group, ncol = 1, return_data = F){
  lapply(sort(unique(gene_group)), function(gg){
    gg_summary <- apply(
      mat[gene_group == gg, , drop = F],
      1, function(x) { tapply(x, cell_group %>% droplevels, mean)}) %>% reshape2::melt()
    colnames(gg_summary) <- c("cluster", "gene", "expression")
    gg_summary$group <- gg

    return(gg_summary)
  }) %>% Reduce(f = rbind) -> ggData

  if (return_data) {
    return(ggData)
  }else{
    ggData %>%
      ggplot(mapping = aes(cluster, expression)) +
      geom_violin(show.legend = F, scale = "width", mapping = aes(fill = cluster)) +
      geom_boxplot(show.legend = F, outlier.shape = NA, fill = "white", width = 0.3) +
      geom_smooth(show.legend = F, mapping = aes(group = group)) +
      facet_wrap(~group, ncol = ncol, scales = "free_y") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
}
