
#' calculate average/median and percentage expression of features (eg. genes) in groups
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled. A matrix is better than data.frame in case of duplicated features.
#' @param features features to be visualized
#' @param group_by a vector including cell annotations
#' @param groups only these selected groups to be visualized
#' @param method using `base::mean` or `stats::median` or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If `TRUE` log_base will be ignored.
#' @param log_base `base::exp(1)` or `2` or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set `0.99` to exclude top 1% cells. (default: 1)
#' @param gene_force Force gene suing zeros, which is missing in data.
#' @param cap_value the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)
#'
#' @return a list of mean and percentage data
#' @export
#'
#' @examples
SummariseDataByGroup <- function(
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
  gene_force = F){

  # check whether features are in data
  gene_diff <- setdiff(features, rownames(data))
  if(length(gene_diff) > 0){
    warning(Message("NOT-IN-DATA: ", gene_diff, ", ", "\n"))
    if(gene_force){
      pseudo_data <- matrix(0, length(gene_diff), ncol = ncol(data),
                            dimnames = list(gene_diff, colnames(data)))
      #data <- as.data.frame(as.matrix(data))
      data <- rbind(data, pseudo_data)
      warning("missing features above will be set to 0\n")
    }
    features <- features[features %in% rownames(data)] # do not use intersect in case of duplicated features
  }

  # subset groups and set group_by as factor and drop group_by levels
  if(is.null(groups) && !is.factor(group_by)) group_by <- factor(group_by)
  # using SubsetDataAndGroup to set group_by as factor when groups is not null
  ret_tmp <- SubsetDataAndGroup(data, group_by, groups, droplevels = T)
  data <- ret_tmp$data; group_by <- ret_tmp$group_by; rm(ret_tmp)

  # do summarise
  meanData <- apply(data[features,,drop = F], 1, function(x){
    base::tapply(x, list(group_by = group_by),
                 SummariseExpr,
                 method = method,
                 by_log = by_log,
                 log_base = log_base,
                 exclude_zero = exclude_zero,
                 outlier_cutoff = outlier_cutoff)
  })

  percData <- apply(data[features,,drop = F], 1, function(x){
    tapply(x, list(group_by = group_by), function(y){
      sum(y > 0)/length(y)
    })
  })

  return(list(mean = meanData, percentage = percData))
}

#' Subset Data and Group
#'
#' @param data expression matrix of feature by cell. Sparse matrix will be densed, so be caution of your memory.
#' @param group_by a vector including cell annotations
#' @param groups subset these selected groups. set `NULL` to use all.
#' @param droplevels drop the unused levels in `group_by`
#'
#' @return a list of data and group_by
#' @export
#'
#' @examples
SubsetDataAndGroup <- function(
  data,
  group_by,
  groups = NULL,
  droplevels = T) {

  groups_ <- setdiff(groups, group_by)
  if(length(groups_) > 0) {
    warning(Message(head = "NOT-IN-group_by: ", items = groups_))
    groups <- setdiff(groups, groups_)
  }

  # data <- as.matrix(data)
  if(!is.null(groups)){
    data <- data[, group_by %in% groups, drop = F]
    group_by <- factor(group_by[group_by %in% groups], levels = groups)
  }

  if(is.factor(group_by) && droplevels){
    group_by <- droplevels(group_by)
  }

  return(list(data = data, group_by = group_by))
}


#' format a text message
#'
#' @param head head of message
#' @param tail tail of message. `\n`` should be used at the end always.
#' @param items join multiple `items` into a single string with `item_sep` as seperator.
#' @param item_sep item seperator
#'
#' @return a string
#' @export
#'
#' @examples
Message <- function(
  head = NULL,
  items = NULL,
  item_sep = ", ",
  tail = "\n"){
  stringr::str_glue(head, stringr::str_c(items, collapse = item_sep), tail)
}

#' check NULL or NA
#'
#' return TRUE if x is NULL or then NULL
#' @param x a value
#'
#' @return a logical value
#' @export
#'
#' @examples
IsNULLorNA <- function(x) {
  if(is.null(x)) {
    return(TRUE)
  }else{
    return(is.na(x)[[1]])
  }
  # is.na(NULL) reseut into logical(0)
}


#' Do summarise (e.g. mean) for a given vector
#'
#' @param x numeric values in log scale
#' @param method base::mean or any other function returned a value
#' @param by_log calculate by log scale ; If `TRUE` log_base will be ignored.
#' @param log_base `base::exp(1)` or `2` or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set `0.99` to exclude top 1% cells. (default: 1)
#'
#' @return a numeric value
#' @export
#'
#' @examples
SummariseExpr <- function(
  x, method = base::mean,
  by_log = F, log_base = base::exp(1),
  exclude_zero = F, outlier_cutoff = 1){
  if(exclude_zero) x <- x[x>0]
  if (length(x) == 0) return(0)
  x <- x[x<=stats::quantile(x, probs = outlier_cutoff)]
  if (length(x) == 0) return(0)
  if(by_log){
    method(x)
  }else{
    logb(method(log_base^x), log_base)
  }
}


#' Check pal in brewer
#'
#' @param pal a string of names of brewer.pal, or a vector of colors to be used to produce n gradient colors.
#' @param n length of returned color vector
#'
#' @return a vector of colors
#' @export
#'
#' @examples
CheckBrewerPal <- function(pal = "YlGn", n = 5){
  if(pal %in% rownames(RColorBrewer::brewer.pal.info)){
    colours <-  c(RColorBrewer::brewer.pal(name = pal, n = n))
  }else{
    colours <- grDevices::colorRampPalette(pal)(n)
  }
  return(colours)
}

#' a helper function to write files for metascape analysis.
#'
#' @param genes a vector of genes (e.g. Symbols)
#' @param group_by a vector of gene annotation (e.g. Cluster). should be equal lenght of genes.
#' @param txt_file text filenamme
#' @param top_n only write top_n genes for each group.
#'
#' @return
#' @export
#'
#' @examples
write_genesToMeatascape <- function(genes, group_by, txt_file, top_n = NULL){
  groups <- sort(unique(group_by))
  txtfile = file(description = txt_file, open = "wt")
  writeLines(text = "#Cluster\tGene", con = txtfile)
  for(i in groups){
    if(is.null(top_n)){
      txtLine <- paste0(i, "\t", paste0(genes[group_by == i], collapse = ","))
    }else if(top_n > 0){
      txtLine <- paste0(i, "\t", paste0(genes[group_by == i][1:round(top_n)], collapse = ","))
    }
    writeLines(text = txtLine, con = txtfile)
  }
  close(txtfile)
}

#' Infer the scale factor used for normalization
#'
#' @param x matrix or data.frame
#' @param log_base normalization log scaled base
#' @param pseudo_count usually 1.
#'
#' @return
#' @export
#'
#' @examples
InferScaleFactor <- function(x, log_base = base::exp(1), pseudo_count = 1){
  unique(round(as.matrix(colSums(log_base^x-pseudo_count))))
}

#' discrete colors imported from scanpy.pl.palettes
#'
#' A list containing the colors used in scanpy
#'
#' @format A list with 4 different palettes of 10, 20, 26, 64 colors.
#' @source \url{https://github.com/theislab/scanpy/issues/387},
#' \url{https://github.com/vega/vega/wiki/Scales#scale-range-literals},
#' \url{https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d},
#' \url{http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html}
"scanpy_colors"
