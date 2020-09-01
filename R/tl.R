#' @include utils.R
#'
NULL


#' calculate inter-similarity or intra-similarity of clusters
#'
#'
#' @param data feature by cell matrix
#' @param group_by vector of cell annotation
#' @param groups only consider subset cells belonging to these groups
#' @param dist_fun function to calculate distance/similarity. see Details.
#' @param metric metric passed to dist_fun
#' @param ... ohter params passed to dist_fun
#' @param transform_to_similarity
#' @param transform_method dist to similarity function, see Details
#' @param subset_by if setted, calc intra-/inter-similarity of groups restricted to each subset
#' @param subsets only consider cells belonging to these subsets
#'
#' @return a slimilarity matrix, usually 0-1
#' @details
#' dist_fun:
#' stats::dist,
#' ClassDiscovery::distanceMatrix (pearson, sqrt pearson, spearman, absolute pearson and so on),
#' philentropy::distance (46 different distances/similarities),
#'
#' transform_method:
#' 1) 1/(1+x)
#' 2) 1-x/max(x)
#' 3) sqrt(1-x/max(x))
#' 4) 1-x
#' 5) 1-abs(x)
#' any meaningful function
#'
#' @source Fig. 4A at https://doi.org/10.1016/j.cell.2019.11.010
#' @export
#'
#' @examples
#'
tl_calcSimilarity <- function(
  data, group_by, groups = NULL,
  subset_by = NULL, subsets = NULL,
  dist_fun = stats::dist, metric = "euclidean",
  transform_to_similarity = T,
  transform_method = function(x) 1(1+x),
  ...){


  if(!IsNULLorNA(subset_by)){
    colnames(data) <- paste0("col", 1:ncol(data)) # keep colnames
    ret1 <- SubsetDataAndGroup(data, group_by, groups, droplevels = T)
    ret2 <- SubsetDataAndGroup(data, subset_by, subsets, droplevels = T)
    ret1_col <- colnames(ret1$data)
    ret2_col <- colnames(ret2$data)
    ret_col <- intersect(ret1_col, ret2_col)
    data <- ret1$data[,ret_col]
    group_by <- droplevels(ret1$group_by[match(ret_col, ret1_col)])
    subset_by <- droplevels(ret2$group_by[match(ret_col, ret2_col)])
  }else{
    ret1 <- SubsetDataAndGroup(data, group_by, groups, droplevels = T)
    data <- ret1$data
    group_by <- ret1$group_by
  }


  distData <- dist_fun(data, metric, ...)
  if(transform_to_similarity){
    distData <- transform_method(distData)
  }

  getBothSimilarity <- function(data, ind_intra, ind_inter){
    return(c(intra=(sum(data[ind_intra, ind_intra]) -	sum(ind_intra))/(sum(ind_intra)*sum(ind_intra) - sum(ind_intra)),
             inter=sum(data[ind_intra, ind_inter]) / (sum(ind_intra) *	sum(ind_inter))))
  }

  if(!is.factor(group_by)){
    group_by <- factor(group_by, sort(unique(group_by)))
  }
  if(!is.factor(subset_by)){
    subset_by <- factor(subset_by, sort(unique(subset_by)))
  }

  if(IsNULLorNA(subset_by)){
    as.data.frame(t(sapply(levels(group_by), function(x){
      c(getBothSimilarity(distData, group_by == x, group_by != x), group = x)
    })))
  }else{
    simData <- as.data.frame(t(
      do.call(cbind, lapply(
        levels(group_by),
        function(y){
          sapply(levels(subset_by), function(x){
            c(getBothSimilarity(pcaMat, group_by == y & subset_by == x, group_by == y & subset_by != x),
              group = y, subset = x)
          })
        }))
    ))
  }
  return(simData)
}

#' a wrapper of slingshot
#'
#' @param object Seurat3 object
#' @param reduction character, umap or others
#'
#' @return Seurat3 object updated with ss_lineages, ss_Pseudotime and reductions$ss
#' @import slingshot
#' @export
#'
#' @examples
#'
tl_RunSlingshot <- function(
  object, reduction = "umap", group_by,
  start.clus, shrink = 0.75, ...){

  clusters <- object@meta.data[,group_by,drop = T]
  sds <- slingshot::slingshot(data = Seurat::Embeddings(object, reduction = "umap"),
                   clusterLabels = clusters, shrink = 0.75,
                   start.clus = start.clus #, end.clus = c("EP5", "EP7")
                   )
  sds_lineages <- slingshot::slingLineages(sds)
  sds_time <- slingshot::slingPseudotime(sds, T)
  sds_curves <- slingshot::slingCurves(sds)

  object$ss_pseudotime <- NA

  for(i in seq_along(sds_lineages)){
    i_ind <- clusters %in% sds_lineages[[i]]

    i_lineage <- paste0("ss_lineage", i)
    object@meta.data[, i_lineage] <- NA
    object@meta.data[i_ind, i_lineage] <- paste0("Lineage",i)

    i_pseudo <- paste0("ss_pseudotime", i)
    object@meta.data[, i_pseudo] <- NA
    object@meta.data[i_ind, i_pseudo] <- sds_time[i_ind, i]

    object$ss_pseudotime[i_ind] <- rowMeans(cbind(sds_time[i_ind, i], object$ss_pseudotime[i_ind]), na.rm = T)

    i_curve <- paste0("ss", i)
    ss_curves <- as.matrix(sds_curves[[i]]$s) + NA
    ss_key <- paste0(toupper(i_curve),"_")
    colnames(ss_curves) <- paste0(ss_key, 1:ncol(ss_curves))

    ss_curves[i_ind,] <- sds_curves[[i]]$s[i_ind,]
    object@reductions[[i_curve]] <- CreateDimReducObject(embeddings = ss_curves, key = ss_key)
  }

  return(object)
}

#' calculate average/median rescale expression of features (eg. genes) for each cell
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled. A matrix is better than data.frame in case of duplicated features.
#' @param features features to be calculated
#' @param method using base::mean or stats::median or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If TRUE log_base will be ignored.
#' @param log_base base::exp(1) or 2 or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set 0.99 to exclude top 1 percent cells. (default: 1)
#' @param gene_force Force gene suing zeros, which is missing in data.
#' @param cap_value the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)
#' @param rescale_max rescale data
#' @param rescale_min rescale data
#' @param rescale_again do rescale again after summarisie
#'
#' @return a list of mean and percentage data
#' @export
#'
#' @examples
#'
tl_calcFeatureScore <- function(
  data,
  features,
  rescale_max = 10, rescale_min = 0,
  method = base::mean,
  by_log = T,
  log_base = base::exp(1),
  exclude_zero = F,
  outlier_cutoff = 1,
  cap_value = NULL,
  gene_force = F,
  rescale_again = T){

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

  # do summarise
  data <- t(apply(data[features, ,drop = F], 1,
                RobustRescale,
                probs = c(0, outlier_cutoff),
                max = rescale_max, min = rescale_min, ground = 0))

  meanData <- apply(data[features,,drop = F], 2,
                    SummariseExpr,
                    method = method,
                    by_log = by_log,
                    log_base = log_base,
                    exclude_zero = exclude_zero,
                    outlier_cutoff = outlier_cutoff)
  if(rescale_again) meanData <- RobustRescale(meanData,
                                              probs = c(0, outlier_cutoff),
                                              max = rescale_max,
                                              min = rescale_min,
                                              ground = 0)

  return(meanData)
}


#' statistics of cross table
#'
#' @param tab a table object
#' @param p.adjust adjust pvalue
#' @param heatmap plot heatmap
#'
#' @return a matrix of pvalues
#' @export
#'
#' @examples
#'
tl_crossTableEnrichment <- function(tab, p.adjust = T,heatmap = T){
  tab.p <- as.matrix(tab + NA)
  tab.d <- dim(tab)
  tab.rs <- as.numeric(rowSums(tab))
  tab.cs <- as.numeric(colSums(tab))
  tab.ss <- sum(tab)
  for(i in 1:tab.d[1]){
    for(j in 1:tab.d[2]){
      tab.c <- tab[i,j]
      tab.res <- fisher.test(matrix(c(tab.c,
                                      tab.rs[i] - tab.c,
                                      tab.cs[j] - tab.c,
                                      tab.ss - tab.rs[i] - tab.cs[j] + tab.c), nrow = 2))
      tab.p[i,j] <- ifelse(tab.res$estimate > tab.res$null.value, tab.res$p.value, -tab.res$p.value)
    }
  }
  tab.p <- as.numeric(tab.p)
  if(p.adjust){
    tab.p <- p.adjust(abs(tab.p), method = "BH") * sign(tab.p)
  }
  tab.p <- matrix(tab.p, nrow = tab.d[1], ncol = tab.d[2])
  dimnames(tab.p) <- dimnames(tab)

  if(heatmap){
    tab.h <- tab
    tab.h[tab == 0] <- NA
    pheatmap::pheatmap(tab.h, cluster_cols = F, cluster_rows = F, display_numbers = tab,
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))(11))
    pheatmap::pheatmap(tab.h, cluster_cols = F, cluster_rows = F,
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))(11), scale = "row")
    tab.ph <- -log10(abs(tab.p)) * sign(tab.p)
    tab.ph[tab.ph < -log10(0.05)] <- NA
    pheatmap::pheatmap(tab.ph, cluster_cols = F, cluster_rows = F,
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(11))
    tab.ph <- log10(abs(tab.p)) * sign(tab.p)
    tab.ph[tab.ph < -log10(0.05)] <- NA
    pheatmap::pheatmap(tab.ph, cluster_cols = F, cluster_rows = F,
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlGn"))(11))
  }
  return(tab.p)
}

