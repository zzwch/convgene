#' @include utils.R
#'
NULL


#' Gene Set Scores
#'
#' @param object Seurat3 object
#' @param genesets a list of one or more sets of features
#' @param method supported methods. It turns the sparse matrix into dense matrix for methods in GSVA.
#' @param ... other params passed to corresponding method
#'
#' @return a list of vectors of score values
#' @export
#'
#' @examples
#'
tl_GeneSetScore <- function(
  object, genesets,
  method = c("pca", "moduleScore", "gsva", "ssgsea", "zscore", "plage"),
  ...){

  genesets <- sapply(genesets, function(features) {
    CheckFeatures(object, features, T, T)
  }, simplify = F)

  if(method == "pca"){
    return(
      sapply(genesets, function(features) {
        object_tmp <- Seurat::RunPCA(object, features = features, ... )
        features_means <- Matrix::colMeans(GetAssayData(object_tmp)[features,])
        return(apply(Embeddings(object_tmp), 2, SignCorrectionByCor, y = features_means))
      }, simplify = F)
    )
  }

  if(method == "moduleScore"){
    object_tmp <- Seurat::AddModuleScore(object, features = genesets, name = "ModuleScore", ... )
    return(object_tmp@meta.data[,grepl("ModuleScore", colnames(object_tmp@meta.data))])
  }

  if(method %in% c( "gsva", "ssgsea", "zscore", "plage")){
    return(GSVA::gsva(as.matrix(GetAssayData(object)), genesets, method = method, ...))
  }
}

#' a robust PCA for Seurat
#'
#' Based on a pre-computed dimensional reduction (typically calculated on a subset of genes and projects this onto the entire dataset (all genes).
#' Then generate a new rPCA reduction.
#'
#' @param object Seurat3 object
#' @param slot use which data to calc rpca
#' @param assay use which assay. NULL to use DefaultAssay(object)
#' @param reduction use which reduction to obtain top features.
#' @param do.score.scale scale robust scores
#' @param topn robust top n genes of loadings
#' @param reduction.name new reduction name
#' @param use_all_genes get top n genes from all genes loadings, not just typically highly variable genes. This will call ProjectDim if inexist.
#'
#' @return updated Seurat3 object with a new pca_score reduction added
#' @import Seurat
#' @export
#'
#' @examples
#'
tl_RunPCAScore <- function(
  object, slot = "data", assay = NULL,
  reduction = "pca",
  use_all_genes = T,
  topn = 30,
  do.score.scale = T,
  reduction.name = "rpca"){

  data <- GetAssayData(object, slot = slot, assay = assay)
  # SetDimReduction(object, "pca_score", "loadingData", object@dr$pca@loadingData)
  # SetDimReduction(object, "pca_score", "loadingData.full", object@dr$pca@loadingData.full)

  loadingData <- Loadings(epc, reduction = reduction, projected = use_all_genes)
  if(min(dim(loadingData)) == 0){
    if(use_all_genes){
      object <- ProjectDim(object, reduction = reduction, assay = assay, do.center = T)
      loadingData <- Loadings(object, reduction = reduction, projected = use_all_genes)
    }else{
      stop("Please RunPCA first. You may also run ProjectDim thereafter.")
    }
  }

  object@reductions[[reduction.name]] <- object@reductions[[reduction]]

  pcaEmbedData <- Embeddings(object, reduction)
  reduction.key <- Key(object)[[reduction]]
  rpcaEmbedData <- matrix(0, nrow = nrow(pcaEmbedData), ncol = 4*ncol(pcaEmbedData))
  rownames(rpcaEmbedData) <- rownames(pcaEmbedData)
  colnames(rpcaEmbedData) <- paste0(rep(paste0(reduction.key, 1:ncol(pcaEmbedData)), each = 4), c("", ".pos",".neg",".score"))
  rpcaEmbedData[,colnames(pcaEmbedData)] <- pcaEmbedData
  for(i in colnames(loadingData)){
    features.pos <- names(sort(loadingData[,i], decreasing = T)[1:topn])
    features.neg <- names(sort(loadingData[,i], decreasing = F)[1:topn])
    features.use <- c(features.pos, features.neg)
    data.scale = data[features.use, ] / apply(data[features.use, ], 1, max)

    score.data <- data.frame(neg = Matrix::colSums(data.scale[features.neg, ]),
                             pos = Matrix::colSums(data.scale[features.pos, ]))
    rownames(score.data) <- colnames(data.scale)
    if(do.score.scale){
      score.data$pos <- score.data$pos / max(score.data$pos)
      score.data$neg <- score.data$neg / max(score.data$neg)
    }
    rpcaEmbedData[,paste0(i,  c(".pos",".neg",".score"))] <- as.matrix(cbind(score.data, score.data$pos - score.data$neg))
  }

  object@reductions[[reduction.name]]@cell.embeddings <- rpcaEmbedData

  return(object)
}



#' clac DEGs among multiple groups
#'
#' @param data feature by cell matrix
#' @param group_by a vector of cell annotations
#' @param groups subset cells belonging to these groups
#' @param features features to calc
#' @param multi_test function to test across multiple groups
#' @param do_pairwise do pairwised test
#' @param pairs a list of pairs to test, NULL for all pairs
#' @param pairwise_test function to test two groups included in pairs.
#' @param ... ohter params passed to pairwise_test.
#' @param p_adjust adjust p value
#'
#' @return a data.frame of gene and pvalues
#' @export
#'
#' @examples
#'
tl_PairwiseDEGs <- function(
  data, group_by, groups = NULL,
  features = NULL, p_adjust = T,
  multi_test = kruskal.test,
  do_pairwise = F,
  pairs = NULL,
  pairwise_test = wilcox.test,
  ...
  ){

  if(IsNULLorNA(features)) features <- rownames(data)

  ret <- SubsetDataAndGroup(data[features,], group_by, groups)
  data <- ret$data; group_by <- ret$group_by

  message("running multi_test")
  multiDEGs <- data.frame(
    row.names = features,
    gene = features,
    multi_pval = apply(data, 1, function(x) {
      multi_test(x, group_by)$p.value
      }))

  if(p_adjust) multiDEGs[['multi_pval_adj']] <- stats::p.adjust(multiDEGs[['multi_pval']])

  if(do_pairwise){
    message("running pairwise_test")

    if(IsNULLorNA(groups)) groups <- sort(unique(group_by))
    groups <- as.character(groups)
    if(IsNULLorNA(pairs)) {
      pairs <- utils::combn(groups, 2, simplify = F)
    }

    for(pair in pairs){
      pair_prefix <- paste0(pair[1],"vs", pair[2])
      message("running pairwise_test ", pair_prefix)
      # pval
      multiDEGs[[paste0(pair_prefix, "_pval")]] <- apply(data, 1, function(x) {
        pairwise_test(x[group_by == pair[1]],
                      x[group_by == pair[2]], ...)$p.value
      })
      # Fold change
      multiDEGs[[paste0(pair_prefix, "_logFC")]] <- apply(data, 1, function(x) {
        ExpMean(x[group_by == pair[1]]) - ExpMean(x[group_by == pair[2]])
      })
      # adjust
      if(p_adjust) multiDEGs[[paste0(pair_prefix, "_pval_adj")]] <- stats::p.adjust(multiDEGs[[paste0(pair_prefix, "_pval")]])
    }
  }

  # reorder columns
  pairs_prefix <- sapply(pairs, function(x) {paste(x, collapse = "vs")})
  columns_order <- c("gene","multi_pval")
  if(do_pairwise) columns_order <- c(columns_order, paste0(pairs_prefix, "_pval"), paste0(pairs_prefix, "_logFC"))
  if(p_adjust) {
    if(do_pairwise){
      columns_order <- c(columns_order, paste0(pairs_prefix, "_pval_adj"))
    }
    columns_order <- c(columns_order, "multi_pval_adj")
  }

  return(multiDEGs[,columns_order])
}

#' calculate inter-similarity or intra-similarity of clusters
#'
#'
#' @param data feature by cell matrix
#' @param group_by vector of cell annotation
#' @param groups only consider subset cells belonging to these groups
#' @param dist_fun function to calculate distance/similarity. see Details.
#' @param metric metric passed to dist_fun
#' @param ... ohter params passed to dist_fun
#' @param transform_to_similarity from dist to similarity
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


  distData <- as.matrix(dist_fun(data, metric, ...))
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


  if(IsNULLorNA(subset_by)){
    simData <- as.data.frame(t(sapply(levels(group_by), function(x){
      c(getBothSimilarity(distData, group_by == x, group_by != x), group = x)
    })))
  }else{
    if(!is.factor(subset_by)){
      subset_by <- factor(subset_by, sort(unique(subset_by)))
    }
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
  sds <- slingshot::slingshot(data = Seurat::Embeddings(object, reduction = reduction),
                   clusterLabels = clusters, shrink = shrink,
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

