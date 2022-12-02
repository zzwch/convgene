#' @include utils.R
#'
NULL


#' Quick add module scores for markers from Atlas of Human Blood Cells
#'
#' @param object Seurat3 object
#' @param topn top n marker genes for each cluster
#' @param ... other parameters passed to `Seurat::AddModuleScore`
#'
#' @return Seurat Object
#' @export
#' @details
#' http://scrna.sklehabc.com/
#' Marker Gene: Information of 43 transcriptional cell clusters and signature genes with AUC score.
#' National Science Review, Volume 8, Issue 3, March 2021, nwaa180, https://doi.org/10.1093/nsr/nwaa180
#'
#' @examples
#'
tl_AddABCScore <- function(object, topn = 50, ...){
  abc_modules <- curated_markers$ABC_human %>%
    group_by(RNA_Cluster) %>%
    top_n(topn, -p_val_adj) %$%
    df2list(Gene, RNA_Cluster, drop = T)

  for (i in names(abc_modules)) {
    object <- Seurat::AddModuleScore(object, abc_modules[i], name = paste0("ABC_",i), ...)
  }

  return(object)
}

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

  method <- match.arg(method)
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
#' @param data.scale.max scale parameter used in calculation of RPCA score. Set 0.99 or 0.95 as you like to exclude outliers.
#' @param score.scale.max scale parameter used in calculation of RPCA score. Set 0.99 or 0.95 as you like to exclude outliers.
#' @param loadings if set, use this loadings to calculate RPCA score
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
  data.scale.max = 1,
  do.score.scale = T,
  score.scale.max = 1,
  reduction.name = "rpca",
  loadings = NULL){

  data <- GetAssayData(object, slot = slot, assay = assay)
  # SetDimReduction(object, "pca_score", "loadingData", object@dr$pca@loadingData)
  # SetDimReduction(object, "pca_score", "loadingData.full", object@dr$pca@loadingData.full)

  if(is.null(loadings)) {
    loadingData <- Loadings(object, reduction = reduction, projected = use_all_genes)
    if(min(dim(loadingData)) == 0){
      if(use_all_genes){
        object <- ProjectDim(object, reduction = reduction, assay = assay, do.center = T)
        loadingData <- Loadings(object, reduction = reduction, projected = use_all_genes)
      }else{
        stop("Please RunPCA first. You may also run ProjectDim thereafter.")
      }
    }
  }else{
    loadingData <- loadings
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
    data.scale = pmin(data[features.use, ] / apply(data[features.use, ], 1, quantile, probs = data.scale.max), 1)
    data.scale[is.na(data.scale)] <- 0

    score.data <- data.frame(neg = Matrix::colSums(data.scale[features.neg, ]),
                             pos = Matrix::colSums(data.scale[features.pos, ]))
    rownames(score.data) <- colnames(data.scale)
    if(do.score.scale){
      score.data$pos <- pmin(score.data$pos / quantile(score.data$pos, probs = score.scale.max), 1)
      score.data$neg <- pmin(score.data$neg / quantile(score.data$neg, probs = score.scale.max), 1)
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


#' calculate DEGs for two groups
#'
#' useful for volcano plot
#'
#' @param data expression matrix
#' @param group.by vector
#' @param group.1 character or vector
#' @param group.2 character or vector
#' @param FC_fun Seurat::ExpMean or base::mean. log2(e) is equal to 1/log(2)
#' @param stat_fun wilcox.test or t.test
#' @param pseudocount avoid NA in calculation of fold change
#'
#' @return data.frame
#' @export
#'
#' @examples
#'
tl_statGene <- function(
  data,
  group.by,
  group.1,
  group.2 = NULL,
  FC_fun = function(x, pseudocount) {Seurat::ExpMean(x + pseudocount)*log2(exp(1))},
  stat_fun = stats::wilcox.test,
  pseudocount = 0){

  ind_1 <- group.by %in% group.1
  ind_2 <- if(is.null(group.2)) !group.by %in% group.1 else group.by %in% group.2

  AVG <- apply(data, 1, function(expr) {c(FC_fun(expr[ind_1], pseudocount), FC_fun(expr[ind_2], pseudocount))}) %>% t %>% as.data.frame()
  FC <- AVG[,1] - AVG[,2]

  P <- apply(data, 1, function(expr) {stat_fun(expr[ind_1], expr[ind_2])$p.value})
  data.frame(
    gene = rownames(data) %>% as.character(),
    avg.1 = AVG[,1],
    avg.2 = AVG[,2],
    log2foldchange = FC,
    pvalue = P,
    p.adjust = p.adjust(P, "fdr")
  )
}
function(x1, x2) {list(avg.1 = Seurat::ExpMean(x1), avg.2 = Seurat::ExpMean(x2), log2foldchange = log2(exp(1))*(Seurat::ExpMean(x1) - Seurat::ExpMean(x2)))}
#' Prepare genelist for GSEA
#'
#' @param data expression matrix
#' @param group.by vector
#' @param group.1 character or vector
#' @param group.2 character or vector
#' @param order.by order method
#' @param stat_fun stat function, used for 'stat_statistic'. Now only t.test is valid.
#'
#' @return geneList
#' @export
#'
#' @examples
#'
tl_OrderGene <- function(
  data,
  group.by,
  group.1,
  group.2 = NULL,
  order.by = c("FC_ExpMean", "FC_mean", "stat_statistic"),
  stat_fun = t.test){

  ind_1 <- group.by %in% group.1
  ind_2 <- if(is.null(group.2)) !group.by %in% group.1 else group.by %in% group.2
  if (order.by == "FC_ExpMean")  ret <- apply(data, 1, function(expr) {Seurat::ExpMean(expr[ind_1]) - Seurat::ExpMean(expr[ind_2])})
  if (order.by == "FC_mean")  ret <- apply(data, 1, function(expr) {mean(expr[ind_1]) - mean(expr[ind_2])})
  if (order.by == "stat_statistic")  ret <- apply(data, 1, function(expr) {stat_fun(expr[ind_1], expr[ind_2])$statistic})

  return(ret %>% sort(decreasing = T))
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
#' @param gene_force Force gene using zeros, if it is not in data.
#' @param cap_value the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)
#' @param rescale_max rescale data
#' @param rescale_min rescale data
#' @param rescale_again do rescale again after summarisie
#'
#' @return a vector of scores
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


#' Plot proportions
#'
#' @param tab table
#' @param proportion_by row or column
#' @param ret_value return proportions or heatmap
#' @inheritParams pheatmap
#' @param ... other params used for pheatmap
#'
#' @return proportions or heatmap
#' @export
#'
#' @examples
#'
tl_crossTable <- function(
  tab, proportion_by = c("column", "row"), ret_value = F,
  color = colorRampPalette(c("white", "orange", "red"))(20),
  cluster_rows = F, cluster_cols = F, display_numbers = T,
  number_format = "%1.1f%%", ...
){
  proportion_by <- match.arg(proportion_by)
  if(proportion_by == "row")
    pmat <- apply(tab, 2, function(x) x/rowSums(tab)) %>% t
  else
    pmat <- apply(tab, 1, function(x) x/colSums(tab))

  pmat %<>% t
  if(ret_value) return(pmat)
  pheatmap::pheatmap(
    pmat * 100,
    color = color, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
    display_numbers = display_numbers, number_format = number_format, ...)
}

#' statistics of cross table
#'
#' @param tab a table object
#' @param p.adjust adjust pvalue
#' @param proportion_by row or column
#' @param pal_1
#' @param pal_2
#' @param pal_3
#' @param pal_4
#' @param heatmap plot heatmap
#'
#' @return a matrix of pvalues
#' @export
#'
#' @examples
#'
tl_crossTableEnrichment <- function(
  tab, p.adjust = T,heatmap = T, proportion_by = c("column", "row"),
  pal_1 = "Purples",
  pal_2 = "Purples",
  pal_3 = "YlOrRd",
  pal_4 = "YlGn"){

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
    proportion_by <- match.arg(proportion_by)
    tab.h <- tab
    tab.h[tab == 0] <- NA
    pheatmap::pheatmap(tab.h, cluster_cols = F, cluster_rows = F, display_numbers = tab, main = "count",
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, pal_1))(11))
    pheatmap::pheatmap(tab.h, cluster_cols = F, cluster_rows = F, main = stringr::str_glue("count, {proportion_by} scaled"),
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, pal_2))(11), scale = proportion_by)

    tab.prop <- tl_crossTable(tab, proportion_by = proportion_by, ret_value = T)

    tab.ph <- -log10(abs(tab.p)) * sign(tab.p)
    tab.ph[tab.ph < -log10(0.05)] <- NA
    pheatmap::pheatmap(tab.ph, cluster_cols = F, cluster_rows = F,
                       display_numbers = formattable::percent(tab.prop),
                       main = "p value, enrichment",
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, pal_3))(11))
    tab.ph <- log10(abs(tab.p)) * sign(tab.p)
    tab.ph[tab.ph < -log10(0.05)] <- NA
    pheatmap::pheatmap(tab.ph, cluster_cols = F, cluster_rows = F,
                       display_numbers = formattable::percent(tab.prop),
                       main = "p value, anti-enrichment",
                       color = colorRampPalette(RColorBrewer::brewer.pal(9, pal_4))(11))
  }
  return(tab.p)
}


#' A wrapper function for clusterprofiler enrichment analysis (compareCluster)
#'
#' @param degs_table deg results from FindAllMarkers
#' @param enrichFun enrich function name in clusterprofiler
#' @param OrgDb "org.Hs.eg.db" or "org.Mm.eg.db according to your organism
#' @param topn
#' @param ont One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#'
#' @return a list of go and simplified go results
#' @export
#'
#' @examples
#'
tl_goCompareCluster <- function(degs_table, topn = 200, enrichFun = "enrichGO",
                               OrgDb = "org.Hs.eg.db", ont = "BP"){
  # clusterprofiler
  message("Enrichment analysis using clusterProfiler...")
  cluster_go <- clusterProfiler::compareCluster(
    data = degs_table %>% group_by(cluster) %>% top_n(topn, -p_val) %>% top_n(topn, avg_logFC),
    geneClusters = gene~cluster, fun = enrichFun,
    OrgDb = OrgDb, keyType = "SYMBOL", ont = ont)
  #message("Simplify enrichment results...")
  #cluster_go_sim <- simplify(cluster_go)

  return(cluster_go)
}


#' A wrapper function for clusterprofiler GO enrichment analysis (enrichGO)
#'
#' @param gene a vector of focused genes
#' @param simplify do simplify or not
#' @param OrgDb org.Hs.eg.db or org.Mm.eg.db, or else
#' @param keyType SYMBOL or ENTREZID, or else
#' @param ont BP, MF, CC
#' @param ... other parameters for enrichGO
#'
#' @return a list of GO and simplified GO results
#' @export
#'
#' @examples
#'
tl_go <- function(gene, simplify = F, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", ...){
  # clusterprofiler
  message("Enrichment analysis using clusterProfiler...")
  cluster_go <- clusterProfiler::enrichGO(gene, OrgDb, keyType, ont, ...)

  if(simplify){
    message("Simplify enrichment results...")
    cluster_go_sim <- clusterProfiler::simplify(cluster_go)
    return(list(GO = cluster_go, GOSIM = cluster_go_sim))
  }else{
    return(cluster_go)
  }
}


#' Using class::knn method in Dimension Reduction space to classify cells with referenced clusters
#'
#' @param object seurat3 object
#' @param reference_group_by a meta.data colname
#' @param reduction dimension reduction
#' @param dims use all dims if NULL, or specify it
#' @param k number of kNN used to predict cluster
#' @param reference_cells cells used as reference
#' @param predict_cells cells to be predicted/classified
#'
#' @return seurat object
#' @importFrom class knn
#' @export
#'
#' @examples
#'
tl_knnClassify <- function(object, reference_group_by, predict_group_by = "knnClassify",
                           reduction = "mnn", dims = NULL, k = 20,
                           reference_cells = NULL, predict_cells = NULL){
  DimSpace <- Embeddings(object, reduction = reduction)
  if(!is.null(dims)) DimSpace <- DimSpace[, dims]

  cell_groups <- object[[reference_group_by]]
  if(is.null(reference_cells)) reference_cells <- rownames(na.omit(cell_groups))
  if(is.null(predict_cells)) predict_cells <- names(na.action(na.omit(cell_groups)))
  res_knn <- class::knn(DimSpace[reference_cells,],
                        DimSpace[predict_cells,],
                        cell_groups[reference_cells, 1],
                        k = k, prob = T)
  cell_groups[predict_cells, 1] <- as.character(res_knn)
  object[[predict_group_by]] <- cell_groups

  return(object)
}


#' Rotate or Flip reduction
#'
#' @param object Seurat object
#' @param reduction reduction. eg. "umap"
#' @param flip_x do x flipping
#' @param flip_y do y flipping
#' @param rotation angle in degree to rotate
#'
#' @return updated Seurat
#' @export
#'
#' @examples
#'
tl_transDim <- function(object, reduction = "umap", flip_x = F, flip_y = F, rotation = 0){
  dr <- Embeddings(object, reduction = reduction)
  if(flip_x) dr[,1] <- -dr[,1]
  if(flip_y) dr[,2] <- -dr[,2]
  rotation <- pi*rotation/180
  trans_matrix <- matrix(c(cos(rotation), sin(rotation), 0,
                           -sin(rotation), cos(rotation), 0,
                           0, 0, 1),
                         nrow = 3, byrow = TRUE)
  object@reductions[[reduction]]@cell.embeddings[,c(1,2)] <- t(trans_matrix %*% t(as.matrix(cbind(dr, 1))))[,c(1,2)]
  return(object)
}


#' Add new Reduction to Seurat object
#'
#' @param object Seurat v3 object
#' @param reduction reduction name
#' @param dim.1 vector of 1st dimension
#' @param dim.2 vector of 2nd dimension
#' @param ... other params passed to CreateDimReducObject()
#'
#' @return updated Seurat object
#' @export
#'
#' @examples
#'
tl_AddReduction <- function(object, reduction, dim.1, dim.2, ...){
  embeddings <- as.matrix(cbind(dim.1, dim.2))
  colnames(embeddings) <- paste(toupper(reduction), 1:2, sep = "_")
  object@reductions[[tolower(reduction)]] <- CreateDimReducObject(embeddings = embeddings, key = paste0(toupper(reduction), "_"), ...)
  return(object)
}

#' Reassign cell cycle pahse and set cc reduction
#'
#' @param object Seurat v3 object
#' @param phase.name name added to meta.data
#' @param s.th,g2m.th,slope threshold of S.Score and G2M.Score used to reassign phase, and slope is used for partition S and G2/M
#' @param set.reduction set NULL if do not want to add cc reduction, Or set a character value for reduction name
#' @param s.scores vector of S.Scores
#' @param g2m.scores vector of G2M.Scores
#'
#' @return updated Seurat object
#' @export
#'
#' @examples
#'
tl_CellCycleAssign <- function(object, phase.name = "Phase4", s.th = 0, g2m.th = 0, set.reduction = "cc", slope = 1, s.scores = "S.Score", g2m.scores = "G2M.Score"){
  s.scores <- object[[s.scores]]
  g2m.scores <- object[[g2m.scores]]

  object[[phase.name]] <- case_when(
    g2m.scores < g2m.th & s.scores < s.th ~ "Q/G0",
    g2m.scores < g2m.th & s.scores >= s.th ~ "G1",
    g2m.scores > g2m.th & (g2m.scores - g2m.th)/(s.scores - s.th) < slope & s.scores >= s.th  ~ "S",
    TRUE ~ "G2/M"
  ) %>% factor(levels = c("Q/G0", "G1", "S", "G2/M"))

  if(!is.null(set.reduction)){
    object <- tl_AddReduction(object = object, reduction = "CC", dim.1 = s.scores, dim.2 = g2m.scores)
  }

  return(object)
}


#' calc Regulon Specificity Score (RSS) or named TFSS
#'
#' @param aucData auc scores matrix derived from SCENIC (e.g., 3.4_regulonAUC.Rds)
#' @param binData binary matrix
#' @param group_by a vector of clusters
#' @param groups clusters to be used
#' @param sig.auc.min,sig.bin.min,sig.rssz.min Regulon with AUCell score > sig.auc.min, Binary Activaty > sig.bin.min and RSSZ > sig.rssz.min were considered to be significant.
#' @param pseudocount used as numerator and denominator constant for calculation of fold change
#'
#' @return list of 1) matrix of group aggregated AUC scores, 2) matrix of Z-score normalized RSS, 3) data.frame of significant regulon
#'
#' @export
#'
#' @examples
#'
tl_regulon_specificity <- function(aucData, binData, group_by, groups = NULL, sig.auc.min = 0, sig.bin.min = 0.1, sig.rssz.min = 1, pseudocount = 0.01){
  if(is.null(groups)) {
    groups <- sort(unique(group_by))
  } else {
    tmp <- SubsetDataAndGroup(aucData, group_by, groups)
    aucData <- tmp$data
    group_by <- tmp$group_by
    groups <- sort(unique(group_by))

    tmp <- SubsetDataAndGroup(binData, group_by, groups)
    binData <- tmp$data

  }

  if(!all(colnames(aucData) == colnames(binData))) stop("data colnames is not identical!")
  group_mat <- lapply(groups, function(i) {as.numeric(group_by == i)}) %>% do.call(what = rbind)
  # apply(aucData, 1, function(i) {
  #   apply(group_mat, 1, function(j) {
  #     philentropy::JSD(rbind(i, j), est.prob = "empirical")
  #   })
  # })
  jsd <- philentropy::JSD(rbind(aucData, group_mat), est.prob = "empirical")[1:nrow(aucData), nrow(aucData)+(1:nrow(group_mat))]
  rss <- 1-sqrt(jsd)
  rownames(rss) <- rownames(aucData)
  colnames(rss) <- groups

  rssz <- t(scale(t(rss)))

  aucg <- apply(aucData, 1, function(x) { tapply(x, group_by, mean)}) %>% t
  bing <- apply(binData, 1, function(x) { tapply(x, group_by, mean)}) %>% t

  myMerge <- function(x,y){merge(x,y)}

  sig <- Reduce(f = myMerge,
                list(
                  aucg %>% as.data.frame() %>% rownames_to_column() %>%
                    pivot_longer(!all_of("rowname"), values_to = "gAUC"),
                  bing %>% as.data.frame() %>% rownames_to_column() %>%
                    pivot_longer(!all_of("rowname"), values_to = "gBinAct"),
                  rss %>% as.data.frame() %>% rownames_to_column() %>%
                    pivot_longer(!all_of("rowname"), values_to = "RSS"),
                  rssz %>% as.data.frame() %>% rownames_to_column() %>%
                    pivot_longer(!all_of("rowname"), values_to = "RSSZ")
                  )
                ) %>%
    mutate(significant = (gAUC > sig.auc.min) & (gBinAct > sig.bin.min) & (RSSZ > sig.rssz.min),
           power = gBinAct*RSS,
           powerz = gBinAct*RSSZ,
           name = factor(name, levels = groups)) %>%
    arrange(name, desc(RSS))
  colnames(sig)[1:2] <- c("regulon", "cluster")

  # differentially active regulons: DAR
  dar <- lapply(groups, function(x) {
    statData <- convgene::tl_statGene(aucData, group.by = group_by, group.1 = x, FC_fun = function(x, p) {log2(pseudocount + mean(x))}, pseudocount = pseudocount)
    statData$cluster = x
    return(statData)
  }) %>% dplyr::bind_rows()

  return(list(
    single = list(AUC = aucData, BIN = binData),
    group = list(AUC = aucg, BIN = bing),
    longer = sig,
    dar = dar)
    )
}


#' categorize gene as expressed or not based on detection rate and non-zero median expression level
#'
#' @param data feature by cell matrix
#' @param group_by cell annotation vector
#'
#' @return data.frame of detection, nzmedian, expressed, cluster
#' @export
#'
#' @examples
#'
tl_geneExpressed <- function(data, group_by = NULL){
  if(IsNULLorNA(group_by)) group_by <- rep("All", ncol(data))

  lapply(sort(unique(group_by)), function(i) {
    det_i <- Matrix::rowSums(data >0)/ncol(data)
    det_i_th <- median(det_i)

    nzmedian_i <- apply(data, 1, function(x) median(x[x>0]))
    nzmedian_i_th <- data %>% median(.[.>0])

    gene_i <- (det_i > det_i_th) & (nzmedian_i > nzmedian_i_th)

    return(data.frame(detection = det_i,
                      nzmedian = nzmedian_i,
                      expressed = gene_i,
                      cluster = i) %>%
             rownames_to_column(var = "gene")
    )
  }) %>% do.call(what = rbind)
}

