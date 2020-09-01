#' Scale data in each group
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.
#' @param group_by a vector including cell annotations
#' @param by_log calculate in log scale or in normalized scale; If TRUE log_base will be ignored.
#' @param log_base base::exp(1) or 2 or other else, depends on your normalization method.
#' @param center a logical value used in base::scale
#' @param scale a logical value used in base::scale
#'
#' @return scaled data
#' @export
#'
#' @examples
#'
pp_centerData <- function(
  data,
  group_by = NULL,
  by_log = F,
  log_base = base::exp(1),
  center = T,
  scale = T){

  # data <- as.matrix(data)
  data <- if(!by_log) log_base^data

  if(IsNULLorNA(group_by)) group_by <- rep("A", ncol(data))

  if(is.factor(group_by)){
    group_by <- droplevels(group_by)
  }else{
    group_by <- factor(group_by)
  }

  for(i in levels(group_by)){
    data[,group_by == i] <- t(base::scale(t(data[,group_by == i]), center = center, scale = scale))
  }

  return(data)
}


#' A Wrapper of harmony::RunHarmony
#'
#' RunHarmony did not support dims.use actually, so I wrote this wrapper function to support that.
#'
#' @param object Seurat object
#' @param dims.use Which PCA dimensions to use for Harmony. By default, use all
#' @param group.by.vars Which variable(s) to remove (character vector).
#' @param ... other parameters to be used in harmony::HarmonyMatrix
#'
#' @return Seurat (version 3) object.
#' @export
#'
#' @examples
#'
pp_runHarmony <- function(object, dims.use, group.by.vars,...){
  embedding <- Seurat::Embeddings(object)[,dims.use]
  harmonyEmbed <- harmony::HarmonyMatrix(embedding, object@meta.data, group.by.vars,
                                FALSE, 0, ...)
  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0("harmony_", seq_len(ncol(harmonyEmbed)))
  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(embeddings = harmonyEmbed,
                                                stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)), key = "harmony")
  })
  object[["harmony"]] <- harmonydata
  if (TRUE) {
    object <- Seurat::ProjectDim(object, reduction = "harmony",
                                 overwrite = TRUE, verbose = FALSE)
  }
  return(object)
}

#' Find variable features across batches
#'
#' Identifies Highly variable features (HVFs) for each batch, then combine them into a single vector of HVFs by count.
#'
#' @param object Seurat object
#' @param split_by a colname of meta.data
#' @param groups only some groups of split_by are considered.
#' @param nfeatures Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'
#' @param ... Arguments passed to Seurat::FindVariableFeatures
#'
#' @return A vector of variable features
#' @export
#'
#' @examples
#'
pp_findVariableFeatures <- function(object, split_by = "ident", groups = NULL, nfeatures = 2000, assay = "RNA", log_var = F, ...){
  object_list <- Seurat::SplitObject(object, split_by)
  if(!is.null(groups)) object_list <- object_list[groups]
  hvgs <- data.frame(row.names = rownames(object))
  for (i in seq_along(object_list)) {
    object_list[[i]] <- Seurat::FindVariableFeatures(object_list[[i]], nfeatures = nfeatures, assay = assay,...)
    var_data <- object_list[[i]]@assays[[assay]]@meta.features
    var_ind <- grep("variance", colnames(var_data), value = T, ignore.case = T)
    hvgs <- cbind(hvgs, var_data[,var_ind[length(var_ind)]])
  }
  if(log_var){
    names(sort(rowMeans(log(hvgs)), decreasing = T)[1:nfeatures])
  }else{
    names(sort(rowMeans(hvgs), decreasing = T)[1:nfeatures])
  }
}
