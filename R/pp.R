#' Scale data in each group
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.
#' @param group_by a vector including cell annotations
#' @param by_log calculate in log scale or in normalized scale; If `TRUE` log_base will be ignored.
#' @param log_base `base::exp(1)` or `2` or other else, depends on your normalization method.
#' @param center a logical value used in `base::scale``
#' @param scale a logical value used in `base::scale``
#'
#' @return scaled data
#' @export
#'
#' @examples
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
