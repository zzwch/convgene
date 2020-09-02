
#' convert genes symbols between organisms
#'
#' using getLDS function of biomaRt to convert human genes to mouse genes or vice versa.
#' now only symbol is supported, which is used more frequently.
#'
#' @param genes genes to be converted
#' @param orgA from which organism, now human, mouse, macaque are supported
#' @param orgB to which organism, now human, mouse, macaque are supported
#' @param verbose report information
#'
#' @return data.frame of genes' mapping
#' @import biomaRt
#' @export
#'
#' @examples
#'
#' @source
#' # genes, homology, orthology
#' http://www.informatics.jax.org/downloads/reports/index.html
#' https://gist.github.com/FloWuenne/f8fc922477df04c1642e9d8945c48d47
#' https://github.com/saeyslab/nichenetr/blob/master/R/supporting_functions.R
#' https://github.com/saeyslab/nichenetr/blob/master/R-raw/developing_package_light.R
#'
ConvertHomoGene <- function(
  genes,
  orgA = "human",
  orgB = "mouse",
  verbose = F
  ){

  orgDict <- list(
    human = c(dataset = "hsapiens_gene_ensembl",
              attr_prefix = "hgnc_"),
    mouse = c(dataset = "mmusculus_gene_ensembl",
              attr_prefix = "mgi_"),
    Macaque = c(dataset = "mmulatta_gene_ensembl",
                attr_prefix = "hgnc_")
  )
  orgA <- orgDict[[orgA]]
  orgB <- orgDict[[orgB]]
  martA = useMart("ensembl", dataset = orgA["dataset"])
  martB = useMart("ensembl", dataset = orgB["dataset"])
  attrA <- paste0(orgA["attr_prefix"], "symbol")
  attrB <- paste0(orgB["attr_prefix"], "symbol")

  geneMap = getLDS(attributes = attrA, filters = attrA,
                   values = genes, mart = martA,
                   attributesL = attrB, martL = martB,
                   uniqueRows=T)
  geneMap <- unique(geneMap)
  if(verbose){
    if(length(unique(geneMap[,1])) == length(unique(genes))){
      message("All genes were translated successfully!")
    }else{
      warning(str_c(length(setdiff(genes, geneMap[,1])),
                    "genes could not be converted. You can use setdiff function to see them."))
    }
  }
  return(geneMap)
}

#' similar to scales::rescale but use quantile 95
#'
#' @param x a vector of numeric
#' @param probs get the actual_min and actual_max. this may remove the outliers effect.
#' @param max rescale actual_max to this max value
#' @param min rescale actual_min to this min value
#' @param ground for gene expression , to set 0 as the actual_min
#'
#' @return a vector
#' @export
#'
#' @examples
#'
RobustRescale <- function(x, probs = c(0.05, 0.95), max = 10, min = 0, ground = NULL){
  prob_values <- quantile(x, probs = probs)
  if(!IsNULLorNA(ground)){
    prob_values[1] <- ground
  }
  x[x > prob_values[2]] <- prob_values[2]
  x[x < prob_values[1]] <- prob_values[1]
  max*(x - min)/(prob_values[2] - prob_values[1])
}


#' update data.frame with new name added
#'
#' @param df data.frame
#' @param x NULL or a value of length one
#' @param default use this value is x is NULL
#' @param name the name to be added to df
#'
#' @return a new df with name added
#' @export
#'
#' @examples
#'
SetDefaultAes <- function(df, x = NULL, default = 1, name = "shape"){
  if(IsNULLorNA(x)){
    df[, name] <- default
  }else{
    df[, name] <- rep(x[[1]], nrow(df))
  }
  return(df)
}


#' calculate average/median and percentage expression of features (eg. genes) in groups
#'
#' @param data expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled. A matrix is better than data.frame in case of duplicated features.
#' @param features features to be visualized
#' @param group_by a vector including cell annotations
#' @param groups only these selected groups to be visualized
#' @param method using base::mean or stats::median or your custom function for each feature vector;
#' @param by_log calculate by log scale ; If TRUE log_base will be ignored.
#' @param log_base base::exp(1) or 2 or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set 0.99 to exclude top 1 percent cells. (default: 1)
#' @param gene_force Force gene suing zeros, which is missing in data.
#' @param cap_value the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)
#' @param do_scale scale the summarised value across groups. Useful to demonstrate the differences among groups.
#'
#' @return a list of summarised and percentage data
#' @export
#'
#' @examples
#'
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
  do_scale = F,
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
  }) # apply always return column for input either each row/column

  if(do_scale) meanData <- apply(meanData, 2, scales::rescale)
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
#' @param groups subset these selected groups. set NULL to use all.
#' @param droplevels drop the unused levels in group_by
#'
#' @return a list of data and group_by (factor)
#' @export
#'
#' @examples
#'
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
#' @param tail tail of message. \\n should be used at the end always.
#' @param items join multiple items into a single string with item_sep as seperator.
#' @param item_sep item seperator
#'
#' @return a string
#' @export
#'
#' @examples
#'
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
#'
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
#' @param by_log calculate by log scale ; If TRUE log_base will be ignored.
#' @param log_base base::exp(1) or 2 or other else, depends on your normalization method.
#' @param exclude_zero exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.
#' @param outlier_cutoff sometimes outliers (several extremely high cells) should be excluded when do summarise. Set 0.99 to exclude top 1 percent cells. (default: 1)
#'
#' @return a numeric value
#' @export
#'
#' @examples
#'
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
#'
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
#'
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
#'
InferScaleFactor <- function(x, log_base = base::exp(1), pseudo_count = 1){
  unique(round(as.matrix(colSums(log_base^x-pseudo_count))))
}

### str_to_gene()

#' convert string vector to list/vector of genes
#'
#' @param x a string with multiple gene symbols seperated by non-numeric/letter/hyphen (Depends on REGEX parameter below); numeric/letter/hyphen charactor will be treated as parts of gene symbols. A vector of strings is also supported, in which case you can choose to return list. If missing of x, the string will be read from system clipboard automately (I like it!).
#' @param to_clipboard If `False` return a charactor vector of gene symbols (default), if `TRUE` return value to clipboard
#' @param as_list If `FALSE` (default) return value is a vector. If `TRUE` reutrn value is structured in list, this is useful when `to_clipboard=T` so that you can paste the return (symbols with quoted) to python.
#' @param print print the return or not. (default: FALSE)
#' @param REGEX a regular expression used for identifying chars of gene symbols. Default setting should work in most cases. You can also customize it yourself.
#' @param tolower 1st prior.
#' @param capitalize 2nd prior. Use both tolower and capitalize to get mouse gene symbol.
#' @param toupper 3nd prior
#'
#' @return A charactor vector of gene symbols or A list of vectors (with names if x named) in case of as_list=T
#' @export
#'
#' @examples
#' # you can copy some text from anywhere (a research article for example) and then str_to_gene() will help to parse it into gene vector.
#' str_to_gene()
#' # custom string of genes
#' str_to_gene(x = "Runx1, Gata4 Gata1;, Dll4, Nkx2-5, NOTAGENE") # note: the function does not check validity of gene symbols
#' str_to_gene(x = c(Hematopoietic ="Runx1 Gata4 GATA1", Ery = "Gypa Gype Ptprc", Endothelial = "Pecam1 Cdh5"), as_list = T)
str_to_gene <- function(
  x, tolower = F, capitalize = F, toupper = F,
  to_clipboard = F, as_list = F,
  print = F, REGEX = "[^a-zA-Z0-9\\-\\.]"){
  if(missing(x)) {
    if(clipr::clipr_available()){
      x <- clipr::read_clip()
    }else{
      stop("system clipboard is not available")
    }
  }
  geneList <-
    stringr::str_split( # split by space
      stringr::str_squish( # squish multiple sapces into one space
        stringr::str_replace_all( # replace all non-symbol-charactor with a space
          x, REGEX, " ")
      ), " ")
  names(geneList) <- names(x) # set names if available

  if(!as_list) {
    geneList <- unlist(geneList)
    if(anyDuplicated(geneList)){
      warning("returned gene list is redundant! you may use unique() to deduplicate it.")
    }

    if(tolower) geneList <- tolower(geneList)
    if(capitalize) geneList <- Hmisc::capitalize(geneList)
    if(toupper) geneList <- toupper(geneList)

  }else{
    dup_ind <- which(sapply(geneList, anyDuplicated) > 0)
    if(length(dup_ind)>0){
      dup_ind <- paste(dup_ind, collapse = ', ')
      message(stringr::str_glue("There are duplicated genes in list of index {dup_ind}"))
    }
    if(tolower) geneList <- lapply(geneList, tolower)
    if(capitalize) geneList <- lapply(geneList, Hmisc::capitalize)
    if(toupper) geneList <- lapply(geneList, toupper)
  }



  if(print) print(geneList)
  if(to_clipboard) suppressWarnings(clipr::write_clip(geneList)) else return(geneList)
}



### str_to_char()

#' convert string into a vector of each char
#'
#' @param x a string or a vector of strings. Note if a vector supplied, it will be collapsed. If you wanna multiple vectors of char returned, you can use sapply() to wrapper the function.
#' @param unique If `TRUE` chars should be uniqued, otherwise as it is.
#' @param sort If `TRUE` chars should be sorted by alphabetic, otherwise as it is.
#' @param ... more parameters to be used in sort(). such as order decreasing = T if you like.
#'
#' @return a vector of chars
#' @export
#'
#' @examples
#' str_to_char(x = c("Here it is"), unique = F, sort = F)
#' str_to_char(x = c("Yes", "or", "Not"), unique = F, sort = F)
#' sapply(c("Yes", "or", "Not"), str_to_char, sort = F, unique = F)
str_to_char <- function(x, unique = F, sort = F, ...){
  chars <- stringr::str_split(
    stringr::str_flatten(unlist(x)),
    "", simplify = F)[[1]]
  if(unique) chars <- unique(chars)
  if(sort) chars <- sort(chars, ...)
  return(chars)
}


#' Write clipboard
#'
#' A wrapper of clipr::write_clip.
#' Write a character vector to the system clipboard
#'
#' @param content An object to be written to the system clipboard.
#' @param quoted quote each character
#' @param ... other params passed to clipr::write_clip
#'
#' @return a string
#' @export
#'
#' @examples
#'
write_clip <- function(content, quoted = T, ...){
  clipr::write_clip(list(content), ...)
}


#' load saved RData
#'
#' A wrapper function of base::load
#' This is usefull when the names of object conflict with variables in .GlobalEnv
#'
#' @param file a (readable binary-mode) connection or a character string giving the name of the file to load (when tilde expansion is done).
#' @param verbose should item names be printed during loading?
#'
#' @return A list of objects created.
#' @export
#'
#' @examples
#'
readRda <- function(file, verbose = F){
  names <- base::load(file, verbose = verbose)
  return(base::mget(names))
}


#' smooth each row of matrix
#'
#' @param mat matrix
#' @param n smooth over n columns
#' @param proportion smooth over proportion*ncol(mat)
#' @param do.scale scale mat after being smoothed
#'
#' @return matrix
#' @export
#'
#' @examples
#'
rowSmooth <- function(mat, n = NULL, proportion = 0.2, do.scale = F){
  if(is.null(n)){
    n <- ceiling(ncol(mat) * proportion)
  }
  if(n == 1) stop("nothing should be smoothed!")
  res <- mat
  for(i in 1:nrow(mat)){
    res[i,] <- vectorSmooth(as.numeric(mat[i,,drop = T]), n)
  }

  if(do.scale){
    res <- t(scale(t(res)))
  }
  return(res)
}

#' smooth a vector
#'
#' @param x a vector
#' @param n smooth over n items
#'
#' @return vector
#' @export
#'
#' @examples
#'
vectorSmooth <- function(x, n){
  sapply(1:length(x), function(i){
    m <- 1:n + i - ceiling(n/2)
    m <- m[m > 0 & m <= length(x)]
    mean(x[m])
  })
}
