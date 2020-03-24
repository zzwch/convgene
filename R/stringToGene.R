str_to_gene <- function(x, to_clipboard = F, as_list = F, print = T){
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
          x, "[^a-zA-Z0-9\\-\\.]", " ")
        ), " ")
  names(geneList) <- names(x) # set names if available
  if(!as_list) {
    geneList <- unlist(geneList)
    if(anyDuplicated(geneList)){
      message("returned gene list is redundant! you may use unique() to deduplicate it.")
    }
  }else{
    dup_ind <- paste(which(sapply(geneList, anyDuplicated) > 0), collapse = ', ')
    message(stringr::str_glue("There are duplicated genes in list of index {dup_ind}"))
  }
  if(print) print(geneList)
  if(to_clipboard) clipr::write_clip(geneList) else return(geneList)
}

str_to_char <- function(x, unique = T, sort = T, ...){
  chars <- stringr::str_split(
    stringr::str_flatten(unlist(x)),
    "", simplify = F)[[1]]
  if(unique) chars <- unique(chars)
  if(sort) chars <- sort(chars, ...)
  return(chars)
}

# vector insertion
insertVector <- function(x, index, value){
  # insert value after index of x
  # if index == 0, value will be insert at the begin of x
  # negative index is also supported
  # index larger than length of x will be warned, but allowed with NAs filled in x.
  if(any(duplicated(index))) stop("index is duplicated")
  index[index < 0] <- length(x)+1+index[index < 0] # in case of negative index
  if(any(duplicated(index))) stop("index is duplicated, beware of negative index")
  if(any(index < 0)) stop("beware of negative index")

  index_ord <- order(index)
  value <- value[seq_along(index)]

  index <- index[index_ord]
  value <- value[index_ord]

  if(index[length(index)] > length(x)) {
    warning("max index is larger than length of x, NA will be introduced")
    x <- c(x, rep(NA, index[length(index)] - length(x)))
  }

  ret <- unlist(lapply(seq_along(index), function(i){
    if(index[i] == 0) return(value[[i]])
    start_ind <- ifelse(i == 1, 0, index[[i-1]]) + 1
    end_ind <- index[[i]]
    return(c(x[start_ind:end_ind], value[[i]]))
  }))

  if(index[length(index)] < length(x)) ret <- c(ret, x[(index[length(index)]+1):length(x)])

  return(ret)
}
