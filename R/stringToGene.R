
str_to_gene <- function(x, doPrint = T){
  if(missing(x)) {
    if(clipr::clipr_available()){
      clipr::
    }
  }
  tmp <- unlist(str_split(str_squish(str_replace_all(x, "[^a-zA-Z0-9\\-\\.]", " ")), " "))
  if(doPrint) print(tmp)
  return(tmp)
}

str_collapse <- function(x){
  #sort(unique(strsplit(paste0(x, sep = "", collapse = ""), "")[[1]]))
  require(stringr)
  unique(str_sort(str_split(str_flatten(x), "", simplify = F)[[1]]))
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
