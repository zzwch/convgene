

#' normalize UMI matrix to TPM
#'
#' @param umi umi matrix
#' @param scale.factor scale factor. default 1e6. recommend the mean of colSums(umi)
#' @param margin Cells in which margin
#'
#' @return a matrix of tpm values
#' @export
#'
#' @examples
#'
umi2tpm <- function(umi, scale.factor= 1e6, margin = 2){
  apply(umi, margin, function(x){
    scale.factor*x/sum(x)
  })
}

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
    macaque = c(dataset = "mmulatta_gene_ensembl",
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


#' Conver MSigDB Human genes to Homologous Organism genes
#'
#' @param files MSigDB files in gmt format
#' @param organism conversion organism
#' @param use.entrez set TRUE if gmt files use entrez gene id but not symbol
#' @param use.biomaRt use BioMart package to do conversion
#' @param hcop.file Not Supportted for now! use HGNC Comparison of Orthology Predictions (HCOP) to do conversion. http://www.genenames.org/cgi-bin/hcop
#' @param try.most Not used
#'
#' @return files
#' @export
#'
#' @examples
#'
homoMSigDB <- function(files = NULL, organism = "mouse", use.entrez = F, use.biomaRt = T, hcop.file = NULL, try.most = T){

  msigdb.files <- files

  genesets <- lapply(msigdb.files, function(x){
    read.gmt(gmtfile = x)
  })
  names(genesets) <- msigdb.files
  allinone <- unique(Reduce(union, lapply(genesets, function(x) unlist(x$gmt))))

  organism <- match.arg(organism)

  homotable <- list(mouse =  c("mmusculus_gene_ensembl", "mgi_symbol"))

  if(use.biomaRt){
    if(is.null(hcop.file)){
      human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      homologous = biomaRt::useMart("ensembl", dataset = homotable[[organism]][1])

      if(use.entrez){
        allinone.mapping <- getLDS(attributes = c("entrezgene"),
                                   filters = "entrezgene", values = allinone, mart = human,
                                   attributesL = c("entrezgene"), martL = homologous)
        colnames(allinone.mapping) <- c("human",organism)
        allinone.mapping[["human"]] <- as.character(allinone.mapping[["human"]])# may be speed up when %in%
        allinone.mapping[[organism]] <- as.character(allinone.mapping[[organism]])#

      }else{

        allinone.mapping <- getLDS(attributes = c("hgnc_symbol"),
                                      filters = "hgnc_symbol", values = allinone, mart = human,
                                      attributesL = c(homotable[[organism]][2]), martL = homologous)
        colnames(allinone.mapping) <- c("human", organism)
      }
    }else{
      stop("Please set hcop.file = NULL when use.biomaRt = TRUE, otherwise confusing!")
    }
  }else{
    if(is.null(hcop.file)){
      stop("Please set hcop.file to the homolog file when use.biomaRt = FALSE, otherwise confusing!\nPlease download mapping file from HGNC Comparison of Orthology Predictions (HCOP). http://www.genenames.org/cgi-bin/hcop")
    }else{
      stop("Sorry, havn't been done! contact lizc07@vip.qq.com")

      ## to be continued
      tmp <- try(hcop <- read.table(file = hcop.file, header = T, sep = "\t", row.names = NULL, quote = ""))
      if(class(tmp) == "try-error"){
        stop("\nCannot read your provided file!\nPlease download mapping file from HGNC Comparison of Orthology Predictions (HCOP). http://www.genenames.org/cgi-bin/hcop")
      }else{
        all.entrez.mapping <- unique(hcop[which(hcop$human_entrez_gene %in% all.entrez), c("human_entrez_gene", "mouse_entrez_gene", "mouse_symbol")])
        all.symbols.mapping <- unique(hcop[which(hcop$human_symbol %in% all.symbols), c("human_symbol", "mouse_entrez_gene", "mouse_symbol")])


      }

    }

  }

  lapply(names(genesets), function(x){
    gmt <- genesets[[x]][["gmt"]]
    desc <- genesets[[x]][["desc"]]
    for(i in names(gmt)){
      gmt[[i]] <- unique(allinone.mapping[[organism]][which(allinone.mapping[["human"]] %in% gmt[[i]])])
    }
    write.gmt(gmt, desc, presuf_file_name(x, suffix = paste0(".",organism, ".", Sys.Date()), index = -1))
  })

  return()
}

#' Read gmt file
#'
#' @param gmtfile file name
#'
#' @return a list of genesets, each with both gene list and description list
#' @export
#'
#' @examples
#'
read.gmt <- function(gmtfile = NULL){
  gmt <- list()
  desc <- list()
  for(i in readLines(con = gmtfile)){
    tmp <- strsplit(i, split = "\t")[[1]]
    gmt[[tmp[1]]] <- tmp[-c(1,2)]
    desc[[tmp[1]]] <- tmp[2]
  }
  return(list(gmt = gmt, desc = desc))
}

#' write gmt file
#'
#' @param gmt.list geneset list
#' @param desc.list description list
#' @param gmtfile file name
#'
#' @return NULL
#' @export
#'
#' @examples
#'
write.gmt <- function(gmt.list = NULL, desc.list = NULL, gmtfile = NULL){
  text <- NULL
  for(i in names(gmt.list)){
    text <- c(text, paste(c(i, desc.list[[i]], gmt.list[[i]]), collapse = "\t"))
  }
  writeLines(text = text, con = gmtfile)
  return()
}

#' Add prefix and/or suffix for file name
#'
#' @param file_name file name to be added
#' @param prefix prefix character
#' @param suffix suffix character
#' @param sep split file name by sep into multiple parts
#' @param index add prefix to the index part. if <= 0, reverse order is used, then 0 is the last one and -1 is the last but one.
#'
#' @return new file name
#' @export
#'
#' @examples
#' presuf_file_name("./database/msigdb.gmt", "homo_")
#' # "./database/homo_msigdb.gmt"
#'
presuf_file_name <- function(file_name, prefix = NULL, suffix = NULL, sep = ".", index = 1){
  dir <- dirname(file_name)
  base <- basename(file_name)
  base.split <- strsplit(base, split = sep, fixed = T)[[1]]
  if(index <= 0){
    index <- length(base.split) + index
  }
  base.split[index] <- paste0(prefix, base.split[index], suffix)
  file.path(dir, paste(base.split, collapse = sep))
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

#' Check Function
#'
#' @param x a vector
#' @param strict if strict, duplicated x will result in stop.
#'
#' @return deduplicated x
#' @export
#'
#' @examples
#'
CheckDuplicate <- function(x, strict = T, rm_dup = T){
  if(anyDuplicated(x)){
    if(strict){
      stop("there are duplicated elements in the vector")
    }
    warning("there are duplicated elements in the vector")
    if(rm_dup){
      return(x[!duplicated(x)])
    }
  }
  return(x)
}



#' Check Feature Duplicate and Not In Data
#'
#' @param object data or seurat3 object
#' @param features features
#' @param rm_dup remove duplicated elements
#' @param rm_notin remove features not in data
#'
#' @return updated features
#' @export
#'
#' @examples
#'
CheckFeatures <- function(object, features, rm_dup = T, rm_notin = T){
  features <- CheckDuplicate(features, strict = F, rm_dup = rm_dup)
  features_ <- setdiff(features, rownames(object))
  if(length(features_) > 0){
    Message("NOT-IN-DATA: ",features_, ", ")
    if(rm_notin){
      features <- features[!features %in% features_]
    }
  }
  return(features)
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


#' Correct sign by correlation
#'
#' @param x a vector to be corrected
#' @param y a vector as reference with equal length of x
#'
#' @return corrected x
#' @export
#'
#' @examples
#'
SignCorrectionByCor <- function(x, y){
  if(cor(x, y ) < 0) return(-x) else return(x)
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


#' write GSEA files
#'
#'
#' @param data feature by cell matrix
#' @param group_by cell annotations
#' @param groups subset groups
#' @param file_prefix gsea file prefix
#' @param add.sys.time add time tag to gsea file name
#'
#' @return produce two files of .txt and .cls
#' @export
#'
#' @examples
#'
write_dataToJavaGSEA <- function(data, group_by, groups = NULL, file_prefix = "gsea", add.sys.time = T){
  ret <- SubsetDataAndGroup(data, group_by, groups)
  data <- ret$data; group_by <- ret$group_by

  time_suffix <- ifelse(add.sys.time, base::Sys.time(), "")

  gsea_data <- cbind(rownames(data), NA, as.matrix(data))
  colnames(gsea_data)[1:2] <- c("Name", "Description")
  utils::write.table(gsea_data,
                     file = paste(file_prefix, time_suffix, "txt",sep = "."),
                     quote = F, sep = "\t", row.names = F)

  file_cluster <- paste(file_prefix, "cluster", time_suffix, "cls",sep = ".")
  write(paste(length(group_by), length(unique(group_by)), 1, sep = " "), file = file_cluster)
  write(paste("#", paste(unique(group_by), collapse = " "), sep = " "), file = file_cluster, append = T)
  write(paste(group_by, collapse = " "), file = file_cluster, append = T)
}

#' write SCENIC data for in-house scenic-pureR
#'
#' @param object seurat3 object
#' @param group_by vector of colnames
#' @param colVars colors vector
#' @param file save data to this file
#' @param slot fetch data
#' @param assay fetch data
#'
#' @return a saved file
#' @export
#'
#' @examples
#'
write_seuratToScenic <- function(object, group_by, colVars, file, slot = "data", assay = "RNA"){
  expr <- as.matrix(GetAssayData(object, slot, assay))
  cellInfo <- object@meta.data[,group_by, drop = F]
  colVars = colVars
  save(expr, cellInfo, colVars, file = file)
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



#' Infer UMI from scale-factor-normalized data
#'
#' Pseudo: the smallest non-zero value (normalized) is corresponding to 1 of UMI.
#'
#' @param data scale-factor-normalized data (normalized by column)
#'
#' @return a matrix of inferred UMI data
#' @export
#'
#' @examples
#'
inferUMI <- function(data, base = exp(1)){

  umi <- apply(data,2, function(x){
    umi_1 <- min(x[x!=0])
    round((base^(x)-1)/(base^(umi_1)-1))
  })

  hist(colSums(umi))
  return(umi)
}


#' export UMI to each file by each Batch
#'
#' this may be used for preparation of files to GEO
#'
#' @param umi umi matrix
#' @param batch a vector for cell batch information
#'
#' @return csv files
#' @export
#'
#' @examples
#'
exportUMIByBatch <- function(umi, batch){
  bs <- unique(as.character(batch))
  for(i in bs){
    print(i)
    umi_ <- umi[,batch == i]
    write.csv(umi_, paste0(i, ".umi.csv"))
  }
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
#' @param n smooth over \code{n} columns
#' @param proportion smooth over \code{proportion*ncol(mat)}
#' @param do.scale scale \code{mat} after being smoothed
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
#' @param n smooth over \code{n} items
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

#' drop unused colors
#'
#' This may be usefull for annotation_colors in pheatmap or pl_heatmap
#'
#' @param colors a named list of color vectors
#' @param df data.frame to check
#' @param drop_col drop which color. NULL is all.
#'
#' @return a list of updated colors
#' @export
#'
#' @examples
#'
dropColors <- function(colors, df, drop_col = NULL){
  sapply(names(colors), function(x) {
    if(x %in% colnames(df)){
      col_ind <- names(colors[[x]]) %in% df[[x]]
      colors[[x]][col_ind]
    }else{
      colors[[x]]
    }
  }, simplify = F)
}

#' convert a color name to RGB character
#'
#' @param col color name
#' @param alpha transparent, 0 to 255
#'
#' @return a character
#' @export
#'
#' @examples
#'
color2RGB <- function(col, alpha = 255){
  grDevices::rgb(t(grDevices::col2rgb(col)), maxColorValue = 255, alpha = alpha)
}


#' convert a list of vectors with different lenght to data.frame
#'
#' @param x a list of vectors
#'
#' @return a data.frame
#' data.frame is easier to view/export to excel
#' @export
#'
#' @examples
#'
list2df <- function(x){
  n.obs <- sapply(x, length)
  seq.max <- seq_len(max(n.obs))
  return(as.data.frame(sapply(x, function(x) {
    if(is.null(x)) x <- NA
    x[seq.max]
  })))
}


#' set calculations
#'
#' @param x set A
#' @param y set B
#'
#' @return a list of diffs, intersect and union
#' @export
#'
#' @examples
#'
setSummary <- function(x, y){
  list(diff_xy = setdiff(x,y),
       diff_yx = setdiff(y, x),
       intersect = intersect(x,y),
       union = union(x,y))
}


#' wrapper for quick transform from human gene to mouse gene
#'
#' @param x vector
#' @param tolower do \code{tolower} first, then \code{Hmisc::capitalize}
#'
#' @return vector
#' @importFrom Hmisc capitalize
#' @export
#'
#' @examples
#'
toCapitablize <- function(x, tolower = T){
  if(tolower)  Hmisc::capitalize(tolower(x)) else Hmisc::capitalize(x)
}


#' Generate Homologous Gene ID Table from MGI's HOM_AllOrganism.rpt file
#'
#' @param MGI_file HOM_AllOrganism.rpt downloaded from http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt
#' @param colnames_keys should not be changed!
#' @param multiple_orgs organisms, in \code{Common Organism Name} column of \code{MGI_file}, to be collapsed or compared,
#' @param multiple_orgs_short short prefix for the organisms in \code{mulitple_orgs}
#' @param centric_org which organism's Symbol should be used as HomoloSymbol
#'
#' @return a data.frame table
#' @import tidyverse
#' @import magrittr
#' @export
#'
#' @examples
#'
getMGIHomoGeneTable <- function(
  MGI_file,
  colnames_keys = c("HomoloGene ID", "Common Organism Name", "Symbol"),
  multiple_orgs = c("human", "mouse, laboratory", "macaque, rhesus"),
  multiple_orgs_short = c("Hs", "Mm", "Mc"),
  centric_org = "human"
){

  orgs <- read.table(
    file = MGI_file,
    header = T, sep = "\t",
    check.names = F,
    quote = "", comment.char = "",
    stringsAsFactors = F)

  hm <- orgs %>%
    filter(`Common Organism Name` %in% multiple_orgs) %>%
    mutate(`Common Organism Name` = factor(`Common Organism Name`, levels = multiple_orgs)) %>%
    dplyr::select(colnames_keys) %>%
    set_colnames(c("homoloID", "organism", "gene")) %>%
    #mutate_at("organism", function(x) str_replace(x, ",.*", "")) %>%
    unique

  hm_bid <- tapply(hm$gene, list(hm$homoloID, hm$organism), c)
  hm_len <- as.data.frame(apply(hm_bid, c(1,2), function(x) length(unlist(x))))

  hm_sym <- hm_len %>%
    rownames_to_column %>%
    mutate(
      centric_symbol = sapply(hm_bid[,centric_org], function(x) ifelse(is.null(unlist(x)[1]), "", unlist(x)[1])),
      hm_prefix = apply(hm_len, 1, function(x){
        if(any(x == 0)) as.character(NA) else paste0(multiple_orgs_short[x>1], collapse = "")
      }),
      homoloSymbol = case_when(
        hm_prefix == "" ~ centric_symbol,
        is.na(hm_prefix) ~ hm_prefix,
        TRUE ~ as.character(stringr::str_glue("{hm_prefix}({centric_symbol})"))
      )) %>%
    column_to_rownames

  table(hm_len$hm_prefix)
  # Hs    Mm  MmHs
  # 16468    60   212    10

  hm_res <- hm %>%
    mutate(homoloSymbol = hm_sym[as.character(homoloID), "homoloSymbol"]) %>%
    dplyr::select(homoloID, homoloSymbol, gene, organism)

  return(hm_res)
}

#' vlookup with multiple keys supported
#'
#' @param df data.frame
#' @param x values to lookup; data.frame of values per row if multiple keys used
#' @param keys colnames of \code{df} as key to search
#' @param select colname(s) to return
#' @param return.key whether to return keys
#'
#' @return vector or data.frame
#' @export
#'
#' @examples
#'
vlookup <- function(df, x, keys = 1, select = NULL, collapse = "_", return.key = FALSE){
  if(is.null(select)) {
    select <- colnames(df)
  }

  if(length(keys) > 1) {
    x <- apply(x, 1, paste, collapse = collapse)
    df$key_collapsed <- apply(df[,keys], 1, paste, collapse = collapse)
    keys <- "key_collapsed"
  }

  ret <- df[match(x, df[,keys]), select, drop = F]
  if(return.key) ret <- cbind(x = x, ret)

  return(ret)
}


#' calculate mean of gene expression values
#'
#' @param x gene expression values in log scale
#' @param log.base base of log scale
#' @param return.log mean returned in log scale or not
#'
#' @return mean value
#' @export
#'
#' @examples
#'
calcLogMean <- function(x, log.base = exp(1), return.log = T){
  m <- mean(log.base^x)
  if(return.log) m <- log(m)/log(log.base)
  return(m)
}

#' calculate Fold Change of gene expression values
#'
#' @param x,y gene expression values in log scale
#' @param log.base base of log scale
#' @param exact usually gene expression values were scaled and plus a pseudo count (eg. 1) and then log transformed. represented as log(UMI/10+1). Fold change would be not exact if do not minus the pseudo count, although pseudo count is usefull when denominator is zero.
#' @param exact.plus pseudo count, usualy be 1.
#' @param return.log fold change value returned in log scale or not
#'
#' @return fold change value of x compared to y
#' @export
#'
#' @examples
#'
calcFoldChange <- function(x, y, log.base = exp(1), exact = T, exact.plus = 1, return.log = F){
  if (!exact) exact.plus <- 0
  fc <- (calcLogMean(x, log.base, return.log = F)-exact.plus)/(calcLogMean(y, log.base, return.log = F) - exact.plus)
  if (return.log) fc <- log(fc)/log(log.base)

  return(fc)
}


#' calculate the probability of collecting all kinds of coupons when collecting n coupons.
#'
#' @param prob probabilitys for each kind of coupon
#' @param k default length of prob
#' @param n.max collected 1:n.max coupons
#' @param B an integer specifying the number of replicates used in the Monte Carlo test.
#' @param seed for reproducibility
#' @param replace sample wi/o replacing
#'
#' @return a vector of probabilities for 1:n.max cllections
#' @export
#'
#' @examples
#' # Given that 6 kinds of coupons, of which each's probability is 1/6,
#' # Let's see the probability of collection all kinds of coupons for having n coupons.
#'
#' pmultinom_sample(rep(1,6)/6)
pmultinom_sample <- function(prob, k = length(prob), n.max = 100, B = 1000, seed = 66, replace = TRUE){
  result <- NULL
  set.seed(seed)
  for (i in seq_len(B)) {
    result <- rbind(result, sample(paste0("k", seq_len(k)), n.max, replace=replace, prob=prob))
  }
  res_p <- NULL
  for (j in 1:n.max) {
    res_p[j] <- sum(apply(result[,1:j,drop = F], 1, function(x) length(unique(x))) == k)/B
  }

  return(res_p)
}


#' unlist by n times
#'
#' @param x list
#' @param depth n depth
#' @param use.names ogical. Should names be preserved?
#'
#' @return NULL or an expression or a vector of an appropriate mode to hold the list components.
#' @export
#'
#' @examples
#'
unlist_depth <- function(x, depth = NA, use.names = T){
  if(is.na(depth) || is.infinite(depth) || is.null(depth)){
    unlist(x, recursive = T, use.names = use.names)
  }else{
    run_fun_n_times(unlist, x, times = depth, recursive = F)
  }
}

#' Run function n times with output as input
#'
#' @param fun function
#' @param x input also output
#' @param times times
#' @param ... other params
#'
#' @return x
#' @export
#'
#' @examples
#'
run_fun_n_times <- function(fun, x, times = 1, ...){
  for(i in seq_len(times)){
    x <- fun(x, ...)
  }
  return(x)
}


#' read table from clipboard
#'
#' @inheritParams utils::read.table
#' @param ... other params passed to \code{\link{utils::read.table}}
#'
#' @return a data.frame
#' @export
#'
#' @examples
#'
readDataFrameFromClipboard <- function(header = T, sep = "\t", check.names = F, ...){
  read.table(file = 'clipboard', header = header, sep = sep, check.names = check.names, ...)
}


#' write data.frame to clipboard
#'
#' @param sep separator
#' @inheritParams utils::write.table
#' @param ... other params passed to \code{\link{utils::write.table}}
#'
#' @return
#' @export
#'
#' @examples
#'
writeDataFrameToClipboard <- function(x, sep = "\t", row.names = F, col.names = F, ...){
  write.table(x = x, file = "clipboard", sep = sep, row.names = row.names, col.names = col.names, ...)
}

#' set levels for a given column
#'
#' @param df data.frame
#' @param x a colname of \code{df}
#' @param levels levels for \code{df[[x]]}
#'
#' @return updated data.frame
#' @export
#'
#' @examples
#'
setLevels <- function(df, x, levels){
  df[[x]] <- factor(df[[x]], levels = levels)
  return(df)
}


#' Show colours
#'
#' A alias function for scales::show_col
#'
#' @inheritParams scales::show_col
#' @return
#' @export
#'
#' @examples
#'
show_colors <- scales::show_col


#' revert vector's levels
#'
#' @param x a vector
#' @param levels relevel x using `levels`
#'
#' @return updated x with levels reverted
#' @export
#'
#' @examples
#'
revLevels <- function(x, levels = NULL){
  if(is.null(levels)) levels <- rev(sort(unique(x)))
  return(factor(x, levels = levels))
}

