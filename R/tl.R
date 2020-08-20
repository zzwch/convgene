#' @include utils.R
#'
NULL

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
