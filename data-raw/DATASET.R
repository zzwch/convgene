## code to prepare `DATASET` dataset goes here

# mSTRT barcodes
strt96barcodes <- read.csv(file = "data-raw/96-8bp-barcode.csv")
strt96barcodes <- paste0("sc", strt96barcodes$BarID)

# genesets
genesets <- list()
gs_mmu <- list(cc.genes = Seurat::cc.genes,
               cc.genes.updated.2019 = Seurat::cc.genes.updated.2019,
               cc.genes.union = mapply(union, Seurat::cc.genes, Seurat::cc.genes.updated.2019)
               )
gs_hsa <- list(cc.genes = Seurat::cc.genes,
               cc.genes.updated.2019 = Seurat::cc.genes.updated.2019,
               cc.genes.union = mapply(union, Seurat::cc.genes, Seurat::cc.genes.updated.2019))
##
sex_genes <- c("XIST", "TSIX", "RPS4Y1", "EIF2S3Y", "DDX3Y")

# scanpy_colors
scanpy_colors <- list(
  vega_10_scanpy = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#7f7f7f', '#b5bd61', '#17becf'),
  default_20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8',
                  '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'),
  default_26 = c('#023fa5', '#7d87b9', '#bec1d4', '#d6bcc0', '#bb7784', '#8e063b', '#4a6fe3', '#8595e1', '#b5bbe3', '#e6afb9',
                 '#e07b91', '#d33f6a', '#11c638', '#8dd593', '#c6dec7', '#ead3c6', '#f0b98d', '#ef9708', '#0fcfc0', '#9cded6',
                 '#d5eae7', '#f3e1eb', '#f6c4e1', '#f79cd4', '#7f7f7f', '#c7c7c7', '#1CE6FF', '#336600'),
  default_64 = c('#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6',
                 '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87', '#5A0007', '#809693', '#FEFFE6', '#1B4400', '#4FC601',
                 '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A', '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA',
                 '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018', '#0AA6D8', '#013349', '#00846F', '#372101',
                 '#FFB500', '#C2FFED', '#A079BF', '#CC0744', '#C0B9B2', '#C2FF99', '#001E09', '#00489C', '#6F0062', '#0CBD66',
                 '#EEC3FF', '#456D75', '#B77B68', '#7A87A1', '#788D66', '#885578', '#FAD09F', '#FF8A9A', '#D157A0', '#BEC459',
                 '#456648', '#0086ED', '#886F4C', '#34362D', '#B4A8BD', '#00A6AA', '#452C2C', '#636375', '#A3C8C9', '#FF913F',
                 '#938A81', '#575329', '#00FECF', '#B05B6F', '#8CD0FF', '#3B9700', '#04F757', '#C8A1A1', '#1E6E00', '#7900D7',
                 '#A77500', '#6367A9', '#A05837', '#6B002C', '#772600', '#D790FF', '#9B9700', '#549E79', '#FFF69F', '#201625',
                 '#72418F', '#BC23FF', '#99ADC0', '#3A2465', '#922329', '#5B4534', '#FDE8DC', '#404E55', '#0089A3', '#CB7E98',
                 '#A4E804', '#324E72', '#6A3A4C'),
  vega_10 = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'),
  vega_20 = c('#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
              '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5')
)

# discret colors
discrete_colors <- list(
  phase = setNames(c("#46ACC8", "#E2D200", "#DD8D29"), c("G0/G1", "S", "G2/M"))
)

# lightening colors
gradient_colors <- list(
  f_gOrRd = function(n) c("grey",RColorBrewer::brewer.pal(n, name = "OrRd")),
  f_g9YlGn = function(n) c("grey90", RColorBrewer::brewer.pal(n, "YlGn")[-1]),
  f_DbGrBr = function(n, interpolate = "linear")
    colorRampPalette(c("darkblue", "royalblue","grey100","tomato", "brown"),
                     interpolate = interpolate, alpha = T)(n),
  f_DbWtRd = function(n, interpolate = "linear")
    colorRampPalette(c("darkblue","white", "red"),
                     interpolate = interpolate, alpha = T)(n),
  f_BuWtRd = function(n, interpolate = "linear")
    colorRampPalette(c("blue","white", "red"),
                     interpolate = interpolate, alpha = T)(n),
  f_BuRd = function(n) rev(RColorBrewer::brewer.pal(n, "RdBu")),
  f_BuYlRd = function(n) rev(RColorBrewer::brewer.pal(n, "RdYlBu")),
  f_Spectral = function(n) rev(RColorBrewer::brewer.pal(n, "Spectral")),
  f_brewer = function(n, name) rev(RColorBrewer::brewer.pal(n, name)),
  gyr = c("grey", "yellow", "red"),
  g9gorb = c("grey90","grey", "orange", "red","brown"),
  gg9gorb = c( "grey96", "grey90","gray","orange", "red","brown")
)

##########################
usethis::use_data(scanpy_colors, discrete_colors, gradient_colors,
                  strt96barcodes, sex_genes,
                  overwrite = TRUE)
