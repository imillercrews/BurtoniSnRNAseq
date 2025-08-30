#### Burtoni snseq seurat analysis
### Neuropeptide network 
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### load libraries ####
##install libraries
# install.packages("tidyverse",
#                  repos = "https://cloud.r-project.org")
# install.packages("Seurat",
#                  repos = "https://cloud.r-project.org")
# install.packages("patchwork",
#                  repos = "https://cloud.r-project.org")
#install.packages('BiocManager')
#BiocManager::install('multtest')
# install.packages('metap',
#                  repos = "https://cloud.r-project.org")
# install.packages('clustree')
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DEsingle")
# install.packages("mathjaxr")
# # install older version of metan that works with R 4.0.3
# install.packages("https://cran.r-project.org/src/contrib/Archive/metan/metan_1.16.0.tar.gz", 
#                  repos=NULL, 
#                  type="source")
# install.packages('ecodist')


#load libraries
library(sp)
library(SeuratObject, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(patchwork)
# library(clustree)
# library(pheatmap)
# library(DEsingle)
library(igraph)
# library(ggridges)
# library(forcats)
# library(metan)
# library(ecodist)
library(tidyverse)
library(ggrepel)
library(ggraph)
library(kSamples)


#### load data ####
### single cell data
# load('burtoni.snseq.combined.sct.all.RData')
load('burtoni.snseq.combined.sct.RData')

### scsorter data
load("burtoni.scsorter.data.scsort.output.RData")

### load souporcell data
burtoni.souporcell.filtered = read.csv('../souporcell/burtoni.souporcell.filtered.csv')

### load sctype.hypo data
burtoni.sctypemarkers.hypo = read.csv('./sctype.hypo/sctypemarkers.hypo.unknown.csv')

### Add metadata
## cell type
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.scsorter.data.scsort.output %>% 
    select(Cell.type,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Cell.type'
)

## genotype
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.souporcell.filtered %>% 
    select(Genotype.id,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Genotype.id'
)

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.sctypemarkers.hypo %>% 
    select(sctypemarkers.hypo,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'sctypemarkers.hypo'
)

## select neurons
burtoni.snseq.combined.sct.all.neurons = burtoni.snseq.combined.sct
#set idents
Idents(object = burtoni.snseq.combined.sct.all.neurons) <- "sctypemarkers.hypo"

#subset to neurons
burtoni.snseq.combined.sct.all.neurons = subset(burtoni.snseq.combined.sct.all.neurons,
                                                idents = c("C7-1: GLU",
                                                           "C7-2: GABA"))

# subset with SCT data
DefaultAssay(burtoni.snseq.combined.sct.all.neurons) = 'SCT'

### neuropeptide list
neuropeptides.list = read_csv("../Gene.lists/neuropeptides.list_orthologs.csv")
# genes
neuropeptides.genes = neuropeptides.list %>% 
  filter(!is.na(Gene.name.nile.tilapia)) %>% 
  pull(Gene.name.nile.tilapia) %>% 
  unique() 


#### functions ####
### create remove_nonexp
## from: https://almeidasilvaf.github.io/BioNERO/reference/remove_nonexp.html 
#' Remove genes that are not expressed based on a user-defined threshold
#'
#' @param exp A gene expression data frame with genes in row names
#' and samples in column names or a `SummarizedExperiment` object.
#' @param method Criterion to filter non-expressed genes out.
#' One of "mean", "median", "percentage", or "allsamples". Default is "median".
#' @param min_exp If method is 'mean', 'median', or 'allsamples',
#' the minimum value for a gene to be considered expressed.
#' If method is 'percentage', the minimum value each gene must have in
#' at least n percent of samples to be considered expressed.
#' @param min_percentage_samples In case the user chooses 'percentage' as method,
#' expressed genes must have expression >= min_exp in at least this percentage.
#' Values must range from 0 to 1.
#'
#' @return Filtered gene expression data frame or `SummarizedExperiment` object.
#' @author Fabricio Almeida-Silva
#' @export
#' @importFrom matrixStats rowMedians
#' @seealso
#'  \code{\link[matrixStats]{rowMedians}}
#'  \code{\link[WGCNA]{goodSamplesGenes}}
#' @rdname remove_nonexp
#' @examples
#' data(zma.se)
#' filt_exp <- remove_nonexp(zma.se, min_exp = 5)
remove_nonexp <- function(exp, method="median", min_exp=1, min_percentage_samples=0.25) {
  # fexp <- handleSE(exp)
  fexp <- exp
  
  if(method == "median") {
    final_exp <- fexp[matrixStats::rowMedians(as.matrix(fexp)) >= min_exp,]
  } else if (method == "mean") {
    final_exp <- fexp[rowMeans(fexp) >= min_exp,]
  } else if (method == "percentage") {
    min_n <- ncol(fexp) * min_percentage_samples
    final_exp <- fexp[rowSums(fexp >= min_exp) >= min_n, ]
  } else if (method == "allsamples") {
    final_exp <- fexp[rowSums(fexp >= min_exp) == ncol(fexp), ]
  } else {
    stop("No method specified. Please, choose a filtering method - mean, median or percentage")
  }
  
  # if(is(exp, "SummarizedExperiment")) {
  #   final_exp <- exp2SE(final_exp, exp)
  # }
  
  return(final_exp)
}


#### add pvalue to mantel test figure
#' Mantel test for a set of correlation matrices
#' @description
#' `r badge('stable')`
#'
#' This function generate a pairwise matrix of plots to compare the similarity
#' of two or more correlation matrices. In the upper diagonal are presented the
#' plots and in the lower diagonal the result of Mantel test based on
#' permutations.
#'
#'
#' @param ... The input matrices. May be an output generated by the function
#'   `lpcor` or a coerced list generated by the function `as.lpcor`
#' @param type The type of correlation if an obect generated by the function
#'   `lpcor` is used. 1 = Linear correlation matrices, or 2 = partial
#'   correlation matrices.
#' @param nrepet The number of permutations. Default is 1000
#' @param names An optional vector of names of the same length of `...` .
#' @param prob The error probability for Mantel test.
#' @param diag Logical argument. If `TRUE`, the Kernel density is shown in
#'   the diagonal of plot.
#' @param export Logical argument. If `TRUE`, then the plot is exported to
#'   the current directory.
#' @param main The title of the plot, set to 'auto'.
#' @param file.type The format of the file if `export = TRUE`.  Set to
#'   `'pdf'`. Other possible values are `*.tiff` using `file.type
#'   = 'tiff'`.
#' @param file.name The name of the plot when exported. Set to `NULL`,
#'   i.e., automatically.
#' @param width The width of the plot, set to `8`.
#' @param height The height of the plot, set to `7`.
#' @param resolution The resolution of the plot if `file.type = 'tiff'` is
#'   used. Set to `300` (300 dpi).
#' @param size.point The size of the points in the plot. Set to `0.5`.
#' @param shape.point The shape of the point, set to ` 19`.
#' @param alpha.point The value for transparency of the points: 1 = full color.
#' @param fill.point The color to fill the points. Valid argument if points are
#'   between 21 and 25.
#' @param col.point The color for the edge of the point, set to `black`.
#' @param minsize The size of the letter that will represent the smallest
#'   correlation coefficient.
#' @param maxsize The size of the letter that will represent the largest
#'   correlation coefficient.
#' @param signcol The colour that indicate significant correlations (based on
#'   the `prob` value.), set to 'green'.
#' @param alpha The value for transparency of the color informed in
#'   `signcol`, when 1 = full color. Set to 0.15.
#' @param diagcol The color in the kernel distribution. Set to 'gray'.
#' @param col.up.panel,col.lw.panel,col.dia.panel The color for the opper, lower
#'   and diagonal pannels. Set to 'gray', 'gray', and 'gray', respectively.
#' @param pan.spacing The space between the pannels. Set to 0.15.
#' @seealso [mantel_test()]
#' @param digits The number of digits to show in the plot.
#' @return An object of class `gg, ggmatrix`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # iris dataset
#' lpc <- iris %>%
#'        group_by(Species) %>%
#'        lpcor() %>%
#'        pairs_mantel(names = c('setosa', 'versicolor', 'virginica'))
#'
#'
#' # mtcars dataset
#' mt_num <- select_numeric_cols(mtcars)
#' lpdata <- as.lpcor(cor(mt_num[1:5]),
#'                    cor(mt_num[1:5]),
#'                    cor(mt_num[2:6]),
#'                    cor(mt_num[4:8])) %>%
#'           pairs_mantel()
#'}
pairs_mantel.pvalue <- function (..., type = 1, nrepet = 1000, names = NULL, prob = 0.05, 
                                 diag = FALSE, export = FALSE, main = "auto", file.type = "pdf", 
                                 file.name = NULL, width = 8, height = 7, resolution = 300, 
                                 size.point = 0.5, shape.point = 19, alpha.point = 1, fill.point = NULL, 
                                 col.point = "black", minsize = 2, maxsize = 3, signcol = "green", 
                                 alpha = 0.15, diagcol = "gray", col.up.panel = "gray", col.lw.panel = "gray", 
                                 col.dia.panel = "gray", pan.spacing = 0.15, digits = 2) 
{
  class <- list(...)
  if (!type %in% c(1:2)) {
    stop("The argument type must be 1 (linear correlation) or 2 (partial correlation).")
  }
  if (sum(lapply(class, function(x) !any(class(x) %in% c("lpcor_group", 
                                                         "lpcor", "mahala_group", "covcor_design", "group_clustering", 
                                                         "clustering") == TRUE)) > 0)) {
    stop("The object must be of the class lpcor. Please use 'as.lpcorr' to convert correlation matrices into the correct format.")
  }
  if (any(class(...) == "lpcor_group")) {
    data <- lapply(...[[2]], function(x) {
      x[["linear.mat"]]
    })
  }
  if (any(class(...) == "group_clustering")) {
    data <- lapply(...[[2]], function(x) {
      x$distance
    })
  }
  if (!any(class(...) %in% c("lpcor_group", "group_clustering"))) {
    data <- lapply(..., function(x) {
      x
    })
  }
  w <- c(21:25)
  if (is.null(fill.point) == TRUE && any(w == shape.point)) {
    stop(call. = FALSE, "If 'shape.point' is a value between 21 and 25, you must provide a color for fill the shape using the argument 'fill.point.'")
  }
  for (i in 1:length(data)) {
    if (i == 1) {
      Dataset <- data.frame(var = as.vector(t(data[[1]])[lower.tri(data[[1]], 
                                                                   diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix 1")
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[1]
      }
    }
    if (i >= 2) {
      Dataset <- mutate(Dataset, var = as.vector(t(data[[i]])[lower.tri(data[[i]], 
                                                                        diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix ", 
                                                                    i)
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[i]
      }
    }
  }
  dim <- nrow(data[[1]])
  my_custom_cor <- function(data, mapping, color = I("black"), 
                            sizeRange = c(minsize, maxsize), ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    sig <- symnum(ct[[3]], corr = FALSE, na = FALSE, cutpoints = c(0, 
                                                                   0.001, 0.01, 0.05, 1), symbols = c("***", "**", 
                                                                                                      "*", ""))
    r <- ct[[1]]
    rt <- paste(format(r, digits = 2)[1],
                format(ct[[3]], digits = 2)[1],
                sep = ' , ')
    cex <- max(sizeRange)
    percent_of_range <- function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }
    GGally::ggally_text(label = as.character(rt), mapping = aes(), 
                        xP = 0.5, yP = 0.5, size = I(percent_of_range(cex * 
                                                                        abs(r), sizeRange)), color = color, ...) + geom_text(aes_string(x = 0.8, 
                                                                                                                                        y = 0.8), label = sig, size = I(cex), color = color, 
                                                                                                                             ...) + theme_classic() + theme(panel.background = element_rect(color = col.lw.panel), 
                                                                                                                                                            axis.line = element_blank(), axis.ticks = element_blank(), 
                                                                                                                                                            axis.text.y = element_blank(), axis.text.x = element_blank())
  }
  my_custom_smooth <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    pval <- ct[[3]]
    p <- ggplot(data = data, mapping = mapping)
    if (is.null(fill.point) == FALSE) {
      p <- p + geom_point(color = I(col.point), fill = fill.point, 
                          shape = shape.point, size = size.point, alpha = alpha.point)
    }
    else {
      p <- p + geom_point(color = I(col.point), shape = shape.point, 
                          size = size.point, alpha = alpha.point)
    }
    p <- p + theme_classic() + theme(panel.background = element_rect(fill = "white", 
                                                                     color = col.up.panel), axis.line = element_blank(), 
                                     axis.ticks = element_blank(), axis.text.y = element_blank(), 
                                     axis.text.x = element_blank())
    if (pval < prob) {
      p <- p + theme(panel.background = element_rect(fill = ggplot2::alpha(signcol, 
                                                                           alpha)))
    }
    p
  }
  ggally_mysmooth <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) + geom_density(fill = alpha(diagcol, 
                                                                       1)) + theme_classic() + theme(panel.background = element_rect(fill = alpha("white", 
                                                                                                                                                  1), color = col.dia.panel))
  }
  if (main == "auto") {
    title <- paste0("Mantel's test with ", nrepet, " resamples")
  }
  else {
    title <- main
  }
  if (diag == TRUE) {
    diag <- list(continuous = ggally_mysmooth)
  }
  else {
    diag <- NULL
  }
  p1 <- GGally::ggpairs(Dataset, title = title, diag = diag, 
                        lower = list(continuous = my_custom_cor), upper = list(continuous = my_custom_smooth), 
                        axisLabels = "none") + theme(panel.spacing = grid::unit(pan.spacing, 
                                                                                "lines"))
  if (export == FALSE) {
    return(p1)
  }
  else if (file.type == "pdf") {
    if (is.null(file.name)) {
      pdf("Pairs of Mantel's test.pdf", width = width, 
          height = height)
    }
    else pdf(paste0(file.name, ".pdf"), width = width, height = height)
    print(p1)
    dev.off()
  }
  if (file.type == "tiff") {
    if (is.null(file.name)) {
      tiff(filename = "Pairs of Mantel's test.tiff", width = width, 
           height = height, units = "in", compression = "lzw", 
           res = resolution)
    }
    else tiff(filename = paste0(file.name, ".tiff"), width = width, 
              height = height, units = "in", compression = "lzw", 
              res = resolution)
    print(p1)
    dev.off()
  }
}

#### add confidence interval to mantel test figure
#' Mantel test for a set of correlation matrices
#' @description
#' `r badge('stable')`
#'
#' This function generate a pairwise matrix of plots to compare the similarity
#' of two or more correlation matrices. In the upper diagonal are presented the
#' plots and in the lower diagonal the result of Mantel test based on
#' permutations.
#'
#'
#' @param ... The input matrices. May be an output generated by the function
#'   `lpcor` or a coerced list generated by the function `as.lpcor`
#' @param type The type of correlation if an obect generated by the function
#'   `lpcor` is used. 1 = Linear correlation matrices, or 2 = partial
#'   correlation matrices.
#' @param nrepet The number of permutations. Default is 1000
#' @param names An optional vector of names of the same length of `...` .
#' @param prob The error probability for Mantel test.
#' @param diag Logical argument. If `TRUE`, the Kernel density is shown in
#'   the diagonal of plot.
#' @param export Logical argument. If `TRUE`, then the plot is exported to
#'   the current directory.
#' @param main The title of the plot, set to 'auto'.
#' @param file.type The format of the file if `export = TRUE`.  Set to
#'   `'pdf'`. Other possible values are `*.tiff` using `file.type
#'   = 'tiff'`.
#' @param file.name The name of the plot when exported. Set to `NULL`,
#'   i.e., automatically.
#' @param width The width of the plot, set to `8`.
#' @param height The height of the plot, set to `7`.
#' @param resolution The resolution of the plot if `file.type = 'tiff'` is
#'   used. Set to `300` (300 dpi).
#' @param size.point The size of the points in the plot. Set to `0.5`.
#' @param shape.point The shape of the point, set to ` 19`.
#' @param alpha.point The value for transparency of the points: 1 = full color.
#' @param fill.point The color to fill the points. Valid argument if points are
#'   between 21 and 25.
#' @param col.point The color for the edge of the point, set to `black`.
#' @param minsize The size of the letter that will represent the smallest
#'   correlation coefficient.
#' @param maxsize The size of the letter that will represent the largest
#'   correlation coefficient.
#' @param signcol The colour that indicate significant correlations (based on
#'   the `prob` value.), set to 'green'.
#' @param alpha The value for transparency of the color informed in
#'   `signcol`, when 1 = full color. Set to 0.15.
#' @param diagcol The color in the kernel distribution. Set to 'gray'.
#' @param col.up.panel,col.lw.panel,col.dia.panel The color for the opper, lower
#'   and diagonal pannels. Set to 'gray', 'gray', and 'gray', respectively.
#' @param pan.spacing The space between the pannels. Set to 0.15.
#' @seealso [mantel_test()]
#' @param digits The number of digits to show in the plot.
#' @return An object of class `gg, ggmatrix`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # iris dataset
#' lpc <- iris %>%
#'        group_by(Species) %>%
#'        lpcor() %>%
#'        pairs_mantel(names = c('setosa', 'versicolor', 'virginica'))
#'
#'
#' # mtcars dataset
#' mt_num <- select_numeric_cols(mtcars)
#' lpdata <- as.lpcor(cor(mt_num[1:5]),
#'                    cor(mt_num[1:5]),
#'                    cor(mt_num[2:6]),
#'                    cor(mt_num[4:8])) %>%
#'           pairs_mantel()
#'}
pairs_mantel.ci <- function (..., type = 1, nrepet = 1000, names = NULL, prob = 0.05, 
                             diag = FALSE, export = FALSE, main = "auto", file.type = "pdf", 
                             file.name = NULL, width = 8, height = 7, resolution = 300, 
                             size.point = 0.5, shape.point = 19, alpha.point = 1, fill.point = NULL, 
                             col.point = "black", minsize = 2, maxsize = 3, signcol = "green", 
                             alpha = 0.15, diagcol = "gray", col.up.panel = "gray", col.lw.panel = "gray", 
                             col.dia.panel = "gray", pan.spacing = 0.15, digits = 2) 
{
  class <- list(...)
  if (!type %in% c(1:2)) {
    stop("The argument type must be 1 (linear correlation) or 2 (partial correlation).")
  }
  if (sum(lapply(class, function(x) !any(class(x) %in% c("lpcor_group", 
                                                         "lpcor", "mahala_group", "covcor_design", "group_clustering", 
                                                         "clustering") == TRUE)) > 0)) {
    stop("The object must be of the class lpcor. Please use 'as.lpcorr' to convert correlation matrices into the correct format.")
  }
  if (any(class(...) == "lpcor_group")) {
    data <- lapply(...[[2]], function(x) {
      x[["linear.mat"]]
    })
  }
  if (any(class(...) == "group_clustering")) {
    data <- lapply(...[[2]], function(x) {
      x$distance
    })
  }
  if (!any(class(...) %in% c("lpcor_group", "group_clustering"))) {
    data <- lapply(..., function(x) {
      x
    })
  }
  w <- c(21:25)
  if (is.null(fill.point) == TRUE && any(w == shape.point)) {
    stop(call. = FALSE, "If 'shape.point' is a value between 21 and 25, you must provide a color for fill the shape using the argument 'fill.point.'")
  }
  for (i in 1:length(data)) {
    if (i == 1) {
      Dataset <- data.frame(var = as.vector(t(data[[1]])[lower.tri(data[[1]], 
                                                                   diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix 1")
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[1]
      }
    }
    if (i >= 2) {
      Dataset <- mutate(Dataset, var = as.vector(t(data[[i]])[lower.tri(data[[i]], 
                                                                        diag = FALSE)]))
      if (is.null(names)) {
        names(Dataset)[which(colnames(Dataset) == "var")] <- paste0("Matrix ", 
                                                                    i)
      }
      else {
        names(Dataset)[which(colnames(Dataset) == "var")] <- names[i]
      }
    }
  }
  dim <- nrow(data[[1]])
  my_custom_cor <- function(data, mapping, color = I("black"), 
                            sizeRange = c(minsize, maxsize), ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    ct.2 <- ecodist::mantel(as.dist(D) ~ as.dist(D2), nboot = nrepet) # use mantel test from ecodist to get confidence interval
    sig <- symnum(ct[[2]], #change to get p-value
                  corr = FALSE, na = FALSE, cutpoints = c(0, 
                                                          0.001, 0.01, 0.05, 1), symbols = c("***", "**", 
                                                                                             "*", ""))
    r <- ct[[1]]
    rt <- paste(format(ct.2[[2]], digits = 2)[1],
                " , ",
                format(ct.2[[1]], digits = 2)[1],
                ' , ',
                format(ct.2[[5]], digits = 2)[1], # add lower confidence interval
                '-',
                format(ct.2[[6]], digits = 2)[1], # add upper confidence interval
                sep = '')
    cex <- max(sizeRange)
    percent_of_range <- function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }
    GGally::ggally_text(label = as.character(rt), mapping = aes(), 
                        xP = 0.5, yP = 0.5, size = I(percent_of_range(cex * 
                                                                        abs(r), sizeRange)), color = color, ...) + geom_text(aes_string(x = 0.8, 
                                                                                                                                        y = 0.8), label = sig, size = I(cex), color = color, 
                                                                                                                             ...) + theme_classic() + theme(panel.background = element_rect(color = col.lw.panel), 
                                                                                                                                                            axis.line = element_blank(), axis.ticks = element_blank(), 
                                                                                                                                                            axis.text.y = element_blank(), axis.text.x = element_blank())
  }
  my_custom_smooth <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    D <- matrix(nrow = dim, ncol = dim)
    D[lower.tri(D, diag = FALSE)] <- x
    D <- make_sym(D, diag = 0)
    D2 <- matrix(nrow = dim, ncol = dim)
    D2[lower.tri(D2, diag = FALSE)] <- y
    D2 <- make_sym(D2, diag = 0)
    ct <- mantel_test(D, D2, nboot = nrepet)
    pval <- ct[[3]]  
    p <- ggplot(data = data, mapping = mapping)
    if (is.null(fill.point) == FALSE) {
      p <- p + geom_point(color = I(col.point), fill = fill.point, 
                          shape = shape.point, size = size.point, alpha = alpha.point)
    }
    else {
      p <- p + geom_point(color = I(col.point), shape = shape.point, 
                          size = size.point, alpha = alpha.point)
    }
    p <- p + theme_classic() + theme(panel.background = element_rect(fill = "white", 
                                                                     color = col.up.panel), axis.line = element_blank(), 
                                     axis.ticks = element_blank(), axis.text.y = element_blank(), 
                                     axis.text.x = element_blank())
    if (pval < prob) {
      p <- p + theme(panel.background = element_rect(fill = ggplot2::alpha(signcol, 
                                                                           alpha)))
    }
    p
  }
  ggally_mysmooth <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) + geom_density(fill = alpha(diagcol, 
                                                                       1)) + theme_classic() + theme(panel.background = element_rect(fill = alpha("white", 
                                                                                                                                                  1), color = col.dia.panel))
  }
  if (main == "auto") {
    title <- paste0("Mantel's test with ", nrepet, " resamples")
  }
  else {
    title <- main
  }
  if (diag == TRUE) {
    diag <- list(continuous = ggally_mysmooth)
  }
  else {
    diag <- NULL
  }
  p1 <- GGally::ggpairs(Dataset, title = title, diag = diag, 
                        lower = list(continuous = my_custom_cor), upper = list(continuous = my_custom_smooth), 
                        axisLabels = "none") + theme(panel.spacing = grid::unit(pan.spacing, 
                                                                                "lines"))
  if (export == FALSE) {
    return(p1)
  }
  else if (file.type == "pdf") {
    if (is.null(file.name)) {
      pdf("Pairs of Mantel's test.pdf", width = width, 
          height = height)
    }
    else pdf(paste0(file.name, ".pdf"), width = width, height = height)
    print(p1)
    dev.off()
  }
  if (file.type == "tiff") {
    if (is.null(file.name)) {
      tiff(filename = "Pairs of Mantel's test.tiff", width = width, 
           height = height, units = "in", compression = "lzw", 
           res = resolution)
    }
    else tiff(filename = paste0(file.name, ".tiff"), width = width, 
              height = height, units = "in", compression = "lzw", 
              res = resolution)
    print(p1)
    dev.off()
  }
}


# #### create networks ####
# ### extract neuropeptide expression and orig.idents from neurons
# ## create expression dataframe of neuropeptides
# burtoni.snseq.combined.sct.neurons.np.expression = burtoni.snseq.combined.sct.all.neurons@assays$SCT@data %>% 
#   as.data.frame() %>% 
#   rownames_to_column('gene') %>% 
#   # mutate(gene = toupper(gene)) %>%
#   filter(gene %in% neuropeptides.genes) %>% 
#   column_to_rownames('gene') %>% 
#   t() %>% 
#   as.data.frame()
# 
# ### get list of sample ids
# burtoni.snseq.combined.sct.neurons.np.expression.samples = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
#   rownames_to_column('Cell.id') %>%  
#   full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#               rownames_to_column("Cell.id") %>% 
#               dplyr::select(c(orig.ident,
#                               Cell.id))) %>% 
#   relocate(orig.ident, .after = Cell.id)
# 
# # remove low expressing genes
# # remove low expressing cells
# burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression %>% 
#   remove_nonexp(method = "percentage",
#                 min_exp = 0.68,
#                 min_percentage_samples = .05)  %>% 
#   t() %>% 
#   remove_nonexp(method = "percentage",
#                 min_exp = 0.68,
#                 min_percentage_samples = .05) 
# 
# # count cells and genes
# ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# # 11523 cells out of 12641
# nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# # 25 genes out of 60
# 
# ###rename to actual gene name
# burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression.filter %>% 
#   t() %>% 
#   as.data.frame() %>% 
#     dplyr::rename(#ADM = ENSONIG00000036130,
#     CHGA = ENSONIG00000000725,
#     #GHRH = ENSONIG00000020252,
#     GNRH1 = ENSONIG00000011023,
#     #KNG1 = ENSONIG00000008378,
#     NMB = ENSONIG00000002537,
#     #NPFF = ENSONIG00000005229,
#     SST = ENSONIG00000033642) %>% 
#   t()
# 
# ### dom males
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom = t(burtoni.snseq.combined.sct.neurons.np.expression.filter) %>% 
#   as.data.frame()
# #create id list
# Male.Dom.cell.ids = burtoni.snseq.combined.sct.neurons.np.expression.samples %>% 
#   filter(orig.ident == 'dom_burtoni_snseq') %>% 
#   pull(Cell.id)
# #filter down cells
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom %>% 
#   filter(rownames(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom) %in% Male.Dom.cell.ids) 
# 
# 
# ## create graph object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g <- graph.adjacency(
#   as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, method="pearson"))),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE)
# 
# # Simplfy the adjacency object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g <- igraph::simplify(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g, remove.multiple=TRUE, remove.loops=TRUE)
# 
# # Colour negative correlation edges as blue
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight<0)]$color <- "negative"
# 
# # Colour positive correlation edges as red
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight>0)]$color <- "positive"
# 
# # Convert edge weights to absolute values
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight <- abs(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight)
# 
# # Remove edges below absolute Pearson correlation 0.01
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g <- delete_edges(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g, 
#                                                                                  E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight<0.01)])
# 
# ### sub males
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub = t(burtoni.snseq.combined.sct.neurons.np.expression.filter) %>% 
#   as.data.frame()
# #create id list
# Male.sub.cell.ids = burtoni.snseq.combined.sct.neurons.np.expression.samples %>% 
#   filter(orig.ident == 'sub_burtoni_snseq') %>% 
#   pull(Cell.id)
# #filter down cells
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub %>% 
#   filter(rownames(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub) %in% Male.sub.cell.ids) 
# 
# 
# ## create graph object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g <- graph.adjacency(
#   as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, method="pearson"))),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE
# )
# 
# # Simplfy the adjacency object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g <- igraph::simplify(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g, remove.multiple=TRUE, remove.loops=TRUE)
# 
# # Colour negative correlation edges as blue
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight<0)]$color <- "negative"
# 
# # Colour positive correlation edges as red
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight>0)]$color <- "positive"
# 
# # Convert edge weights to absolute values
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight <- abs(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight)
# 
# # Remove edges below absolute Pearson correlation 0.01
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g <- delete_edges(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g, 
#                                                                                    E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight<0.01)])
# 
# 
# 
# # Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# # Multiply scaled vales by a factor of 10
# scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
# vSizes.all <- (scale01(apply(burtoni.snseq.combined.sct.neurons.np.expression.filter, 1, mean)) + 1.0) 
# 
# # create size list for each sample 
# vSizes.dom = apply(t(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom), 1, mean) + 1.0
# vSizes.sub = apply(t(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub), 1, mean) + 1.0
# 
# # graph sizes
# as.data.frame(vSizes.all) %>% 
#   rownames_to_column('Gene') %>% 
#   full_join(as.data.frame(vSizes.dom) %>% 
#               rownames_to_column('Gene')) %>% 
#   full_join(as.data.frame(vSizes.sub) %>% 
#               rownames_to_column('Gene')) %>% 
#   pivot_longer(cols = c(vSizes.all,
#                         vSizes.dom,
#                         vSizes.sub),
#                names_to = 'Sample') %>% 
#   ggplot(aes(x = reorder(Gene,
#                          -value),
#              y = value,
#              color = Sample)) +
#   geom_point() +
#   theme_classic() +
#   ggtitle('Neuropeptide relative size comparison') +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   xlab('')
# ggsave('neuropeptides/comparison/network/neuropeptides relative count comparison.png')

# #### graph networks ####
# ### get weight values for each sample
# ## max
# dom.male.network.weight.max = E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight %>% 
#   max()
# sub.male.network.weight.max = E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight %>% 
#   max()
# ## min
# dom.male.network.weight.min = E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g)$weight %>% 
#   min()
# sub.male.network.weight.min = E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$weight %>% 
#   min()
# 
# ## make data frame
# network.weight.df = data.frame(max.weight = c(dom.male.network.weight.max,sub.male.network.weight.max),
#                                min.weight = c(dom.male.network.weight.min,sub.male.network.weight.min),
#                                sample = c('dom.male',
#                                           'sub.male'))
# 
# # add scale value per sample
# network.weight.df = network.weight.df %>% 
#   mutate(Group.max = max(max.weight),
#          Group.min = max(min.weight)) %>% 
#   mutate(scalingmax = 10*max.weight/Group.max,
#          scalingmin = 1*min.weight/Group.min)
# # use to set parameters per sample
# # scale_edge_width from scalingmin to scalingmax
# # scale_edge_alpha from 0 to scalingmax/10
# network.weight.df
# 
# ## males
# # dom
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g %>%
#   ggraph(layout = 'linear',
#          circular = TRUE,
#          sort.by = vSizes.all
#   ) +
#   geom_edge_parallel(aes(width = weight ,
#                          alpha = weight ,
#                          start_cap = label_rect(node1.name),
#                          end_cap = label_rect(node2.name),
#                          edge_colour = color)) +
#   geom_node_label(aes(label = name,
#                       size = vSizes.dom
#   )) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5),
#         title =element_text(size=20, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank(),
#         legend.position = 'none') +
#   ggtitle('Dominant Male') +
#   scale_edge_color_manual(values = c('blue','red')) +
#   scale_edge_width(range = c(.993,10))+
#   scale_edge_alpha(range = c(0,1)) 
# ggsave('neuropeptides/comparison/network/neuropeptides correlation network male dom.pdf',
#        height = 5.25,
#        width = 5.25)
# 
# 
# # sub
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g %>%
#   ggraph(layout = 'linear',
#          circular = TRUE,
#          sort.by = vSizes.all
#   ) +
#   geom_edge_parallel(aes(width = weight ,
#                          alpha = weight ,
#                          start_cap = label_rect(node1.name),
#                          end_cap = label_rect(node2.name),
#                          edge_colour = color)) +
#   geom_node_label(aes(label = name,
#                       size = vSizes.sub
#   )) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5),
#         title =element_text(size=20, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank(),
#         legend.position = 'none'
#   ) +
#   ggtitle('Subordinate Male') +
#   scale_edge_color_manual(values = c('blue','red')) +
#   scale_edge_width(range = c(1,3.6))+
#   scale_edge_alpha(range = c(0,.36)) 
# ggsave('neuropeptides/comparison/network/neuropeptides correlation network male sub.pdf',
#        height = 5.25,
#        width = 5.25)
# 
# 
# #### testing covariance matrices ####
# ### create covariance matrix per sample
# library(HDtest)
# library(parallel)
# 
# ### run hd test
# males.hd = testCov(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, 
#                    burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, 
#                    method = "HD",
#                    J = 1000, 
#                    alpha = 0.05, 
#                    n.core = 12)
# 
# ## run clx test
# males.CLX = testCov(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, 
#                     burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, 
#                     method = "CLX",
#                     J = 1000, 
#                     alpha = 0.05, 
#                     n.core = 12)
# 
# 
# 
# #### create dataframe
# results.testcov = data.frame(comparison = c('Male'),
#                              HD.p.val = c(males.hd[[2]]),
#                              HD.stat = c(males.hd[[1]]),
#                              CLX.p.val = c(males.CLX[[2]]),
#                              CLX.stat = c(males.CLX[[1]])
# )
# 
# results.testcov
# # comparison HD.p.val  HD.stat CLX.p.val CLX.stat
# # Statistic       Male        0 54.56802         0 2977.669
# 
# # results.testcov.short = data.frame(comparison = c('Male'),
# #                                    HD.stat = round(c(males.hd[[1]])),
# #                                    CLX.stat = round(c(males.CLX[[1]]))
# # )
# 
# ### mantel test 
# mantel.np.ci = as.lpcor(as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, method="pearson"))),
#          as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, method="pearson")))) %>%
#   pairs_mantel.ci()
# 
# # rename labels
# mantel.np.ci[["xAxisLabels"]]  = c("Male.Dom",
#                                                    "Male.Sub")
# mantel.np.ci[["yAxisLabels"]]  = c("Male.Dom",
#                                                    "Male.Sub")
# mantel.np.ci[["title"]] = paste0( mantel.np.ci[["title"]],
#                                                   ' filter 50 ci')
# mantel.np.ci
# ggsave('neuropeptides/comparison/network/Mantel test.png',
#        height = 5.25,
#        width = 5.25)
# 
# ## Mantel test: P-value = 0.003, r = 0.18, 95% CI for r = 0.07-0.25 
# 
# #### clustering coefficient/ transitivity ####
# ## global
# #dom = 0.5159744
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g %>% 
#   transitivity(type = 'global')
# #sub = 0.5771842
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g %>% 
#   transitivity(type = 'global')
# 
# 
# transitivity.local = data.frame(dom.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g %>% 
#                     transitivity(type = 'local'),
#                   dom.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g %>% 
#                     transitivity(type = 'weighted'),
#                   vertices = V(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$name) %>% 
#   full_join(data.frame(sub.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g %>% 
#                          transitivity(type = 'local'),
#                        sub.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g %>% 
#                          transitivity(type = 'weighted'),
#                        vertices = V(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g)$name))
# 
# 
# 
# ## graph
# # local
# transitivity.local %>% 
#   ggplot(aes(x = dom.local,
#              y = sub.local,
#              label = vertices)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_label() +
#   theme_classic() +
#   xlim(0.25,.75) +
#   ylim(0.25,.75)
# 
# 
# # weighted
# transitivity.local %>% 
#   ggplot(aes(x = dom.weight,
#              y = sub.weight,
#              label = vertices)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_label() +
#   theme_classic()  +
#   xlim(0.25,.75) +
#   ylim(0.25,.75)
# 
# 
# 
# #### clustering coefficient/transitivity positive vs negative correlation ####
# ### dom
# ## create graph object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans <- graph.adjacency(
#   as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, method="pearson"))),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE
# )
# 
# # Simplfy the adjacency object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans <- igraph::simplify(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans, remove.multiple=TRUE, remove.loops=TRUE)
# 
# # Colour negative correlation edges as blue
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)$weight<0)]$color <- "negative"
# 
# # Colour positive correlation edges as red
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)$weight>0)]$color <- "positive"
# 
# ## Remove edges
# # keep positive pearson correlation
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos <- delete_edges(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans, 
#                                                                                              E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)$weight<0.01)])
# 
# # keep positive pearson correlation
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.neg <- delete_edges(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans, 
#                                                                                              E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans)$weight>0.01)])
# 
# 
# # Convert edge weights to absolute values
# # positive
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos)$weight <- abs(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos)$weight)
# # negative
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos)$weight <- abs(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos)$weight)
# 
# 
# ### sub
# ## create graph object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans <- graph.adjacency(
#   as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, method="pearson"))),
#   mode="undirected",
#   weighted=TRUE,
#   diag=FALSE
# )
# 
# # Simplfy the adjacency object
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans <- igraph::simplify(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans, remove.multiple=TRUE, remove.loops=TRUE)
# 
# # Colour negative correlation edges as blue
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)$weight<0)]$color <- "negative"
# 
# # Colour positive correlation edges as red
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)$weight>0)]$color <- "positive"
# 
# ## Remove edges
# # keep positive pearson correlation
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos <- delete_edges(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans, 
#                                                                                              E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)$weight<0.01)])
# 
# # keep positive pearson correlation
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.neg <- delete_edges(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans, 
#                                                                                              E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)[which(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans)$weight>0.01)])
# 
# 
# # Convert edge weights to absolute values
# # positive
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos)$weight <- abs(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos)$weight)
# # negative
# E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos)$weight <- abs(E(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos)$weight)
# 
# 
# 
# 
# 
# ### compare dom and sub
# ### global
# #dom.pos = 0.4233871
# #sub.pos = 0.2011173
# transitivity.global.pos = data.frame(dom.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos %>% 
#                                          transitivity(type = 'global'), 
#                                      sub.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos %>% 
#                                          transitivity(type = 'global'), 
#                                      dom.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos %>% 
#                                        transitivity(type = 'global'), 
#                                      sub.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos %>% 
#                                        transitivity(type = 'global'), 
#                                      vertices = c('global'))
# 
# #dom.neg = 0.6306569
# #sub.neg = 0.8456376
# transitivity.global.neg = data.frame(dom.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.neg %>% 
#                                        transitivity(type = 'global'), 
#                                      sub.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.neg %>% 
#                                        transitivity(type = 'global'), 
#                                      dom.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.neg %>% 
#                                        transitivity(type = 'global'), 
#                                      sub.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.neg %>% 
#                                        transitivity(type = 'global'), 
#                                      vertices = c('global'))
# 
# ### local
# ## positive
# # replace NAN with 0
# # add global
# transitivity.local.pos = data.frame(dom.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos %>% 
#                                       transitivity(type = 'local'),
#                                     dom.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos %>% 
#                                       transitivity(type = 'weighted'),
#                                     vertices = V(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.pos)$name) %>% 
#   full_join(data.frame(sub.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos %>% 
#                          transitivity(type = 'local'),
#                        sub.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos %>% 
#                          transitivity(type = 'weighted'),
#                        vertices = V(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.pos)$name)) %>% 
#   full_join(transitivity.global.pos) %>% 
#   mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) %>% 
#   mutate(color = ifelse(vertices == 'global',
#                         'global',
#                         'local'))
#   
# 
# ## graph
# # local
# transitivity.local.pos %>% 
#   ggplot(aes(x = dom.local,
#              y = sub.local,
#              label = vertices,
#              color = color)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_label() +
#   theme_classic() +
#   xlim(0,1) +
#   ylim(0,1) +
#   ggtitle('Local transitivity positive') +
#   scale_color_manual(values = c('red', 'black')) +
#   theme(legend.position = 'null')
# ggsave('neuropeptides/comparison/network/Transitivity positive.png',
#        height = 5.25,
#        width = 5.25)
# 
# 
# # weighted
# transitivity.local.pos %>% 
#   ggplot(aes(x = dom.weight,
#              y = sub.weight,
#              label = vertices,
#              color = color)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_label() +
#   theme_classic()  +
#   xlim(0,1) +
#   ylim(0,1)+
#   ggtitle('Local transitivity positive weighted') +
#   scale_color_manual(values = c('red', 'black'))+
#   theme(legend.position = 'null')
# ggsave('neuropeptides/comparison/network/Transitivity positive weighted.png',
#        height = 5.25,
#        width = 5.25)
# 
# ## negative
# # replace NAN with 0
# # add global
# transitivity.local.neg = data.frame(dom.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.neg %>% 
#                                       transitivity(type = 'local'),
#                                     dom.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.neg %>% 
#                                       transitivity(type = 'weighted'),
#                                     vertices = V(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom.g.trans.neg)$name) %>% 
#   full_join(data.frame(sub.local = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.neg %>% 
#                          transitivity(type = 'local'),
#                        sub.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.neg %>% 
#                          transitivity(type = 'weighted'),
#                        vertices = V(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub.g.trans.neg)$name)) %>% 
#   full_join(transitivity.global.neg) %>% 
#   mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) %>% 
#   mutate(color = ifelse(vertices == 'global',
#                         'global',
#                         'local'))
# 
# 
# ## graph
# # local
# transitivity.local.neg %>% 
#   ggplot(aes(x = dom.local,
#              y = sub.local,
#              label = vertices,
#              color = color)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_label() +
#   theme_classic() +
#   xlim(0,1) +
#   ylim(0,1) +
#   ggtitle('Local transitivity negative') +
#   scale_color_manual(values = c('red', 'black')) +
#   theme(legend.position  = 'null')
# ggsave('neuropeptides/comparison/network/Transitivity negative.png',
#        height = 5.25,
#        width = 5.25)
# 
# 
# # weighted
# transitivity.local.neg %>% 
#   ggplot(aes(x = dom.weight,
#              y = sub.weight,
#              label = vertices,
#              color = color)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_label() +
#   theme_classic()  +
#   xlim(0,2) +
#   ylim(0,2)+
#   ggtitle('Local transitivity negative weighted') +
#   scale_color_manual(values = c('red', 'black'))+
#   theme(legend.position = 'null')
# ggsave('neuropeptides/comparison/network/Transitivity negative weighted.png',
#        height = 5.25,
#        width = 5.25)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### compare sample genes ####
# ### combine doms and subs into single data frame
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all = full_join(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom %>% 
#                                                                                rownames_to_column("ID") %>% 
#                                                                                mutate(Status = 'DOM'),
#                                                                              burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub %>% 
#                                                                                rownames_to_column("ID") %>% 
#                                                                                mutate(Status = 'SUB'))
# 
# ### graph
# ## OXT vs AVP
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all %>% 
#   ggplot(aes(x = oxt,
#              y = avp,
#              color = Status,
#              group = Status)) +
#   geom_jitter(height = 0.05,
#               width = 0.05) +
#   geom_smooth(method = 'lm') +
#   theme_bw()
# ggsave('neuropeptides/comparison/network/examples/oxt vs avp expression.png')
#        
# ## OXT vs NXPH4
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all %>% 
#   ggplot(aes(x = oxt,
#              y = NXPH4,
#              color = Status,
#              group = Status)) +
#   geom_jitter(height = 0.05,
#               width = 0.05) +
#   geom_smooth(method = 'lm') +
#   theme_bw()
# ggsave('neuropeptides/comparison/network/examples/oxt vs NXPH4 expression.png')
# 
# ## OXT vs nxph1
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all %>% 
#   ggplot(aes(x = oxt,
#              y = nxph1,
#              color = Status,
#              group = Status)) +
#   geom_jitter(height = 0.05,
#               width = 0.05) +
#   geom_smooth(method = 'lm') +
#   theme_bw()
# ggsave('neuropeptides/comparison/network/examples/oxt vs nxph1 expression.png')
# 
# ### co-expressed cells only
# ## OXT vs AVP
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all %>% 
#   filter(oxt> 0,
#          avp >0 ) %>% 
#   ggplot(aes(x = oxt,
#              y = avp,
#              color = Status,
#              group = Status)) +
#   geom_jitter(height = 0.05,
#               width = 0.05) +
#   geom_smooth(method = 'lm') +
#   theme_bw()
# ggsave('neuropeptides/comparison/network/examples/oxt vs avp coexpression.png')
# 
# ## OXT vs NXPH4
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all %>% 
#   filter(oxt> 0,
#          NXPH4 >0 ) %>% 
#   ggplot(aes(x = oxt,
#              y = NXPH4,
#              color = Status,
#              group = Status)) +
#   geom_jitter(height = 0.05,
#               width = 0.05) +
#   geom_smooth(method = 'lm') +
#   theme_bw()
# ggsave('neuropeptides/comparison/network/examples/oxt vs NXPH4 coexpression.png')
# 
# ## OXT vs nxph1
# burtoni.snseq.combined.sct.neurons.np.expression.filter.male.all %>% 
#   filter(oxt> 0,
#          nxph1 >0 ) %>% 
#   ggplot(aes(x = oxt,
#              y = nxph1,
#              color = Status,
#              group = Status)) +
#   geom_jitter(height = 0.05,
#               width = 0.05) +
#   geom_smooth(method = 'lm') +
#   theme_bw()
# ggsave('neuropeptides/comparison/network/examples/oxt vs nxph1 coexpression.png')
# 
# 
# #### presence/absence matrix ####
# ### extract neuropeptide expression and orig.idents from neurons
# ## create expression dataframe of neuropeptides
# burtoni.snseq.combined.sct.neurons.np.expression = burtoni.snseq.combined.sct.all.neurons@assays$SCT@data %>% 
#   as.data.frame() %>% 
#   rownames_to_column('gene') %>% 
#   # mutate(gene = toupper(gene)) %>%
#   filter(gene %in% neuropeptides.genes) %>% 
#   column_to_rownames('gene') %>% 
#   t() %>% 
#   as.data.frame()
# 
# 
# # remove low expressing genes
# # remove low expressing cells
# burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression %>% 
#   remove_nonexp(method = "percentage",
#                 min_exp = 0.68,
#                 min_percentage_samples = .05)  %>% 
#   t() %>% 
#   remove_nonexp(method = "percentage",
#                 min_exp = 0.68,
#                 min_percentage_samples = .05) 
# 
# # count cells and genes
# ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# # 11523 cells out of 12641
# nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# # 25 genes out of 60
# 
# ###rename to actual gene name
# burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression.filter %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   dplyr::rename(#ADM = ENSONIG00000036130,
#     CHGA = ENSONIG00000000725,
#     #GHRH = ENSONIG00000020252,
#     GNRH1 = ENSONIG00000011023,
#     #KNG1 = ENSONIG00000008378,
#     NMB = ENSONIG00000002537,
#     #NPFF = ENSONIG00000005229,
#     SST = ENSONIG00000033642) %>% 
#   t()
# 
# ### get list of sample ids
# burtoni.snseq.combined.sct.neurons.np.expression.samples = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
#   rownames_to_column('Cell.id') %>%  
#   full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#               rownames_to_column("Cell.id") %>% 
#               dplyr::select(c(orig.ident,
#                               Cell.id))) %>% 
#   relocate(orig.ident, .after = Cell.id)
# 
# ### get list of genotype ids
# burtoni.snseq.combined.sct.neurons.np.expression.genotypes = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
#   rownames_to_column('Cell.id') %>%  
#   full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#               rownames_to_column("Cell.id") %>% 
#               dplyr::select(c(Genotype.id,
#                               Cell.id))) %>% 
#   relocate(Genotype.id, .after = Cell.id)
# 
# #### convert to presence/absence matrix 
# # create new matrix
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj = burtoni.snseq.combined.sct.neurons.np.expression.filter
# # set all non-zero values to zero
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj > 0] <- 1 
# 
# # # change threshold to 4
# # burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj < 1.37] <- 0
# # burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj > 0] <- 1 
# 
# # count cells and genes
# ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj)
# # 11523 cells out of 12641
# nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj)
# # 25 genes out of 60
# 
# ### extract overlap of cell types for each genotype
# # make genotype ID list
# Genotype.id.list = burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#   pull(Genotype.id) %>% 
#   unique()
# 
# # create expression matrix dataframe
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df = t(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj) %>% 
#   as.data.frame()
# 
# 
# # create dummy dataframe
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = data.frame(from = as.character(),
#                                                                                     to = as.character(),
#                                                                                     weight = as.numeric(),
#                                                                                     Genotype.id = as.character())
# 
# ## run across genotype IDs
# for (i in Genotype.id.list) {
#   #create id list
#   tmp.names = burtoni.snseq.combined.sct.neurons.np.expression.genotypes %>% 
#     filter(Genotype.id == i) %>% 
#     pull(Cell.id)
#   # number of cells 
#   tmp.names.length = length(tmp.names)
#   #filter down cells
#   tmp = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df %>% 
#     filter(rownames(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df) %in% tmp.names) 
#   
#   ### create overlap square matrix of genes 
#   # calculate the overlap 
#   tmp.mat <- crossprod(as.matrix(tmp)) 
# 
#   ## create graph object
#   # undirected 
#   tmp.mat.g <- graph.adjacency(
#     tmp.mat,
#     mode="undirected",
#     weighted=TRUE,
#     diag=FALSE
#   )
#   
#   # Simplfy the adjacency object
#   ### needed?
#   tmp.mat.g <- igraph::simplify(tmp.mat.g, 
#                                 remove.multiple=TRUE, 
#                                 remove.loops=TRUE)
#   
#   ## extract edge list
#   burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
#     full_join(
#       get.data.frame(tmp.mat.g) %>% 
#     mutate(Genotype.id = i,
#            cell.count = tmp.names.length,
#            weight.scaled = 100*weight/cell.count)
#     )
# }
# 
# # add origin identity to genotype
# # add comparison ID
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
#   full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#               select(Genotype.id,
#                      orig.ident) %>% 
#               distinct()) %>% 
#   mutate(From.To = paste(from,
#                          to,
#                          sep = '.'))
# 
# ## create null dataframe
# # each genotype with each gene overlap
# weights.null = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
#   expand(Genotype.id,
#          From.To) %>% 
#   full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#               select(Genotype.id,
#                      orig.ident) %>% 
#               distinct())
# 
# # add null data
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
#   full_join(weights.null) %>% 
#   mutate(weight.scaled = ifelse(is.na(weight.scaled),
#                                 0,
#                                 weight.scaled)) 
# 
# ## run t.test across each neuropeptide pair
# burtoni.snseq.neuropeptide.network.sig = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
#   group_by(From.To) %>%
#   summarise(p.value = t.test(weight.scaled ~ orig.ident, paired = TRUE)$p.value) %>%
#   ungroup()
# 
# # FDR correction
# burtoni.snseq.neuropeptide.network.sig= burtoni.snseq.neuropeptide.network.sig %>% 
#   mutate(p.value.adj = p.adjust(p.value,
#                                 method = 'fdr'))
# # Create vertices for network 
# burtoni.snseq.neuropeptide.network.sig = burtoni.snseq.neuropeptide.network.sig %>% 
#   separate_wider_delim(From.To,
#                        delim = '.',
#                        names = c('From',
#                                  'To')) %>% 
#   select(-c(p.value))
# 
# ## calculate difference dataframe
# burtoni.snseq.neuropeptide.network.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
#   group_by(From.To,
#            orig.ident) %>%
#   summarise(avg.weight.scaled = mean(weight.scaled)) %>%
#   ungroup() %>% 
#   mutate(avg.weight.scaled = ifelse(orig.ident == 'sub_burtoni_snseq',
#                                     -1*avg.weight.scaled,
#                                     avg.weight.scaled)) %>% 
#   group_by(From.To) %>% 
#   summarise(avg.weight.scaled.diff = sum(avg.weight.scaled)) %>% 
#   ungroup()
# # Create vertices for network 
# burtoni.snseq.neuropeptide.network.weight = burtoni.snseq.neuropeptide.network.weight %>% 
#   separate_wider_delim(From.To,
#                        delim = '.',
#                        names = c('From',
#                                  'To')) 
# # combine dataframes
# burtoni.snseq.neuropeptide.network = burtoni.snseq.neuropeptide.network.weight %>% 
#   full_join(burtoni.snseq.neuropeptide.network.sig)
# 
# # create igraph
# burtoni.snseq.neuropeptide.network.g <- graph_from_data_frame(d = burtoni.snseq.neuropeptide.network,
#                                 directed = FALSE)
# 
# # # Remove edges below adjusted pvalue 0.05
# burtoni.snseq.neuropeptide.network.g <- delete_edges(burtoni.snseq.neuropeptide.network.g,
#                                                      E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$p.value.adj>0.05)])
# 
# # create adjacency matrix
# burtoni.snseq.neuropeptide.network.g.mat = as.matrix(burtoni.snseq.neuropeptide.network.g, 
#                 matrix.type = "adjacency",
#                 sparse = FALSE,
#                 attr = "avg.weight.scaled.diff")
# 
# 
# # Colour negative correlation edges as blue
# E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff<0)]$color <- "Sub.biased"
# 
# # Colour positive correlation edges as red
# E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff>0)]$color <- "Dom.biased"
# 
# # Convert edge weights to absolute values
# E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff <- abs(E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff)
# 
# # # Remove edges below adjusted pvalue 0.05
# burtoni.snseq.neuropeptide.network.g <- delete_edges(burtoni.snseq.neuropeptide.network.g,
#                        E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$p.value.adj>0.05)])
# 
# # get list of most expressed genes
# scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
# vSizes.all.gene <- (scale01(apply(burtoni.snseq.combined.sct.neurons.np.expression.filter, 1, mean )) + 1.0) 
# 
# # add scale data names to igraph 
# burtoni.snseq.neuropeptide.network.g <- permute(burtoni.snseq.neuropeptide.network.g, 
#                    match(V(burtoni.snseq.neuropeptide.network.g)$name, 
#                                   names(sort(vSizes.all.gene))))
# # create names variable
# V(burtoni.snseq.neuropeptide.network.g)$node_label <- names(V(burtoni.snseq.neuropeptide.network.g))
# # create size variable
# V(burtoni.snseq.neuropeptide.network.g)$node_size = sort(vSizes.all.gene)
# 
# ## graph
# burtoni.snseq.neuropeptide.network.g %>%
#   ggraph(layout = 'linear',
#          circular = TRUE
#   ) +
#   geom_edge_link(aes(width = avg.weight.scaled.diff ,
#                          alpha = -(p.value.adj),
#                          edge_colour = color)) +
#   geom_node_label(aes(label = node_label,
#                       size = node_size
#   )) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5),
#         title =element_text(size=20, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank(),
#         # legend.position = 'none'
#         ) +
#   ggtitle('Dominant vs Subordinate neuropeptide co-occurance') +
#   scale_edge_color_manual(values = c('red','blue')) 
#   # scale_edge_width(range = c(.993,10))+
#   # scale_edge_alpha(range = c(0,1)) 
# ggsave('neuropeptides/comparison/network/neuropeptides significant presence network.png',
#        height = 10,
#        width = 10)
# 
# # paper
# burtoni.snseq.neuropeptide.network.g %>%
#   ggraph(layout = 'linear',
#          circular = TRUE
#   ) +
#   geom_edge_link(aes(width = avg.weight.scaled.diff ,
#                      alpha = -(p.value.adj),
#                      edge_colour = color)) +
#   geom_node_label(aes(label = node_label,
#                       size = node_size
#   )) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5),
#         title =element_text(size=5, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank(),
#         legend.position = 'none'
#   ) +
#   # ggtitle('Dominant vs Subordinate neuropeptide co-occurance') +
#   scale_edge_color_manual(values = c('#4e499e','#60bb46'))
# # scale_edge_width(range = c(.993,10))+
# # scale_edge_alpha(range = c(0,1)) 
# ggsave('neuropeptides/comparison/network/neuropeptides significant presence network presentation.png',
#        height = 3,
#        width = 3)
# 
# ### graph heatmap
# # setup matrix for graphing
# # order by number of cells
# burtoni.snseq.neuropeptide.network.g.mat = burtoni.snseq.neuropeptide.network.g.mat[names(sort(-vSizes.all.gene)), names(sort(-vSizes.all.gene))]
# # set diagnol to NA
# diag(burtoni.snseq.neuropeptide.network.g.mat) <- NA
# 
# ## create pheatmap
# library(pheatmap)
# 
# burtoni.snseq.neuropeptide.network.g.mat %>% 
#   pheatmap(na_col = 'grey',
#            cluster_cols = F,
#            cluster_rows = F,
#            breaks = c(-.3,-0.01,0.01,5,10,20,80),
#            color = c('blue',
#                      'white',
#                      'red',
#                      'red1',
#                      'red2',
#                      'red4'),
#            angle_col = 45,
#            filename = 'neuropeptides/comparison/network/neuropeptides significant presence network heatmap.png',
#            height = 10,
#            width = 10)
# 
# 
# ## calculate global transitivity 
# # 0.4918
# burtoni.snseq.neuropeptide.network.g %>%
#   transitivity(type = 'global')
# 
# # ###graph dom and sub average networks
# # ##dom
# # burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.dom = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
# #   filter(orig.ident == 'dom_burtoni_snseq') %>% 
# #   group_by(from,
# #            to,
# #            From.To,
# #            orig.ident) %>%
# #   summarise(avg.weight.scaled = mean(weight.scaled)) %>%
# #   ungroup()
# # 
# # ##sub
# # burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.sub = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
# #   filter(orig.ident == 'sub_burtoni_snseq') %>% 
# #   group_by(from,
# #            to,
# #            From.To,
# #            orig.ident) %>%
# #   summarise(avg.weight.scaled = mean(weight.scaled)) %>%
# #   ungroup()
# 
# 
# #### clustering coefficient/ transitivity ####
# ### create networks for each genotype
# # Create vertices for network 
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
#   separate_wider_delim(From.To,
#                        delim = '.',
#                        names = c('From',
#                                  'To')) %>% 
#   dplyr::select(c(From,
#                   To,
#                   weight.scaled,
#                   orig.ident,
#                   Genotype.id))
# 
# # graph distribution
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>% 
#   ggplot(aes(weight.scaled,
#              fill = orig.ident)) +
#   geom_histogram(binwidth = 1) +
#   theme_classic() +
#   facet_grid(Genotype.id ~ .)
# ggsave('neuropeptides/comparison/network/Histogram genotypes weight scaled.png',
#        height = 5.25,
#        width = 5.25)
# 
# 
# #### remove zero weights?
# # burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>%
# #   mutate(weight.scaled = ifelse(weight.scaled < 1,
# #                                 0,
# #                                 weight.scaled))
# # 
# # burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>%
# #   filter(weight.scaled != 0)
# 
# ### calculate local and global transitivity across samples
# ## create dummy dataframe
# transitivty.df = data.frame(Genotype.id = as.character(),
#                                    transitivity = as.numeric(),
#                                    vertice = as.character())
#   
# 
# ## run igraph and calculate transitivity for each sample
# # sample list
# # Genotype.id.list = unique(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all$Genotype.id)
# 
# # i = Genotype.id.list[4]
# for (i in Genotype.id.list) {
#   ## create igraph
#   tmp <- graph_from_data_frame(d = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>% 
#                                  filter(Genotype.id == i),
#                                                                 directed = FALSE)
#   
#   # # Remove edges with less than 1% 
#   tmp <- delete_edges(tmp,
#                       E(tmp)[which(E(tmp)$weight.scaled<1)])
# 
#   
#   ## calculate global transitivity 
#   tmp.global = tmp %>% 
#     transitivity(type = 'global')
#   
#   # add to dataframe 
#   transitivty.df = bind_rows(transitivty.df,
#               data.frame(Genotype.id = i,
#                      transitivity = tmp.global,
#                      vertice = 'global'))
#   
#   ## calculate local transitivity 
#   tmp.local = tmp %>%
#     transitivity(type = 'local')
# 
#   # tmp.local = tmp %>%
#   #   transitivity(type = 'weighted')
#   
#   # add to dataframe 
#   transitivty.df = bind_rows(transitivty.df,
#                                    data.frame(Genotype.id = i,
#                                               transitivity = tmp.local,
#                                               vertice = V(tmp)$name))
# }
# 
# ## set NA to zero
# transitivty.df  = transitivty.df %>% 
#   mutate(transitivity = ifelse(is.na(transitivity),
#                                0,
#                                transitivity),
#          color = ifelse(vertice == 'global',
#                         'global',
#                         'local'))
# 
# ## add orig.ident
# transitivty.df = transitivty.df %>% 
#   full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
#               select(Genotype.id,
#                      orig.ident) %>% 
#               distinct()) 
# 
# ## histogram
# # transitivty.df %>% 
# #   ggplot(aes(x = transitivity,
# #              fill  = orig.ident,
# #              color = Genotype.id)) +
# #   geom_histogram() +
# #   facet_grid(.~color) +
# #   theme_classic() 
# 
# ## convert dataframe to wide format for graphing
# # average across status
# transitivty.df.wide = transitivty.df %>% 
#   group_by(orig.ident,
#            vertice) %>% 
#   summarise(transitivity.avg = mean(transitivity)) %>%
#   ungroup() %>% 
#   pivot_wider(id_cols = 'vertice',
#               names_from = 'orig.ident',
#               values_from = 'transitivity.avg') %>% 
#   mutate(color = ifelse(vertice == 'global',
#                         'global',
#                         'local')) %>% 
#   mutate(dom_burtoni_snseq = ifelse(is.na(dom_burtoni_snseq),
#                                     0,
#                                     dom_burtoni_snseq),
#          sub_burtoni_snseq = ifelse(is.na(sub_burtoni_snseq),
#                                     0,
#                                     sub_burtoni_snseq))
# 
# 
# ## graph
# transitivty.df.wide %>% 
#   ggplot(aes(x = dom_burtoni_snseq,
#              y = sub_burtoni_snseq,
#              label = vertice,
#              color = color)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_point() +
#   geom_label_repel(max.overlaps = 20) +
#   theme_classic() +
#   xlim(0,1) +
#   ylim(0,1) +
#   ggtitle('transitivity') +
#   scale_color_manual(values = c('red', 'black')) +
#   theme(legend.position = 'null')
# ggsave('neuropeptides/comparison/network/Transitivity genotypes average.png',
#        height = 5.25,
#        width = 5.25)
# 
# 
# 
# #### testing covariance matrices ####
# ### create covariance matrix per sample
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights
# 
# library(HDtest)
# library(parallel)
# 
# ### run hd test
# males.hd = testCov(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, 
#                    burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, 
#                    method = "HD",
#                    J = 1000, 
#                    alpha = 0.05, 
#                    n.core = 12)
# 
# ## run clx test
# males.CLX = testCov(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, 
#                     burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, 
#                     method = "CLX",
#                     J = 1000, 
#                     alpha = 0.05, 
#                     n.core = 12)
# 
# 
# 
# #### create dataframe
# results.testcov = data.frame(comparison = c('Male'),
#                              HD.p.val = c(males.hd[[2]]),
#                              HD.stat = c(males.hd[[1]]),
#                              CLX.p.val = c(males.CLX[[2]]),
#                              CLX.stat = c(males.CLX[[1]])
# )
# 
# results.testcov
# # comparison HD.p.val  HD.stat CLX.p.val CLX.stat
# # Statistic       Male        0 54.56802         0 2977.669
# 
# # results.testcov.short = data.frame(comparison = c('Male'),
# #                                    HD.stat = round(c(males.hd[[1]])),
# #                                    CLX.stat = round(c(males.CLX[[1]]))
# # )
# 
# ### mantel test 
# mantel.np.ci = as.lpcor(as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.dom, method="pearson"))),
#                         as.matrix(as.dist(cor(burtoni.snseq.combined.sct.neurons.np.expression.filter.male.sub, method="pearson")))) %>%
#   pairs_mantel.ci()
# 
# # rename labels
# mantel.np.ci[["xAxisLabels"]]  = c("Male.Dom",
#                                    "Male.Sub")
# mantel.np.ci[["yAxisLabels"]]  = c("Male.Dom",
#                                    "Male.Sub")
# mantel.np.ci[["title"]] = paste0( mantel.np.ci[["title"]],
#                                   ' filter 50 ci')
# mantel.np.ci
# ggsave('neuropeptides/comparison/network/Mantel test.png',
#        height = 5.25,
#        width = 5.25)
# 
# ## Mantel test: P-value = 0.003, r = 0.18, 95% CI for r = 0.07-0.25 
# 
# 
# #### compare GnRH1, prl, and pomc ####
# ### use feature plot
# burtoni.snseq.combined.sct.all.neurons %>% 
#   FeaturePlot(features = c('ENSONIG00000011023',
#                            'prl',
#                            'pomc'),
#               max.cutoff = 1,
#               pt.size	= 1,
#               order= T) 
# ggsave('neuropeptides/comparison/UMAP Gnrh1 prl pomc.png',
#        height = 20,
#        width = 20)
# 
# burtoni.snseq.combined.sct.all.neurons %>% 
#   FeaturePlot(features = c('ENSONIG00000033642',
#                            'cort'),
#               max.cutoff = 1,
#               pt.size	= 1,
#               order= T) 
# ggsave('neuropeptides/comparison/UMAP SST cort.png',
#        height = 10,
#        width = 20)
# 
# 
# 
# 
# 
# 
# 

# #### neuropeptides comparison across social status ####
# ## set idents
# Idents(object = burtoni.snseq.combined.sct.all.neurons) = "orig.ident"
# ### create matrix of neuropeptide expression level 
# ## use neuropeptides.list
# ## reduce columns to neuropeptides
# ## seperate into dom and sub
# # create presence absence matrix
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = GetAssayData(subset(burtoni.snseq.combined.sct.all.neurons,
#                                                                                                            idents = 'dom_burtoni_snseq'),
#                                                                                                     assay = 'SCT',
#                                                                                                     slot = 'data') %>% 
#   as.data.frame() %>% 
#   filter(rownames(burtoni.snseq.combined.sct.all.neurons@assays$SCT@data) %in% neuropeptides.genes) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   select(any_of(neuropeptides.genes)) %>% 
#   mutate(across(where(is.numeric), 
#                 function(x) ifelse(x < 1, 0, x))) %>% 
#   mutate(across(where(is.numeric), 
#                 function(x) ifelse(x >= 1, 1, x))) 
# 
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = GetAssayData(subset(burtoni.snseq.combined.sct.all.neurons,
#                                                                                                            idents = 'sub_burtoni_snseq'),
#                                                                                                     assay = 'SCT',
#                                                                                                     slot = 'data') %>% 
#   as.data.frame() %>% 
#   filter(rownames(burtoni.snseq.combined.sct.all.neurons@assays$SCT@data) %in% neuropeptides.genes) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   select(any_of(neuropeptides.genes)) %>% 
#   mutate(across(where(is.numeric), 
#                 function(x) ifelse(x < 1, 0, x))) %>% 
#   mutate(across(where(is.numeric), 
#                 function(x) ifelse(x >= 1, 1, x))) 
# 
# ##rename to actual gene name
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
#   dplyr::rename(#ADM = ENSONIG00000036130,
#     CHGA = ENSONIG00000000725,
#     GHRH = ENSONIG00000020252,
#     GNRH1 = ENSONIG00000011023,
#     #KNG1 = ENSONIG00000008378,
#     NMB = ENSONIG00000002537,
#     NPFF = ENSONIG00000005229,
#     SST = ENSONIG00000033642)
# 
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
#   dplyr::rename(#ADM = ENSONIG00000036130,
#     CHGA = ENSONIG00000000725,
#     GHRH = ENSONIG00000020252,
#     GNRH1 = ENSONIG00000011023,
#     #KNG1 = ENSONIG00000008378,
#     NMB = ENSONIG00000002537,
#     NPFF = ENSONIG00000005229,
#     SST = ENSONIG00000033642)
# 
# #how many neurons have neuropeptides
# #dom has 6454 neurons and 6061 expressing neuropeptide
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   filter(Total.neuropeptides != 0) %>% 
#   nrow()
# # total neurons
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
#   nrow()
# 
# 
# #sub has 6187 neurons and 3608 expressing neuropeptide
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   filter(Total.neuropeptides != 0) %>% 
#   nrow()
# 
# # total neurons
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
#   nrow()
# 
# 
# #graph together
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>%
#   select(Total.neuropeptides) %>%
#   mutate(Status = 'dom') %>%
#   rbind(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#           mutate(Total.neuropeptides = rowSums(.)) %>%
#           select(Total.neuropeptides) %>%
#           mutate(Status = 'sub')) %>%
#   mutate(count = 1) %>% 
#   group_by(Status,
#            Total.neuropeptides) %>%
#   summarise(count.total = sum(count)) %>% 
#   group_by(Status) %>%
#   mutate(Total = sum(count.total)) %>% 
#   ungroup() %>%
#   mutate(Percent = 100*count.total/Total) %>%
#   ggplot(aes(y = Percent,
#              x = Total.neuropeptides,
#              group = Status,
#              color = Status)) +
#   geom_line(size = 4) +
#   geom_point(size = 10) +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(0,11,1)) +
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         title =element_text(size=20, face='bold')) +
#   scale_color_manual(values = c("#4e499e",
#                                 "#60bb46"))+
#   xlab('Neuropeptides per neuron')
# ggsave('neuropeptides/comparison_new/social.status/Neuropeptides per neuron across status.png',
#        width = 10,
#        height = 10)
# 
# 

#### neuropeptides across social status ####
### extract neuropeptide expression and orig.idents from neurons
## create expression dataframe of neuropeptides
burtoni.snseq.combined.sct.neurons.np.expression = burtoni.snseq.combined.sct.all.neurons@assays$SCT@data %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  # mutate(gene = toupper(gene)) %>%
  filter(gene %in% neuropeptides.genes) %>% 
  column_to_rownames('gene') %>% 
  t() %>% 
  as.data.frame()


# remove low expressing genes
# remove low expressing cells
burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression #%>% 
  # remove_nonexp(method = "percentage",
  #               min_exp = 0.68,
  #               min_percentage_samples = .05)  %>% 
  # t() %>% 
  # remove_nonexp(method = "percentage",
  #               min_exp = 0.68,
  #               min_percentage_samples = .05) 

# count cells and genes
ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# 11523 cells out of 12641
nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# 25 genes out of 60

###rename to actual gene name
burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression.filter %>% 
  # t() %>% 
  as.data.frame() %>% 
  dplyr::rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642) %>% 
  t()

### get list of sample ids
burtoni.snseq.combined.sct.neurons.np.expression.samples = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
  rownames_to_column('Cell.id') %>%  
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              rownames_to_column("Cell.id") %>% 
              dplyr::select(c(orig.ident,
                              Cell.id))) %>% 
  relocate(orig.ident, .after = Cell.id)

### get list of genotype ids
burtoni.snseq.combined.sct.neurons.np.expression.genotypes = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
  rownames_to_column('Cell.id') %>%  
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              rownames_to_column("Cell.id") %>% 
              dplyr::select(c(Genotype.id,
                              Cell.id))) %>% 
  relocate(Genotype.id, .after = Cell.id)

#### convert to presence/absence matrix 
# create new matrix
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj = burtoni.snseq.combined.sct.neurons.np.expression.filter
# set all non-zero values to one
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj > 0] <- 1

### extract overlap of cell types for each genotype
# make genotype ID list
Genotype.id.list = burtoni.snseq.combined.sct.all.neurons@meta.data %>%
  pull(Genotype.id) %>% 
  unique()

# create expression matrix dataframe
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df = t(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj) %>% 
  as.data.frame()


# create dummy dataframe
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total = data.frame(Total.neuropeptides = as.numeric(),
                                                                                         Genotype.id = as.character())

## run across genotype IDs
for (i in Genotype.id.list) {
  #create id list
  tmp.names = burtoni.snseq.combined.sct.neurons.np.expression.genotypes %>% 
    filter(Genotype.id == i) %>% 
    pull(Cell.id)
  # number of cells 
  tmp.names.length = length(tmp.names)
  #filter down cells
  tmp = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df %>% 
    filter(rownames(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df) %in% tmp.names) 
  
  ## calculate total count
  tmp = tmp %>%
    mutate(Total.neuropeptides = rowSums(.)) %>%
    dplyr::select(Total.neuropeptides)  %>% 
    mutate(Genotype.id = i)
 
  ## extract edge list
  burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total %>% 
    full_join(tmp)
}

# summarise
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total %>%
  mutate(count = 1) %>%
  group_by(Genotype.id,
           Total.neuropeptides) %>%
  summarise(count.total = sum(count)) %>%
  group_by(Genotype.id) %>%
  mutate(Total = sum(count.total)) %>%
  ungroup() %>%
  mutate(Percent = 100*count.total/Total)

# add orig.ident
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count %>% 
  separate_wider_delim(Genotype.id,
                              delim = '.',
                              names = c('Status',
                                        'Genotype'),
                              cols_remove = F)


### stats
## compare distributions
# dom
sink('neuropeptides/comparison_new/social.status/dom.ad.txt')
kSamples::ad.test(count.total ~ Genotype,
                           data = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count %>% 
                             filter(Status == 'Dom'))
sink()
# sub
sink('neuropeptides/comparison_new/social.status/sub.ad.txt')
kSamples::ad.test(count.total ~ Genotype,
                           data = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count %>% 
                             filter(Status == 'Sub'))
sink()
# all
sink('neuropeptides/comparison_new/social.status/all.ad.txt')
kSamples::ad.test(count.total ~ Genotype.id,
                           data = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count)
sink()

## compare number of neuropeptides
# by status
sink('neuropeptides/comparison_new/social.status/status.glm.txt')
summary(glm(Total.neuropeptides ~ Status, 
            family=poisson,
            data = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total%>% 
              separate_wider_delim(Genotype.id,
                                   delim = '.',
                                   names = c('Status',
                                             'Genotype'),
                                   cols_remove = F)))
sink()

# within doms
sink('neuropeptides/comparison_new/social.status/dom.glm.txt')
summary(glm(Total.neuropeptides ~ Genotype.id, 
            family=poisson,
            data = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total%>% 
              separate_wider_delim(Genotype.id,
                                   delim = '.',
                                   names = c('Status',
                                             'Genotype'),
                                   cols_remove = F) %>% 
              filter(Status == 'Dom')))
sink()

# within sub
sink('neuropeptides/comparison_new/social.status/sub.glm.txt')
summary(glm(Total.neuropeptides ~ Genotype.id, 
            family=poisson,
            data = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total%>% 
              separate_wider_delim(Genotype.id,
                                   delim = '.',
                                   names = c('Status',
                                             'Genotype'),
                                   cols_remove = F) %>% 
              filter(Status == 'Sub')))
sink()



### graph together
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.status.total.count %>% 
  group_by(Status,
           Total.neuropeptides) %>% 
  mutate(Percent.median = median(Percent)) %>% 
  ggplot() +
  geom_line(aes(y = Percent,
                x = Total.neuropeptides,
                group = Genotype.id,
                color = Status),
            size = 1) +
  geom_point(aes(y = Percent,
                 x = Total.neuropeptides,
                 group = Status,
                 fill = Status),
             size = 5,
             shape = 21) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,15,1)) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title =element_text(size=20, face='bold'),
        legend.position = "none") +
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))+
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  xlab('Neuropeptides per neuron')
ggsave('neuropeptides/comparison_new/social.status/Neuropeptides per neuron across status error.png',
       width = 10,
       height = 10)


#### presence/absence matrix ####
### extract neuropeptide expression and orig.idents from neurons
## create expression dataframe of neuropeptides
burtoni.snseq.combined.sct.neurons.np.expression = burtoni.snseq.combined.sct.all.neurons@assays$SCT@data %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  # mutate(gene = toupper(gene)) %>%
  filter(gene %in% neuropeptides.genes) %>% 
  column_to_rownames('gene') %>% 
  t() %>% 
  as.data.frame()


# remove low expressing genes
# remove low expressing cells
burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression %>% 
  remove_nonexp(method = "percentage",
                min_exp = 0.68,
                min_percentage_samples = .05)  %>% 
  t() %>% 
  remove_nonexp(method = "percentage",
                min_exp = 0.68,
                min_percentage_samples = .05) 

# count cells and genes
ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# 11523 cells out of 12641
nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# 25 genes out of 60

###rename to actual gene name
burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression.filter %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    #GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    #NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642) %>% 
  t()

### get list of sample ids
burtoni.snseq.combined.sct.neurons.np.expression.samples = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
  rownames_to_column('Cell.id') %>%  
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              rownames_to_column("Cell.id") %>% 
              dplyr::select(c(orig.ident,
                              Cell.id))) %>% 
  relocate(orig.ident, .after = Cell.id)

### get list of genotype ids
burtoni.snseq.combined.sct.neurons.np.expression.genotypes = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
  rownames_to_column('Cell.id') %>%  
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              rownames_to_column("Cell.id") %>% 
              dplyr::select(c(Genotype.id,
                              Cell.id))) %>% 
  relocate(Genotype.id, .after = Cell.id)

#### convert to presence/absence matrix 
# create new matrix
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj = burtoni.snseq.combined.sct.neurons.np.expression.filter
# set all non-zero values to zero
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj > 0] <- 1 

# # change threshold to 4
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj < 1.37] <- 0
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj > 0] <- 1 

# count cells and genes
ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj)
# 11523 cells out of 12641
nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj)
# 25 genes out of 60

### extract overlap of cell types for each genotype
# make genotype ID list
Genotype.id.list = burtoni.snseq.combined.sct.all.neurons@meta.data %>%
  pull(Genotype.id) %>% 
  unique()

# create expression matrix dataframe
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df = t(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj) %>% 
  as.data.frame()


# create dummy dataframe
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = data.frame(from = as.character(),
                                                                                    to = as.character(),
                                                                                    weight = as.numeric(),
                                                                                    Genotype.id = as.character())

## run across genotype IDs
for (i in Genotype.id.list) {
  #create id list
  tmp.names = burtoni.snseq.combined.sct.neurons.np.expression.genotypes %>% 
    filter(Genotype.id == i) %>% 
    pull(Cell.id)
  # number of cells 
  tmp.names.length = length(tmp.names)
  #filter down cells
  tmp = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df %>% 
    filter(rownames(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df) %in% tmp.names) 
  
  ### create overlap square matrix of genes 
  # calculate the overlap 
  tmp.mat <- crossprod(as.matrix(tmp)) 
  
  ## create graph object
  # undirected 
  tmp.mat.g <- graph.adjacency(
    tmp.mat,
    mode="undirected",
    weighted=TRUE,
    diag=FALSE
  )
  
  # Simplfy the adjacency object
  ### needed?
  tmp.mat.g <- igraph::simplify(tmp.mat.g, 
                                remove.multiple=TRUE, 
                                remove.loops=TRUE)
  
  ## extract edge list
  burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
    full_join(
      get.data.frame(tmp.mat.g) %>% 
        mutate(Genotype.id = i,
               cell.count = tmp.names.length,
               weight.scaled = 100*weight/cell.count)
    )
}

# add origin identity to genotype
# add comparison ID
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              select(Genotype.id,
                     orig.ident) %>% 
              distinct()) %>% 
  mutate(From.To = paste(from,
                         to,
                         sep = '.'))

## create null dataframe
# each genotype with each gene overlap
weights.null = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
  expand(Genotype.id,
         From.To) %>% 
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              select(Genotype.id,
                     orig.ident) %>% 
              distinct())

# add null data
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
  full_join(weights.null) %>% 
  mutate(weight.scaled = ifelse(is.na(weight.scaled),
                                0,
                                weight.scaled)) 

## run t.test across each neuropeptide pair
burtoni.snseq.neuropeptide.network.sig = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
  group_by(From.To) %>%
  summarise(p.value = t.test(weight.scaled ~ orig.ident, paired = TRUE)$p.value) %>%
  ungroup()

# FDR correction
burtoni.snseq.neuropeptide.network.sig= burtoni.snseq.neuropeptide.network.sig %>% 
  mutate(p.value.adj = p.adjust(p.value,
                                method = 'fdr'))
# Create vertices for network 
burtoni.snseq.neuropeptide.network.sig = burtoni.snseq.neuropeptide.network.sig %>% 
  separate_wider_delim(From.To,
                       delim = '.',
                       names = c('From',
                                 'To')) %>% 
  select(-c(p.value))

## calculate difference dataframe
burtoni.snseq.neuropeptide.network.weight = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
  group_by(From.To,
           orig.ident) %>%
  summarise(avg.weight.scaled = mean(weight.scaled)) %>%
  ungroup() %>% 
  mutate(avg.weight.scaled = ifelse(orig.ident == 'sub_burtoni_snseq',
                                    -1*avg.weight.scaled,
                                    avg.weight.scaled)) %>% 
  group_by(From.To) %>% 
  summarise(avg.weight.scaled.diff = sum(avg.weight.scaled)) %>% 
  ungroup()
# Create vertices for network 
burtoni.snseq.neuropeptide.network.weight = burtoni.snseq.neuropeptide.network.weight %>% 
  separate_wider_delim(From.To,
                       delim = '.',
                       names = c('From',
                                 'To')) 
# combine dataframes
burtoni.snseq.neuropeptide.network = burtoni.snseq.neuropeptide.network.weight %>% 
  full_join(burtoni.snseq.neuropeptide.network.sig)

# create igraph
burtoni.snseq.neuropeptide.network.g <- graph_from_data_frame(d = burtoni.snseq.neuropeptide.network,
                                                              directed = FALSE)

# # Remove edges below adjusted pvalue 0.05
burtoni.snseq.neuropeptide.network.g <- delete_edges(burtoni.snseq.neuropeptide.network.g,
                                                     E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$p.value.adj>0.05)])

# create adjacency matrix
burtoni.snseq.neuropeptide.network.g.mat = as.matrix(burtoni.snseq.neuropeptide.network.g, 
                                                     matrix.type = "adjacency",
                                                     sparse = FALSE,
                                                     attr = "avg.weight.scaled.diff")


# Colour negative correlation edges as blue
E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff<0)]$color <- "Sub.biased"

# Colour positive correlation edges as red
E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff>0)]$color <- "Dom.biased"

# Convert edge weights to absolute values
E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff <- abs(E(burtoni.snseq.neuropeptide.network.g)$avg.weight.scaled.diff)

# # Remove edges below adjusted pvalue 0.05
burtoni.snseq.neuropeptide.network.g <- delete_edges(burtoni.snseq.neuropeptide.network.g,
                                                     E(burtoni.snseq.neuropeptide.network.g)[which(E(burtoni.snseq.neuropeptide.network.g)$p.value.adj>0.05)])

# get list of most expressed genes
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes.all.gene <- (scale01(apply(burtoni.snseq.combined.sct.neurons.np.expression.filter, 1, mean )) + 1.0) 

# add scale data names to igraph 
burtoni.snseq.neuropeptide.network.g <- permute(burtoni.snseq.neuropeptide.network.g, 
                                                match(V(burtoni.snseq.neuropeptide.network.g)$name, 
                                                      names(sort(vSizes.all.gene))))
# create names variable
V(burtoni.snseq.neuropeptide.network.g)$node_label <- names(V(burtoni.snseq.neuropeptide.network.g))
# create size variable
V(burtoni.snseq.neuropeptide.network.g)$node_size = sort(vSizes.all.gene)

# Put Sub edges on top in new variable
# Colour positive correlation edges as red
E(burtoni.snseq.neuropeptide.network.g)[color != "Sub.biased"]$avg.weight.scaled.diff2 <- NA
E(burtoni.snseq.neuropeptide.network.g)[color != "Dom.biased"]$avg.weight.scaled.diff2 <- E(burtoni.snseq.neuropeptide.network.g)[color != "Dom.biased"]$avg.weight.scaled.diff

## graph
# burtoni.snseq.neuropeptide.network.g %>%
#   ggraph(layout = 'linear',
#          circular = TRUE
#   ) +
#   geom_edge_link(aes(width = avg.weight.scaled.diff ,
#                      alpha = -(p.value.adj),
#                      edge_colour = color)) +
#   geom_node_label(aes(label = node_label,
#                       size = node_size
#   )) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5),
#         title =element_text(size=20, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank(),
#         # legend.position = 'none'
#   ) +
#   ggtitle('Dominant vs Subordinate neuropeptide co-occurance') +
#   scale_edge_color_manual(values = c('red','blue')) 
# # scale_edge_width(range = c(.993,10))+
# # scale_edge_alpha(range = c(0,1)) 
# ggsave('neuropeptides/comparison_new//network/neuropeptides significant presence network.png',
#        height = 10,
#        width = 10)
# 
# # presentation
# burtoni.snseq.neuropeptide.network.g %>%
#   ggraph(layout = 'linear',
#          circular = TRUE
#   ) +
#   geom_edge_link(aes(width = avg.weight.scaled.diff ,
#                      alpha = -(p.value.adj),
#                      edge_colour = color)) +
#   geom_node_label(aes(label = node_label,
#                       size = node_size
#   )) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5),
#         title =element_text(size=5, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank(),
#         legend.position = 'none'
#   ) +
#   # ggtitle('Dominant vs Subordinate neuropeptide co-occurance') +
#   scale_edge_color_manual(values = c('#4e499e','#60bb46'))
# # scale_edge_width(range = c(.993,10))+
# # scale_edge_alpha(range = c(0,1)) 
# ggsave('neuropeptides/comparison_new/network/neuropeptides significant presence network presentation.png',
#        height = 3,
#        width = 3)

# paper
p = burtoni.snseq.neuropeptide.network.g %>%
  ggraph(layout = 'linear',
         circular = TRUE
  ) +
  geom_edge_link(aes(width = avg.weight.scaled.diff ,
                     # alpha = -log(p.value.adj),
                     edge_colour = color
                     )) +
  geom_edge_link(aes(width = avg.weight.scaled.diff2 ,
                     # alpha = -log(p.value.adj),
                     edge_colour = color
  )) +
  geom_node_point(aes(size = node_size),
                  fill = 'white',
                  shape = 21) +
  # geom_node_label(aes(label = node_label),
  #                 repel = T) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=4, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank(),
        legend.position = 'none'
  ) +
  # ggtitle('Dominant vs Subordinate neuropeptide co-occurance') +
  scale_edge_color_manual(values = c('#4e499e','#60bb46'))
# scale_edge_width(range = c(.993,10))+
# scale_edge_alpha(range = c(0,1)) 

p + 
  geom_node_label(aes(label = node_label),
                 nudge_x = p$data$x * .1, 
                 nudge_y = p$data$y * .1)

ggsave('neuropeptides/comparison_new/network/neuropeptides significant presence network paper.png',
       height = 6,
       width = 6)

### graph heatmap
# setup matrix for graphing
# order by number of cells
burtoni.snseq.neuropeptide.network.g.mat = burtoni.snseq.neuropeptide.network.g.mat[names(sort(-vSizes.all.gene)), names(sort(-vSizes.all.gene))]
# set diagonal to NA
diag(burtoni.snseq.neuropeptide.network.g.mat) <- NA

## create pheatmap
library(pheatmap)

burtoni.snseq.neuropeptide.network.g.mat %>% 
  pheatmap(na_col = 'grey',
           cluster_cols = F,
           cluster_rows = F,
           breaks = c(-.3,-0.01,0.01,80),
           color = c('#60bb46',
                     'white',
                     '#4e499e'),
           angle_col = 45,
           legend = F,
           filename = 'neuropeptides/comparison_new/network/neuropeptides significant presence network heatmap paper.png',
           height = 6.5,
           width = 6.5)

burtoni.snseq.neuropeptide.network.g.mat %>% 
  pheatmap(na_col = 'grey',
           cluster_cols = F,
           cluster_rows = F,
           breaks = c(-.3,-0.01,0.01,5,10,20,80),
           color = c('blue',
                     'white',
                     'red',
                     'red1',
                     'red2',
                     'red4'),
           angle_col = 45,
           filename = 'neuropeptides/comparison_new/network/neuropeptides significant presence network heatmap.png',
           height = 10,
           width = 10)


## calculate global transitivity 
# 0.4918
burtoni.snseq.neuropeptide.network.g %>%
  transitivity(type = 'global')

# ###graph dom and sub average networks
# ##dom
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.dom = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
#   filter(orig.ident == 'dom_burtoni_snseq') %>% 
#   group_by(from,
#            to,
#            From.To,
#            orig.ident) %>%
#   summarise(avg.weight.scaled = mean(weight.scaled)) %>%
#   ungroup()
# 
# ##sub
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.sub = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>%
#   filter(orig.ident == 'sub_burtoni_snseq') %>% 
#   group_by(from,
#            to,
#            From.To,
#            orig.ident) %>%
#   summarise(avg.weight.scaled = mean(weight.scaled)) %>%
#   ungroup()


#### clustering coefficient/ transitivity ####
### create networks for each genotype
# Create vertices for network 
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights %>% 
  separate_wider_delim(From.To,
                       delim = '.',
                       names = c('From',
                                 'To')) %>% 
  dplyr::select(c(From,
                  To,
                  weight.scaled,
                  orig.ident,
                  Genotype.id))

# graph distribution
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>% 
  ggplot(aes(weight.scaled,
             fill = orig.ident)) +
  geom_histogram(binwidth = 1) +
  theme_classic() +
  facet_grid(Genotype.id ~ .)
ggsave('neuropeptides/comparison_new/network/Histogram genotypes weight scaled.png',
       height = 5.25,
       width = 5.25)


#### remove zero weights?
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>%
#   mutate(weight.scaled = ifelse(weight.scaled < 1,
#                                 0,
#                                 weight.scaled))
# 
# burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>%
#   filter(weight.scaled != 0)

### calculate local and global transitivity across samples
## create dummy dataframe
transitivty.df = data.frame(Genotype.id = as.character(),
                            transitivity = as.numeric(),
                            vertice = as.character())


## run igraph and calculate transitivity for each sample
# sample list
# Genotype.id.list = unique(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all$Genotype.id)

# i = Genotype.id.list[4]
for (i in Genotype.id.list) {
  ## create igraph
  tmp <- graph_from_data_frame(d = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df.weights.all %>% 
                                 filter(Genotype.id == i),
                               directed = FALSE)
  
  # # Remove edges with less than 1% 
  tmp <- delete_edges(tmp,
                      E(tmp)[which(E(tmp)$weight.scaled<1)])
  
  
  ## calculate global transitivity 
  tmp.global = tmp %>% 
    transitivity(type = 'global')
  
  # add to dataframe 
  transitivty.df = bind_rows(transitivty.df,
                             data.frame(Genotype.id = i,
                                        transitivity = tmp.global,
                                        vertice = 'global'))
  
  ## calculate local transitivity 
  tmp.local = tmp %>%
    transitivity(type = 'local')
  
  # tmp.local = tmp %>%
  #   transitivity(type = 'weighted')
  
  # add to dataframe 
  transitivty.df = bind_rows(transitivty.df,
                             data.frame(Genotype.id = i,
                                        transitivity = tmp.local,
                                        vertice = V(tmp)$name))
}

## set NA to zero
transitivty.df  = transitivty.df %>% 
  mutate(transitivity = ifelse(is.na(transitivity),
                               0,
                               transitivity),
         color = ifelse(vertice == 'global',
                        'global',
                        'local'))

## add orig.ident
transitivty.df = transitivty.df %>% 
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              select(Genotype.id,
                     orig.ident) %>% 
              distinct()) 

## histogram
# transitivty.df %>% 
#   ggplot(aes(x = transitivity,
#              fill  = orig.ident,
#              color = Genotype.id)) +
#   geom_histogram() +
#   facet_grid(.~color) +
#   theme_classic() 

## convert dataframe to wide format for graphing
# average across status
transitivty.df.wide = transitivty.df %>% 
  group_by(orig.ident,
           vertice) %>% 
  summarise(transitivity.avg = mean(transitivity)) %>%
  ungroup() %>% 
  pivot_wider(id_cols = 'vertice',
              names_from = 'orig.ident',
              values_from = 'transitivity.avg') %>% 
  mutate(color = ifelse(vertice == 'global',
                        'global',
                        'local')) %>% 
  mutate(dom_burtoni_snseq = ifelse(is.na(dom_burtoni_snseq),
                                    0,
                                    dom_burtoni_snseq),
         sub_burtoni_snseq = ifelse(is.na(sub_burtoni_snseq),
                                    0,
                                    sub_burtoni_snseq))


## graph
# transitivty.df.wide %>% 
#   ggplot(aes(x = dom_burtoni_snseq,
#              y = sub_burtoni_snseq,
#              label = vertice,
#              color = color)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_point() +
#   geom_label_repel(max.overlaps = 30) +
#   theme_classic() +
#   xlim(0,1) +
#   ylim(0,1) +
#   ggtitle('transitivity') +
#   scale_color_manual(values = c('red', 'black')) +
#   theme(legend.position = 'null')
# ggsave('neuropeptides/comparison_new/network/Transitivity genotypes average.png',
#        height = 5.25,
#        width = 5.25)

# paper
transitivty.df.wide %>% 
  ggplot(aes(x = dom_burtoni_snseq,
             y = sub_burtoni_snseq,
             label = vertice,
             color = color)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point() +
  geom_label_repel(max.overlaps = 30) +
  theme_classic() +
  xlim(0,1) +
  ylim(0,1) +
  # ggtitle('transitivity') +
  scale_color_manual(values = c('red', 'black')) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        title =element_text(size=10, face='bold'),
        legend.position = "none") +
  xlab('Dominant transitivity score')+
  ylab('Subordinate transitivity score')
ggsave('neuropeptides/comparison_new/network/Transitivity genotypes average paper.png',
       height = 5.25,
       width = 5.25)



#### neuropeptide count statistics ####
### compare number of neuropeptide expressing cells across samples
### GLM binomial
library(multcomp)
library(emmeans)

## create expression dataframe of neuropeptides
burtoni.snseq.combined.sct.neurons.np.expression = burtoni.snseq.combined.sct.all.neurons@assays$SCT@data %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  # mutate(gene = toupper(gene)) %>%
  filter(gene %in% neuropeptides.genes) %>% 
  column_to_rownames('gene') %>% 
  t() %>% 
  as.data.frame()


# remove low expressing genes
# remove low expressing cells
burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression #%>% 
# remove_nonexp(method = "percentage",
#               min_exp = 0.68,
#               min_percentage_samples = .05)  %>% 
# t() %>% 
# remove_nonexp(method = "percentage",
#               min_exp = 0.68,
#               min_percentage_samples = .05) 

# count cells and genes
ncol(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# 11523 cells out of 12641
nrow(burtoni.snseq.combined.sct.neurons.np.expression.filter)
# 60 genes out of 60

###rename to actual gene name
burtoni.snseq.combined.sct.neurons.np.expression.filter = burtoni.snseq.combined.sct.neurons.np.expression.filter %>% 
  # t() %>% 
  as.data.frame() %>% 
  dplyr::rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642) %>% 
  t()

### get list of sample ids
burtoni.snseq.combined.sct.neurons.np.expression.samples = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
  rownames_to_column('Cell.id') %>%  
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              rownames_to_column("Cell.id") %>% 
              dplyr::select(c(orig.ident,
                              Cell.id))) %>% 
  relocate(orig.ident, .after = Cell.id)

### get list of genotype ids
burtoni.snseq.combined.sct.neurons.np.expression.genotypes = burtoni.snseq.combined.sct.neurons.np.expression  %>% 
  rownames_to_column('Cell.id') %>%  
  full_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>%
              rownames_to_column("Cell.id") %>% 
              dplyr::select(c(Genotype.id,
                              Cell.id))) %>% 
  relocate(Genotype.id, .after = Cell.id)

#### convert to presence/absence matrix 
# create new matrix
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj = burtoni.snseq.combined.sct.neurons.np.expression.filter
# set all non-zero values to one
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj[burtoni.snseq.combined.sct.neurons.np.expression.filter.adj > 0] <- 1

### extract overlap of cell types for each genotype
# make genotype ID list
Genotype.id.list = burtoni.snseq.combined.sct.all.neurons@meta.data %>%
  pull(Genotype.id) %>% 
  unique()

# create expression matrix dataframe
burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df = t(burtoni.snseq.combined.sct.neurons.np.expression.filter.adj) %>% 
  as.data.frame()


## create neuropeptides per nuclei per count matrix
neuropeptides.per.neuron.df.geno.count = burtoni.snseq.combined.sct.neurons.np.expression.filter.adj.df %>% 
  as.data.frame() %>% 
  rownames_to_column('Cell.id') %>% 
  left_join(burtoni.snseq.combined.sct.all.neurons@meta.data %>% 
              rownames_to_column('Cell.id') %>% 
              dplyr::select(orig.ident,
                            Cell.id, 
                            Genotype.id)) %>% 
  dplyr::select(-c(Cell.id)) %>% 
  group_by(Genotype.id) %>% 
  mutate(keep = 1,
         total.count = sum(keep)) %>% 
  dplyr::select(-c(keep)) %>% 
  ungroup() %>% 
  pivot_longer(cols = -c(orig.ident,
                         Genotype.id,
                         total.count),
               names_to = 'neuropeptide',
               values_to = 'presence') %>% 
  group_by(orig.ident,
           Genotype.id,
           neuropeptide,
           total.count) %>% 
  summarise(neuropeptide.nuclei.count = sum(presence)) %>% 
  ungroup() %>% 
  separate_wider_delim(Genotype.id,
                       delim = '.',
                       names = c('Status',
                                 'Genotype'),
                       cols_remove = F) %>% 
  mutate(percent = 100*(neuropeptide.nuclei.count/total.count))

### run glm on every cluster
## use with proportion data
# create empty data frame
neuropeptide.genes.glm = data.frame(matrix(ncol = 6,
                                           # ncol = 7,
                                           nrow = 0))

#provide column names
colnames(neuropeptide.genes.glm) <- c("contrast",
                                      "Estimate",
                                      "Std.error",
                                      "z.value",
                                      "adj.pvalue",
                                      # "z.ratio",
                                      "Genes")
## loop through each cluster
for (i in unique(neuropeptides.per.neuron.df.geno.count$neuropeptide)) {
  ## run interaction binomial GLM on sex and status
  tmp = glm(cbind(neuropeptide.nuclei.count, Others) ~ Status, 
            data = neuropeptides.per.neuron.df.geno.count %>% 
              filter(neuropeptide == i) %>% 
              mutate(Others = total.count - neuropeptide.nuclei.count,
                     Status = as.factor(Status)), 
            family = binomial)
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  
  # save results
  tmp.res.df = data.frame(Estimate = tmp.res$test$coefficients,
                          Std.error = tmp.res$test$sigma,
                          z.value = tmp.res$test$tstat,
                          adj.pvalue = tmp.res$test$pvalues
                          # z.ratio = NA
                          ) %>% 
    rownames_to_column('contrast') %>% 
    mutate(Genes = i)
  
  # 
  # ## run all pairwise comparison
  # # Tukey needed for all pairwise comparisons
  # tmp2 = glm(cbind(neuropeptide.nuclei.count, Others) ~ orig.ident, 
  #            data = neuropeptides.per.neuron.df.geno.count %>% 
  #              filter(neuropeptide == i) %>% 
  #              mutate(Others = total.count - neuropeptide.nuclei.count,
  #                     Status = as.factor(Status)), 
  #            family = binomial)
  # 
  # # use emmeans to get pairwise differences
  # emm = emmeans(tmp2,
  #               ~ "orig.ident")
  # emm.res = pairs(emm)
  # 
  # 
  # # save results
  # emm.res.df = data.frame(Estimate = summary(emm.res)$estimate,
  #                         Std.error = summary(emm.res)$SE,
  #                         z.ratio = summary(emm.res)$z.ratio,
  #                         adj.pvalue = summary(emm.res)$p.value,
  #                         contrast = summary(emm.res)$contrast,
  #                         z.value = NA) %>%
  #   mutate(Genes = i)
  # 
  ## combine in data frame
  neuropeptide.genes.glm = neuropeptide.genes.glm %>% 
    # rbind(emm.res.df) %>% 
    rbind(tmp.res.df) 
}


## FDR correction for all comparisons
neuropeptide.genes.glm.fdr = neuropeptide.genes.glm %>%
  # filter(is.na(z.ratio)) %>%
  filter(contrast != "(Intercept)") %>%
  group_by(Genes) %>% 
  mutate(p.adjust.fdr = p.adjust(adj.pvalue,
                                 method = 'fdr',
                                 n = 60))

# combine with data frame
neuropeptide.genes.glm = neuropeptide.genes.glm %>%
  left_join(neuropeptide.genes.glm.fdr)



## create rounded p.value to make it easier to read
neuropeptide.genes.glm = neuropeptide.genes.glm %>%
  mutate(adj.pvalue.round = ifelse(is.na(p.adjust.fdr),
                                   round(adj.pvalue, digits = 4),
                                   round(p.adjust.fdr, digits = 4)))

## save glm table
write_csv(neuropeptide.genes.glm,
          file = 'neuropeptides/comparison_new/social.status/neuropeptide.genes.glm.csv')

# neuropeptide.genes.glm = read.csv('neuropeptides/comparison_new/social.status/neuropeptide.genes.glm.csv')

### graph
## get average percentage of neuropeptide nuclei state
neuropeptides.per.neuron.df.geno.count.percent = neuropeptides.per.neuron.df.geno.count %>% 
  group_by(Status,
           neuropeptide) %>% 
  summarise(percent.avg = mean(percent))

# get list of significant genes for males
neuropeptide.genes.glm.genes = neuropeptide.genes.glm %>% 
  filter(contrast == 'StatusSub') %>% 
  filter(adj.pvalue.round <= 0.05) %>% 
  mutate(direction = ifelse(Estimate < 0,
                                 "Dom.Bias",
                                 "Sub.Bias")) %>% 
  dplyr::select(direction,
                Genes) %>% 
  dplyr::rename(neuropeptide = Genes)

## graph
# label with significance 
neuropeptides.per.neuron.df.geno.count.percent.fig = neuropeptides.per.neuron.df.geno.count.percent %>% 
  pivot_wider(id_cols = c('neuropeptide'),
              names_from = 'Status',
              values_from = 'percent.avg') %>% 
  full_join(neuropeptide.genes.glm.genes) %>% 
  mutate(label = ifelse(is.na(direction),
                        '',
                        neuropeptide),
         direction = ifelse(is.na(direction),
                                 'no bias',
                            direction)) 

# graph
neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
  ggplot(aes(x = Dom,
             y = Sub,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(fill = direction),
             size = 10,
             shape = 21) +
  # geom_label_repel(max.overlaps = 30,
  #                  size =7.5,
  #                  position = ggpp::position_nudge_to(x = c(15, 60)),
  #                  box.padding = 0.5)+
  geom_label_repel(data = neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
                     filter(direction == 'Sub.Bias'), 
                   aes(label = label),
                   direction = "y",  
                   # nudge_x = -0.5,      
                   nudge_y = 1,        
                   # xlim = c(0, 25),
                   ylim = c(5, NA),
                   seed = 1,
                   size =7.5,
                   max.overlaps = 100) +
  geom_label_repel(data = neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
                     filter(direction != 'Sub.Bias'), 
                   aes(label = label),
                   direction = "both",  
                   # nudge_x = 1,
                   # nudge_y = -0.5,
                   xlim = c(15, NA),
                   # ylim = c(5, NA),
                   seed = 1,
                   max.overlaps = 100,
                   # box.padding = 1,
                   size =7.5,
                   segment.colour = "black") +
  theme_classic() +
  scale_size(guide = 'none')+
  ylab('Subordinate neuropeptide nuclei (%)') +
  xlab('Dominant neuropeptide nuclei (%)')+
  # ggtitle('Male Neuropeptide Cell Proportion') +
  xlim(0,100)+
  ylim(0,100)+ 
  scale_fill_manual(
    labels = c("Dominant bias",
               "Subordinate bias",
               "No bias"),
    values = c('Dom.Bias' = "#4e499e",
               'Sub.Bias' = "#60bb46",
               'No bias' = "grey")) +
  theme(legend.position = c(0.25,
                            0.9),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_blank())+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title =element_text(size=20, face='bold')) 
ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status paper glm.png', 
       height = 10,
       width = 10)

# graph no labels
neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
  ggplot(aes(x = Dom,
             y = Sub,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(fill = direction),
             size = 10,
             shape = 21) +
  # geom_label_repel(max.overlaps = 30,
  #                  size =7.5,
  #                  position = ggpp::position_nudge_to(x = c(15, 60)),
  #                  box.padding = 0.5)+
  # geom_label_repel(data = neuropeptides.per.neuron.df.geno.count.percent.fig %>%
  #                    filter(direction == 'Sub.Bias'),
  #                  aes(label = label),
  #                  direction = "y",
  #                  # nudge_x = -0.5,
  #                  nudge_y = 1,
  #                  # xlim = c(0, 25),
  #                  ylim = c(5, NA),
  #                  seed = 1,
  #                  size =7.5,
  #                  max.overlaps = 100) +
  geom_label_repel(data = neuropeptides.per.neuron.df.geno.count.percent.fig %>%
                     filter(direction != 'Sub.Bias') %>% 
                     filter(label %in% c('GNRH1',
                                         'oxt',
                                         'prl')),
                   aes(label = label),
                   direction = "both",
                   # nudge_x = 1,
                   # nudge_y = -0.5,
                   xlim = c(15, NA),
                   # ylim = c(5, NA),
                   seed = 1,
                   max.overlaps = 100,
                   # box.padding = 1,
                   size =7.5,
                   segment.colour = "black") +
  theme_classic() +
  scale_size(guide = 'none')+
  ylab('Subordinate neuropeptide nuclei (%)') +
  xlab('Dominant neuropeptide nuclei (%)')+
  # ggtitle('Male Neuropeptide Cell Proportion') +
  xlim(0,100)+
  ylim(0,100)+ 
  scale_fill_manual(
    labels = c("Dominant bias",
               "Subordinate bias",
               "No bias"),
    values = c('Dom.Bias' = "#4e499e",
               'Sub.Bias' = "#60bb46",
               'No bias' = "grey")) +
  theme(legend.position = c(0.25,
                            0.9),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_blank())+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title =element_text(size=20, face='bold')) 
ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status paper glm subset label.png', 
       height = 10,
       width = 10)

# graph no labels
neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
  ggplot(aes(x = Dom,
             y = Sub,
             label = label)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(fill = direction),
             size = 10,
             shape = 21) +
  # geom_label_repel(max.overlaps = 30,
  #                  size =7.5,
  #                  position = ggpp::position_nudge_to(x = c(15, 60)),
  #                  box.padding = 0.5)+
  # geom_label_repel(data = neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
  #                    filter(direction == 'Sub.Bias'), 
  #                  aes(label = label),
  #                  direction = "y",  
  #                  # nudge_x = -0.5,      
  #                  nudge_y = 1,        
  #                  # xlim = c(0, 25),
  #                  ylim = c(5, NA),
  #                  seed = 1,
  #                  size =7.5,
  #                  max.overlaps = 100) +
  # geom_label_repel(data = neuropeptides.per.neuron.df.geno.count.percent.fig %>% 
  #                    filter(direction != 'Sub.Bias'), 
  #                  aes(label = label),
  #                  direction = "both",  
  #                  # nudge_x = 1,
  #                  # nudge_y = -0.5,
  #                  xlim = c(15, NA),
  #                  # ylim = c(5, NA),
  #                  seed = 1,
  #                  max.overlaps = 100,
  #                  # box.padding = 1,
  #                  size =7.5,
  #                  segment.colour = "black") +
  theme_classic() +
  scale_size(guide = 'none')+
  ylab('Subordinate neuropeptide nuclei (%)') +
  xlab('Dominant neuropeptide nuclei (%)')+
  # ggtitle('Male Neuropeptide Cell Proportion') +
  xlim(0,100)+
  ylim(0,100)+ 
  scale_fill_manual(
    labels = c("Dominant bias",
               "Subordinate bias",
               "No bias"),
    values = c('Dom.Bias' = "#4e499e",
               'Sub.Bias' = "#60bb46",
               'No bias' = "grey")) +
  theme(legend.position = c(0.25,
                            0.9),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_blank())+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title =element_text(size=20, face='bold')) 
ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status paper glm no label.png', 
       height = 10,
       width = 10)
 