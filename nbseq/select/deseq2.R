cat(R.version$version.string, "\n")
args <- commandArgs(TRUE)

# path to the metadata file, in tsv format. rows are samples, columns are variables.
# all columns must have headers, e.g. no implied index
col_data_path <- args[[1]]
# col_data_path <- opts$col_data

# path to the count data, in tsv format, with samples in rows, and ASV hashes
# in columns
count_data_path <- args[[2]]
# count_data_path <- opts$count_data

# name of the column from the qiime2 metadata file to use to identify samples;
# values of this column should match row names in count_data_path
identifier_name <- args[[3]]
# identifier_name <- opts$identifier_name

# formula in R syntax for the design matrix to fit, starting with ~
design_formula <- args[[4]]
# design_formula <- opts$design_formula

# coefficients to examine
coef_names <- args[[5]]
# coef_names <- opts$coef_names

# path to the output file(s)
out_data_path <- args[[6]]
# out_data_path <- opts$out_path

# load.library('optparse')
# load.library('biom')
# load.library('metagenomeSeq', bioconductor=TRUE)
#
# make option list and parse command line
# option_list <- list(
#     make_option(c("--col_data"), type="character",
#         help="Path to the qiime2 metadata file, in TSV format [required]."),
#     make_option(c("--count_data"), type="character",
#         help="Path to the count data, in TSV format [required]."),
#     make_option(c("--identifier_name"), type="character",
#         help="name of the column from the metadata file to use to identify samples; values in this column should match rownames in the count data table.", default="Specimen_ID"),
#     make_option(c("--design_formula"), type="character",
#         help="formula in R syntax for the design matrix to fit, starting with ~ [required]."),
#     make_option(c("--coef_names"), type="character",
#         help="comma-separated list of names of the coefficients for which to print comparisons; should be members resultsNames (possible values will be printed) [required]."),

#     make_option(c("-i", "--input_path"), type="character",
#         help="Input otu table [required]."),
#     make_option(c("-o", "--out_path"), type="character", default='.',
#         help="Output directory [default %default]"),
#     make_option(c("-z", "--DESeq_negatives_to_zero"), type="character", default=NULL,
#         help="set the negatives that result from DESeq transformation to zero")
# )
# opts <- parse_args(OptionParser(option_list=option_list), args=args)
#
# # Error checking
# if(is.null(opts$input_path)) stop('Please supply an otu table.')

# identifier_name <- "Specimen_ID"
# design_formula <- " ~ run + Age_at_visit + CF_Y_N + Age_at_visit*CF_Y_N"
# col_data_path <- "/Users/caseygrun/Scratch/mapfile-qiime2.txt"
# count_data_path <- "/Users/caseygrun/Scratch/output/tmp/1456b-table.csv"
# coef_name <- "CF_Y_N_Y_vs_N"

read_biom_hd5_matrix <- function(biom_file) {
	x = rhdf5::h5read(biom_file,"/",read.attributes = TRUE)
	library(Matrix)

	generate_matrix <- function(x){

	    nrow = length(x$observation$ids)
	    ncol = length(x$sample$ids)


	    # x$sample$matrix stores data in compressed sparse column (CSC) format. csc_matrix in scipy, dgCMatrix in R Matrix
	    # x$observation$matrix stores data in compressed sparse row (CSR) format. csr_matrix in scipy, dgRMatrix in R Matrix

	    # for details, see
	    # class? CsparseMatrix
	    # class? dgCMatrix
	    dx = new("dgCMatrix",
	             # sample/matrix/indices      : <int32> A (nnz,) dataset containing the row indices (e.g., maps into observation/ids)

	             # i: Object of class "integer" of length nnzero (number of non-zero elements).
	             # These are the 0-based row numbers for each non-zero element in the matrix,
	             # i.e., i must be in 0:(nrow(.)-1).
	             i=as.integer(x$sample$matrix$indices),

	             # sample/matrix/indptr       : <int32> A (N+1,) dataset containing the compressed column offsets

	             # p: integer vector for providing pointers, one for each column, to the initial
	             # (zero-based) index of elements in the column. .@p is of length ncol(.) + 1,
	             # with p[1] == 0 and p[length(p)] == nnzero, such that in fact, diff(.@p) are
	             # the number of non-zero elements for each column. In other words, m@p[1:ncol(m)]
	             # contains the indices of those elements in m@x that are the first elements in the
	             # respective column of m.
	             p=as.integer(x$sample$matrix$indptr),

	             # sample/matrix/data         : <float64> A (nnz,) dataset containing the actual matrix data
	             # x: Object of class "numeric" - the non-zero elements of the matrix.
	             x=as.numeric(x$sample$matrix$data),

	             Dim=c(nrow, ncol),
	             Dimnames = list(x$observation$ids, x$sample$ids)

	            )
	    dx
	}
	generate_matrix(x)
}

read_feature_table_hd5 <- function(biom_file) {
	return(t(read_biom_hd5_matrix(biom_file)))
}



suppressMessages(library("DESeq2"))
suppressMessages(library("BiocParallel"))

cat("DESeq2 R package version:", as.character(packageVersion("DESeq2")), "\n")
cat("BiocParallel R package version:", as.character(packageVersion("BiocParallel")), "\n")

design = as.formula(design_formula)

colData <- read.csv(col_data_path,sep="\t",na.strings="NA")

# drop rows of colData where one of the variables in the design formula is na
design_vars <- all.vars(design)
narows <- apply(colData[,design_vars], 1, function(x) any(is.na(x)))
colData <- colData[ !narows, ]

# attach the data
row.names(colData) <- colData[,identifier_name]

# default qiime2 outputs table with samples in rows, ASVs in columns; we need
# the transpose of this (ASVs in rows, samples in columns)
countData <- read.csv(count_data_path,row.names=1,,check.names=FALSE)
countData <- t(countData)

countData <- read_feature_table_hd5(count_data_path)
countData <- as.matrix(countData)


# drop rows of countData with na's
narows <- apply(countData, 1, function(x) any(is.na(x)))
countData <- countData[ !narows, ]

# only look at every n samples to speed up processing
# every_howmany_rows = 5
# countData <- countData[, seq(1, ncol(countData), every_howmany_rows)]


# it is critical that the order of the columns in the countData exactly
# matches the order of the rows in the colData, so choose samples which are
# represented in both colData and countData, then so sort the rows of colData
# (samples) and the columns of colData (samples) into that same order.

# find all samples which are represented in the colData and countData
sample_names = intersect(rownames(colData),colnames(countData))

cat("Excluding samples in count data which have no metadata:\n")
setdiff(colnames(countData),sample_names)

cat("Excluding samples in metadata which have no count data:\n")
setdiff(rownames(colData),sample_names)

# cat("Samples to analyze:\n")
# sample_names

# colData <- colData[colnames(countData),]
# countData <- countData[,row.names(colData)]
colData <- colData[sample_names,]
countData <- countData[,sample_names]

# check to make sure the order of identifiers is the same between the countData
# and the colData
all(colData[,identifier_name] %in% colnames(countData))
all(row.names(colData) == colnames(countData))

cat("Structure of colData to be used:\n")
str(colData[,design_vars])

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# build DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = as.matrix(countData),
                              colData = colData,
                              design = design) #~ Cohort + GAE + Cohort*GAE)

# use a custom geometric mean function that is zero-resistant
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# estimate the size factors using the custom geometric means
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

# filter the DESeq dataset to remove low-power samples
power_threshold = 10
keep <- rowSums(counts(dds, normalized=TRUE) >= power_threshold) >= power_threshold
dds <- dds[keep,]

write("Running DESeq:", stderr())
dds <- DESeq(dds, fitType="local", parallel=TRUE, BPPARAM=MulticoreParam(2))

cat("Available contrasts for results (resultsNames): \n")
resultsNames(dds)
cat("\n")
write("", stderr())


cat("Evaluating the following coefficients: \n")
coef_names <- strsplit(coef_names,",")[[1]]
coef_names

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# this block is the problem; if placed inside the loop, it does not work
# (res prints nothing)

# coef_name <- "CF_Y_N_Y_vs_N" # coef_names[[i]]
# cat("Coefficient: ", coef_name, "\n")
# coef_name
#
# cat("Results before shrinkage: \n") #"%s",coef_name)
# res <- results(dds,name=coef_name)
# res
#
# # res <- results(dds,name="COMPARISON")
# write("", stderr())
# cat("\n")
#
# lfc.type <- "ashr"
#
# if (lfc.type == "ashr") {
#     cat("Results after shrinkage (ashr): \n") #",coef_name)
#     res <- lfcShrink(dds, coef=coef_name, type="ashr")
# } else if (lfc.type == "apeglm") {
#     cat("Results after shrinkage (apeglm): \n") #",coef_name)
#     res <- lfcShrink(dds, coef=coef_name, type="apeglm")
# }
#
# res[order(res$pvalue),]
# write.csv(res[order(res$log2FoldChange),],cat(out_data_path,"-",coef_name))


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


for (i in 1:length(coef_names)) {

    coef_name <- coef_names[[i]]
    cat("Coefficient: ", coef_name, "\n")
    print(coef_name)

    cat("Results before shrinkage: \n") #"%s",coef_name)
    res <- results(dds,name=coef_name)
    print(res)

    write("", stderr())
    cat("\n")

    lfc.type <- "ashr"

    if (lfc.type == "ashr") {
        cat("Results after shrinkage (ashr): \n") #",coef_name)
        res <- lfcShrink(dds, coef=coef_name, type="ashr")
    } else if (lfc.type == "apeglm") {
        cat("Results after shrinkage (apeglm): \n") #",coef_name)
        res <- lfcShrink(dds, coef=coef_name, type="apeglm")
    }

    # out_data <- res[order(res$log2FoldChange),]
    out_data <- res[order(res$padj),]
    # print(out_data)
    out_path_coef <- paste(out_data_path,"-",coef_name,".csv",sep="")
    cat("Writing to CSV file:\n")
    print(out_path_coef)
    write.csv(out_data,out_path_coef)

    cat("\n")
}
