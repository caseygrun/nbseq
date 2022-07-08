library(Matrix)
library(scran)
library(scuttle)
library(scater)
library(tidyverse)

norm_scran <- function(x, subset.row=NULL) {
	if (verbose) { message("Estimating size factors...") }

	# calculate initial estimate of library size factors per sample
	sf <- scuttle::librarySizeFactors(x, subset_row=subset.row)
	y <- scater::normalizeCounts(x, size_factors=sf, subset_row=subset.row)

	if (verbose) { message("Modeling gene variance structure...") }

	# estimate total, technical and biological components to per-gene variance.
	fit <- scran::modelGeneVar(y)
	# alternative fitting which assumes technical variation follows a poisson distribution
	fit_poisson <- scran::modelGeneVarByPoisson(y)

	# select the top 10% or 500 most highly-variable genes (whichever is higher)
	top_hvgs <- scran::getTopHVGs(fit, n=500, prop=0.1)

	if (verbose) { message("Using ", length(top_hvgs), " most highly-variable genes") }

	# the above is normally performed internally by scran::quickCluster. We replicated it
	# here 1) to get the intermediate info, and 2) because scran::getDenoisedPCs complains
	# that identical(rownames(y), rownames(fit)) == FALSE, even though
	# all(rownames(y) == rownames(fit)) == TRUE.
	# anyway, we manually replicate that call here, which identifies the number of PCs
	# such that the sum of variance of remaining PCs is less than the sum of technical
	# variance (noise). That is, find the minimum number of PCs that would be needed to
	# capture the biological variance. Pass that to `d` which will use this metric
	# for its clustering.
	rownames(y) <- rownames(fit)
	yh <- scran::getDenoisedPCs(y, technical=fit, subset.row=top_hvgs)$components
	d <- ncol(yh)

	if (verbose) { message("Using ", d, " principal components for clustering") }

	# assign samples into clusters
	clusters <- scran::quickCluster(ft, d=d)

	if (verbose) { message("Calculating sum factors by deconvolution...") }

	# calculate size factors for normalization by the deconvolution method (Lun, Bach, and Marioni 2016 <https://doi.org/10.1186/s13059-016-0947-7>)
	size_factors_deconvoluted <- scran::calculateSumFactors(ft, cluster=clusters)


	if (verbose) { message("Normalizing counts...") }

	# logNormCounts requires a SingleCellExperiment, but calls this method underneath
	xh = scater::normalizeCounts(x, size.factors=size_factors_deconvoluted)


	obs <- data.frame(
		size_factor=sf,
		size_factor_deconvoluted=size_factors_deconvoluted,
		cluster=clusters
	)

	var <- dplyr::rename(as.data.frame(fit), var=total, var_biological=bio, var_technical=tech, p_value=p.value)

	return(list(
		'x' = xh,
		'obs' = obs,
		'var' = var
	))
}

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

write_feature_table_sparse <- function(ft, dir, samples = NULL, features = NULL) {
    dir.create(dir, recursive=TRUE, showWarnings = FALSE)
    writeMM(ft, file = file.path(dir, 'table.mtx'))
    if (is.null(samples)) samples <- rownames(ft)
    write.table(samples, file.path(dir, 'samples.tsv'), sep="\t", row.names = FALSE, col.names = FALSE)

    if (is.null(features)) features <- colnames(ft)
    write.table(features, file.path(dir, 'features.tsv'), sep="\t", row.names = FALSE, col.names = FALSE)
}


write_norm_scran <- function(out, path_x, path_obs, path_var) {
	write_feature_table_sparse(ft=t(out$x), dir=path_x)
	write_csv(tibble::rownames_to_column(out$obs, 'ID'), path_obs)
	write_csv(tibble::rownames_to_column(out$var, 'feature'), path_var)
}


args <- commandArgs(TRUE)

input_ft_path <- args[[1]]
output_x <- args[[2]]
output_obs <- args[[3]]
output_var <- args[[4]]
verbose <- as.logical(as.integer(args[[5]]))

if (verbose) { message("Reading input...") }
ft <- read_biom_hd5_matrix(input_ft_path)
out <- norm_scran(ft)
if (verbose) { message("Writing output...") }
write_norm_scran(out, output_x, output_obs, output_var)
