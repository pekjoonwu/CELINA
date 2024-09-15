#####################################################################
# Package: Celina
# Version: 1.0.0
# Date : 2024-7-11
######################################################################

##########################################################
#       Celina construct object functions and preprocessing related functions   #
##########################################################
#' @title Each Celina object has a number of slots which store information. Key slots to access
#' are listed below.
#'
#' @slot celltype_mat The raw celltype indicator or proportion matrix. 
#' Rows are celltypes, columns are spots/cells.
#' @slot gene_expression_mat Original expression matrix. By default we use 
#' normalizeCounts normalization in Scater R package.
#' @slot genes_list Different tested gene lists for cell types.
#' @slot project Name of the project (for record keeping).
#' @slot covariates The covariates in experiments (if any covariates is included).
#' @slot location Cell/spot spatial coordinates to compute the kernel matrix.
#' @slot kernelmat The kernel matrix for spatial relationship between locations.
#' @slot kernelmat_approx_U The calculated approximated kernel matrix for spatial relationship between locations.
#' @slot bandwidth The bandwidth in Gaussian kernel, users can also specify their preferred bandwidth.
#' @slot result The summary of cell type specific SE genes results for different cell types.
#' @slot approximation logical value to indicate which type of algorithms to be used in the analysis
#' @export

setClass("Celina", slots = list(
  celltype_mat = "ANY",
  gene_expression_mat = "ANY",
  genes_list = "ANY",
  project = "character",
  covariates = "ANY",
  location = "matrix", 
  kernelmat = "ANY",
  kernelmat_approx_U = "ANY",
  bandwidth = "numeric",
  result = "ANY",
  approximation = "logical"
))


#' @title Create the Celina object.
#' @param celltype_mat The raw celltype indicator or proportion matrix. Rows are celltypes, columns are spots/cells.
#' @param gene_expression_mat Original expression matrix. By default we use 
#' normalizeCounts normalization in Scater R package.
#' @param location Spatial location matrix (matrix), the dimension is n x d, n is the number of locations, d is dimensin of spatial coordinates, 
#' e.g. d = 2 for locations on 2D space. The rownames of locations and the colnames of count matrix should be matched.
#' @param covariates The covariates in experiments (matrix, if any covariates included), n x q, 
#' n is the number of locations, q is the number of covariates. 
#' The rownames of covariates and the rownames of locations should be matched.
#' @param project Name of the project (User defined).
#' @return Returns Celina object.
#' 
#' @export
Create_Celina_Object <- function(celltype_mat = NULL, 
                              gene_expression_mat = NULL, 
                              location = NULL,
                              covariates = NULL,
                              project = "Celina"){
  
  ## check dimension
  if (ncol(celltype_mat) != nrow(location) | ncol(gene_expression_mat) != nrow(location)) {
    stop ("The number of columns in two matrices, and the number of locations should be consistent.")
  } # end fi

  ## check names
  if (any(colnames(celltype_mat) != rownames(location)) | any(colnames(gene_expression_mat) != rownames(location))) {
    stop ("The cells/spots names should match between celltype matrix, gene expression matrix and location matrix")
  }# end fi
  
  ## inheriting
  object <- methods::new(
    Class = "Celina",
    celltype_mat = celltype_mat,
    gene_expression_mat = gene_expression_mat,
    location = location,
    project = project
  )
  
  if (!is.null(covariates)) {
    ## check data order should consistent
    if (!identical(rownames(covariates), rownames(location))) {
      stop("The row names of covariates and row names of location should be should be matched (covariates -- n locations x k covariates; location -- n locations x d dimension).")
    }
    q <- dim(covariates)[2]
    n_covariates <- dim(covariates)[1]
    ## remove the intercept if added by user, later intercept will add automatically
    if (length(unique(covariates[, 1])) == 1){
      covariates <- covariates[, -1]
      q <- q - 1
    }
    object@covariates <- as.matrix(covariates, n_covariates, q)
  } else {
    object@covariates <- NULL
  }
  
  object@result <- list()
  object@location <- scale(object@location)
  
  ## Return created object
  return(object)
}


#' @title Preprocess the gene expression matrix - normalize the counts and filter genes
#' @param object Celina object
#' @param cell_types_to_test a vector of cell types to be used for testing
#' @param scRNA_count a gene expression matrix of reference scRNA-seq data with gene names as rownames and cells as columns
#' @param sc_cell_type_labels a vector of cell type labels for each cell in \code{scRNA_count}
#' @param threshold a number for filtering the genes
#' @return Returns Celina object.
#' 
#' @export
preprocess_input <- function(object, cell_types_to_test,  
                             scRNA_count, sc_cell_type_labels,
                             threshold = 5e-5) {
  
  ## Pick gene list 
  genes_list_all <- get_cell_type_gene_list(spatial_count = object@gene_expression_mat, 
                          cell_types_to_test = cell_types_to_test, 
                          celltype_proportion = t(object@celltype_mat), 
                          scRNA_count = scRNA_count, 
                          sc_cell_type_labels = sc_cell_type_labels, threshold = threshold)
  object@genes_list <- genes_list_all
  
  ## Normalize the data and subset the gene expression matrix on the subset of genes
  all_genes <- unique(unlist(genes_list_all))
  normalized_counts <- scater::normalizeCounts(object@gene_expression_mat)
  filter_index <- rownames(normalized_counts) %in% all_genes
  normalized_counts <- normalized_counts[filter_index, ]
  
  # make sure gene expression and cell type prop matrix have matched locations
  common_colnames <- intersect(colnames(object@celltype_mat), colnames(normalized_counts))
  object@gene_expression_mat <- normalized_counts[, match(common_colnames, colnames(normalized_counts))]
  object@celltype_mat <- object@celltype_mat[, match(common_colnames, colnames(object@celltype_mat))]
  object@location <- object@location[match(common_colnames, rownames(object@location)),]
  
  ## Return object with preprocess gene expression matrix and 
  ## gene lists for testing for each cell type
  return(object)
}

#' Filter the genes from the count matrix based on the threshold
#' 
#' @noRd
filter_genes <- function(counts, threshold = 5e-5) {
  cat(paste0('filtering genes based on threshold = ', threshold), "\n")
  
  ## Calculate the size factor and normalize
  size_factor <- colSums(counts)
  adjusted_counts <- sweep(counts, 2, size_factor, '/')
  gene_means <- rowMeans(adjusted_counts)
  names(gene_means) <- rownames(counts)
  
  ## Filter the lowly expressed genes and mt genes
  gene_list <- names(which(gene_means > threshold))
  if (length(grep("mt-", gene_list, ignore.case = T)) > 0) {
    gene_list <- gene_list[-grep("mt-", gene_list, ignore.case = T)]
  }
  
  return(gene_list)
}

#' Filtering out the genes that are lowly expressed acorss cells or the cells that are a few number of genes expressed
#' @param counts Raw gene expression counts, p x n matrix
#' @param percentage_cell The percentage of cells that are expressed
#' @param min_total_counts Minimum counts for each cell
#' @export
FilteringGenesCells <- function(counts, percentage_cell = 0.1, min_total_counts = 10){
  idx <- which(apply(counts ,1 ,function(x) {
    sum(x != 0)}) >= (floor(percentage_cell * ncol(counts))))
  ## filtering out genes
  counts <- counts[idx,]
  ## filtering out cells
  counts <- counts[, which(apply(counts, 2, sum) > min_total_counts)]
  return(counts)
}



#' @title Obtain a list of genes to test different cell types
#'
#' @param spatial_count a gene expression matrix of spatial transcriptomics data with gene names as rownames and cells as columns
#' @param cell_types_to_test a vector of cell types to be used for testing
#' @param celltype_proportion a matrix of cell type proportions for each location in \code{spatial_count}, with locations as rows and cell types as columns
#' @param scRNA_count a gene expression matrix of reference scRNA-seq data with gene names as rownames and cells as columns
#' @param sc_cell_type_labels a vector of cell type labels for each cell in \code{scRNA_count}
#' @param threshold a number for filtering the genes
#' @return cell_type_means, a data_frame of genes by cell types for mean normalized expression
#' @export
get_cell_type_gene_list <- function(spatial_count, cell_types_to_test, celltype_proportion, 
                                    scRNA_count, sc_cell_type_labels, threshold) {
  ## Get a initial gene list and some preprocessing
  spatial_count <- as.matrix(spatial_count)
  scRNA_count <- as.matrix(scRNA_count)
  gene_list_tot <- filter_genes(counts = spatial_count, threshold = threshold)
  gene_list_tot_use <- intersect(gene_list_tot, rownames(scRNA_count))
  scRNA_count_use <- scRNA_count[gene_list_tot_use, ]
  celltype_proportion <- celltype_proportion[, cell_types_to_test]
  proportions <- colMeans(celltype_proportion)/sum(colMeans(celltype_proportion))
  
  ## Calculate scRNA mean expression for all cell types and renormalize to match spatial data
  scRNA_norm <- get_scRNA_info(scRNA_count_use, sc_cell_type_labels, cell_type_names = cell_types_to_test)
  scRNA_renormalized <- get_renormalized_scRNA(spatial_count = spatial_count, 
                                              cell_type_means = scRNA_norm, 
                                              gene_list = gene_list_tot_use,
                                              proportions = proportions)
  
  ## Filter spatial locations
  location_names <- rownames(celltype_proportion)
  res_filter <- filter_locations(location_names, cell_types_to_test, celltype_proportion, thresh = 0.8)
  location_names <- res_filter$location_names
  celltype_proportion <- res_filter$celltype_proportion
  spatial_nUMI <- colSums(spatial_count[, location_names])
  
  ## Get gene list for each cell type
  cat("Get gene list for each cell type.\n")
  filtered_genelist <- list()
  for(cell_type in cell_types_to_test) {
    tmp_gene_list <- get_gene_list(celltype_proportion = celltype_proportion, 
                               location_names = location_names, cell_type = cell_type,
                               spatial_nUMI = spatial_nUMI,
                              gene_list_type = gene_list_tot, scRNA_renormalized = scRNA_renormalized,
                              cell_types_present = cell_types_to_test)
    filtered_genelist[[cell_type]] <- tmp_gene_list
  }
  
  return(filtered_genelist)
}


#' Get tested gene list for a cell type
#' 
#' @noRd
get_gene_list <- function(celltype_proportion, location_names, 
                          cell_type, spatial_nUMI, gene_list_type, 
                          scRNA_renormalized, cell_types_present) {
  ## Normalize cell type proportion and
  celltype_proportion <- sweep(celltype_proportion, 1, apply(celltype_proportion, 1, sum), '/')
  total_proportions <- colSums(celltype_proportion[location_names, , drop = FALSE])[cell_type]
  
  ## Threshold refer to spacexr packages
  UMI_list <- spatial_nUMI[names(which(celltype_proportion[, cell_type] >= .99))]
  if (length(UMI_list) < 10) {
    UMI_list <- spatial_nUMI[names(which(celltype_proportion[, cell_type] >= .8))]
  }
  
  if (length(UMI_list) < 10) {
    UMI_list <- spatial_nUMI[names(which(celltype_proportion[, cell_type] >= .5))]
  }
  
  if (length(UMI_list) < 10) {
    UMI_list <- spatial_nUMI[names(which(celltype_proportion[, cell_type] >= .01))]
  }
  UMI_m <- median(UMI_list)
  expr_thresh <-  15 / (total_proportions * UMI_m)
  
  ## Pick gene list for cell type
  gene_list_type_intersect <- intersect(gene_list_type, rownames(scRNA_renormalized))
  gene_list_type <- setdiff(gene_list_type_intersect, 
                            gene_list_type_intersect[which(scRNA_renormalized[gene_list_type_intersect, cell_type] < expr_thresh)])
  cell_type_means <- scRNA_renormalized[gene_list_type, cell_types_present]
  if (dim(celltype_proportion)[2] > 1) {
    cell_prop <- sweep(cell_type_means, 1, apply(cell_type_means, 1, max), '/')
    gene_list_type <- gene_list_type[which(cell_prop[gene_list_type, cell_type] > 0.5)]
  }
  
  ## Return results
  return(gene_list_type)
}

#' Renormalizes cell_type_means to have same average as the spatial count data
#' Output is the renormalized gene by cell types matrix representing cell type specific mean
#' 
#' @noRd
get_renormalized_scRNA <- function(spatial_count, cell_type_means, gene_list, proportions) {
  gene_list_use <- intersect(gene_list, rownames(spatial_count))
  gene_list_use <- intersect(gene_list_use, rownames(cell_type_means))
  
  bulk_vec <- rowSums(spatial_count[gene_list_use, ])
  weight_avg <- rowSums(sweep(cell_type_means[gene_list_use, ], 2, proportions / sum(proportions), '*'))
  target_means_genes <- bulk_vec[gene_list_use]/sum(spatial_count[gene_list_use, ])
  cell_type_means_renormalized <- sweep(cell_type_means[gene_list_use, ], 1, 
                                        weight_avg / target_means_genes, '/')
  
  return(cell_type_means_renormalized)
}

#' Filter the spatial locations which have summation of tested cell type proportions less than a threshold
#' 
#' @noRd
filter_locations <- function(location_names, cell_types, celltype_proportion, thresh = 0.9999) {
  location_names <- location_names[(rowSums(celltype_proportion[location_names, cell_types, drop = FALSE]) >= thresh)]
  celltype_proportion <- celltype_proportion[location_names, cell_types, drop = FALSE]
  
  return (list(location_names = location_names, celltype_proportion = celltype_proportion))
}


#' Computes averaged normalized expression (summing to 1) for all cells within a cell type
#' Output is the gene by cell types matrix representing cell type specific mean
#' 
#' @noRd
get_scRNA_info <- function(scRNA_count, sc_cell_type_labels, cell_type_names = NULL) {
  ## Sanity check
  sc_cell_type_labels <- as.factor(sc_cell_type_labels)
  unique_celltype <- unique(as.character(sc_cell_type_labels))
  if (is.null(cell_type_names)) {
    cell_type_names <- levels(sc_cell_type_labels)
  } else {
    if(length(setdiff(cell_type_names, unique_celltype)) > 0){
      stop("cell_type_names must be the a subset of sc_cell_type_labels")
    }
  }
  
  nUMI <- colSums(scRNA_count)
  n_cell_types <- length(cell_type_names)
  
  get_gene_mean <- function(cell_type) {
    cell_type_data <- scRNA_count[, sc_cell_type_labels == cell_type]
    cell_type_umi <- nUMI[sc_cell_type_labels == cell_type]
    normData <- sweep(cell_type_data, 2, cell_type_umi, `/`)
    return(rowSums(normData)/dim(normData)[2])
  }
  
  cell_type <- cell_type_names[1]
  cell_type_means <- data.frame(get_gene_mean(cell_type))
  colnames(cell_type_means)[1] <- cell_type
  for (cell_type in cell_type_names[2:length(cell_type_names)]) {
    cell_type_means[cell_type] <- get_gene_mean(cell_type)
  }
  return(cell_type_means)
}



##########################################################
#       Celina calculate kernel matrices related functions             #
##########################################################
#' @title Calculate the kernel matrix.
#' @param object The Celina object.
#' @param approximation Option to use the default 11 kernels combined results or 1 kernel approximated results
#' True is to proceed with approximated algorithm
#' @param bandwidth.set.by.user User could select their own bandwidth (a numeric value) 
#' if the recommended bandwidth doesn't work in their dataset.
#' @param bandwidthtype The type of bandwidth to be used in Gaussian kernel, 
#' "SJ" for Sheather & Jones (1991) method (usually used in small size datasets), 
#' "Silverman" for Silverman's ‘rule of thumb’ method (1986) (usually used in large size datasets).
#' @param sparseKernel Select "TURE" if the user wants to use a sparse kernel matrix or "FALSE" if not. 
#' It is recommended to choose sparseKernel = "TRUE" when sample size is large.
#' @param sparseKernel_tol When sparseKernel = TRUE, 
#' the cut-off value when building sparse kernel matrix, 
#' any element in the kernel matrix greater than sparseKernel_tol will be kept, 
#' otherwise will be set to 0 to save memory.
#' @param sparseKernel_ncore When sparseKernel=TRUE, 
#' the number of CPU cores to build sparse kernel matrix.
#' @return Returns Celina object.
#' 
#' @export
Calculate_Kernel <- function(object, 
                             approximation = FALSE,
                             bandwidth.set.by.user = NULL,
                             bandwidthtype = "Silverman", 
                             sparseKernel = FALSE,
                             sparseKernel_tol = 1e-10,
                             sparseKernel_ncore = 1) {
  
  if (approximation == F) {
    ## Default - construct 11 kernels
    kernel_res <- Celina_kernel(object@location)
    object@kernelmat <- kernel_res
    object@approximation <- FALSE
    
    return(object)
  } else {
    ## Alternative - 1 kernels with approximation 
    object@approximation <- TRUE
    sparseKernel <- TRUE
    ## Select the bandwidth for the approximation kernel
    if (is.null(bandwidth.set.by.user) & is.null(bandwidthtype)) {
      stop("Please specify a bandwidth for kernel matrix, or specify a bandwidthtype and 
           provide a normalized expression (gene by location).")
    } else if(!is.null(bandwidth.set.by.user)) {
      ## User defined bandwidth
      object@bandwidth <- bandwidth.set.by.user
      cat(paste("## The bandwidth is: ", object@bandwidth, " \n"))
    } else if(!is.null(bandwidthtype)) {
      ## Silverman selected bandwidth
      object@bandwidth <- bandwidth_select(object@gene_expression_mat, method = bandwidthtype)
      cat(paste("## The bandwidth is: ", object@bandwidth, " \n"))
    }
    
    ## Calculate the approximated kernel matrix
    if (sparseKernel == FALSE) {
      cat(paste("## Calculating kernel matrix\n"))
      object@kernelmat <- kernel_build(kerneltype = "gaussian", 
                                      location = object@location, bandwidth = object@bandwidth)
      
      ## SVD on the original gaussian kernel and extract the first two eigen values
      svd_results <- RSpectra::eigs_sym(object@kernelmat, k = 2)
      newkernel <- svd_results$vectors %*% diag(svd_results$values) %*% t(svd_results$vectors)
      ## Normalize the kernel
      k_norm_diag <- 1/sqrt(diag(newkernel))
      newkernel <- newkernel * (k_norm_diag %*% t(k_norm_diag))
      object@kernelmat_approx_U <- as.matrix(svd_results$vectors %*% sqrt(diag(svd_results$values)) * k_norm_diag)
    } else if (sparseKernel == TRUE) {
      cat(paste("## Calculating sparse kernel matrix\n"))
      object@kernelmat <- kernel_build_sparse(kerneltype = "gaussian", 
                                              location = object@location, 
                                              bandwidth = object@bandwidth,
                                              tol = sparseKernel_tol, 
                                              ncores = sparseKernel_ncore)
      ## SVD on the original gaussian kernel and extract the first two eigen values
      svd_results <- RSpectra::eigs_sym(object@kernelmat, k = 2)
      newkernel <- svd_results$vectors %*% diag(svd_results$values) %*% t(svd_results$vectors)
      ## Normalize the kernel
      k_norm_diag <- 1/sqrt(diag(newkernel))
      newkernel <- newkernel * (k_norm_diag %*% t(k_norm_diag))
      object@kernelmat_approx_U <- as.matrix(svd_results$vectors %*% sqrt(diag(svd_results$values)) * k_norm_diag)
    }
    cat(paste("## Finished calculating kernel matrix.\n"))
    return(object)
  }  
}


#' @title Build kernel matrix.
#' @description This function calculates kernel matrix from spatial locations.
#' @param kerneltype The type of kernel to be used, either "gaussian", 
#' or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @return The kernel matrix for spatial relationship between locations.
kernel_build <- function (kerneltype = "gaussian", location, bandwidth)
{
  if (kerneltype == "gaussian") {
    K <- exp(-1 * as.matrix(stats::dist(location) ^ 2)/bandwidth)
  }
  else if (kerneltype == "cauchy") {
    K <- 1/(1 + 1 * as.matrix(stats::dist(location) ^ 2)/bandwidth)
  }
  else if (kerneltype == "quadratic") {
    ED2 <- 1 * as.matrix(stats::dist(location) ^ 2)
    K <- 1 - ED2/(ED2 + (bandwidth))
  }
  return(K)
}


#' @title Build sparse kernel matrix.
#' @description This function calculates kernel matrix.
#' @param kerneltype The type of kernel to be used, either "gaussian", 
#' or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @param tol A numeric value of cut-off value when building sparse kernel matrix.
#' @param ncores A integer value of number of CPU cores to use when building sparse kernel matrix.
#' @return The sparse kernel matrix for spatial relationship between locations.
kernel_build_sparse <- function(kerneltype, location, bandwidth, tol, ncores) {
  if (kerneltype == "gaussian") {
    fx_gaussian <- function(i){
      line_i <- rep(0, dim(location)[1])
      line_i[i] <- 1
      line_i[-i] <- exp(-(pdist::pdist(location[i, ], location[-i, ])@dist^2)/(bandwidth))
      ind_i <- which(line_i >= tol)
      return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)), "val_i" = line_i[ind_i]))
    }
    ## Aggregate the sparse matrix
    results <- parallel::mclapply(1:dim(location)[1], fx_gaussian, mc.cores = ncores)
    tib <- tidyr::tibble(results) %>% tidyr::unnest_wider(results)
    K_sparse <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]), 
                                     x = unlist(tib[[3]]), dims = c(dim(location)[1],
                                                                    dim(location)[1]))
  } else if (kerneltype == "cauchy") {
    fx_cauchy <- function(i){
      line_i <- rep(0, dim(location)[1])
      line_i[i] <- 1
      line_i[-i] <- 1/(1 + (pdist::pdist(location[i, ],location[-i, ])@dist^2)/(bandwidth))
      ind_i <- which(line_i >= tol)
      return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)), "val_i" = line_i[ind_i]))
    }
    ## Aggregate the sparse matrix
    results <- parallel::mclapply(1:dim(location)[1], fx_cauchy, mc.cores = ncores)
    tib <- tidyr::tibble(results)  %>%  tidyr::unnest_wider(results)
    K_sparse <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]), 
                                     x = unlist(tib[[3]]), dims = c(dim(location)[1],
                                                                    dim(location)[1]))
  } else if (kerneltype == "quadratic") {
    fx_quadratic <- function(i){
      line_i <- rep(0, dim(location)[1])
      line_i[i] <- 1
      ED2 <- pdist::pdist(location[i, ],location[-i, ])@dist^2
      line_i[-i] <- 1 - ED2/(ED2 + (bandwidth))
      ind_i<- which(line_i >= tol)
      return(list("ind_i" = ind_i, "ind_j" = rep(i, length(ind_i)),
                  "val_i" = line_i[ind_i]))
    }
    
    results <- parallel::mclapply(1:dim(location)[1], fx_quadratic, mc.cores = ncores)
    tib <- tidyr::tibble(results)  %>%  tidyr::unnest_wider(results)
    K_sparse <- Matrix::sparseMatrix(i = unlist(tib[[1]]), j = unlist(tib[[2]]), 
                             x = unlist(tib[[3]]), dims = c(dim(location)[1],
                                                            dim(location)[1]))
  }
  return(K_sparse)
}

#' @title Select bandwidth in Gaussian kernel.
#' @description This function selects bandwidth in Gaussian kernel.
#' @param expr A m gene by n location matrix of normalized gene expression matrix.
#' @param method The method used in bandwidth selection, "SJ" usually for small sample size data, 
#' "Silverman" usually for large sample size data.
#' @return A numeric value of calculated bandwidth.
bandwidth_select <- function (expr, method = "Silverman") {
  N <- dim(expr)[2]
  if (method == "SJ") {
    bw_SJ <- c()
    for (i in 1:dim(expr)[1]) {
      tryCatch ({
        bw_SJ[i] <- stats::bw.SJ(expr[i, ], method = "dpi")
      }, error = function(e){cat("Gene", i, " :", conditionMessage(e), "\n")})
    }
    bw <- stats::median(na.omit(bw_SJ))
  } else if (method == "Silverman") {
    bw_Silverman <- c()
    for (i in 1:dim(expr)[1]) {
      tryCatch({
        bw_Silverman[i] <- stats::bw.nrd0(expr[i, ])
      }, error = function(e){cat("Gene", i, " :", conditionMessage(e), "\n")})
    }
    bw <- median(na.omit(bw_Silverman))
  }
  return(bw)
}


#' @title Construct the kernel matrix
#' @description This function construct the kernel matrices in a list.
#' @param location Cell/spot spatial coordinates to compute the kernel matrix.
#' @param basis_number Number of spline functions to construct non parametric covariance matrix
#' @return A list with 11 pre-calculated kernels
Celina_kernel <- function(location, basis_number = 4) {
  ## Calculate the Euclidean distance matrix
  ED <- as.matrix(dist(location[ , 1:2]))
  
  kernels_list <- list()
  lrang <- quantile(abs(location), probs = seq(0.2, 1, by = 0.2))
  for (ikernel in c(1:5)) {
    ## Gaussian kernel
    kernel_mat <- exp(-ED ^ 2/(2 * lrang[ikernel] ^ 2))
    kernels_list <- c(kernels_list, list(kernel_mat))
    rm(kernel_mat)
    
    ## matern kernels
    kernel_mat <- fields::Matern(d = ED, range = 1, alpha = lrang[ikernel], smoothness = 1.5)
    kernels_list <- c(kernels_list, list(kernel_mat))
    rm(kernel_mat)
  }# end for ikernel
  
  ## Add in the spline basis kernels
  center_coords <- sweep(location, 2, apply(location, 2, mean), '-')
  center_coords <- as.data.frame(center_coords / sd(as.matrix(center_coords)))
  colnames(center_coords) <- c("x", "y")
  sm <- mgcv::smoothCon(mgcv::s(x, y, 
    k = basis_number, fx = T, bs = 'tp'), data = center_coords)[[1]]
  mm <- as.matrix(data.frame(sm$X))
  mm_inv <- solve(crossprod(mm, mm))
  spline_kernels <- list(mm %*% mm_inv %*% t(mm))
  
  ## Final output kernels
  output_kernels <- c(kernels_list, spline_kernels)
  return(output_kernels)
}


##########################################################
#       Celina main model related functions             #
##########################################################
#' @title Functions for testing the variance term using multiple kernels
#' @param object The SVC object.
#' @param kernel_mat A list of optional kernel matrix to do the testing, 
#' otherwise used previous calculated kernels
#' @param celltype_to_test Cell types selected for testing
#' @param num_cores The number of cores used in calculation
#' @return A CELINA object with calculated testing results recorded
#' @export
Testing_interaction_all <- function(object, kernel_mat = NULL, 
                                                  celltype_to_test = NULL,
                                                  num_cores) {
  ## Replace the kernel matrices by the user provided kernel_mat
  if (!is.null(kernel_mat)) {
    object@kernelmat <- kernel_mat
  }
  
  object@result <- list()
  
  ## Subset the cell type proportion matrix for remaining cell types
  celltype_mat <- t(object@celltype_mat[names(object@genes_list), ])
  celltype_mat <- sweep(celltype_mat, 1, rowSums(celltype_mat), "/")
  object@celltype_mat <- t(celltype_mat)

  ## Remove the cells from the cell types that have been filtered out for single cell resolution
  
  if (is.null(celltype_to_test)) {
    celltype_to_test <- rownames(object@celltype_mat)
  }
  
  ## Regress cell type information for each target cell type
  for (each_celltype in celltype_to_test) {
    ## Extract the gene expression matrix using selected gene list
    target_normalized_counts <- object@gene_expression_mat[object@genes_list[[each_celltype]], , drop = F]
    
    ## Construct new gene expression matrix after regressing out the remaining cell type information
    target_normalized_counts <- t(apply(target_normalized_counts, 1, function(x) {
      tmp_lm <- lm(x ~ t(object@celltype_mat[-which(rownames(object@celltype_mat) == each_celltype), , drop = F]))
      return(tmp_lm$residuals)
    }))
    
    ## Construct the combination pairs
    gene_names <- rownames(target_normalized_counts)
    combinations <- expand.grid(each_celltype, gene_names, stringsAsFactors = F)

    ## Run different algorithms based on the approximation
    if (object@approximation == FALSE) {
      ## Run the default 11 kernels algorithm
      pvalues_results <- pbmcapply::pbmclapply(1:nrow(combinations), 
                           function(i) {Testing_interaction_multi_kernels(target_normalized_counts[combinations[i, 2], ],
                                                                          object@celltype_mat[combinations[i, 1], ],
                                                                          covariates = object@covariates, kernel_mat = object@kernelmat)},  
                           mc.cores = num_cores)
      
      pvalues_results <- as.data.frame(do.call(rbind, pvalues_results))
      rownames(pvalues_results) <- gene_names
      object@result[[each_celltype]] <- pvalues_results
      rm(pvalues_results)
    } else {
      ## Run the alternative 1 kernel approximation algorithm
      ## Run the Testing_interaction paralleled for each pair
      ## Add in the covariates version
      pvalues_results <- pbmcapply::pbmclapply(1:nrow(combinations), 
                                               function(i) {Testing_interaction(target_normalized_counts[combinations[i, 2], ],
                                                                                object@celltype_mat[combinations[i, 1], ],
                                                                                covariates = object@covariates, kernelmat_approx_U = object@kernelmat_approx_U)},  
                                               mc.cores = num_cores)
      
      pvalues_results <- as.data.frame(do.call(rbind, pvalues_results))
      rownames(pvalues_results) <- gene_names
      object@result[[each_celltype]] <- pvalues_results
      rm(pvalues_results)
    }
  }

  return(object)
}


#' @title Test heterogeneity term using null model for one pair. Used by Testing_interaction_all
#' @param Y Gene expression vector
#' @param X celltype proportion for testing 
#' @param covariates covariates for the null model
#' @param kernelmat_approx_U The calculated approximated kernel matrix for spatial relationship between locations.
#' @export

Testing_interaction <- function(Y, X, covariates, kernelmat_approx_U) {
  ##------------------------
  ## Fit NULL model using REML + AI algorithm
  ##------------------------
  num_spots <- length(Y)
  if (is.null(covariates)) {
    model_mat <- data.frame(y = Y, intercept = rep(1, length(Y)), x = X)
    model0 <- stats::glm(formula = as.numeric(y) ~ . - 1, data = model_mat, family = "gaussian")
  } else{## fit the model with covariates
    model_mat <- cbind(data.frame(y = Y, intercept = rep(1, length(Y))), x = X,
                       covariates)
    model0 <- stats::glm(formula = as.numeric(y) ~ . - 1, data = model_mat, family = "gaussian")
  }# end fi

  ## fit the model using Celina.fit
  if (class(model0)[1] != "try-error"){
    model1 <- try(Celina.fit(model0, ct_info = X, 
                             covariates = as.matrix(model_mat[, -1])))
  } else{
    class(model1) <- "try-error"
  }# end fi
  
  ## Summarize the results
  beta <- model1$coefficients
  var_beta <- model1$cov
  sigma_r2 <- model1$theta[2]
  sigma_e2 <- model1$theta[3]
  converged   <- model1$converged
  
  XUk <- c(X) * kernelmat_approx_U 
  y <- Y
  
  ##-------------------------
  ## calculate P matrix under the null model
  ##-------------------------
  H <- sigma_r2 * (X ^ 2) + sigma_e2
  Hinv <- 1/H
  W <- stats::model.matrix(model0)
  HinvW <- apply(W, 2, function(i){return(i * Hinv)})
  WtHinvW <- t(W) %*% HinvW
  WtHinvW_inv <- solve(WtHinvW)
  P <- diag(Hinv) - HinvW %*% (WtHinvW_inv %*% t(HinvW))
  
  ##-------------------------
  ## calculate test statistics using only two eigen values
  ##-------------------------
  Eigenkb <- eigen(t(XUk) %*% XUk)
  Ekb <- Eigenkb$values
  U <- sqrt(1/2) * (H ^ 0.5) * (P %*% (XUk %*% Eigenkb$vectors))
  Q <- U %*% t(U)
  
  statistics <- as.numeric(crossprod(y * (H ^ -0.5), Q) %*% (y * (H ^ -0.5)))
  e <- eigen(t(U) %*% U)
  lambda <- sort(e$values, decreasing = TRUE)
  
  r1 <- CompQuadForm::davies(statistics, lambda, h = rep(1, length(lambda)),
                             delta = rep(0, length(lambda)), sigma = 0, lim = 10000)
  pvalue_b <- r1$Qq
  if (pvalue_b <= 0) {
    pvalue_b <- CompQuadForm::liu(statistics, lambda)
  }
  
  return(pvalue_b)
}



#' @title Test heterogeneity term using null model for one pair. Used by Testing_interaction_all
#' @param Y Gene expression vector
#' @param X celltype proportion for testing 
#' @param covariates covariates for the null model
#' @param kernel_mat A list of kernel matrix to do the testing
#' either by 'MM' - matching moments or by 'davies' - davies method
#' @export

Testing_interaction_multi_kernels <- function(Y, X, covariates, kernel_mat) {
  ##------------------------
  ## Fit NULL model using REML + AI algorithm
  ##------------------------
  num_spots <- length(Y)
  if (is.null(covariates)) {
    model_mat <- data.frame(y = Y, intercept = rep(1, length(Y)), x = X)
    model0 <- stats::glm(formula = as.numeric(y) ~ . - 1, data = model_mat, family = "gaussian")
  } else{## fit the model with covariates
    model_mat <- cbind(data.frame(y = Y, intercept = rep(1, length(Y))), x = X,
                       covariates)
    model0 <- stats::glm(formula = as.numeric(y) ~ . - 1, data = model_mat, family = "gaussian")
  }# end fi
  
  # kernel_r <- diag(X ^ 2)
  # tmp_relatedness <- list(kernel_r, diag(num_spots))
  ## fit the model using Celina
  if (class(model0)[1] != "try-error") {
    model1 <- try(Celina.fit(model0, ct_info = X, 
                             covariates = as.matrix(model_mat[, -1])))
  } else{
    class(model1) <- "try-error"
  }
  
  ## Summarize the results
  beta <- model1$coefficients
  var_beta <- model1$cov
  sigma_r2 <- model1$theta[2]
  sigma_e2 <- model1$theta[3]
  converged   <- model1$converged
  
  ## Testing each kernel based on the null model
  test_results <- Celina_test(model1, kernel_mat = kernel_mat)

  ## Return the results
  return(test_results)
}





#' @title Testing multiple kernel matrices
#' 
#' @param model1 Celina fitted null models with parameters stored
#' @param kernel_mat A list to store the pre-defined kernel matrix, 
#' @param check_positive Check the kernel matrix is positive or not
Celina_test <- function(model1,
                        kernel_mat = NULL, 
                        check_positive = TRUE) {
  
  ## Check the kernel matrix first
  if (!is.null(kernel_mat) & !is.list(kernel_mat)) { 
    stop("kernel_mat must be a list consisted of different spatial kernel matrices")
  }
  
  res_pval <- NULL
  
  ## H inverse 
  sigma_r2 <- model1$theta[2]
  sigma_e2 <- model1$theta[3]
  H <- sigma_r2 * (model1$X ^ 2) + sigma_e2
  Hinv <- 1/(H + 1e-5)
  HinvW <- apply(model1$W, 2, function(i){return(i * Hinv)})
  WtHinvW <- t(model1$W) %*% HinvW
  WtHinvW_inv <- solve(WtHinvW)
  HinvW_WtHinvW_inv <- (HinvW %*% WtHinvW_inv)
  ## P matrix
  P <- - HinvW_WtHinvW_inv %*% t(HinvW)
  diag(P) <- diag(P) + Hinv
  precal_list <- list(HinvW_WtHinvW_inv, P)
  names(precal_list) <- c("P_right", "P")
  
  for (ikernel in 1:11) {
    test_results <- Celina_test_each_kernel(model1, kernel_mat_each = kernel_mat[[ikernel]],
                                            precal_list = precal_list)
    res_pval <- cbind(res_pval, test_results)
  }
  colnames(res_pval) <- c(paste0("Gaussian", 1:5), paste0("Matern", 1:5), "Spline")
  
  combined_pvalue <- ACAT(res_pval)
  pvals_results <- cbind(res_pval, combined_pvalue)
  colnames(pvals_results)[ncol(pvals_results)] <- "CombinedPvals"
  ## return the results
  return(pvals_results)
}



#' @title Testing one kernel matrix to identify spatial pattern
#' 
#' @param model1 Celina fitted null models with parameters
#' @param kernel_mat_each The single kernel matrix to be tested
#' @param check_positive Check the kernel matrix is positive or not
#' @param precal_list Precalculated identities list to save computational time
Celina_test_each_kernel <- function(model1, kernel_mat_each, 
                                    precal_list = NULL, check_positive = FALSE) {
  ## Check the input kernel matrix to see if its positive definite
  if (check_positive == T) {
    kernel_mat_each <- SysMatEigen(kernel_mat_each)$kernel_mat_each
  }
  
  ## Calculate XKX
  XKX <- (c(model1$X) %*% t(c(model1$X))) * kernel_mat_each
  ## Calculate the quantities based on fitting model and kernel matrix
  teststats_results <- TestStatfast(yin = model1$Y, Pyin = model1$Py, 
                                    Win = model1$W, Xin = model1$X, XKXin = XKX, 
                                    tauin = model1$theta, 
                                    Prightin = precal_list$P_right)
  I_raw <- teststats_results$I_raw
  E <- teststats_results$E
  
  ## calculate the scale parameter and degree of freedom
  k_raw <- I_raw/(2 * E)
  df_raw <- (2 * E ^ 2)/(I_raw)
  pvs_raw <- stats::pchisq(teststats_results$S0/k_raw, df_raw, lower.tail = F)
  
  ## Return the pvalues
  return(pvs_raw)
}

##########################################################
#         Celina optimization related functions                     #
##########################################################
#' Fitting function using the AI algorithm
#' @noRd
Celina.fit <- function(model0, ct_info = NULL, covariates = NULL,
                       maxiter = 500, tol = 1e-5, verbose = FALSE) {
  fixtau.old <- rep(0, 3)
  # to use average information method to fit alternative model
  model1 <- Celina.AI(model0, ct_info = ct_info, covariates = covariates, 
                      maxiter = maxiter, tol = tol, 
                      verbose = verbose)
  fixtau.new <- 1 * (model1$theta < 1.01 * tol)
  
  while (any(fixtau.new != fixtau.old)) {
    fixtau.old <- fixtau.new
    # to use average information method to fit alternative model
    model1 <- Celina.AI(model0, ct_info = ct_info, covariates = covariates, 
                        fixtau = fixtau.old, 
                        maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new <- 1 * (model1$theta < 1.01 * tol)
  } # end while
  return(model1)
} ## end func


#' Fitting function using the AI algorithm - main
#' @noRd
Celina.AI <- function(model0, ct_info, covariates = NULL, tau = rep(0, 3), 
                      fixtau = rep(0, 3), maxiter = 500, 
                      tol = 1e-5, verbose = FALSE) {
  
  ## Extracting the input from model0
  Y <- model0$y 
  num_idv <- length(Y)
  X <- stats::model.matrix(model0) ## Create design matrix(including the intercept)
  alpha <- model0$coef ## Coefficients initial estimates
  
  tau[1] <- 1
  fixtau[1] <- 1 ## 0 stands for to be estimated, 1 for not estimating
  
  numK <- 2 ## Total number of parameters for estimation
  idxtau <- which(fixtau == 0) ## Index for the parameters waiting to be estimated
  numK2 <- sum(fixtau == 0) ## Total number of parameters waiting to be estimated
  
  ## Initializa the variance parameters
  if (numK2 > 0) {
    tau[fixtau == 0] <- rep(min(0.9, stats::var(Y)/(numK + 1)), numK2)
  } 
  
  ## main loop for parameter estimation
  for (iter in seq_len(maxiter)) {
    alpha0 <- alpha
    tau0 <- tau
    model1 <- LmmAICovarites(Yin = Y, Xin = ct_info, 
                             Win = covariates, numVarin = 2, 
                             tauin = tau, fixtauin = fixtau, tol)
    
    tau <- as.numeric(model1$tau)
    alpha <- as.numeric(model1$alpha)
    
    if (2 * max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {
      break
    }
    if (max(tau) > tol^ (-2)) {
      iter <- maxiter
      break
    } 
  } 
  
  converged <- ifelse(iter < maxiter, TRUE, FALSE)
  res <- Y - X %*% alpha
  cov_beta <- model1$cov
  return(list(theta = tau, coefficients = alpha, residuals = res, Y = Y, W = X,
              X = ct_info, converged = converged, cov_beta = cov_beta, Py = model1$Py))
} # end function



##########################################################
#         Celina cauchy combination pvalue calculation function                     #
##########################################################
ACAT <- function(Pvals, Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    Pvals[which(Pvals==0)] <- 5.55e-17
  }
  if (is.one){
    Pvals[which((1-Pvals)<1e-3)] <- 0.99
  }
  
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }
  
  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-stats::pcauchy(cct.stat)
  }
  return(pval)
}