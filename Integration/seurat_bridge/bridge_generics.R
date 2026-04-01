
#' Leverage Score Calculation
#'
#' This function computes the leverage scores for a given object
#' It uses the concept of sketching and random projections. The function provides an approximation 
#' to the leverage scores using a scalable method suitable for large matrices.
#'
#' @param object A matrix-like object
#' @param ... Arguments passed to other methods
#' 
#' @references Clarkson, K. L. & Woodruff, D. P.
#' Low-rank approximation and regression in input sparsity time.
#' JACM 63, 1–45 (2017). \doi{10.1145/3019134};
#'
#' @export
#'
LeverageScore <- function(object, ...) {
  UseMethod(generic = 'LeverageScore', object = object)
}

#' Normalize Raw Data
#'
#' @param data Matrix with the raw count data
#' @param scale.factor Scale the data; default is \code{1e4}
#' @param margin Margin to normalize over
#' @param verbose Print progress
#'
#' @return A matrix with the normalized and log-transformed data
#'
#' @template param-dotsm
#'
#' @param object A Seurat object
#' @param ... Arguments passed to
#' \code{\link[RSpectra:eigs_sym]{RSpectra::eigs_sym}}
#'
#' @return Returns Seurat object with the Graph laplacian eigenvector
#' calculation stored in the reductions slot
#'
#' @rdname RunGraphLaplacian
#' @export RunGraphLaplacian
#'

RunGraphLaplacian <- function(object, ...) {
  UseMethod(generic = 'RunGraphLaplacian', object = object)
}


