#' @export
#' @concept integration
#'
NNtoGraph <- function(
  nn.object,
  col.cells = NULL,
  weighted = FALSE
) {
  select_nn <- Indices(object = nn.object)
  col.cells <- col.cells %||% Cells(x = nn.object)
  ncol.nn <- length(x = col.cells)
  k.nn <- ncol(x = select_nn)
  j <- as.numeric(x = t(x = select_nn))
  i <- ((1:length(x = j)) - 1) %/% k.nn + 1
  if (weighted) {
    select_nn_dist <- Distances(object = nn.object)
    dist.element <- as.numeric(x = t(x = select_nn_dist))
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = dist.element,
      dims = c(nrow(x = select_nn), ncol.nn)
    )
  } else {
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(nrow(x = select_nn), ncol.nn)
    )
  }
  rownames(x = nn.matrix) <- Cells(x = nn.object)
  colnames(x = nn.matrix) <- col.cells
  nn.matrix <- as.Graph(x = nn.matrix)
  return(nn.matrix)
}


# Find Anchor directly from assay
#
#
# @return Returns a TranserAnchor or Integration set
FindAssayAnchor <- function(
  object.list,
  reference = NULL,
  anchor.type = c("Transfer", "Integration"),
  assay = "Bridge",
  slot = "data",
  reduction =  NULL,
  k.anchor = 20,
  k.score = 50,
  verbose = TRUE
) {
  anchor.type <- match.arg(arg = anchor.type)
  reduction.name <- reduction %||% paste0(assay, ".reduc")
  if ( is.null(x = reduction) || !reduction %in% Reductions(object.list[[1]])) {
    object.list <- lapply(object.list, function(x) {
      if (is.null(reduction)) {
        x[[reduction.name]] <- CreateDimReducObject(
          embeddings = t(GetAssayData(
            object = x,
            layer = slot,
            assay = assay
          )),
          key = "L_",
          assay = assay
        )
      }
    DefaultAssay(x) <- assay
    x <- DietSeurat(x, assays = assay, dimreducs = reduction.name)
    return(x)
    }
  )
}
    object.both <- merge(object.list[[1]], object.list[[2]], merge.dr = reduction.name)
    objects.ncell <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
    offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
    if (verbose) {
      message("Finding ",  anchor.type," anchors from assay ", assay)
    }
    anchors <- FindAnchors(object.pair = object.both,
                           assay = c(DefaultAssay(object.both), DefaultAssay(object.both)),
                           slot = 'data',
                           cells1 = colnames(object.list[[1]]),
                           cells2 = colnames(object.list[[2]]),
                           internal.neighbors = NULL,
                           reduction = reduction.name,
                           k.anchor = k.anchor,
                           k.score = k.score,
                           dims = 1:ncol(object.both[[reduction.name]]),
                           k.filter = NA,
                           verbose = verbose
    )
    inte.anchors <- anchors
    inte.anchors[, 1] <- inte.anchors[, 1] + offsets[1]
    inte.anchors[, 2] <- inte.anchors[, 2] + offsets[2]
    # determine all anchors
    inte.anchors <- rbind(inte.anchors, inte.anchors[, c(2, 1, 3)])
    inte.anchors <- AddDatasetID(
      anchor.df = inte.anchors,
      offsets = offsets,
      obj.lengths = objects.ncell
      )
    command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
    anchor.features <- rownames(object.both)
    if (anchor.type == "Integration") {
      anchor.set <- new(Class = "IntegrationAnchorSet",
                        object.list = object.list,
                        reference.objects = reference %||% seq_along(object.list),
                        anchors = inte.anchors,
                        weight.reduction = object.both[[reduction.name]],
                        offsets = offsets,
                        anchor.features = anchor.features,
                        command = command
      )
    } else if (anchor.type == "Transfer") {
      reference.index <- reference
      reference <- object.list[[reference.index]]
      query  <- object.list[[setdiff(c(1,2), reference.index)]]
      query <- RenameCells(
        object = query,
        new.names = paste0(Cells(x = query), "_", "query")
      )
      reference <- RenameCells(
        object = reference,
        new.names = paste0(Cells(x = reference), "_", "reference")
      )
      combined.ob <- suppressWarnings(expr = merge(
        x = reference,
        y = query,
        merge.dr = reduction.name
      ))
      anchor.set <- new(
        Class = "TransferAnchorSet",
        object.list = list(combined.ob),
        reference.cells = colnames(x = reference),
        query.cells = colnames(x = query),
        anchors = anchors,
        anchor.features = anchor.features,
        command = command
      )
    }
    return(anchor.set)
}


#' Construct a dictionary representation for each unimodal dataset
#'
#'
#' @param object.list A list of Seurat objects
#' @param bridge.object A multi-omic bridge Seurat which is used as the basis to
#' represent unimodal datasets
#' @param object.reduction A list of dimensional reductions from object.list used
#' to be reconstructed by bridge.object
#' @param bridge.reduction A list of dimensional reductions from bridge.object used
#' to reconstruct object.reduction
#' @param laplacian.reduction Name of bridge graph laplacian dimensional reduction
#' @param laplacian.dims Dimensions used for bridge graph laplacian dimensional reduction
#' @param bridge.assay.name Assay name used for bridge object reconstruction value (default is 'Bridge')
#' @param return.all.assays Whether to return all assays in the object.list.
#' Only bridge assay is returned by default.
#' @param l2.norm Whether to l2 normalize the dictionary representation
#' @param verbose Print messages and progress
#'
#' @importFrom MASS ginv
#' @return Returns a object list in which each object has a bridge cell derived assay
#'
#' @export
#' @concept integration
#'
BridgeCellsRepresentation <- function(object.list,
                                      bridge.object,
                                      object.reduction,
                                      bridge.reduction,
                                      laplacian.reduction = 'lap',
                                      laplacian.dims = 1:50,
                                      bridge.assay.name = "Bridge",
                                      return.all.assays = FALSE,
                                      l2.norm = TRUE,
                                      verbose = TRUE
) {
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  if (verbose) {
    message("Constructing Bridge-cells representation")
  }
  single.object = FALSE
  if (length(x = object.list) == 1 &
      inherits(x = object.list, what = 'Seurat')
  ) {
    object.list <- list(object.list)
    single.object = TRUE
  }
  dims.list <- list()
  for (i in 1:length(object.reduction)) {
   ref.dims <- list(
    object= Misc(object.list[[i]][[object.reduction[[i]]]], slot = 'ref.dims'),
    bridge = Misc( bridge.object[[bridge.reduction[[i]]]], slot = 'ref.dims')
   )
   all.dims <- list(
     object = 1:ncol(object.list[[i]][[object.reduction[[i]]]]),
     bridge = 1:ncol( bridge.object[[bridge.reduction[[i]] ]])
     )
   projected.dims.index <- which(sapply(ref.dims, function(x) !is.null(x)))
   if (length(projected.dims.index) == 0) {
     warning('No reference dims found in the dimensional reduction,',
             ' all dims in the dimensional reduction will be used.')
     if (all.dims[[1]] == all.dims[[2]]) {
       dims.list[[i]]  <- all.dims
     } else {
       stop( 'The number of dimensions in the object.list ',
             object.reduction[[i]],
             ' (', length(all.dims[[1]]), ') ',
       ' and the number of dimensions in the bridge object ',
       bridge.reduction[[i]],
       ' (', length(all.dims[[2]]), ') ',
       ' is different.')
     }
   } else {
     reference.dims.index <- setdiff(c(1:2), projected.dims.index)
     dims.list[[i]] <- list()
     dims.list[[i]][[reference.dims.index]] <- ref.dims[[projected.dims.index ]]
     dims.list[[i]][[projected.dims.index]] <- all.dims[[projected.dims.index]]
     names(dims.list[[i]]) <- c('object', 'bridge')
   }
    }
  object.list <- my.lapply(
    X = 1:length(x = object.list),
    FUN = function(x) {
      SA.inv <- ginv(
        X = Embeddings(
          object = bridge.object,
          reduction = bridge.reduction[[x]]
        )[ ,dims.list[[x]]$bridge]
      )
        if (!is.null(laplacian.reduction)) {
          lap.vector <- Embeddings(bridge.object[[laplacian.reduction]])[,laplacian.dims]
          X <- Embeddings(
            object = object.list[[x]],
            reduction = object.reduction[[x]]
          )[, dims.list[[x]]$object] %*% (SA.inv %*% lap.vector)
        } else {
          X <- Embeddings(
            object = object.list[[x]],
            reduction = object.reduction[[x]]
          )[,  dims.list[[x]]$object] %*% SA.inv
          colnames(X) <- Cells(bridge.object)
        }
      if (l2.norm) {
        X <- L2Norm(mat = X, MARGIN = 1)
      }
      colnames(x = X) <- paste0('bridge_',  colnames(x = X))
      suppressWarnings(
        object.list[[x]][[bridge.assay.name]] <- CreateAssayObject(data = t(X))
        )
      object.list[[x]][[bridge.assay.name]]@misc$SA.inv <- SA.inv
      DefaultAssay(object.list[[x]]) <- bridge.assay.name
      VariableFeatures(object = object.list[[x]]) <- rownames(object.list[[x]])
      return (object.list[[x]])
    }
  )
  if (!return.all.assays) {
    object.list <- my.lapply(
      X = object.list,
      FUN = function(x) {
        x <- DietSeurat(object = x, assays = bridge.assay.name, scale.data = TRUE)
        return(x)
      }
    )
  }
  if (single.object) {
    object.list <- object.list[[1]]
  }
  return(object.list)
}

#' Find bridge anchors between two unimodal datasets
#'
#' First, bridge object is used to reconstruct two single-modality profiles and
#' then project those cells into bridage graph laplacian space.
#' Next, find a set of anchors between two single-modality objects. These
#' anchors can later be used to integrate embeddings or transfer data from the reference to
#' query object using the \code{\link{MapQuery}} object.
#'
#'  \itemize{
#'  \item{ Bridge cells reconstruction
#'  }
#'   \item{ Find anchors between objects. It can be either IntegrationAnchors or TransferAnchor.
#'  }
#' }
#'
#' @inheritParams BridgeCellsRepresentation
#' @param anchor.type The type of anchors. Can
#' be one of:
#' \itemize{
#'   \item{Integration: Generate IntegrationAnchors for integration}
#'   \item{Transfer: Generate TransferAnchors for transfering data}
#' }
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration or transfer data.
#' @param reduction Dimensional reduction to perform when finding anchors. Can
#' be one of:
#' \itemize{
#'   \item{cca: Canonical correlation analysis}
#'   \item{direct: Use assay data as a dimensional reduction}
#' }
#' @param reference.bridge.stored If refernece has stored the bridge dictionary representation
#' @param k.anchor How many neighbors (k) to use when picking anchors
#' @param k.score How many neighbors (k) to use when scoring anchors
#' @param verbose Print messages and progress
#' @param ... Additional parameters passed to \code{FindIntegrationAnchors} or
#' \code{FindTransferAnchors}
#'
#'
#' @return Returns an \code{\link{AnchorSet}} object that can be used as input to
#' \code{\link{IntegrateEmbeddings}}.or \code{\link{MapQuery}}
#'
#' @keywords internal
#'
FindBridgeAnchor <- function(object.list,
                             bridge.object,
                             object.reduction,
                             bridge.reduction,
                             anchor.type = c("Transfer", "Integration"),
                             reference = NULL,
                             laplacian.reduction = "lap",
                             laplacian.dims = 1:50,
                             reduction = c("direct", "cca"),
                             bridge.assay.name = "Bridge",
                             reference.bridge.stored = FALSE,
                             k.anchor = 20,
                             k.score = 50,
                             verbose = TRUE,
                             ...
                             ) {
  anchor.type <- match.arg(arg = anchor.type)
  reduction <- match.arg(arg = reduction)
  if (!is.null(laplacian.reduction)) {
    bridge.method <- "bridge graph"
  } else {
    bridge.method <- "bridge cells"
  }
  if (verbose) {
    message("Finding ", anchor.type," anchors")
    switch(
      EXPR = bridge.method,
      "bridge graph" = {
        message('Transform cells to bridge graph laplacian space')
      },
      "bridge cells" = {
        message('Transform cells to bridge cells space')
      }
    )
  }
  reference <- reference %||% c(1)
  query <- setdiff(c(1,2), reference)
  if (anchor.type == "Transfer") {
    stored.bridge.weights <- FALSE
    # check weight matrix
    if (is.null(bridge.object@tools$MapQuery)) {
      warning("No weights stored between reference and bridge obejcts.",
           "Please set store.weights to TRUE in MapQuery")
    } else if (is.null(object.list[[query]]@tools$MapQuery)) {
      warning("No weights stored between query and bridge obejcts.",
           "Please set store.weights to TRUE in MapQuery")
    } else {
      stored.bridge.weights <- TRUE
    }
  }
  if (reference.bridge.stored) {
    object.list[[query]] <- BridgeCellsRepresentation(
      object.list = object.list[[query]] ,
      bridge.object = bridge.object,
      object.reduction = object.reduction[[query]] ,
      bridge.reduction = bridge.reduction[[query]] ,
      bridge.assay.name = bridge.assay.name,
      laplacian.reduction = laplacian.reduction,
      laplacian.dims = laplacian.dims,
      verbose = verbose
    )
  } else {
    object.list <- BridgeCellsRepresentation(
      object.list = object.list ,
      bridge.object = bridge.object,
      object.reduction = object.reduction,
      bridge.reduction = bridge.reduction,
      bridge.assay.name = bridge.assay.name,
      laplacian.reduction = laplacian.reduction,
      laplacian.dims = laplacian.dims,
      verbose = verbose
    )
  }
  if (reduction == "direct") {
    anchor <- FindAssayAnchor(
      object.list = object.list ,
      reference = reference,
      slot = "data",
      anchor.type = anchor.type,
      assay = bridge.assay.name,
      k.anchor = k.anchor,
      k.score = k.score,
      verbose = verbose
    )
  } else if (reduction == "cca") {
    # set data slot to scale.data slot
    object.list <- lapply(
      X = object.list,
      FUN = function(x) {
     x <- SetAssayData(
       object = x,
       layer = "scale.data",
       new.data = as.matrix(
         x = GetAssayData(object = x, layer = "data")
         ))
     return(x)
     }
     )
    anchor <- switch(EXPR = anchor.type,
                     "Integration" = {
                       anchor <- FindIntegrationAnchors(
                         object.list = object.list,
                         k.filter = NA,
                         reference = reference,
                         reduction = "cca",
                         scale = FALSE,
                         k.anchor = k.anchor,
                         k.score = k.score,
                         verbose = verbose,
                         ...)
                       object.merge <- merge(x = object.list[[1]],
                                             y = object.list[2:length(object.list)]
                                             )
                       slot(
                         object = anchor,
                         name = "weight.reduction"
                         ) <- CreateDimReducObject(
                           embeddings = t(GetAssayData(
                             object = object.merge,
                             layer = 'data'
                           )),
                           key = "L_",
                           assay = bridge.assay.name
                         )
                       anchor
                     },
                     "Transfer" = {
                       anchor <-  FindTransferAnchors(
                         reference = object.list[[reference]],
                         query = object.list[[query]],
                         reduction = "cca",
                         scale = FALSE,
                         k.filter = NA,
                         k.anchor = k.anchor,
                         k.score = k.score,
                         verbose = verbose,
                         ...
                       )
                     }
    )
  }
  if (anchor.type == "Transfer") {
    if (stored.bridge.weights) {
      slot( object = anchor,name = "weight.reduction"
      )@misc$bridge.sets <- list(
        bridge.weights =   slot(object = bridge.object,
                                name = "tools"
        )$MapQuery_PrepareBridgeReference$weights.matrix,
        bridge.ref_anchor =  slot(object = bridge.object,
                                  name = "tools"
        )$MapQuery_PrepareBridgeReference$anchor[,1],
        query.weights =  slot(object = object.list[[query]],
                              name = "tools"
        )$MapQuery$weights.matrix,
        query.ref_anchor =  slot(object = object.list[[query]],
                                 name = "tools"
        )$MapQuery$anchor[,1]
      )
    }
  }
  slot(object = anchor, name = "command") <- LogSeuratCommand(
    object = object.list[[1]],
    return.command = TRUE
    )
  return(anchor)
}


# Helper function to transfer labels based on neighbors object
# @param nn.object  the query neighbors object
# @param reference.object the reference seurat object
# @param group.by  A vector of variables to group cells by
# @param weight.matrix A reference x query cell weight matrix
# @return Returns a list for predicted labels, prediction score and matrix
#' @importFrom Matrix sparseMatrix
#' @importFrom fastDummies dummy_cols
#' @importFrom Matrix rowMeans t
#'
TransferLablesNN <- function(
  nn.object = NULL,
  weight.matrix = NULL,
  reference.labels
){
  reference.labels.matrix <- CreateCategoryMatrix(labels = as.character(reference.labels))
  if (!is.null(x = weight.matrix) & !is.null(x = nn.object)) {
    warning('both nn.object and weight matrix are set. Only weight matrix is used for label transfer')
  }
  if (is.null(x = weight.matrix)) {
    select_nn <- Indices(nn.object)
    k.nn <- ncol(select_nn)
    j <- as.numeric(x = t(x = select_nn ))
    i <- ((1:length(x = j)) - 1) %/% k.nn + 1
    nn.matrix <- sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(nrow(select_nn), nrow(reference.labels.matrix))
    )
    rownames(nn.matrix) <- Cells(nn.object)
  } else if (nrow(weight.matrix) == nrow(reference.labels.matrix)) {
    nn.matrix <- t(weight.matrix)
    k.nn <- 1
  } else if (ncol(weight.matrix) == nrow(reference.labels.matrix)) {
    nn.matrix <- weight.matrix
    k.nn <- 1
  } else {
    stop('wrong weights matrix input')
  }
  query.label.mat <- nn.matrix %*% reference.labels.matrix
  query.label.mat <- query.label.mat/k.nn
  prediction.max <- apply(X = query.label.mat, MARGIN = 1, FUN = which.max)

  query.label <- colnames(x = query.label.mat)[prediction.max]
  query.label.score <- apply(X = query.label.mat, MARGIN = 1, FUN = max)
  names(query.label) <- names(query.label.score) <- rownames(query.label.mat)
  if (is.factor(reference.labels)) {
    levels(query.label) <- levels(reference.labels)
  }
  output.list <- list(labels = query.label,
                      scores = query.label.score,
                      prediction.mat = query.label.mat
                      )
  return(output.list)
}

# transfer continuous value based on neighbors
#
TransferExpressionNN<- function(
  nn.object,
  reference.object,
  var.name = NULL
) {
  nn.matrix <- NNtoGraph(nn.object = nn.object,
                         col.cells = Cells(reference.object)
                         )
  reference.exp.matrix <- FetchData(object = reference.object, vars = var.name)
  # remove NA
  reference.exp.matrix <- reference.exp.matrix[complete.cases(reference.exp.matrix), ,drop= F]
  nn.matrix <- nn.matrix[, rownames(reference.exp.matrix)]

  # remove NO neighbor query
  nn.sum <- RowSumSparse(mat = nn.matrix)
  nn.matrix <- nn.matrix[nn.sum > 2, ]
  nn.sum <- nn.sum[nn.sum>2]

  # transfer data
  reference.exp.matrix <- as.matrix(reference.exp.matrix)
  query.exp.mat <- nn.matrix %*% reference.exp.matrix
  query.exp.mat <- sweep(x = query.exp.mat, MARGIN = 1, STATS = nn.sum, FUN = "/")

  # set output for all query cells
  query.exp.all <- data.frame(row.names = Cells(nn.object))
  query.exp.all[rownames(query.exp.mat),1] <- query.exp.mat[,1]
  colnames(query.exp.all) <- var.name
  return(query.exp.all)
}


#' @param reduction.name dimensional reduction name, lap by default
#' @param graph The name of graph
#' @rdname RunGraphLaplacian
#' @concept dimensional_reduction
#' @export
#' @method RunGraphLaplacian Seurat
#'
RunGraphLaplacian.Seurat <- function(
  object,
  graph,
  reduction.name = "lap",
  reduction.key ="LAP_",
  n = 50,
  verbose = TRUE,
  ...
) {
  lap_dir <- RunGraphLaplacian(object = object[[graph]],
                               n = n,
                               reduction.key = reduction.key ,
                               verbose = verbose,
                               ...
                               )
  object[[reduction.name]] <- lap_dir
  return(object)
}



#' @param n Total Number of Eigenvectors to compute and store (50 by default)
#' @param reduction.key dimensional reduction key, specifies the string before
#' the number for the dimension names. LAP by default
#' @param verbose Print message and process
#' @param ... Arguments passed to eigs_sym
#'
#'
#' @concept dimensional_reduction
#' @rdname RunGraphLaplacian
#' @export
#'
#' @importFrom Matrix diag t rowSums
#' @importFrom RSpectra eigs_sym
RunGraphLaplacian.default <- function(object,
                                      n = 50,
                                      reduction.key ="LAP_",
                                      verbose = TRUE,
                                      ...
) {
 if (!all(
   slot(object = t(x = object), name = "x") == slot(object = object, name = "x")
   )) {
   stop("Input graph is not symmetric")
 }
  if (verbose) {
    message("Generating normalized laplacian graph")
  }
  D_half <- sqrt(x = rowSums(x = object))
  L <- -1 * (t(object / D_half) / D_half)
  diag(L) <- 1 + diag(L)
  if (verbose) {
    message("Performing eigendecomposition of the normalized laplacian graph")
  }
  L_eigen <- eigs_sym(L, k = n + 1, which = "SM", ...)
  #delete the first eigen vector
  new_order <- n:1
  lap_output <- list(eigen_vector = Re(L_eigen$vectors[, new_order]),
                         eigen_value = L_eigen$values[new_order]
  )
  rownames(lap_output$eigen_vector) <- colnames(object)
  colnames(lap_output$eigen_vector) <- paste0(reduction.key, 1:n )
  lap_dir <- CreateDimReducObject(embeddings = lap_output$eigen_vector,
                       key = reduction.key,
                       assay = DefaultAssay(object),
                       stdev = lap_output$eigen_value
  )
  return(lap_dir)
}


# Check if the var.name already existed in the meta.data
#
CheckMetaVarName <- function(object, var.name) {
  if (var.name %in% colnames(x = object[[]])) {
    var.name.exist <- var.name
    var.name <- rev(
      x = make.unique(
        names = c(colnames(object[[]]), var.name.exist)
        )
      )[1]
    warning(var.name.exist, " is already existed in the meta.data. ",
            var.name, " will store leverage score value")
  }
  return(var.name)
}



# Run hnsw to find neighbors
#
# @param data Data to build the index with
# @param query A set of data to be queried against data
# @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
# "hamming"
# @param k Number of neighbors
# @param ef_construction  A larger value means a better quality index, but increases build time.
# @param ef Higher values lead to improved recall at the expense of longer search time.
# @param n_threads Maximum number of threads to use.
# @param index optional index object, will be recomputed if not provided
#' @importFrom RcppHNSW hnsw_build hnsw_search
#
HnswNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    k,
                    ef_construction = 200,
                    ef = 10,
                    index = NULL,
                    n_threads = 0
) {
  idx <- index %||% hnsw_build(
    X = data,
    distance = metric,
    ef = ef_construction,
    n_threads = n_threads
    )
  nn <- hnsw_search(
    X = query,
    ann = idx,
    k = k,
    ef = ef,
    n_threads = n_threads
    )
  names(nn) <- c("nn.idx", "nn.dists")
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}


# Calculate reference index from the integrated object
#
IntegrationReferenceIndex <- function(object) {
  if (is.null(object@tools$Integration@sample.tree)) {
    reference.index <- object@commands$FindIntegrationAnchors$reference
    if (length(x = reference.index) > 1) {
      stop('the number of the reference is bigger than 1')
    }
  } else {
    reference.index <- SampleIntegrationOrder(tree = object@tools$Integration@sample.tree)[1]
  }
  return(reference.index)
}


# Calculate mean and sd
#
SparseMeanSd <- function(object,
                         assay = NULL,
                         slot = 'data',
                         features = NULL,
                         eps = 1e-8
){
  assay <- assay%||% DefaultAssay(object)
  features <- features %||% rownames(object[[assay]])
  assay <- assay %||% DefaultAssay(object = object)
  mat <- GetAssayData(object = object[[assay]], layer = slot)[features,]
  if (class(mat)[1] !='dgCMatrix'){
    stop('Matrix is not sparse')
  }
  mat.mean <-  RowMeanSparse(mat)
  mat.sd <-  sqrt(RowVarSparse(mat))
  names(mat.mean) <- names(mat.sd) <- rownames(mat)
  mat.sd <- MinMax(data = mat.sd, min = eps, max = max(mat.sd))
  output <- list(mean = mat.mean, sd = mat.sd)
  return(output)
}



# Run PCA on sparse matrix
#
#' @importFrom Matrix t
#' @importFrom rlang exec
#' @importFrom irlba irlba
#
#
RunPCA_Sparse <- function(
  object,
  features = NULL,
  reduction.key = "PCsp_",
  reduction.name = "pca.sparse",
  npcs = 50,
  do.scale = TRUE,
  verbose = TRUE
) {
  features <- features %||% VariableFeatures(object)
  data <- GetAssayData(object = object, layer = "data")[features,]
  n <- npcs
  args <- list(A = t(data), nv = n)
  args$center <- RowMeanSparse(data)
  feature.var <- RowVarSparse(data)
  args$totalvar <- sum(feature.var)
  if (do.scale) {
    args$scale <- sqrt(feature.var)
    args$scale <- MinMax(args$scale, min = 1e-8, max = max(args$scale))
  } else {
    args$scale <- FALSE
  }
  if (verbose) {
    message("Running PCA")
  }
  pca.irlba <- exec(.fn = irlba, !!!args)
  sdev <- pca.irlba$d/sqrt(max(1, ncol(data) - 1))
  feture.loadings <- pca.irlba$v
  rownames(feture.loadings) <- rownames(data)
  embeddings <- sweep(x = pca.irlba$u, MARGIN = 2, STATS = pca.irlba$d, FUN = "*")
  rownames(embeddings) <- colnames(data)
  colnames(feture.loadings) <- colnames(embeddings) <- paste0(reduction.key, 1:npcs)
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = embeddings,
    loadings = feture.loadings,
    stdev = sdev,
    key = reduction.key,
    assay = DefaultAssay(object),
    misc = list(d = pca.irlba$d)
  )
  return(object)
}

# Smoothing labels based on the clusters
# @param labels the original labels
# @param clusters the clusters that are used to smooth labels
#
SmoothLabels <- function(labels, clusters) {
  cluster.set <- unique(clusters)
  smooth.labels <- labels
  for (c in cluster.set) {
    cell.c <- which(clusters == c)
    smooth.labels[cell.c] <- names(sort(table(labels[cell.c]), decreasing = T)[1])
  }
  return(smooth.labels)
}



#' Project query data to reference dimensional reduction
#'
#' @param query Query object
#' @param reference Reference object
#' @param mode Projection mode name for projection
#'  \itemize{
#' \item{pcaproject: PCA projection}
#' \item{lsiproject: LSI projection}
#' }
#' @param reference.reduction Name of dimensional reduction in the reference object
#' @param combine Determine if query and reference objects are combined
#' @param query.assay Assay used for query object
#' @param reference.assay Assay used for reference object
#' @param features Features used for projection
#' @param do.scale Determine if scale expression matrix in the pcaproject mode
#' @param reduction.name dimensional reduction name, reference.reduction is used by default
#' @param reduction.key dimensional reduction key, the key in reference.reduction
#' is used by default
#' @param verbose Print progress and message
#'
#' @return Returns a query-only or query-reference combined seurat object
#'
#' @export
#' @concept integration
#'
ProjectDimReduc <- function(query,
                            reference,
                            mode = c('pcaproject', 'lsiproject'),
                            reference.reduction,
                            combine = FALSE,
                            query.assay = NULL,
                            reference.assay = NULL,
                            features = NULL,
                            do.scale = TRUE,
                            reduction.name = NULL,
                            reduction.key= NULL,
                            verbose = TRUE
) {
  query.assay <- query.assay %||% DefaultAssay(object = query)
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  DefaultAssay(object = query) <- query.assay
  DefaultAssay(object = reference) <- reference.assay
  reduction.name <- reduction.name %||% reference.reduction
  reduction.key <- reduction.key %||% Key(object = reference[[reference.reduction]])
  if (reduction.name %in% Reductions(object = query)) {
    warning(reduction.name,
            ' already exists in the query object. It will be overwritten.'
    )
  }
  features <- features %||% rownames(x = Loadings(object = reference[[reference.reduction]]))
  features <- intersect(x = features, y = rownames(x = query))
  if (mode == 'lsiproject') {
    if (verbose) {
      message('LSI projection to ', reference.reduction)
    }
    projected.embeddings <- ProjectSVD(
      reduction = reference[[reference.reduction]],
      data = GetAssayData(object = query, assay = query.assay, layer = "data"),
      mode = "lsi",
      do.center = FALSE,
      do.scale = FALSE,
      features = features,
      use.original.stats = FALSE,
      verbose = verbose
    )
  } else if (mode == 'pcaproject') {
    if (inherits(query[[query.assay]], what = 'SCTAssay')) {
      if (verbose) {
        message('PCA projection to ', reference.reduction, ' in SCT assay')
      }
      query <- suppressWarnings(
        expr = GetResidual(object = query,
                           assay = query.assay,
                           features = features,
                           verbose = FALSE)
      )
      query.mat <- GetAssayData(object = query, layer = 'scale.data')[features,]

      projected.embeddings <- t(
        crossprod(x = Loadings(
          object = reference[[reference.reduction]])[features, ],
          y = query.mat
        )
      )
    } else {
      if (verbose) {
        message('PCA projection to ', reference.reduction)
      }
      projected.embeddings <- ProjectCellEmbeddings(
        reference = reference,
        reduction = reference.reduction,
        query = query,
        scale = do.scale,
        dims = 1:ncol(reference[[reference.reduction]]),
        feature.mean = NULL,
        verbose = verbose
      )
    }
  }
  query[[reduction.name]] <- CreateDimReducObject(
    embeddings = projected.embeddings,
    loadings = Loadings(reference[[reference.reduction]])[features,],
    assay = query.assay,
    key = reduction.key,
    misc = Misc(reference[[reference.reduction]])
  )
  if (combine) {
    query <- DietSeurat(object = query,
                        dimreducs = reduction.name,
                        features = features,
                        assays = query.assay
    )
    reference <- DietSeurat(object = reference,
                            dimreducs = reference.reduction,
                            features = features,
                            assays = reference.assay)
    suppressWarnings(
      combine.obj <- merge(query, reference,
                           merge.dr = c(reduction.name, reference.reduction)
      )
    )
    Idents(combine.obj) <- c(rep(x = 'query', times = ncol(query)),
                            rep(x = 'reference', times = ncol(reference))
                            )
    return(combine.obj)
  } else {
    return(query)
  }
}


#' Prepare the bridge and reference datasets
#'
#' Preprocess the multi-omic bridge and unimodal reference datasets into
#' an extended reference.
#' This function performs the following three steps:
#' 1. Performs within-modality harmonization between bridge and reference
#' 2. Performs dimensional reduction on the SNN graph of bridge datasets via
#' Laplacian Eigendecomposition
#' 3. Constructs a bridge dictionary representation for unimodal reference cells
#'
#' @param reference A reference Seurat object
#' @param bridge A multi-omic bridge Seurat object
#' @param reference.reduction Name of dimensional reduction of the reference object (default is 'pca')
#' @param reference.dims Number of dimensions used for the reference.reduction (default is 50)
#' @param normalization.method Name of normalization method used: LogNormalize
#' or SCT
#' @param reference.assay Assay name for reference (default is \code{\link{DefaultAssay}})
#' @param bridge.ref.assay Assay name for bridge used for reference mapping. RNA by default
#' @param bridge.query.assay Assay name for bridge used for query mapping. ATAC by default
#' @param supervised.reduction Type of supervised dimensional reduction to be performed
#' for integrating the bridge and query.
#' Options are:
#' \itemize{
#'    \item{slsi: Perform supervised LSI as the dimensional reduction for
#'    the bridge-query integration}
#'    \item{spca: Perform supervised PCA as the dimensional reduction for
#'    the bridge-query integration}
#'    \item{NULL: no supervised dimensional reduction will be calculated.
#'    bridge.query.reduction is used for the bridge-query integration}
#' }
#' @param bridge.query.reduction Name of dimensions used for the bridge-query harmonization.
#' 'bridge.query.reduction' and 'supervised.reduction' cannot be NULL together.
#' @param bridge.query.features Features used for bridge query dimensional reduction
#' (default is NULL which uses VariableFeatures from the bridge object)
#' @param laplacian.reduction.name Name of dimensional reduction name of graph laplacian eigenspace (default is 'lap')
#' @param laplacian.reduction.key Dimensional reduction key (default is 'lap_')
#' @param laplacian.reduction.dims Number of dimensions used for graph laplacian eigenspace (default is 50)
#' @param verbose Print progress and message (default is TRUE)
#'
#' @return Returns a \code{BridgeReferenceSet} that can be used as input to
#'  \code{\link{FindBridgeTransferAnchors}}.
#' The parameters used are stored in the \code{BridgeReferenceSet} as well
#'
#' @export
#' @concept integration
#'
PrepareBridgeReference <- function (
  reference,
  bridge,
  reference.reduction = 'pca',
  reference.dims = 1:50,
  normalization.method = c('SCT', 'LogNormalize'),
  reference.assay = NULL,
  bridge.ref.assay = 'RNA',
  bridge.query.assay = 'ATAC',
  supervised.reduction = c('slsi', 'spca', NULL),
  bridge.query.reduction = NULL,
  bridge.query.features = NULL,
  laplacian.reduction.name = 'lap',
  laplacian.reduction.key = 'lap_',
  laplacian.reduction.dims = 1:50,
  verbose = TRUE
) {
  ## checking
  if (!is.null(supervised.reduction)) {
  supervised.reduction <- match.arg(arg = supervised.reduction)
  }
  if (!is.null(x = bridge.query.reduction) & !is.null(x = supervised.reduction)) {
    stop('bridge.query.reduction and supervised.reduction can only set one.',
         'If you want to set bridge.query.reduction, supervised.reduction should set to NULL')
  }
  if (is.null(x = bridge.query.reduction) & is.null(x = supervised.reduction)) {
    stop('Both bridge.query.reduction and supervised.reduction are NULL. One of them needs to be set')
  }
  bridge.query.features <- bridge.query.features %||%
    VariableFeatures(object = bridge[[bridge.query.assay]])
  if (length(x = bridge.query.features) == 0) {
    stop('bridge object ', bridge.query.assay,
         ' assay has no variable genes and bridge.query.features has no input')
  }
  # modality harmonization
  reference.assay <- reference.assay %||% DefaultAssay(reference)
  DefaultAssay(reference) <- reference.assay
  DefaultAssay(bridge) <- bridge.ref.assay
  ref.anchor  <- FindTransferAnchors(
    reference =  reference,
    reference.reduction = reference.reduction,
    normalization.method = normalization.method,
    dims = reference.dims,
    query = bridge,
    recompute.residuals = TRUE,
    features = rownames(reference[[reference.reduction]]@feature.loadings),
    k.filter = NA,
    verbose = verbose
  )
  bridge <- MapQuery(anchorset = ref.anchor,
                     reference = reference,
                     query = bridge,
                     store.weights = TRUE,
                     verbose = verbose
  )
  bridge.ref.reduction <- paste0('ref.', reference.reduction)
  bridge <- FindNeighbors(object = bridge,
                          reduction = bridge.ref.reduction,
                          dims = 1:ncol(x = bridge[[bridge.ref.reduction]]),
                          return.neighbor = FALSE,
                          graph.name = c('bridge.ref.nn', 'bridge.ref.snn'),
                          prune.SNN = 0)
  bridge <- RunGraphLaplacian(object = bridge,
                              graph = "bridge.ref.snn",
                              reduction.name = laplacian.reduction.name,
                              reduction.key = laplacian.reduction.key,
                              verbose = verbose)
  DefaultAssay(object = bridge) <- bridge.query.assay
  if (!is.null(supervised.reduction)) {
    bridge <- switch(EXPR = supervised.reduction,
                     'slsi' = {
                       bridge.reduc <- RunSLSI(object = bridge,
                                               features = VariableFeatures(bridge),
                                               graph = 'bridge.ref.nn',
                                               assay = bridge.query.assay
                       )
                       bridge.reduc
                     },
                     'spca' = {
                       bridge.reduc <- RunSPCA(object = bridge,
                                               features = VariableFeatures(bridge),
                                               graph = 'bridge.ref.snn',
                                               assay = bridge.query.assay
                       )
                       bridge.reduc
                     }
    )
  }
  # bridge representation
  reference.bridge <- BridgeCellsRepresentation(
    object.list =  reference,
    bridge.object = bridge,
    object.reduction = c(reference.reduction),
    bridge.reduction =  c(bridge.ref.reduction),
    laplacian.reduction = laplacian.reduction.name,
    laplacian.dims = laplacian.reduction.dims
  )
  reference[['Bridge']] <- reference.bridge[['Bridge']]
  reference <- merge(x = reference, y = bridge, merge.dr = NA)
  reference@tools$MapQuery_PrepareBridgeReference <- bridge@tools$MapQuery
  command <- LogSeuratCommand(object = reference, return.command = TRUE)
  slot(object = command, name = "params")$bridge.query.features <- NULL
  command.name <- slot(object = command, name = "name")
  reference[[command.name]] <- command
  return(reference)
}


#' Find bridge anchors between query and extended bridge-reference
#'
#' Find a set of anchors between unimodal query and the other unimodal reference
#' using a pre-computed \code{\link{BridgeReferenceSet}}.
#' This function performs three steps:
#' 1. Harmonize the bridge and query cells in the bridge query reduction space
#' 2. Construct the bridge dictionary representations for query cells
#' 3. Find a set of anchors between query and reference in the bridge graph laplacian eigenspace
#' These anchors can later be used to integrate embeddings or transfer data from the reference to
#' query object using the \code{\link{MapQuery}} object.

#' @param extended.reference BridgeReferenceSet object generated from
#'  \code{\link{PrepareBridgeReference}}
#' @param query A query Seurat object
#' @param query.assay Assay name for query-bridge integration
#' @param scale Determine if scale the query data for projection
#' @param dims Number of dimensions for query-bridge integration
#' @param reduction Dimensional reduction to perform when finding anchors.
#' Options are:
#' \itemize{
#'    \item{pcaproject: Project the PCA from the bridge onto the query. We
#'    recommend using PCA when bridge and query datasets are from scRNA-seq}
#'    \item{lsiproject: Project the LSI from the bridge onto the query. We
#'    recommend using LSI when bridge and query datasets are from scATAC-seq or scCUT&TAG data.
#'    This requires that LSI or supervised LSI has been computed for the bridge dataset, and the
#'    same features (eg, peaks or genome bins) are present in both the bridge
#'    and query.
#' }
#' }
#' @param bridge.reduction Dimensional reduction to perform when finding anchors. Can
#' be one of:
#' \itemize{
#'   \item{cca: Canonical correlation analysis}
#'   \item{direct: Use assay data as a dimensional reduction}
#' }
#' @param verbose Print messages and progress
#'
#' @return Returns an \code{AnchorSet} object that can be used as input to
#' \code{\link{TransferData}}, \code{\link{IntegrateEmbeddings}} and
#' \code{\link{MapQuery}}.
#'
#' @export
#' @concept integration
#'
FindBridgeTransferAnchors <- function(
  extended.reference,
  query,
  query.assay = NULL,
  dims = 1:30,
  scale = FALSE,
  reduction = c('lsiproject', 'pcaproject'),
  bridge.reduction = c('direct', 'cca'),
  verbose = TRUE
) {
  bridge.reduction <- match.arg(arg = bridge.reduction)
  reduction <-  match.arg(arg = reduction)
  query.assay <- query.assay %||% DefaultAssay(query)
  DefaultAssay(query) <- query.assay
  command.name <- grep(pattern = 'PrepareBridgeReference',
       x = names(slot(object = extended.reference, name = 'commands')),
       value = TRUE)
  params <- Command(object = extended.reference, command = command.name)
  bridge.query.assay <- params$bridge.query.assay
  bridge.query.reduction <- params$bridge.query.reduction %||% params$supervised.reduction
  reference.reduction <- params$reference.reduction
  bridge.ref.reduction <-  paste0('ref.', reference.reduction)
  DefaultAssay(extended.reference) <- bridge.query.assay
  extended.reference.bridge <- DietSeurat(
    object = extended.reference,
    assays = bridge.query.assay,
    dimreducs = c(bridge.ref.reduction, bridge.query.reduction, params$laplacian.reduction.name)
    )
    query.anchor <- FindTransferAnchors(
      reference = extended.reference.bridge,
      reference.reduction = bridge.query.reduction,
      dims = dims,
      query = query,
      reduction = reduction,
      scale = scale,
      features = rownames(Loadings(extended.reference[[bridge.query.reduction]])),
      k.filter = NA,
      verbose = verbose
    )

  query <- MapQuery(anchorset =  query.anchor,
                    reference = extended.reference.bridge,
                    query = query,
                    store.weights = TRUE
  )
  DefaultAssay(extended.reference) <- 'Bridge'
  bridge_anchor  <- FindBridgeAnchor(
    object.list = list(DietSeurat(object = extended.reference, assays = 'Bridge'), query),
    bridge.object = extended.reference.bridge,
    object.reduction = c(reference.reduction, paste0('ref.', bridge.query.reduction)),
    bridge.reduction = c(bridge.ref.reduction, bridge.query.reduction),
    anchor.type = "Transfer",
    reduction = bridge.reduction,
    reference.bridge.stored = TRUE,
    verbose = verbose
  )
  return(bridge_anchor)
}



#' Find integration bridge anchors between query and extended bridge-reference
#'
#' Find a set of anchors between unimodal query and the other unimodal reference
#' using a pre-computed \code{\link{BridgeReferenceSet}}.
#' These integration anchors can later be used to integrate query and reference
#' using the \code{\link{IntegrateEmbeddings}} object.
#'
#' @inheritParams FindBridgeTransferAnchors
#' @param integration.reduction Dimensional reduction to perform when finding anchors
#' between query and reference.
#' Options are:
#' \itemize{
#'    \item{direct: find anchors directly on the bridge representation space}
#'    \item{cca: perform cca on the on the bridge representation space and then find anchors
#' }
#' }
#'
#' @return Returns an \code{AnchorSet} object that can be used as input to
#' \code{\link{IntegrateEmbeddings}}.
#'
#' @export
#' @concept integration
#'
FindBridgeIntegrationAnchors <- function(
  extended.reference,
  query,
  query.assay = NULL,
  dims = 1:30,
  scale = FALSE,
  reduction = c('lsiproject', 'pcaproject'),
  integration.reduction = c('direct', 'cca'),
  verbose = TRUE
) {
  reduction <-  match.arg(arg = reduction)
  integration.reduction <-  match.arg(arg = integration.reduction)
  query.assay <- query.assay %||% DefaultAssay(query)
  DefaultAssay(query) <- query.assay
  command.name <- grep(pattern = 'PrepareBridgeReference',
                       x = names(slot(object = extended.reference, name = 'commands')),
                       value = TRUE)
  params <- Command(object = extended.reference, command = command.name)
  bridge.query.assay <- params$bridge.query.assay
  bridge.query.reduction <- params$bridge.query.reduction %||% params$supervised.reduction
  reference.reduction <- params$reference.reduction
  bridge.ref.reduction <- paste0( 'ref.', params$bridge.ref.reduction)
  DefaultAssay(extended.reference) <- bridge.query.assay

  extended.reference.bridge <- DietSeurat(
    object = extended.reference,
    assays = bridge.query.assay,
    dimreducs = c(bridge.query.reduction, bridge.ref.reduction, params$laplacian.reduction.name)
    )

  query.anchor <- FindTransferAnchors(
    reference = extended.reference.bridge,
    reference.reduction = bridge.query.reduction,
    dims = dims,
    query = query,
    reduction = reduction,
    scale = scale,
    features = rownames(Loadings(extended.reference.bridge[[bridge.query.reduction]])),
    k.filter = NA,
    verbose = verbose
  )
  query <- MapQuery(anchorset =  query.anchor,
                    reference = extended.reference.bridge,
                    query = query,
                    store.weights = TRUE
  )
  DefaultAssay(extended.reference) <- 'Bridge'
  bridge_anchor  <- FindBridgeAnchor(
    object.list = list(DietSeurat(object = extended.reference, assays = 'Bridge'), query),
    bridge.object = extended.reference.bridge,
    reduction = integration.reduction,
    object.reduction = c(reference.reduction, paste0('ref.', bridge.query.reduction)),
    bridge.reduction = c(bridge.ref.reduction, bridge.query.reduction),
    anchor.type = "Integration",
    reference.bridge.stored = TRUE,
    verbose = verbose
  )
  return(bridge_anchor)
}


#' Perform integration on the joint PCA cell embeddings.
#'
#' This is a convenience wrapper function around the following three functions
#' that are often run together when perform integration.
#' \code{\link{FindIntegrationAnchors}}, \code{\link{RunPCA}},
#' \code{\link{IntegrateEmbeddings}}.
#'
#' @inheritParams FindIntegrationAnchors
#' @param new.reduction.name Name of integrated dimensional reduction
#' @param npcs Total Number of PCs to compute and store (50 by default)
#' @param findintegrationanchors.args A named list of additional arguments to
#' \code{\link{FindIntegrationAnchors}}
#' @param verbose Print messages and progress
#'
#' @importFrom rlang exec
#' @return Returns a Seurat object with integrated dimensional reduction
#' @export
#' @concept integration
#'
FastRPCAIntegration <- function(
  object.list,
  reference = NULL,
  anchor.features = 2000,
  k.anchor = 20,
  dims = 1:30,
  scale = TRUE,
  normalization.method = c("LogNormalize", "SCT"),
  new.reduction.name = 'integrated_dr',
  npcs = 50,
  findintegrationanchors.args = list(),
  verbose = TRUE
) {
  npcs <- max(npcs, dims)
  my.lapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pblapply,
    no = future_lapply
  )
  reduction <- 'rpca'
  if (is.numeric(x = anchor.features)) {
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = anchor.features,
      verbose = FALSE
    )
  }
  if (normalization.method == 'SCT') {
    scale <- FALSE
    object.list <- PrepSCTIntegration(object.list = object.list,
                                      anchor.features = anchor.features
    )
  }

  if (verbose) {
    message('Performing PCA for each object')
  }
  object.list <- my.lapply(X = object.list,
                           FUN = function(x) {
                             if (normalization.method != 'SCT') {
                               x <- ScaleData(x, features = anchor.features, do.scale = scale, verbose = FALSE)
                             }
                             x <- RunPCA(x, features = anchor.features, verbose = FALSE, npcs = npcs)
                             return(x)
                           }
  )
  fia.allarguments <- c(list(
    object.list = object.list,
    reference = reference,
    anchor.features = anchor.features,
    reduction = reduction,
    normalization.method = normalization.method,
    scale = scale,
    k.anchor = k.anchor,
    dims = dims,
    verbose = verbose
    ), findintegrationanchors.args
  )
  anchor <- exec("FindIntegrationAnchors",!!!fia.allarguments)
  object_merged <- merge(x = object.list[[1]],
                         y = object.list[2:length(object.list)]

  )

  anchor.feature <- slot(object = anchor, name = 'anchor.features')
  if (normalization.method != 'SCT') {
    object_merged <- ScaleData(object = object_merged,
                               features = anchor.feature,
                               do.scale = scale,
                               verbose = FALSE
    )
  }
  object_merged <- RunPCA(object_merged,
                          features = anchor.feature,
                          verbose = FALSE,
                          npcs = npcs

  )

  temp <- object_merged[["pca"]]
  object_merged <- IntegrateEmbeddings(
    anchorset = anchor,
    reductions = object_merged[['pca']],
    new.reduction.name = new.reduction.name,
    verbose = verbose)
  object_merged[['pca']] <- temp
  VariableFeatures(object = object_merged) <- anchor.feature
  return(object_merged)

}


#' Transfer embeddings from sketched cells to the full data
#'
#' @param atom.data Atom data
#' @param atom.cells Atom cells
#' @param orig.data Original data
#' @param embeddings Embeddings of atom cells
#' @param sketch.matrix Sketch matrix
#'
#' @importFrom MASS ginv
#' @importFrom Matrix t
#'
#' @export
#' @concept integration
#'
UnSketchEmbeddings <- function(
  atom.data,
  atom.cells = NULL,
  orig.data,
  embeddings,
  sketch.matrix = NULL
) {
  if(!all(rownames(atom.data) == rownames(orig.data))) {
    stop('features in atom.data and orig.data are not identical')
  } else {
    features = rownames(atom.data)
  }
  atom.cells <- atom.cells %||% colnames(x = atom.data)
  if (inherits(x = orig.data, what = 'DelayedMatrix') ) {
    stop("PseudobulkExpression does not support DelayedMatrix objects")
  } else if(inherits(x = orig.data, what = 'IterableMatrix')) {
    matrix.prod.function <- crossprod_BPCells
  } else {
    matrix.prod.function <- crossprod
  }
  sketch.matrix <- sketch.matrix %||% as.sparse(diag(length(features)))
  atom.data <- atom.data[, atom.cells]
  embeddings <- embeddings[atom.cells,]
  exp.mat <- as.matrix(x = t(x = atom.data) %*% sketch.matrix)
  sketch.transform <- ginv(X = exp.mat) %*% embeddings
  emb <- matrix.prod.function(
    x = as.matrix(sketch.matrix %*% sketch.transform),
    y = orig.data
  )
  emb <- as.matrix(x = emb)
  return(emb)
}

FeatureSketch <- function(features, ratio = 0.8, seed = 123) {
  sketch.R <- t(x = CountSketch(
    nsketch = round(x = ratio *  length(x = features)),
    ncells = length(x = features),
    seed = seed)
  )
  return(sketch.R)
}
