#' @slot reference The Reference object only containing bridge representation assay
#' @slot params A list of parameters used in the PrepareBridgeReference
#' @slot command Store log of parameters that were used
#'
#' @name BridgeReferenceSet-class
#' @rdname BridgeReferenceSet-class
#' @concept objects
#' @exportClass BridgeReferenceSet
#'
BridgeReferenceSet <- setClass(
  Class = "BridgeReferenceSet",
  slots = list(
    bridge = "ANY",
    reference = "ANY",
    params = "list",
    command = "ANY"
  )
)

#' The IntegrationData Class
#'
#' The IntegrationData object is an intermediate storage container used internally throughout the
#' integration procedure to hold bits of data that are useful downstream.
#'
#' @slot neighbors List of neighborhood information for cells (outputs of \code{RANN::nn2})
#' @slot weights Anchor weight matrix
#' @slot integration.matrix Integration matrix
#' @slot anchors Anchor matrix
#' @slot offsets The offsets used to enable cell look up in downstream functions
#' @slot objects.ncell Number of cells in each object in the object.list
#' @slot sample.tree Sample tree used for ordering multi-dataset integration
#'
#' @name IntegrationData-class
#' @rdname IntegrationData-class
#' @concept objects
#' @exportClass IntegrationData
#'
IntegrationData <- setClass(
  Class = "IntegrationData",
  slots = list(
    neighbors = "ANY",
)

setMethod(
  f = 'show',
  signature = 'BridgeReferenceSet',
  definition = function(object) {
    cat(
      'A BridgeReferenceSet object has a bridge object with ',
      ncol(slot(object = object, name = 'bridge')),
      'cells and a reference object with ',
      ncol(slot(object = object, name = 'reference')),
      'cells. \n','The bridge query reduction is ',
      slot(object = object, name = 'params')$bridge.query.reduction %||%
        slot(object = object, name = 'params')$supervised.reduction,
   "\n This can be used as input to FindBridgeTransferAnchors and FindBridgeIntegrationAnchors")
  }
)

setMethod(
  f = 'show',
  signature = 'SCTModel',
  definition = function(object) {
    cat(
