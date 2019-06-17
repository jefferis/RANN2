#' WANN objects encapsulating points/kd trees allowing repeated searches
#'
#' @details
#'
#' WANN objects will primarily be useful if you make repeated queries. You can
#' also delay building the k-d tree if you are not sure if the object in
#' question will be searched; the tree will automatically be built when it is
#' queried. You can also explicitly control when the tree is built or deleted
#' (for memory management). The tree is wrapped in an R reference class (R5)
#' object which imposes a significant performance penalty for building small
#' trees (< ~ 1000 points). Summarising, WANN objects provide experimental
#' functionality to:
#'
#' \itemize{
#'
#' \item  keep ANN points in memory to avoid repeated copying
#'
#' \item  keep the ANN k-d tree in memory to avoid repeated building
#'
#' \item  separate building the k-d tree from allocating the points
#'
#' \item  permit very fast self queries
#'
#' \item  permit queries of the points from one ANN tree against a second tree
#'
#' }
#'
#' @description Methods for WANN objects include
#'
#'   \itemize{
#'
#'   \item \code{build_tree}
#'
#'   \item \code{delete_tree}
#'
#'   \item \code{getPoints}
#'
#'   \item \code{querySelf}
#'
#'   \item \code{objectPointer}
#'
#'   \item \code{objectPointer}
#'
#'   }
#' @name WANN-class
#' @aliases WANN
#' @exportClass WANN
#' @export WANN
#' @examples
#' p1=kcpoints[[1]]
#' w1=WANN(p1)
#' # 2 neighbours since the first neighbour will be a self match of each point
#' w1sq=w1$querySelf(2,0)
#' w1sq$nn.dists[,2]
#' w1$query(p1, 2, 0)
#' p2=kcpoints[[2]]
#' w1$query(p2, 1, 0)
#' # equivalent to
#' nn2(p1, p2, k=1)
#'
#' # Explicitly build or delete k-d tree
#' w1$build_tree()
#' w1$delete_tree()
WANN <- setRcppClass("WANN")
