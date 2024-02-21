#' get Infinium manifest GRanges
#'
#' Note that some unaligned probes are not included.
#' For full manifest, please visit
#' \url{http://zwdzwd.github.io/InfiniumAnnotation}
#'
#' @param platform Mammal40, MM285, EPIC, and HM450
#' @param genome hg38, mm10, ... will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @return GRanges
#' @examples
#' sesameDataCache("Mammal40.address")
#' res <- sesameData_getManifestGRanges("Mammal40")
#' @export
sesameData_getManifestGRanges <- function(
    platform, genome = NULL) {

    if ("GenomicRanges" %in% class(genome)) { # if is a GRanges object already
        return(genome)
    }

    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)
    ## only one genome is supported natively.
    addr <- sesameDataGet(sprintf("%s.address", platform))
    if (genome %in% names(addr)) {
        return(addr[[genome]])
    } else {
        return(NULL)
    }
}

