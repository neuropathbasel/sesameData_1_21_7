#' infer platform from Probe_IDs
#'
#' @param Probe_IDs probe IDs
#' @param silent suppress message
#' @return a platform code
#' @examples
#' sesameDataCache("probeIDSignature")
#' inferPlatformFromProbeIDs(c("cg14620903","cg22464003"))
#' @export
inferPlatformFromProbeIDs <- function(Probe_IDs, silent = FALSE) {
    sig <- sesameDataGet("probeIDSignature")
    cnts <- vapply(sig, function(x) sum(Probe_IDs %in% x), integer(1))
    if(sum(cnts == max(cnts)) > 1) {
        stop("Ambiguous platform. Please provide platform explicitly.") }

    platform <- names(which.max(cnts))
    if (!silent) {
        message("Platform set to: ", platform)
    }
    platform
}

#' Check platform code
#'
#' Note: custome platforms lead to error here.
#' 
#' @param platform input platform
#' @param probes probes by which the platform may be guessed
#' @return platform code
#' @examples
#' sesameData_check_platform("HM450")
#' @export
sesameData_check_platform <- function(platform = NULL, probes = NULL) {
    if (is.null(platform)) {
        if (is.null(probes)) {
            platform <- "EPIC" # random guess
        } else {
            platform <- inferPlatformFromProbeIDs(probes)
        }
    }
    stopifnot(platform %in% c(
        "MSA", "EPICv2", "EPIC", "HM27", "HM450", "MM285", "Mammal40"))
    platform
}

