
## fall back data retrieval in case ExperimentHub is down
.sesameDataGet_fallback <- function(title) {
    u1 <- sprintf('%s/sesameData/%s.rda', alt_base, title)
    if (valid_url(u1)) {
        sesameDataGet_assignEnv(title, get(load(url(u1))))
        TRUE
    } else {
        warning(sprintf("Resource %s cannot be retrieved.", title))
        FALSE
    }
    TRUE
}

.sesameDataGet <- function(title) {

    if ((!is.null(options("SESAMEDATA_USE_ALT")[[1]])) &&
        options("SESAMEDATA_USE_ALT")[[1]]) {
        if (!exists(title, envir=cacheEnv, inherits=FALSE)) {
            .sesameDataGet_fallback(title)
        }
        return(get(title, envir=cacheEnv, inherits=FALSE))
    } else {
        eh_id <- df_master$EHID[match(title, df_master$Title)]
        if (eh_id %in% c("TBD", "NA")) { stopAndCache(title); }
        stopifnot(length(eh_id) == 1)
        
        if (exists(eh_id, envir=cacheEnv, inherits=FALSE)) {
            return(get(eh_id, envir=cacheEnv, inherits=FALSE))
        }

        if (!file.exists(getExperimentHubOption("CACHE"))) {
            stopAndCache(title) }
        tryCatch({
            eh <- query(ExperimentHub(localHub=TRUE), 'sesameData')
        }, error = function(cond) { stopAndCache(title); })
        if (!(eh_id %in% names(eh))) {
            stopAndCache(title); }
        sesameDataGet_assignEnv(eh_id, eh[[eh_id]])
        return(get(eh_id, envir=cacheEnv, inherits=FALSE))
    }
}

#' Get SeSAMe data
#'
#' @param title title of the data
#' @param verbose whether to output ExperimentHub message
#' @return data object
#' @import ExperimentHub
#' @import AnnotationHub
#' @importFrom stringr str_replace
#' @examples
#'
#' sesameDataCache("EPIC.1.SigDF")
#' EPIC.1.SigDF <- sesameDataGet('EPIC.1.SigDF')
#' @export
sesameDataGet <- function(title, verbose = FALSE) {

    title <- str_replace(title, "MMB", "MM285") # fix potential code discrepancy
    
    if (verbose) {
        .sesameDataGet(title)
    } else {
        suppressMessages(
            log <- capture.output(obj <- .sesameDataGet(title)))
        obj
    }
}

#' List all SeSAMe data
#'
#' @param filter keyword to filter title, optional
#' @param full whether to display all columns
#' @return all titles from SeSAMe Data
#' @examples
#' sesameDataList("KYCG")
#' @export
sesameDataList <- function(filter = NULL, full = FALSE) {
    df <- df_master
    if (!full) {
        df <- df[,c("EHID","Title")]
    }
    if (!is.null(filter)) {
        df <- df[grep(filter, df$Title),]
    }
    df
}

#' Whether sesameData has
#'
#' @param data_titles data titles to check
#' @return a boolean vector the same length as data_titles
#' @examples
#' sesameDataHas(c("EPIC.address","EPIC.address.Nonexist"))
#' @export
sesameDataHas <- function(data_titles) {
    data_titles %in% sesameDataList()$Title
}

