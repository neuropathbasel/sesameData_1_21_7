sesameDataRecache <- function(eh_id) {
    stopifnot(startsWith(eh_id, "EH"))
    tryCatch({
        sesameDataCache0(eh_id)
    }, error=function(cond) {
        stopAndCache(eh_id)
    })
    query(ExperimentHub(localHub=TRUE), 'sesameData')
}

stopAndCache <- function(title) {
    stop(sprintf('
| File %s either not found or needs to be cached to be
| used in sesame.
| Please make sure you have updated ExperimentHub and try
| > sesameDataCache("%s")
| or download all data
| > sesameDataCache()
| to retrieve and cache needed sesame data.', title, title))
}

sesameDataCache0 <- function(eh_ids) {
    ## load meta data
    message(sprintf("Metadata (N=%d):\n", length(eh_ids)))
    suppressMessages(log <- capture.output(
        eh <- query(ExperimentHub(), "sesameData")[eh_ids]))
    
    ## load actual data
    tmp2 <- lapply(seq_along(eh), function(i) {
        message(sprintf(
            "(%d/%d) %s:\n", i, length(eh_ids), eh_ids[i]))
        suppressMessages(log <- capture.output(cache(eh[i])))
    })
}

#' Cache SeSAMe data
#'
#' @import ExperimentHub
#' @import AnnotationHub
#' @param data_titles data to cache, if not given will cache all
#' @return TRUE
#' @examples
#' sesameDataCache("genomeInfo.hg38")
#' @export
sesameDataCache <- function(data_titles = NULL) {
    setExperimentHubOption(arg="MAX_DOWNLOADS", 100)
    if (is.null(data_titles)) { # cache all
        eh_ids <- unique(df_master$EHID)
    } else if (all(data_titles %in% df_master$Title)) {
        eh_ids <- df_master$EHID[df_master$Title %in% data_titles]
    } else {
        data_titles <- data_titles[!(data_titles %in% df_master$Title)]
        stop(sprintf("%s is missing from sesameData.", data_titles[1]))
    }
    eh_ids <- eh_ids[eh_ids != "TBD"]
    eh_ids <- eh_ids[!is.na(eh_ids)]
    stopifnot(length(eh_ids) > 0)
    suppressMessages(try({
        eh_ids <- eh_ids[!(eh_ids %in% names(ExperimentHub(localHub=TRUE)))]
    }, silent = TRUE))
    if (length(eh_ids) == 0) return(invisible(TRUE));
    tryCatch({
        sesameDataCache0(eh_ids)
    }, error = function(cond) {
        message("ExperimentHub Caching fails:")
        message(cond)
        return(invisible(FALSE))
    })
    invisible(TRUE)
}

#' Cache all SeSAMe data
#' 
#' @param data_titles data to cache, if not given will cache all
#' @return TRUE
#' @import ExperimentHub
#' @import AnnotationHub
#' @examples
#' sesameDataCache("genomeInfo.hg38")
#' @export
sesameDataCacheAll <- sesameDataCache
