#' Master data frame for all object to cache
#'
#' This is an internal object which will be updated on every new release
#' library(ExperimentHub)
#' eh <- query(ExperimentHub(localHub=FALSE), c("sesameData", "v1.13.1"))
#' data.frame(name=eh$title, eh=names(eh))
#'
#' Cache location is default to
#' /Users/zhouw3/Library/Caches/org.R-project.R/R/ExperimentHub/
#'
#' @name df_master
#' @docType data
#' @return master sheet of sesameData objects
NULL

cacheEnv <- new.env()
alt_base <- "https://zhouserver.research.chop.edu"

#' Check whether the title exists in cacheEnv
#'
#' @param title the title to check
#' @return the data associated with the title or NULL if title doesn't exist
sesameDataGet_checkEnv <- function(title) {
    if (exists(title, envir = cacheEnv, inherits = FALSE)) {
        return(get(title, envir = cacheEnv, inherits = FALSE))
    } else {
        return(NULL)
    }
}

sesameDataGet_assignEnv <- function(title, data) {
    assign(title, data, envir = cacheEnv)
}

#' Empty cache environment to free memory
#'
#' When this function is called sesameDataGet will
#' retrieve all data from disk again instead of using the in-memory
#' cache, i.e., sesameData:::cacheEnv.
#'
#' Note this is different from sesameDataClearHub which empties the
#' ExperimentHub on disk.
#'
#' @return gc() output
#' @examples
#' sesameDataGet_resetEnv()
#' @export
sesameDataGet_resetEnv <- function() {
    rm(list=ls(envir=cacheEnv), envir=cacheEnv)
    gc()
}
