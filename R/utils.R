
#' Extend a GRanges
#'
#' source: https://support.bioconductor.org/p/78652/
#' @param gr a GenomicRanges::GRanges
#' @param upstream distance to expand upstream
#' @param downstream distance to expand downstream
#' @importFrom GenomicRanges trim
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @return a GenomicRanges::GRanges
extend <- function(gr, upstream=0, downstream=0) {
    if (any(strand(gr) == "*")) {
        warning("'*' ranges were treated as '+'") }

    on_plus <- strand(gr) == "+" | strand(gr) == "*"
    ranges(gr) <- IRanges(
        start(gr) - ifelse(on_plus, upstream, downstream),
        end(gr) + ifelse(on_plus, downstream, upstream))
    trim(gr)
}

valid_url <- function(url_in,t=2){
    con <- url(url_in)
    check <- suppressWarnings(try(open.connection(
        con,open="rt",timeout=t),silent=TRUE)[1])
    suppressWarnings(try(close.connection(con),silent=TRUE))
    ifelse(is.null(check),TRUE,FALSE)
}
