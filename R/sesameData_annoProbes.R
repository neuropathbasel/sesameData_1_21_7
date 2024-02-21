#' Annotate Probes by Probe ID
#'
#' Columns in the manifests will be added to the annotation.
#' Please note that if unfound, the annotation will be NA.
#' The probe will always be kept in the output.
#'
#' @param Probe_IDs a character vector of probe IDs
#' @param regs a GenomicRanges::GRanges object against which
#' probes will be annotated, default to genes if not given
#' @param collapse whether to collapse multiple regs into one
#' @param chooseOne choose an arbitrary annotation if multiple exist
#' default to FALSE. which concatenates all with ","
#' @param sep the delimiter for collapsing
#' @param column which column in regs to annotate, if not given
#' return all overlapping probes
#' @param return_ov_probes if TRUE will return overlapping
#' probes in a GRanges object.
#' @param return_ov_features if TRUE will return overlapping
#' features in a GRanges object.
#' @param out_name column header of the annotation, use column if not given
#' @param platform EPIC, MM285 etc. will infer from Probe_IDs if not given
#' @param genome hg38, mm10, ... will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @param silent suppress messages
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges mcols
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @return a GRanges with annotated column
#' If a probe has no overlap with regs, it will be included in the results
#' with NA. But if a probe is not included in the manifest (due to
#' mappability), it won't be included in the results.
#' 
#' @examples
#' library(GenomicRanges)
#' sesameDataCache(c(
#'     "genomeInfo.mm10", "MM285.address",
#'     "genomeInfo.hg38", "Mammal40.address"))
#' 
#' regs = sesameData_getTxnGRanges("mm10")
#' Probe_IDs = names(sesameData_getManifestGRanges("MM285"))
#' anno = sesameData_annoProbes(Probe_IDs, promoters(regs), column="gene_name")
#'
#' ## get all genes associated with a probe set
#' genes = sesameData_getTxnGRanges("hg38", merge2gene = TRUE)
#' anno = sesameData_annoProbes(
#'     c("cg14620903","cg22464003"), genes, return_ov_features=TRUE)
#' @export
sesameData_annoProbes <- function(Probe_IDs, regs = NULL,
    collapse = TRUE, chooseOne = FALSE, column = NULL, sep=",",
    return_ov_probes = FALSE, return_ov_features = FALSE,
    out_name = NULL, platform = NULL, genome = NULL, silent = FALSE) {

    stopifnot(is.character(Probe_IDs))
    if(is.null(platform)) {
        platform <- inferPlatformFromProbeIDs(Probe_IDs, silent = silent)
    }

    if (is.null(regs)) { # default to annotate genes
        if (is.null(genome)) {
            genome <- sesameData_check_genome(NULL, platform) }
        regs <- sesameData_getTxnGRanges(genome)
        if (is.null(column)) column <- "gene_name"
    }

    gr <- sesameData_getManifestGRanges(platform, genome = genome)
    in_mft <- Probe_IDs %in% names(gr)
    if (sum(!in_mft) > 0) { warning(sprintf(
        "%d probes out of manifest were excluded.", sum(!in_mft))) }
    probes <- gr[Probe_IDs[in_mft]]
    if (length(probes) == 0) { return(probes); } # empty GRanges
    hits <- findOverlaps(probes, regs, ignore.strand = TRUE)

    if (return_ov_probes) {
        return(probes[unique(queryHits(hits))])
    } else if (return_ov_features) {
        return(regs[unique(subjectHits(hits))])
    }

    if (is.null(column)) {
        label <- names(regs[subjectHits(hits)])
        if (is.null(out_name)) { out_name <- "anno" }
    } else {
        stopifnot(column %in% colnames(mcols(regs)))
        label <- mcols(regs[subjectHits(hits)])[[column]]
        if (is.null(out_name)) { out_name <- column }
    }
    if (collapse) {
        if (chooseOne) {
            pid2label <- vapply(split(label, queryHits(hits)),
                function(x) x[1], character(1))
        } else {
            pid2label <- vapply(split(label, queryHits(hits)),
                function(x) paste0(unique(x), collapse=sep), character(1))
        }
        mcols(probes)[[out_name]] <- NA
        mcols(probes[as.integer(names(pid2label))])[[out_name]] <- pid2label
    } else {
        unfound <- probes[!(seq_along(probes) %in% queryHits(hits))]
        if (length(unfound) > 0) {
            mcols(unfound)[[out_name]] <- NA
        }
        found <- probes[queryHits(hits)]
        if (length(found) > 0) {
            mcols(found)[[out_name]] <- label
        }
        probes <- sort(c(unfound, found))
    }
    probes
}
