#' Get probes by genomic region
#'
#' The function takes a genomic coordinate and output the a vector of probes
#' on the specified platform that falls in the given genomic region.
#'
#' @param regs GRanges
#' @param chrm chromosome, when given regs are ignored
#' @param chrm_to_exclude chromosome to exclude.
#' @param beg begin, 1 if omitted
#' @param end end, chromosome end if omitted
#' @param platform EPICv2, EPIC, HM450, ...
#' @param genome hg38, mm10, ... will infer if not given.
#' For additional mapping, download the GRanges object from
#' http://zwdzwd.github.io/InfiniumAnnotation
#' and provide the following argument
#' ..., genome = sesameAnno_buildManifestGRanges("downloaded_file"),...
#' to this function.
#' @return GRanges of selected probes
#' @importMethodsFrom IRanges subsetByOverlaps
#' @examples
#'
#' ## download needed data
#' sesameDataCache(c("Mammal40.address", "genomeInfo.hg38"))
#' 
#' ## get probes in a region
#' library(GenomicRanges)
#' probes = sesameData_getProbesByRegion(
#'     GRanges('chr5', IRanges(135313937, 135419936)), platform = 'Mammal40')
#'
#' ## get all probes on chromosome 5
#' probes = sesameData_getProbesByRegion(chrm = "chr5", platform = "Mammal40")
#'
#' ## get all probes on chromosome X
#' probes = sesameData_getProbesByRegion(chrm = 'chrX', platform = "Mammal40")
#'
#' ## get all probes on both chromosome X and Y
#' probes = sesameData_getProbesByRegion(
#'     chrm = c('chrX', 'chrY'), platform = "Mammal40")
#'
#' ## get all autosomal probes
#' probes = sesameData_getProbesByRegion(
#'     chrm_to_exclude = c("chrX", "chrY"), platform = "Mammal40")
#' 
#' @export
sesameData_getProbesByRegion <- function(
    regs, chrm = NULL, beg = 1, end = -1,
    platform = NULL, chrm_to_exclude = NULL, genome = NULL) {

    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)
    probes <- sesameData_getManifestGRanges(platform, genome=genome)

    if (!is.null(chrm) || !is.null(chrm_to_exclude)) {
        if (beg == 1 && end < 0) {
            if (is.null(chrm_to_exclude)) {
                return(probes[as.character(
                    GenomicRanges::seqnames(probes)) %in% chrm])
            } else {
                return(probes[!(as.character(
                    GenomicRanges::seqnames(probes)) %in% chrm_to_exclude)])
            }
        }
        if (end < 0) {
            seqLength <- sesameData_getGenomeInfo(genome)$seqLength
            stopifnot(chrm %in% names(seqLength))
            end <- seqLength[chrm]
        }
        message(sprintf('Retrieve probes from %s:%d-%d.\n', chrm, beg, end))
        regs <- GenomicRanges::GRanges(chrm, IRanges::IRanges(beg, end))
    }
    subsetByOverlaps(probes, regs)
}

#' Get Probes by Genes or Gene Promoters
#'
#' Get probes mapped to a gene. All transcripts for the gene are considered.
#' The function takes a gene name as appears in UCSC RefGene database. The
#' platform and reference genome build can be changed with `platform` and
#' `genome` options. The function returns a vector of probes that falls
#' into the given gene.
#'
#' @param gene_name gene name, if NULL return all genes
#' @param platform EPIC or HM450
#' @param promoter if TRUE, use TSS instead of the whole gene
#' @param upstream number of bases to expand upstream of target gene
#' @param downstream number of bases to expand downstream of target gene
#' @param genome hg38 or hg19
#' @return GRanges containing probes that fall into the given gene
#' @examples
#'
#' ## download needed data
#' sesameDataCache(c("Mammal40.address", "genomeInfo.hg38"))
#' 
#' ## get all probes overlapping with DNMT3A
#' probes <- sesameData_getProbesByGene(
#'     'DNMT3A', "Mammal40", upstream=500, downstream=500)
#'
#' ## get the promoter-associated probes
#' probes <- sesameData_getProbesByGene('DNMT3A', "Mammal40", promoter = TRUE)
#' 
#' @export
sesameData_getProbesByGene <- function(
    gene_name = NULL, platform = NULL, promoter = FALSE,
    upstream = 1500, downstream = 1500, genome = NULL) {

    platform <- sesameData_check_platform(platform)
    genome <- sesameData_check_genome(genome, platform)
    
    txns <- sesameData_getTxnGRanges(genome)
    if (!is.null(gene_name)) {
        txns <- txns[txns$gene_name == gene_name]
    }
    stopifnot(length(txns) > 0)

    if (promoter) {
        reg <- promoters(
            txns, upstream = upstream, downstream = downstream)
    } else {
        up <- ifelse(as.vector(GenomicRanges::strand(
            txns)) == '-', downstream, upstream)
        dw <- ifelse(as.vector(GenomicRanges::strand(
            txns)) == '-', upstream, downstream)
        reg <- GRanges(
            as.vector(GenomicRanges::seqnames(txns)),
            IRanges(
                GenomicRanges::start(txns) - up,
                GenomicRanges::end(txns) + dw))
    }
    
    sesameData_getProbesByRegion(reg, platform = platform, genome = genome)
}

