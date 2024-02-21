#' Find genome assembly version(s) supported for a platform
#'
#' @param genome mm10, hg38, ..., or NULL
#' @param platform HM27, HM450, EPIC, EPICv2, MSA, ...
#' @return genome as string
#' @examples
#' sesameData_check_genome(NULL, "Mammal40")
#' @export
sesameData_check_genome <- function(genome, platform) {
    if ("GRanges" %in% class(genome)) { return(genome); }
    platform <- sesameData_check_platform(platform)
    supported_genomes <- c("hg19", "hg38", "mm10")
    default_genome <- c(
        HM27 = "hg38", HM450 = "hg38",
        EPIC = "hg38", EPICv2 = "hg38", MSA = "hg38",
        Mammal40 = "hg38", MM285 = "mm10")
    if (is.null(genome)) { genome <- default_genome[platform] }
    stopifnot(!is.null(genome))
    stopifnot(genome %in% supported_genomes)
    genome
}

#' Get genome info files
#' 
#' @param genome hg38, mm10, or GRanges with a metadata(genome)[["genome"]]
#' @return a list of genome info files
#' @examples
#' sesameDataCache("genomeInfo.hg38")
#' res <- sesameData_getGenomeInfo("hg38")
#' @export
sesameData_getGenomeInfo <- function(genome) {
    if ("GenomicRanges" %in% class(genome)) {
        genome <- attr(genome, "genome")
    }
    stopifnot(is.character(genome))
    key <- paste0('genomeInfo.', genome)
    if (!sesameDataHas(key)) {
        stop(sprintf("%s genome info is not found in Bioconductor.
Please go to http://zwdzwd.github.io/InfiniumAnnotation
for additional genome info files.
", genome))
    }
    sesameDataGet(key)
}

read_GENCODE_gtf <- function(x) {
    ## https://www.gencodegenes.org/pages/data_format.html
    download.file(x, sprintf("%s/gtf.gz", tempdir()), mode="wb")
    gtf <- read_tsv(sprintf("%s/gtf.gz", tempdir()), comment="#",
        col_names = c("chrm", "source", "feature_type", "start", "end",
            "score_not_used", "strand", "cds_phase", "additional"),
        col_types=cols(
            start=col_integer(),
            end=col_integer(),
            .default=col_character()))
    gtf
}

read_GENCODE_gtf_transcript <- function(gtf) {
    g <- gtf[gtf$feature_type == "transcript", ]
    g$transcript_id <- str_match(g$additional, 'transcript_id "([^"]*)"')[,2]
    g$transcript_name <- str_match(
        g$additional, 'transcript_name "([^"]*)"')[,2]
    g$transcript_type <- str_match(
        g$additional, 'transcript_type "([^"]*)"')[,2]
    g$gene_id <- str_match(g$additional, 'gene_id "([^"]*)"')[,2]
    g$gene_name <- str_match(g$additional, 'gene_name "([^"]*)"')[,2]
    g$gene_type <- str_match(g$additional, 'gene_type "([^"]*)"')[,2]
    ## there is also transcript_name, is it useful? not included
    g$level <- str_match(g$additional, 'level ([^;]*)')[,2]
    ## gene_status and transcript_status are obsolete after 25 and M11
    g
}

read_GENCODE_gtf_exon <- function(gtf) {
    g <- gtf[gtf$feature_type == "exon", ]
    g$transcript_id <- str_match(g$additional, 'transcript_id "([^"]*)"')[,2]
    g$exon_id <- str_match(g$additional, 'exon_id "([^"]*)"')[,2]
    g$exon_number <- str_match(g$additional, 'exon_number ([^;]*)')[,2]
    g
}

guess_chrmorder <- function(chrms) {
    chrms1 <- chrms[!(chrms %in% c("chrX","chrY","chrM"))]
    paste0("chr",c(as.character(seq_len(max(as.integer(str_replace(
        sort(unique(chrms1)), "chr", "")), na.rm=TRUE))), c("X","Y","M")))
}

#' build GENCODE gtf
#'
#' @param x GENCODE ftp url
#' @return GRangesList
#' @importFrom readr read_tsv
#' @importFrom readr cols
#' @importFrom readr col_integer
#' @importFrom readr col_character
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom IRanges IRanges
#' @import stringr
#' @import GenomicRanges
build_GENCODE_gtf <- function(x) {

    gtf <- read_GENCODE_gtf(x)

    ## transcript
    g1 <- read_GENCODE_gtf_transcript(gtf)
    stopifnot(length(g1$transcript_id) == length(unique(g1$transcript_id)))

    ## exon
    g2 <- read_GENCODE_gtf_exon(gtf)
    chrms <- guess_chrmorder(g2$chrm)
    gr <- GRanges(seqnames = g2$chrm, ranges=IRanges(g2$start, g2$end),
        strand = g2$strand, seqinfo=Seqinfo(chrms))
    mcols(gr)$exon_number <- as.integer(g2$exon_number)
    names(gr) <- g2$exon_id
    grl <- GRangesList(split(gr, g2$transcript_id)) # slow
    stopifnot(length(grl) == length(g1$transcript_id))
    stopifnot(all(sort(names(grl)) == sort(g1$transcript_id)))

    ## CDS
    g3 <- gtf[gtf$feature_type == "CDS", ]
    g3$transcript_id <- str_match(g3$additional, 'transcript_id "([^"]*)"')[,2]
    tid2start <- vapply(split(g3$start, g3$transcript_id), min, numeric(1))
    tid2end <- vapply(split(g3$end, g3$transcript_id), max, numeric(1))
    g1$cdsStart <- tid2start[g1$transcript_id]
    g1$cdsEnd <- tid2end[g1$transcript_id]
    
    ## put together
    g1 <- g1[order(factor(g1$chrm, levels=chrms), g1$start),]
    grl <- grl[g1$transcript_id]
    mcl <- g1[match(names(grl), g1$transcript_id), c(
        "chrm", "start", "end", "strand",
        "transcript_id", "transcript_type", "transcript_name",
        "gene_name", "gene_id", "gene_type", "source", "level",
        "cdsStart", "cdsEnd")]
    colnames(mcl)[2] <- "transcript_start"
    colnames(mcl)[3] <- "transcript_end"
    colnames(mcl)[4] <- "transcript_strand"
    mcols(grl) <- mcl
    grl
}

#' convert transcript GRanges to gene GRanges
#'
#' @param txns GRanges object
#' @return a GRanges object
#' @examples
#' txns <- sesameData_getTxnGRanges("mm10")
#' genes <- sesameData_txnToGeneGRanges(txns)
#' 
#' @export
sesameData_txnToGeneGRanges <- function(txns) {
    gene_ids <- unique(txns$gene_id)
    gene2starts <- split(start(txns), txns$gene_id)[gene_ids]
    gene2ends <- split(end(txns), txns$gene_id)[gene_ids]
    genes <- GRanges(seqnames = seqnames(txns)[match(gene_ids, txns$gene_id)],
        IRanges(
            vapply(gene2starts, min, integer(1)),
            vapply(gene2ends, max, integer(1))),
        strand = strand(txns)[match(gene_ids, txns$gene_id)])

    names(genes) <- gene_ids
    mcols(genes)$gene_name <- txns$gene_name[match(names(genes), txns$gene_id)]
    mcols(genes)$gene_type <- txns$gene_type[match(names(genes), txns$gene_id)]
    sort(genes, ignore.strand=TRUE)
}

#' convert GRangesList to transcript GRanges
#'
#' @param genome hg38, mm10, ...
#' @param grl GRangesList object
#' @param merge2gene merge transcript to genes
#' @return a GRanges object
#' @examples
#'
#' ## all mm10 transcripts
#' txns <- sesameData_getTxnGRanges("mm10")
#' 
#' ## verified protein-coding transcripts
#' txns[(txns$transcript_type == "protein_coding" & txns$level <= 2)]
#'
#' ## merged to genes
#' sesameData_getTxnGRanges("mm10", merge2gene = TRUE)
#' 
#' @export
sesameData_getTxnGRanges <- function(
    genome = NULL, grl = NULL, merge2gene = FALSE) {
    
    if (is.null(grl)) {
        genome <- sesameData_check_genome(genome, NULL)
        grl <- sesameData_getGenomeInfo(genome)$txns
    }
    mcl <- mcols(grl)
    gr <- GRanges(
        seqnames = mcl$chrm, ranges = IRanges(
            mcl$transcript_start, mcl$transcript_end),
        strand = mcl$transcript_strand, seqinfo = seqinfo(grl))
    names(gr) <- mcl$transcript_id
    mcols(gr) <- mcl[,colnames(mcl)[!(colnames(mcl) %in% c(
        "chrm","transcript_id", "transcript_start",
        "transcript_end","transcript_strand"))]]

    gr <- sort(gr, ignore.strand = TRUE)
    if (merge2gene) {
        gr <- sesameData_txnToGeneGRanges(gr)
    }
    gr
}

