###############################
######### inst/extdata ########
###############################

## TCGA Ovarian tumor HM27 IDATs
## barcode: TCGA-13-0799-01A-01D-0359-05
## 4207113116_A_Grn.idat
## 4207113116_A_Red.idat

## TCGA Lymphoblastoid cell line control included in Ovarian cancer
## barcode: TCGA-07-0227-20A-01D-0402-05
## 4207113116_B_Grn.idat
## 4207113116_B_Red.idat

## EPIC.sset.LNCaP.Rep1.rds, EPIC.sset.LNCaP.Rep1.chr4.rds
## EPIC.ssets.5normals.rds, EPIC.ssets.5normals.chr4.rds, EPIC.seg.LNCaP.Rep1.rds
## source: GSE86829
library(sesame)
ssets <- readIDATsFromDir('~/tools/sesame/extdata/GSE86829_EPIC_evaluation')
names(ssets) <- gsub('\ ', '.', samplenames[match(substr(names(ssets), 1, 10), samplenames$V1), 'V2'])
EPIC.sset.LNCaP.Rep1 <- ssets[['LNCaP.Rep1']]
EPIC.sset.LNCaP.Rep1.chr4 <- EPIC.sset.LNCaP.Rep1[names(probes[seqnames(probes)=='chr4'])]
EPIC.ssets.5normals <- lapply(ssets[grepl('PrEC', names(ssets)) | grepl('NAF', names(ssets))], function(s) dyeBiasCorr(noob(s)))
EPIC.ssets.5normals.chr4 <- lapply(EPIC.ssets.5normals, function(s) s[names(probes[seqnames(probes)=='chr4'])])
EPIC.seg.LNCaP.Rep1 <- cnSegmentation(EPIC.sset.LNCaP.Rep1.chr4, ssets.normal.chr4)

## HM450.betas.10TCGAnormalPAAD.rds
## HM450.betas.76matchedTCGAchr20.rds
## HM450.sampleinfo.76matchedTCGAchr20.rds
## PAAD    TCGA-FZ-5919-11A-02D-1744-05    6055424077_R04C02
## PAAD    TCGA-FZ-5920-11A-01D-1744-05    6055424077_R05C01
## PAAD    TCGA-FZ-5922-11A-01D-1744-05    6055424077_R02C01
## PAAD    TCGA-FZ-5923-11A-01D-1744-05    6055424077_R03C02
## PAAD    TCGA-FZ-5924-11A-01D-1744-05    6055424077_R04C01
## PAAD    TCGA-FZ-5926-11A-01D-1744-05    6055424077_R06C01
## PAAD    TCGA-H6-8124-11A-01D-2399-05    8795194038_R04C02
## PAAD    TCGA-H6-A45N-11A-12D-A26Q-05    8795194086_R01C02
## PAAD    TCGA-HV-A5A3-11A-11D-A26Q-05    8795194086_R05C02
## PAAD    TCGA-YB-A89D-11A-11D-A368-05    9482801040_R03C01
merged.mapping <- read.table('450k/merged_mapping', sep='\t', stringsAsFactors = F, col.names = c('cancer','barcode','idat'))
merged.mapping.normal <- merged.mapping[substr(merged.mapping$barcode, 14,14) == '1',]
normalPAAD <- merged.mapping.normal[merged.mapping.normal$cancer == 'PAAD',]
betasPAAD <- apply(normalPAAD, 1, function(x) {
    getBetas(dyeBiasCorr(noob(readIDATs(paste0('450k/IDATs/', x[['idat']], '.rda'))[[1]])))
})
colnames(betasPAAD) <- normalPAAD$barcode

## HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05.rds
## A PAAD tumor
## PAAD    TCGA-2L-AAQA-01A-21D-A38H-05    9806233125_R03C02
HM450.betas.TCGA-2L-AAQA-01A-21D-A38H-05 <- getBetas(dyeBiasCorr(noob(readIDATs('450k/IDATs/9806233125_R03C02.rda'))[[1]]))

## HM450.ssets.10normals.rds, 10 TCGA BLCA adjacent normal samples
## "TCGA-BL-A13J-11A-13D-A10W-05" "TCGA-BT-A20J-11A-11D-A14Z-05"
## "TCGA-BT-A20N-11A-11D-A14Z-05" "TCGA-BT-A20P-11A-11D-A14Z-05"
## "TCGA-BT-A20R-11A-11D-A16P-05" "TCGA-BT-A20U-11A-11D-A14Z-05"
## "TCGA-BT-A20V-11A-11D-A14Z-05" "TCGA-BT-A20W-11A-11D-A14Z-05"
## "TCGA-BT-A20X-11A-12D-A16P-05" "TCGA-BT-A2LA-11A-11D-A18G-05"
ms <- readIDATs(merged.mapping.normal$idat, base.dir='450k/IDATs', raw=T, mc=T)
message('Processing')
ssets <- mclapply(names(dms), function(nm) {
  dm <- dms[[nm]]
  sset <- chipAddressToSignal(dm)
  sset
})
names(ssets) <- idat2barcode[names(dms)]
HM450.ssets.10normal <- ssets[1:10]


###############
### data/ #####
###############

##### cytoband ####
library(biovizBase)
## on R 3.4, not working on R 3.5, should have the package
## import biovizBase when this is fixed
cytoband.hg19 <- as.data.frame(getIdeogram('hg19', cytoband = TRUE))
cytoband.hg38 <- as.data.frame(getIdeogram('hg38', cytoband = TRUE))

library(rtracklayer)
sess <- browserSession('UCSC')
genome(sess) <- 'hg19'
cytoBand <- getTable(ucscTableQuery(sess, table='cytoBand'))
saveRDS(cytoBand, 'sesame/latest/cytoBand.hg19.rds')

genome(sess) <- 'hg38'
cytoBand <- getTable(ucscTableQuery(sess, table='cytoBand'))
saveRDS(cytoBand, 'sesame/latest/cytoBand.hg38.rds')

##### ordering #####
a <- readRDS(gzcon(url('http://zwdzwd.io/InfiniumAnnotation/current/EPIC/EPIC.hg19.manifest.rds')))
EPIC.ordering <- data.frame(Probe_ID=names(a), M=a$address_B, U=a$address_A, DESIGN=a$designType, COLOR_CHANNEL=a$channel, col=as.factor(ifelse(a$designType=='I', substr(a$channel,1,1), NA)))

a <- readRDS(gzcon(url('http://zwdzwd.io/InfiniumAnnotation/current/hm450/hm450.hg19.manifest.rds')))
HM450.ordering <- data.frame(Probe_ID=names(a), M=a$address_B, U=a$address_A, DESIGN=a$designType, COLOR_CHANNEL=a$channel, col=as.factor(ifelse(a$designType=='I', substr(a$channel,1,1), NA)))

a <- readRDS(gzcon(url('http://zwdzwd.io/InfiniumAnnotation/current/hm27/hm27.hg19.manifest.rds')))
HM27.ordering <- data.frame(Probe_ID=names(a), M=a$address_B, U=a$address_A, DESIGN=a$designType, COLOR_CHANNEL=a$channel, col=as.factor(ifelse(a$designType=='I', substr(a$channel,1,1), NA)))

##### controls ######
manifest <- read.csv('EPIC.official.manifest/MethylationEPIC_v-1-0_B4_noheader.csv', header=T, stringsAsFactors=F, row.names=1)
ct <- which(rownames(manifest) == '[Controls]')
manifest.controls <- manifest[(ct+1):nrow(manifest),]
EPIC.controls <- with(manifest.controls, data.frame(Address=rownames(manifest.controls), Type=Name, Color_Channel=AddressA_ID, Name=AlleleA_ProbeSeq))

##### non-pseudoautosomal chrY probes ######
load('/secondary/projects/shen/projects/2016_07_01_dyebias/allnormal.MPU.chrY.rda')
load('/secondary/projects/shen/projects/2016_07_01_dyebias/allnormal.pvaloob.chrY.rda')
samplemean <-  apply(allnormal.MPU.chrY, 2, mean)
femalesamples <- names(samplemean[samplemean<3000])
malesamples <- names(samplemean)[!(names(samplemean) %in% femalesamples)]
female.chrY.numdetect <- apply(allnormal.pvaloob.chrY[,femalesamples], 1, function(x) {sum(x<0.05)})
head(sort(female.chrY.numdetect, decreasing=T), n=100)
hm450.female.clean.chrY.probes <- names(female.chrY.numdetect)[female.chrY.numdetect <= 5]
head(hm450.female.clean.chrY.probes)
length(hm450.female.clean.chrY.probes) # 337 probes
save(hm450.female.clean.chrY.probes, file='~/tools/sesame/sesame/data/hm450.female.clean.chrY.probes.rda')
EPIC.female.clean.chrY.probes <- female.clean.chrY.probes[female.clean.chrY.probes %in% names(EPIC.manifest.hg38)]
length(EPIC.female.clean.chrY.probes) # 314 probes
save(EPIC.female.clean.chrY.probes, file='~/tools/sesame/sesame/data/EPIC.female.clean.chrY.probes.rda')

##### X-linked probes #####
load('/secondary/projects/shen/projects/2016_07_01_dyebias/allnormal.betas.chrX.rda')
head(allnormal.betas.chrX[,femalesamples])
head(allnormal.betas.chrX[,malesamples])
allnormal.betas.chrX <- allnormal.betas.chrX[grep('cg',rownames(allnormal.betas.chrX)),]
female.chrX.mediumCnt <- apply(allnormal.betas.chrX[,femalesamples], 1, function(x) sum(x>0.3 & x<0.7, na.rm = T))
hm450.female.xlinked.chrX.probes <- names(female.chrX.mediumCnt[female.chrX.mediumCnt > 300])
length(hm450.female.xlinked.chrX.probes) # 3797 probes
save(hm450.female.xlinked.chrX.probes, file='~/tools/sesame/sesame/data/hm450.female.xlinked.chrX.probes.rda')
EPIC.female.xlinked.chrX.probes <- female.xlinked.chrX.probes[female.xlinked.chrX.probes %in% names(EPIC.manifest.hg38)]
length(EPIC.female.xlinked.chrX.probes) # 3433 probes
save(EPIC.female.xlinked.chrX.probes, file='~/tools/sesame/sesame/data/EPIC.female.xlinked.chrX.probes.rda')

###### probe2chr #######
EPIC.hg19.probe2chr <- setNames(with(manifest.ordering, as.factor(ifelse(CHR=="", as.character(rschr), as.character(paste0('chr',CHR))))),rownames(manifest.ordering))
more.missing <- which(is.na(EPIC.hg19.probe2chr))
load('data/hm450.hg19.probe2chr.rda')
EPIC.hg19.probe2chr[more.missing] <- hm450.hg19.probe2chr[names(more.missing)]
any(is.na(EPIC.hg19.probe2chr)) # no NA
save(EPIC.hg19.probe2chr, file='data/EPIC.hg19.probe2chr.rda')

####### probe mapping #######
load('GR.InfiniumMethylation/20160711//EPIC/EPIC.manifest.hg38.rda')
mani <- EPIC.manifest.hg38
mani <- mani[!mani$MASK.mapping]
probedf <- as.data.frame(mani)
probedf$seqnames <- as.character(probedf$seqnames)
probedf <- probedf[probedf$seqnames != '*',]
EPIC.mapped.probes.hg38 <- GRanges(as.character(probedf$seqnames), IRanges(probedf$start, probedf$start), seqinfo=hg38.chrominfo$seqinfo)
names(EPIC.mapped.probes.hg38) <- rownames(probedf)

load('GR.InfiniumMethylation/20160711//EPIC/EPIC.manifest.rda')
mani <- EPIC.manifest
mani <- mani[!mani$MASK.mapping]
probedf <- as.data.frame(mani)
probedf$seqnames <- as.character(probedf$seqnames)
probedf <- probedf[probedf$seqnames != '*',]
EPIC.mapped.probes.hg19 <- GRanges(as.character(probedf$seqnames), IRanges(probedf$start, probedf$start), seqinfo=hg19.chrominfo$seqinfo)
names(EPIC.mapped.probes.hg19) <- rownames(probedf)

load('GR.InfiniumMethylation/20160711//hm450/hm450.manifest.hg38.rda')
mani <- hm450.manifest.hg38
mani <- mani[!mani$MASK.mapping]
probedf <- as.data.frame(mani)
probedf$seqnames <- as.character(probedf$seqnames)
probedf <- probedf[probedf$seqnames != '*',]
hm450.mapped.probes.hg38 <- GRanges(as.character(probedf$seqnames), IRanges(probedf$start, probedf$start), seqinfo=hg38.chrominfo$seqinfo)
names(hm450.mapped.probes.hg38) <- rownames(probedf)

load('GR.InfiniumMethylation/20160711//hm450/hm450.manifest.rda')
mani <- hm450.manifest
mani <- mani[!mani$MASK.mapping]
probedf <- as.data.frame(mani)
probedf$seqnames <- as.character(probedf$seqnames)
probedf <- probedf[probedf$seqnames != '*',]
hm450.mapped.probes.hg19 <- GRanges(as.character(probedf$seqnames), IRanges(probedf$start, probedf$start), seqinfo=hg19.chrominfo$seqinfo)
names(hm450.mapped.probes.hg19) <- rownames(probedf)

####### mask #######
EPIC.manifest <- readRDS(gzcon(url('http://zwdzwd.io/InfiniumAnnotation/current/EPIC/EPIC.hg19.manifest.rds')))
EPIC.mask <- names(EPIC.manifest[EPIC.manifest$MASK.general])

hm450.manifest <- readRDS(gzcon(url('http://zwdzwd.io/InfiniumAnnotation/current/hm450/hm450.hg19.manifest.rds')))
hm450.mask <- names(hm450.manifest[hm450.manifest$MASK.general])

hm27.manifest <- readRDS(gzcon(url('http://zwdzwd.io/InfiniumAnnotation/current/hm27/hm27.hg19.manifest.rds')))
hm27.mask <- names(hm27.manifest[hm27.manifest$MASK.general])

######## C/T extension probes ########
load('GR.InfiniumMethylation/latest//EPIC/EPIC.manifest.hg38.rda')
EPIC.typeI.extC <- names(EPIC.manifest.hg38[(!EPIC.manifest.hg38$MASK.general) & EPIC.manifest.hg38$designType=='I' & EPIC.manifest.hg38$nextBaseRef=='C'])
EPIC.typeI.extT <- names(EPIC.manifest.hg38[(!EPIC.manifest.hg38$MASK.general) & EPIC.manifest.hg38$designType=='I' & EPIC.manifest.hg38$nextBaseRef=='T'])
load('GR.InfiniumMethylation/latest//hm450/hm450.manifest.hg38.rda')
hm450.typeI.extC <- names(hm450.manifest.hg38[(!hm450.manifest.hg38$MASK.general) & hm450.manifest.hg38$designType=='I' & hm450.manifest.hg38$nextBaseRef=='C'])
hm450.typeI.extT <- names(hm450.manifest.hg38[(!hm450.manifest.hg38$MASK.general) & hm450.manifest.hg38$designType=='I' & hm450.manifest.hg38$nextBaseRef=='T'])

######## ethnicity, probes and models #######
load('GR.InfiniumMethylation/20160711//EPIC/EPIC.manifest.rda')
load('GR.InfiniumMethylation/20160711//hm450//hm450.manifest.rda')
rsprobes <- intersect(names(EPIC.manifest)[grep('rs', names(EPIC.manifest))], names(hm450.manifest)[grep('rs', names(hm450.manifest))])
rsprobes <- sort(rsprobes)
ccsprobes <- colnames(SNPswitched.maf001)
ccsprobes <- sort(ccsprobes)
samples <- intersect(rownames(SNPswitched.maf001), colnames(rs.betas))
samples <- samples[!is.na(samplerace[samples])]

df.train <- cbind(t(rs.betas[rsprobes,samples]), SNPswitched.maf001[samples, ccsprobes])
df.train[is.na(df.train)] <- 0.5
library(randomForest)
fit <- randomForest(x=df.train, y=as.factor(samplerace[rownames(df.train)]), importance=TRUE, ntree=200)
load('450k/signals.dyebias/3999492009_R01C02.rda')
load('450k/signals/3999492009_R01C02.rda')
sset[rsprobes]
b <- sset[ccsprobes]$toBetaTypeIbySum(na.mask=FALSE)

a <- getBetas(sset[rsprobes], quality.mask = F, nondetection.mask=F)
b <- getBetasTypeIbySumAlleles(sset[ccsprobes], quality.mask = F, nondetection.mask = F)
ab <- c(a,b)
predict(fit, ab)

ethnicity.ccs.probes <- ccsprobes
ethnicity.rs.probes <- rsprobes
save(ethnicity.ccs.probes, file='~/tools/sesame/sesame/data/ethnicity.ccs.probes.rda')
save(ethnicity.rs.probes, file='~/tools/sesame/sesame/data/ethnicity.rs.probes.rda')
ethnicity.model <- fit
save(ethnicity.model, file='~/tools/sesame/sesame/data/ethnicity.model.rda')

######## chrominfo #########
library(rtracklayer)
sess <- browserSession('UCSC')
## set/get genome build
genome(sess) <- 'hg19'
tableQuery <- ucscTableQuery(sess)
## list all table names:
## tableNames(tableQuery)
seqnames <- paste0('chr',c(1:22,'X','Y','M'))
chromInfo <- getTable(ucscTableQuery(sess, table='chromInfo'))
chromInfo <- chromInfo[chromInfo$chrom %in% seqnames,]
head(chromInfo)
seqlengths <- chromInfo$size[match(seqnames, chromInfo$chrom)]
seqinfo <- Seqinfo(seqnames, seqlengths)
gap <- getTable(ucscTableQuery(sess, table='gap'))
gap <- gap[gap$chrom %in% seqnames,]
cytoBand <- getTable(ucscTableQuery(sess, table='cytoBand'))
cytoBand <- cytoBand[cytoBand$chrom %in% seqnames,]
gap <- sort(GRanges(as.vector(gap$chrom), IRanges(gap$chromStart+1, gap$chromEnd), seqinfo=seqinfo))
hg19.chrominfo <- list(gap=gap, seqinfo=seqinfo)
save(hg19.chrominfo, file='~/tools/sesame/sesame/data/hg19.chrominfo.rda', compress='xz')

genome(sess) <- 'hg38'
tableQuery <- ucscTableQuery(sess)
seqnames <- paste0('chr',c(1:22,'X','Y','M'))
chromInfo <- getTable(ucscTableQuery(sess, table='chromInfo'))
chromInfo <- chromInfo[chromInfo$chrom %in% seqnames,]
head(chromInfo)
seqlengths <- chromInfo$size[match(seqnames, chromInfo$chrom)]
seqinfo <- Seqinfo(seqnames, seqlengths)
gap <- getTable(ucscTableQuery(sess, table='gap'))
gap <- gap[gap$chrom %in% seqnames,]
head(gap)
cytoBand <- getTable(ucscTableQuery(sess, table='cytoBand'))
cytoBand <- cytoBand[cytoBand$chrom %in% seqnames,]
head(cytoBand)
gap <- sort(GRanges(as.vector(gap$chrom), IRanges(gap$chromStart+1, gap$chromEnd), seqinfo=seqinfo))
hg38.chrominfo <- list(gap=gap, seqinfo=seqinfo)
save(hg38.chrominfo, file='~/tools/sesame/sesame/data/hg38.chrominfo.rda', compress='xz')

######### transcripts ##########
library(GenomicRanges)
refgene <- read.table('/references/hg19/annotations/UCSC_refGene_hg19_20161101.tsv', stringsAsFactors = F, check.names = F, header=T, comment.char = '')
colnames(refgene)
txns <- do.call(GRangesList, apply(refgene, 1, function(row) {
  g <- GRanges(row[['chrom']],
        ranges=IRanges(sapply(unlist(strsplit(row[['exonStarts']],',')), as.integer),
          sapply(unlist(strsplit(row[['exonEnds']],',')), as.integer)),
        strand=row[['strand']])
  mcols(g)$cdsStart <- row[['cdsStart']]
  mcols(g)$cdsEnd <- row[['cdsEnd']]
  sort(g)
}))

txn2gene <- split(refgene$name2, refgene$name)
gene2txn <- split(refgene$name, refgene$name2)
names(txns) <- refgene$name

saveRDS(txns, file='sesame/latest/UCSC.refGene.txns.hg19.rds')
saveRDS(txn2gene, file='sesame/latest/UCSC.refGene.txn2gene.hg19.rds')
saveRDS(gene2txn, file='sesame/latest/UCSC.refGene.gene2txn.hg19.rds')

library(GenomicRanges)
refgene <- read.table('/references/hg38/annotations/UCSC_refGene_hg38_20161101.tsv', stringsAsFactors = F, check.names = F, header=T, comment.char = '')
colnames(refgene)
txns <- do.call(GRangesList, apply(refgene, 1, function(row) {
  g <- GRanges(row[['chrom']],
        ranges=IRanges(sapply(unlist(strsplit(row[['exonStarts']],',')), as.integer),
          sapply(unlist(strsplit(row[['exonEnds']],',')), as.integer)),
        strand=row[['strand']])
  mcols(g)$cdsStart <- as.integer(row[['cdsStart']])
  mcols(g)$cdsEnd <- row[['cdsEnd']])
  sort(g)
}))

txn2gene <- split(refgene$name2, refgene$name)
gene2txn <- split(refgene$name, refgene$name2)
names(txns) <- refgene$name
saveRDS(txns, file='sesame/latest/UCSC.refGene.txns.hg38.rds')
saveRDS(txn2gene, file='sesame/latest/UCSC.refGene.txn2gene.hg38.rds')
saveRDS(gene2txn, file='sesame/latest/UCSC.refGene.gene2txn.hg38.rds')

######### age #########
a <- read.csv('Rtutorial/AdditionalFile3.csv',header=T)
saveRDS(a, file='/secondary/projects/shen/projects/2016_12_06_sesame_home/Horvath353.rds')

######### cell type reference ########
load('GR.InfiniumMethylation/20160711/hm450/hm450.manifest.hg38.rda')
load('GR.InfiniumMethylation/20160711/EPIC/EPIC.manifest.hg38.rda')

discretize <- function(m) {
  md <- as.numeric(cut(m, c(-1,0.3,0.7,1)))
  md <- matrix(md, ncol=ncol(m))
  dimnames(md) <- dimnames(m)
  md
}

construct.reference <- function(m) {
  md <- discretize(m)
  nsamples <- ncol(m)
  cnt <- apply(md, 1, function(x) sum(x==3))
  g <- setNames(rep(NA, nrow(m)), rownames(m))
  g[names(na.omit(cnt[cnt==ncol(m)]))] <- 1
  cnt <- apply(md, 1, function(x) sum(x==1))
  g[names(na.omit(cnt[cnt==ncol(m)]))] <- 0
  g
}

liftEPIC <- function(g) {
  aa <- match(names(g), EPIC.probes)
  ge <- EPIC.sig.na
  ge[na.omit(aa)] <- g[!is.na(aa)]
  ge
}

## relaxed construction, tune uniform.frac to set
## the minimum uniformity in setting the state of a probe
construct.reference.relaxed <- function(m, uniform.frac=0.8) {
  md <- discretize(m)
  nsamples <- ncol(m)
  cnt <- apply(md, 1, function(x) sum(x==3))
  g <- setNames(rep(NA, nrow(m)), rownames(m))
  g[names(na.omit(cnt[cnt > uniform.frac*ncol(m)]))] <- 1
  cnt <- apply(md, 1, function(x) sum(x==1))
  g[names(na.omit(cnt[cnt > uniform.frac*ncol(m)]))] <- 0
  g
}

construction.pipeline <- function(betas, cls, clsname) {
  b <- betas[,cls]
  attr(b, 'tissue') <- clsname
  saveRDS(b, file=sprintf('curation_betas/hm450.%s.rds', clsname))
  saveRDS(liftEPICb(b), file=sprintf('curation_betas/EPIC.%s.rds', clsname))
  g <- construct.reference(b)
  message(sum(is.na(g)))
  attr(g, 'tissue') <- 'granulocytes'
  saveRDS(g, file=sprintf('curation/hm450.%s.rds', clsname))
  saveRDS(liftEPIC(g), file=sprintf('curation/EPIC.%s.rds', clsname))
}

construction.pipeline(betas, c('Granulocytes_5684819001_R01C02', 'Granulocytes_5684819001_R02C02', 'Granulocytes_5684819001_R03C02'), 'granulocytes')
construction.pipeline(betas, c('CD4+_T-cells_5727920027_R04C02', 'CD4+_T-cells_5727920027_R05C02', 'CD4+_T-cells_5727920027_R06C02'), 'CD4T')
construction.pipeline(betas, c('CD8+_T-cells_5727920033_R01C01', 'CD8+_T-cells_5727920033_R02C01', 'CD8+_T-cells_5727920033_R03C01'), 'CD8T')
construction.pipeline(betas, c('CD14+_Monocytes_5727920033_R01C02', 'CD14+_Monocytes_5727920033_R02C02', 'CD14+_Monocytes_5727920033_R03C02'), 'CD14Monocytes')
construction.pipeline(betas, c('CD56+_NK-cells_5727920033_R04C02', 'CD56+_NK-cells_5727920033_R05C02', 'CD56+_NK-cells_5727920033_R06C02'), 'CD56NK')
construction.pipeline(betas, c('CD19+_B-cells_5727920033_R04C01', 'CD19+_B-cells_5727920033_R05C01', 'CD19+_B-cells_5727920033_R06C01'), 'CD19B')
construction.pipeline(betas, c('Eosinophils_5727920038_R01C02', 'Eosinophils_5727920038_R02C02', 'Eosinophils_5727920038_R03C02', 'Eosinophils_5727920038_R04C02', 'Eosinophils_5727920038_R05C02', 'Eosinophils_5727920038_R06C02'), 'eosinophil')
construction.pipeline(betas, c('Neutrophils_5727920038_R01C01', 'Neutrophils_5727920038_R02C01', 'Neutrophils_5727920038_R03C01', 'Neutrophils_5727920038_R04C01', 'Neutrophils_5727920038_R05C01', 'Neutrophils_5727920038_R06C01'), 'neutrophil')
construction.pipeline(betas, c('Hair_PT1', 'Hair_PT2', 'Hair_PT3', 'Hair_PT4', 'Hair_PT5'), 'hair') # 168536
construction.pipeline(betas, c('Muscle_IT9', 'Muscle_IT13', 'Muscle_IT10', 'Muscle_IT15', 'Muscle_IT16', 'Muscle_IT17'), 'muscle') # 193000
construction.pipeline(betas, c('SC_Fat_IT9', 'SC_Fat_IT13', 'SC_Fat_IT10', 'SC_Fat_IT15', 'SC_Fat_IT16', 'SC_Fat_IT17'), 'scFat') # 201708
construction.pipeline(betas, c('Omentum_IT9', 'Omentum_IT13', 'Omentum_IT10', 'Omentum_IT15', 'Omentum_IT16', 'Omentum_IT17'), 'omentum') # 207868
construction.pipeline(betas, c('Saliva_PT1', 'Saliva_PT2', 'Saliva_PT3', 'Saliva_PT4', 'Saliva_PT5'), 'saliva') # 161780
construction.pipeline(betas, c('Liver_IT9', 'Liver_IT13', 'Liver_IT10', 'Liver_IT15', 'Liver_IT17'), 'liver') # 192919
construction.pipeline(betas, c('Spleen_IT9', 'Spleen_IT10', 'Spleen_IT16'), 'spleen') # 182984
construction.pipeline(betas, c('Buccal_PT1', 'Buccal_PT2', 'Buccal_PT3', 'Buccal_PT4', 'Buccal_PT5'), 'buccal') # 185727

load('dataset/skin.rda')
colnames(betas)
b <- betas
attr(b, 'tissue') <- 'skin'
saveRDS(b, file='curation_betas/hm450.skin.rds')
saveRDS(liftEPICb(b), file='curation_betas/EPIC.skin.rds')
g <- construct.reference.relaxed(b)
sum(is.na(g)) # 204920 a bit high, need to relax
attr(g, 'tissue') <- 'skin'
saveRDS(g, file='curation/hm450.skin.rds')
saveRDS(liftEPIC(g), file='curation/EPIC.skin.rds')
