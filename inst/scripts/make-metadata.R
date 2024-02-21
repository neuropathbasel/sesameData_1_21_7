library(readxl)
metadata <- as.data.frame(read_excel('~/Dropbox/sesame/sesameData/inst/extdata/metadata0.xlsx', sheet='metadata0'))
colnames(metadata) <- c(
    'Title', 'Description', 'BiocVersion', 'Genome', 'SourceType', 'SourceUrl',
    'SourceVersion', 'Species', 'TaxonomyId', 'Coordinate_1_based', 'DataProvider',
    'Maintainer', 'RDataClass', 'DispatchClass', 'RDataPath', 'Tags', 'Notes')

metadata$BiocVersion <- as.numeric(metadata$BiocVersion)
metadata$TaxonomyId <- as.numeric(metadata$TaxonomyId)
metadata$Coordinate_1_based <- as.logical(metadata$Coordinate_1_based)

write.csv(metadata, file='~/tools/sesame/sesameData/inst/extdata/metadata.csv', row.names = FALSE)

