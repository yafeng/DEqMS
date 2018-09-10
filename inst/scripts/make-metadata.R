## validate with `AnnotationHub::readMetadataFromCsv("miRProteomicsData")`
## (above pkg directory)

main.data <- data.frame(
    Title = c("microRNA treated A431 cell proteomics data"),
    Description = c("High resolution isoelectric focusing LC-MS/MS data for A431 cells treated by three different microRNA mimics"),
    RDataPath = c("miRProteomicsData/PXD004163.rds"),
    BiocVersion="3.7",
    SourceType="text",
    SourceUrl="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/06/PXD004163/Yan_miR_PSM_table.flatmzidtsv.txt",
    Species="Homo Sapiens",
    TaxonomyId="9606",
    DataProvider="Swedish BioMS infrastructure",
    Maintainer="Yafeng Zhu <yafeng.zhu@ki.se>",
    RDataClass="character",
    DispatchClass="Rds",
    stringsAsFactors = FALSE
)

write.csv(file="metadata.csv", main.data, row.names=FALSE)
