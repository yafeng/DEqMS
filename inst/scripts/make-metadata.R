## validate with `AnnotationHub::readMetadataFromCsv("DEqMS")`
## (above pkg directory)

main.data <- data.frame(
    Title = c("microRNA treated A431 cell proteomics data"),
    Description = c("High resolution isoelectric focusing LC-MS/MS data for A431 cells treated by three different microRNA mimics"),
    RDataPath = c("DEqMS/PXD004163.rds"),
    BiocVersion="3.8",
    SourceType="text",
    Genome = "hg38",
    SourceVersion = "b3046142103a2033da2137086b5afd36",
    Coordinate_1_based = TRUE,
    SourceUrl="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/06/PXD004163/Yan_miR_PSM_table.flatmzidtsv.txt",
    Species="Homo Sapiens",
    TaxonomyId="9606",
    DataProvider="Swedish BioMS infrastructure",
    Maintainer="Yafeng Zhu <yafeng.zhu@ki.se>",
    RDataClass="data.frame",
    DispatchClass="Rds",
    stringsAsFactors = FALSE
)

write.csv(file="metadata.csv", main.data, row.names=FALSE)
