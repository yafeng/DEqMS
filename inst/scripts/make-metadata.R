## validate with `AnnotationHub::readMetadataFromCsv("DEqMS")`
## (above pkg directory)

main.data <- data.frame(
    Title = c("microRNA treated A431 cell proteomics data",
              "microRNA treated A431 cell proteomics data-ProteinTable",
              "MaxQuant LFQ benchmark dataset",
              "DIA quantification UPS1 spike-in dataset"),
    Description = c("High resolution isoelectric focusing LC-MS/MS data for A431 cells treated by three different microRNA mimics",
                    "TMT10plex labelled data or A431 cells treated by three different microRNA mimics",
                    "E.coli proteome spiked into Hela proteome at fixed ratio",
                    "48 human protein standards spiked into yeast proteome at different concentrations in triplicates"),
    RDataPath = c("DEqMS/PXD004163.rds",
                  "DEqMS/miR_Proteins.rda",
                  "DEqMS/LFQ_benchmark.rda",
                  "DEqMS/UPS1-spike-in.rda"),
    BiocVersion=c("3.8","3.16","3.16","3.16"),
    SourceType="TXT",
    Genome = "hg38",
    SourceVersion = c("v1","v2","v2","v2"),
    Coordinate_1_based = TRUE,
    SourceUrl=c("https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_PSM_table.flatmzidtsv.txt",
                "https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt",
                "https://ftp.ebi.ac.uk/pride-archive/2014/09/PXD000279/proteomebenchmark.zip",
                "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087597/quant/RD139_Software_OutputRScript/Spectronaut_Mixed_Fasta.txt&forceDownload=true"),
    Species=c("Homo sapiens","Homo sapiens","Homo sapiens","Saccharomyces cerevisiae"),
    TaxonomyId=c("9606","9606","9606","4932"),
    DataProvider="ProteomeXchange",
    Maintainer="Yafeng Zhu <yafeng.zhu@outlook.com>",
    RDataClass="data.frame",
    DispatchClass=c("Rds","Rda","Rda","Rda"),
    stringsAsFactors = FALSE
)

write.csv(file="metadata_v2.csv", main.data, row.names=FALSE)
