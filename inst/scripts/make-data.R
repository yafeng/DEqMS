# Downloading the data.

url1 <- "https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_PSM_table.flatmzidtsv.txt"
download.file(url1, destfile = "./miR_PSMtable.txt",method = "auto")

url2 <- "https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt"
download.file(url2, destfile = "./miR_Proteintable.txt",method = "auto")


df.psm = read.table("miR_PSMtable.txt",stringsAsFactors = FALSE,header = TRUE,
                    quote = "", comment.char = "",sep = "\t")

df.prot = read.table("miR_Proteintable.txt",stringsAsFactors = FALSE,
                     header = TRUE, quote = "", comment.char = "",sep = "\t")

colnames(df.psm)

# filter at 1% peptide FDR
df.psm.sub = df.psm[df.psm$peptide.q.value<0.01,c(12,18,31:40)]

# filter at 1% protein FDR
df.psm.filter = df.psm.sub[df.psm.sub$HGNC.Gene.Symbol %in% df.prot$Protein.accession,]

# remove PSMs with missing values
df.psm.filter[df.psm.filter==0] <-NA 
df.psm.filter = na.omit(df.psm.filter)

colnames(df.psm.filter )[2]  = "gene"

saveRDS(df.psm.filter, "PXD004163.rds")
save(df.prot,file = "miR_Proteins.rda")

url3 <- "https://ftp.ebi.ac.uk/pride-archive/2014/09/PXD000279/proteomebenchmark.zip"
download.file(url3, destfile = "./PXD000279.zip",method = "auto")
unzip("PXD000279.zip")
df.prot3 = read.table("proteinGroups.txt",header=T,sep="\t",stringsAsFactors = F,
                     comment.char = "",quote ="")
save(df.prot3,file = "LFQ_benchmark.rda")


url4 <- "ftp://massive.ucsd.edu/MSV000087597//quant/RD139_Software_Exports/Spectronaut_Mixed_Fasta_Export.csv"
download.file(url2, destfile = "Spectronaut_Mixed_Fasta_Export.csv",method = "auto")

df.dia = read.csv("Spectronaut_Mixed_Fasta.txt",row.names=1,sep = "\t")
save(df.dia,file = "UPS1-spike-in.rda")

