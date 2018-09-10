# Downloading the data.

url1 <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/06/PXD004163/Yan_miR_PSM_table.flatmzidtsv.txt"
download.file(url1, destfile = "./miR_PSMtable.txt",method = "auto")

url2 <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt"
download.file(url2, destfile = "./miR_Proteintable.txt",method = "auto")


df.psm = read.table("miR_PSMtable.txt",stringsAsFactors = F,header = T,
                    quote = "", comment.char = "",sep = "\t")

df.prot = read.table("miR_Proteintable.txt",stringsAsFactors = F,header = T,
                     quote = "", comment.char = "",sep = "\t")

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