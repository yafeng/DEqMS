# DEqMS
DEqMS is a tool for quantitative proteomic analysis, developed by Yafeng Zhu @ Karolinska Institutet. Manuscript in preparation.

## Installation
```{r}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DEqMS", version = "3.8")# current bioconductor release

# To install from latest development branch
BiocManager::install("DEqMS", version = "3.9")

```
## Introduction
DEqMS works on top of Limma. However, Limma assumes same prior variance for all genes, the function `spectraCounteBayes` in DEqMS package  is able to correct the biase of prior variance estimate for genes identified with different number of PSMs/peptides. It works in a similar way to the intensity-based hierarchical Bayes method (Maureen A. Sartor et al BMC Bioinformatics 2006).
Outputs of `spectraCounteBayes`:

object is augmented form of "fit" object from `eBayes` in Limma, with the additions being:

`sca.t`     - Spectra Count Adjusted posterior t-value

`sca.p`     - Spectra Count Adjusted posterior p-value

`sca.dfprior` - estimated prior degrees of freedom

`sca.priorvar`- estimated prior variance

`sca.postvar` - estimated posterior variance

`model` - fitted model

`fit.method` - the method used to fit the model 

## analyze TMT labelled dataset
### 1. load R packages
```{r}
library(DEqMS)
```
### Read PSM input data
Since the input data used in DEqMS is PSM or peptide level data, it is highly recommended to filter them based protein level 1% FDR. (Grouping PSMs or peptides usually generate larger list of protiens)
The first two columns in input table should be peptide sequence and protein/gene names, intensity values for different samples start from 3rd columns. It is important the input file is arranged in this way.

Here we analyzed a published protemoics dataset (TMT10plex labelling) in which A431 cells (human epidermoid carcinoma cell line) were treated with three different miRNA mimics (Zhou Y et al. Oncogene 2016). [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/27477696)

```{r}
download.file("http://lehtiolab.se/Supplementary_Files/PXD004163.zip",
destfile = "./PXD004163.zip",method = "auto")
unzip("./PXD004163.zip")
dat.psm = readRDS("PXD004163.rds")
dat.psm[dat.psm == 0] <- NA # convert 0 to NA
dat.psm = na.omit(dat.psm) # remove rows with NAs

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])  # log2 transformation

psm.count.table = as.data.frame(table(dat.psm$gene)) # generate PSM count table
rownames(psm.count.table)=psm.count.table$Var1
```
### 3. Generate sample annotation table.
```{r}
cond = c("ctrl","miR191","miR372","miR519","ctrl",
"miR372","miR519","ctrl","miR191","miR372")

sampleTable <- data.frame(
row.names = colnames(dat.psm)[3:12],
cond = as.factor(cond)
)
```

### 4. Summarization and Normalization
Choose one of the following functions to summarize peptide data to protein level. (Recommend median sweeping method)

 `group_col` is the column number you want to group by, set 2 if genes/proteins are in second column.
`ref_col`  is the columns where reference samples are.

1. use median sweeping method. [D'Angelo G et al JPR 2017](https://www.ncbi.nlm.nih.gov/pubmed/28745510) , [Herbrich SM et al JPR 2013](https://www.ncbi.nlm.nih.gov/pubmed/23270375)
```{r}
# medianSweeping does equal median normalization in the end
data.gene.nm = medianSweeping(dat.psm.log,group_col = 2)
```

2. calculate relative ratio using control/reference channels as denominator and then summarize to protein level by the median of all PSMs/Peptides.
```{r}
dat.gene = medianSummary(dat.psm.log,group_col = 2, ref_col=c(3,7,10))
dat.gene.nm = equalMedianNormalization(dat.gene)
```

3. summarize using Tukey's median polish procedure
```{r}
dat.gene = medpolishSummary(dat.psm.log,group_col = 2)
dat.gene.nm = equalMedianNormalization(dat.gene)
```

4. use Factor Analysis for Robust Microarray Summarization (FARMS)
see [Hochreiter S et al Bioinformatic 2007](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/22/8/943), [Zhang B et al MCP 2017](https://www.ncbi.nlm.nih.gov/pubmed/28302922)
```{r}
# input is psm raw intensity, not log transformed values

dat.gene = farmsSummary(dat.psm,group_col = 2)
dat.gene.nm = equalMedianNormalization(dat.gene)
```

### 5. Differential gene expression analysis
Use the dataframe `dat.gene.nm` from Step 4 Summarization and Normalization as the input.

Since the genes in this data frame are aggregated from the PSM table, which has no FDR contorl at protein/gene level.
it is recommended that the data frame `dat.gene.nm` is filtered according to a gene ID list with 1% FDR.

```{r}
##skip this if you just want to follow this tutorial
genelist = read.table("example_genetable_0.01fdr.txt",header=T,sep="\t")
data.gene.nm = dat.gene.nm[rownames(dat.gene.nm) %in% genelist$gene,]
```
contitnue to DEqMS analysis
```{r}
gene.matrix = as.matrix(dat.gene.nm)
design = model.matrix(~cond,sampleTable)

fit1 <- eBayes(lmFit(gene.matrix,design))
fit1$count <- psm.count.table[rownames(fit1$coefficients),2]  # add an attribute containing PSM/peptide count for each gene

##check the values in the vector fit1$count
##if min(fit1$count) return NA or 0, you should troubleshoot the error before you continue
min(fit1$count)
head(fit1$count)

fit2 = spectraCounteBayes(fit1)
```
### 6. plot the fitted prior variance
Check fitted relation between piror variance and peptide/PSMs count works as expected. It should look similar to the plot below. Red curve is fitted value for prior variance, y is log pooled variances calculated for each gene.
```{r}
VarianceBoxplot(fit2,n=30, main="TMT10 dataset PXD004163", xlab="PSM count) # this function only available in latest devel branch
```

![My image](https://github.com/yafeng/DEqMS/blob/master/image/PXD004163.png)

### 7. Output the results
```{r}
DEqMS.results = outputResult(fit2,coef_col=3)
write.table(DEqMS.results, "DEqMS.analysis.out.txt", quote=F,sep="\t",row.names = F)
head(DEqMS.results,n=5)
```
| logFC        | AveExpr      | t            | P.Value  | adj.P.Val   | B           | gene    | PSMcount | sca.t        | sca.P.Value | sca.adj.pval |
|--------------|--------------|--------------|----------|-------------|-------------|---------|----------|--------------|-------------|--------------|
| -1.192424423 | -0.093472789 | -18.24327059 | 3.33E-08 | 0.000156241 | 8.895452085 | ANKRD52 | 17       | -19.01031999 | 4.21E-10    | 3.42E-06     |
| -1.177468714 | -0.051976035 | -16.52824778 | 7.66E-08 | 0.000235167 | 8.286987366 | CROT    | 21       | -17.7744113  | 8.96E-10    | 3.42E-06     |
| -1.241322465 | 0.072630242  | -18.20504475 | 3.39E-08 | 0.000156241 | 8.882911735 | TGFBR2  | 8        | -17.43066408 | 1.12E-09    | 3.42E-06     |
| -0.78072293  | 0.007763848  | -13.13285085 | 5.22E-07 | 0.00096131  | 6.742716671 | PDCD4   | 40       | -14.34371369 | 9.75E-09    | 2.24E-05     |
| -0.7976368   | -0.000657979 | -14.41697245 | 2.41E-07 | 0.000553767 | 7.388343266 | PHLPP2  | 8        | -12.7927884  | 3.42E-08    | 6.08E-05     |

Column `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B` are original values generated from Limma.
Last three columns `sca.t`, `sca.P.Value` and `sca.adj.pval` are values produced from `spectraCounteBayes`, which takes into account the number of quantified spectra/peptides.

## analyze label free dataset
### 1. load R packages
```{r}
library(DEqMS)
```
### 2. Read input data and experimental design.
Here we analyze a published label-free dataset in which they did quantitative proteomic analysis to detect proteome changes in FOXP3-overexpressed gastric cancer (GC) cells. (Pan D. et al 2017 Sci Rep) [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/29089565). The data was searched by MaxQuant Software and the output file "peptides.txt" was used here. (The column "Leading razor protein" in peptides.txt table was extracted as Protein column here).
```{r}
pepTable = readRDS(system.file("int/extdata/","PXD007725.rds",package = "DEqMS"))
```

### 3. Filter peptides based on missing values (DEqMS requires minimum two observations in each condition)
```{r}
pepTable[pepTable==0] <- NA
pepTable$cond1_na_count  = apply(pepTable,1, function(x) sum(is.na(x[3:7])))
pepTable$cond2_na_count  = apply(pepTable,1, function(x) sum(is.na(x[8:12])))

#require missing values no more than 3 in each condition
df.pep.filter =pepTable[pepTable$cond1_na_count<3 & pepTable$cond2_na_count <3,1:12]
```

In our tests, imputing methods have negatively affects the statistical accuracy. Therefore, we don't impute missing values here.

### 4.  calculate ratio using control samples and then summarize to protein level by the median of all PSMs/Peptides.
```{r}
df.pep.log = df.pep.filter
df.pep.log[,3:12] = log2(df.pep.log[,3:12])

protein.df = medianSummary(df.pep.log,group_col = 2,ref_col =8:12)
protein.df.nm = equalMedianNormalization(protein.df)
```
### 5. Differential expression analysis
```{r}
protein.matrix = as.matrix(protein.df.nm)

pep.count.table = as.data.frame(table(protein.df.nm$Protein))
rownames(pep.count.table) = pep.count.table$Var1

cond = as.factor(exp_design$condition)

design = model.matrix(~0+cond) # fitting without intercept
colnames(design) = c("AF","ANC")

fit1 = lmFit(protein.matrix,design = design)
cont <- makeContrasts(AF-ANC, levels = design)
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count.table[rownames(fit3$coefficients),2]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, you should troubleshoot the error before you continue
min(fit3$count)
head(fit3$count)

fit4 = spectraCounteBayes(fit3)
```

### 6. plot the fitted prior variance
```{r}
VarianceBoxplot(fit4, n=20, main = "Label-free dataset PXD0007725",xlab="peptide count")
```
![My image](https://github.com/yafeng/DEqMS/blob/master/image/PXD007725.png)

### 7. Output the results
```{r}
AF.results = outputResult(fit4,coef_col = 1)
write.table(AF.results,"AF.DEqMS.results.txt",sep = "\t",row.names = F,quote=F)
```

## Analyze multiple TMT sets, start from gene tabe
```{r}
df.gene = read_xlsx("BC_cellline_gene.xlsx") # this gene table contains gene ratios from six TMT10plex sets
df.gene = df.gene[df.gene$`Cross-set picked q-value`<0.01,] # filter genes by 1% FDR

df.ratio = as.data.frame(df.gene[,3:56]) # extract quant values of 54 samples from column 3 to 56
rownames(df.ratio) = df.gene$gene

df.log = log2(df.ratio) # ratios are log2 transformed
df.log.narm = na.omit(df.log) # remove genes with missing values

library(matrixStats)
df.count = as.data.frame(df.gene[,57:116]) ## extract number of quantified PSMs for in total 60 TMT channels
rownames(df.count) = df.gene$gene
df.count = na.omit(df.count)

#use median number of PSM count across sets
df.count$median = round(rowMedians(as.matrix(df.count[,1:60])))
gene.matrix = as.matrix(df.log.narm)

#read in annotation table
df.annot = read.table("BC_cellline_annotation.txt",header = T,sep = "\t")
```
Specificaly make comparisons between two cell lines, eg LCC2 vs MCF7
```{r}
design = model.matrix(~0+Cell.line,df.annot) # no intercept
colnames(design) = gsub("Cell.line","",colnames(design))
fit1 = lmFit(gene.matrix,design = design)
cont <- makeContrasts(LCC2-T47D, levels = design)
fit2 = eBayes(contrasts.fit(fit1,contrasts = cont))

fit2$count = df.count[rownames(fit2$coefficients),]$median
fit3 = spectraCounteBayes(fit2)

DEqMS.result = outputResult(fit3)
```
