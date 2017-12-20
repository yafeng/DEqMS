# DEqMS
DEqMS is a tool for quantitative proteomic analysis, developped by Yafeng Zhu @ Karolinska Institutet. Manuscript in preparation.

## Installation
git clone https://github.com/yafeng/DEqMS

or click green button (clone or download) choose Download ZIP, and unzip it.

## Introduction
DEqMS works on top of Limma. However, Limma assumes same prior variance for all genes, the function `spectra.count.eBayes` in DEqMS package  is able to correct the biase of prior variance estimate for genes identified with different number of PSMs/peptides. It works in a similar way to the intensity-based hierarchical Bayes method (Maureen A. Sartor et al BMC Bioinformatics 2006).
Outputs of `spectra.count.eBayes`:

object is augmented form of "fit" object from `eBayes` in Limma, with the additions being:

`sca.t`     - Spectra Count Adjusted posterior t-value

`sca.p`     - Spectra Count Adjusted posterior p-value

`sca.dfprior` - estimated prior degrees of freedom

`sca.priorvar`- estimated prior variance

`sca.postvar` - estimated posterior variance

`loess.model` - fitted model

## analyze TMT labelled dataset
### 1. load R packages
```{r}
source("DEqMS.R")
library(matrixStats)
library(plyr)
library(limma)
```

### 2. Read input data and generate count table.
The first two columns in input table should be peptide sequence and protein/gene names, intensity values for different samples start from 3rd columns. It is important the input file is arranged in this way.


Here we analyzed a published protemoics dataset (TMT10plex labelling) in which A431 cells (human epidermoid carcinoma cell line) were treated with three different miRNA mimics (Zhou Y et al. Oncogene 2016). [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/27477696)

```{r}
dat.psm = readRDS("./data/PXD004163.rds")

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])  # log2 transformation

count.table = as.data.frame(table(dat.psm$gene)) # generate count table
rownames(count.table)=count.table$Var1
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
Here we show how to use different functions to summarize peptide data to protein level.

 `group_col` is the column number you want to group by, set 2 if genes/proteins are in second column.
`ref_col`  is the columns where reference samples are.

1. use median sweeping method. [D'Angelo G et al JPR 2017](https://www.ncbi.nlm.nih.gov/pubmed/28745510) , [Herbrich SM et al JPR 2013](https://www.ncbi.nlm.nih.gov/pubmed/23270375)
```{r}
# median.sweep does equal median normalization for you automatically
data.gene.nm = median.sweep(dat.psm.log,group_col = 2)
```

2. calculate ratio using control samples and then summarize to protein level by the median of all PSMs/Peptides.
```{r}
dat.gene = median.summary(dat.psm.log,group_col = 2, ref_col=c(3,7,10))
dat.gene.nm = equal.median.normalization(dat.gene)
```

3. sumarize using Tukey's median polish procedure
```{r}
dat.gene = medpolish.summary(dat.psm.log,group_col = 2)
dat.gene.nm = equal.median.normalization(dat.gene)
```

4. use Factor Analysis for Robust Microarray Summarization (FARMS)
see [Hochreiter S et al Bioinformatic 2007](http://bioinformatics.oxfordjournals.org/cgi/content/abstract/22/8/943), [Zhang B et al MCP 2017](https://www.ncbi.nlm.nih.gov/pubmed/28302922)
```{r}
# input is psm raw intensity, not log transformed values

dat.gene = farms.summary(dat.psm,group_col = 2)
dat.gene.nm = equal.median.normalization(dat.gene)
```

### 5. Differential expression analysis
```{r}
gene.matrix = as.matrix(dat.gene.nm)
design = model.matrix(~cond,sampleTable)

fit1 <- eBayes(lmFit(gene.matrix,design))
fit1$count <- count.table[names(fit1$sigma),]$Freq  # add PSM/peptide count values
fit2 = spectra.count.eBayes(fit1,coef_col=3) # two arguements, a fit object from eBayes() output, and the column number of coefficients
```
### 6. plot the fitted prior variance
Check fitted relation between piror variance and peptide/PSMs count works as expected. It should look similar to the plot below. Red curve is fitted value for prior variance, y is log pooled variances calculated for each gene.
```{r}
plot.fit.curve(fit2,title="TMT10 dataset PXD004163", xlab="PSM count",type = "boxplot")
```

![My image](https://github.com/yafeng/DEqMS/blob/master/image/PXD004163.png)

### 7. Output the results
```{r}
sca.results = output_result(fit2,coef_col=3)
write.table(sca.results, "DEqMS.analysis.out.txt", quote=F,sep="\t",row.names = F)
```
## analyze label free dataset
### 1. load R packages
```{r}
source("DEqMS.R")
library(matrixStats)
library(plyr)
library(limma)
```
### 2. Read input data and experimental design.
Here we analyze a published label-free dataset in which they did quantitative proteomic analysis to detect proteome changes in FOXP3-overexpressed gastric cancer (GC) cells. (Pan D. et al 2017 Sci Rep) [Pubmed](https://www.ncbi.nlm.nih.gov/pubmed/29089565)
```{r}
pepTable = readRDS("./data/PXD007725.rds")
exp_design = read.table("./data/PXD007725_design.txt",header = T,sep = "\t",stringsAsFactors = F)
```

### 3. impute missing data using [DEP package](https://www.bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html)
```{r}
library(DEP)
library(dplyr)

pepTable[pepTable==0] <- NA # convert zero to NAs

colnames(pepTable)[1:2] = c("name","ID")
data_se <- make_se(pepTable, columns = 3:12, expdesign = exp_design)
plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 2)
plot_frequency(data_filt)
plot_numbers(data_filt)
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)
plot_detect(data_filt)

data_imp <- impute(data_norm, fun = "QRILC") # left-censored imputation method
plot_imputation(data_norm, data_imp)
```

### 4.  calculate ratio using control samples and then summarize to protein level by the median of all PSMs/Peptides.
```{r}
imp.matrix = assays(data_imp)[[1]]
imp.df  = as.data.frame(imp.matrix)
imp.df$Sequence = rownames(imp.df)
imp.df$Protein = df.peptide[imp.df$Sequence,]$Leading.razor.protein
imp.df = imp.df[,c(11,12,1:10)]


protein.df = median.summary(imp.df,group_col = 2,ref_col =8:12 )
protein.df.nm = equal.median.normalization(protein.df)
```
### 5. Differential expression analysis
```{r}
protein.matrix = as.matrix(protein.df.nm)

pep.count.table = as.data.frame(table(imp.df$Protein))
rownames(pep.count.table) = pep.count.table$Var1

cond = as.factor(exp_design$condition)

design = model.matrix(~0+cond) # fitting without intercept
colnames(design) = c("AF","ANC")

fit1 = lmFit(protein.matrix,design = design)
cont <- makeContrasts(AF-ANC, levels = design)
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count = pep.count.table[names(fit3$sigma),2]
fit4 = spectra.count.eBayes(fit3,coef_col = 1)
```

### 6. plot the fitted prior variance
```{r}
plot.loess.fit(fit4,title = "Label-free dataset PXD0007725",type = "boxplot")
```
### 7. Output the results
```{r}
AF.results = output_result(fit4,coef_col = 1)
write.table(AF.results,"AF.DEqMS.results.txt",sep = "\t",row.names = F,quote=F)
```
## Package vignette
more functioanlities are in HTML vignette.  [Go to HTML vignette](https://yafeng.github.io/DEqMS/index.html)




