\name{spectraCounteBayes}
\alias{spectraCounteBayes}
\title{Peptide/Spectra Count Based Empirical Bayes Statistics for Differential 
Expression}
\usage{
spectraCounteBayes(fit, fit.method="loess", coef_col)
}
\arguments{
\item{fit}{an list object produced by Limma eBayes function, it should have one
additional attribute named count, which stored the peptide or PSM count 
quantified for the gene in label-free or isobaric labelled data}

\item{fit.method}{the method used to fit prior variance against the number of 
peptides. Two available methods: "loess" and "nls". default "loess".}

\item{coef_col}{an integer vector indicating the column(s) of fit$coefficients 
for which the function is to be performed. if not specified, all contrasts are 
calculated}

}
\value{
a list object with the additional attributes being:
sca.t - Spectra Count Adjusted posterior t-value
sca.p - Spectra Count Adjusted posterior p-value
sca.dfprior - Spectra Count Adjusted prior degrees of freedom
sca.priorvar- Spectra Count Adjusted estimated prior variance
sca.postvar - Spectra Count Adjusted posterior variance
loess.model - fitted loess model
}
\description{
This function is to calculate peptide/PSM count adjusted t-statistics, p-values.
}
\details{
This function adjusts the T-statistics and p-values for quantitative MS 
proteomics experiment according to the number of peptides/PSMs used for 
quantification. The method is similar in nature to intensity-based Bayes 
method (Maureen A. Sartor et al BMC Bioinformatics 2006).
}
\examples{
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])

dat.gene.nm = medianSweeping(dat.psm.log,group_col = 2)
    
psm.count.table = as.data.frame(table(dat.psm$gene)) # generate PSM count table
rownames(psm.count.table)=psm.count.table$Var1
    
cond = c("ctrl","miR191","miR372","miR519","ctrl",
"miR372","miR519","ctrl","miR191","miR372")

sampleTable <- data.frame(
row.names = colnames(dat.psm)[3:12],
cond = as.factor(cond)
)
    
gene.matrix = as.matrix(dat.gene.nm)
design = model.matrix(~cond,sampleTable)

fit1 <- eBayes(lmFit(gene.matrix,design))
# add PSM count for each gene
fit1$count <- psm.count.table[rownames(fit1$coefficients),2]  

fit2 = spectraCounteBayes(fit1)
}
\author{
Yafeng Zhu
}
