\name{outputResult}
\alias{outputResult}
\title{output the DEqMS analysis results in a data frame}
\usage{
outputResult(fit, coef_col=1)
}
\arguments{
\item{fit}{an list object produced by spectraCounteBayes function}

\item{coef_col}{is an integer indicating the column of fit$coefficients for 
which corresponding t-statistics and p-values are extracted in the output}

}
\value{
a data frame object with the last three columns being:
sca.t - Peptide or Spectra Count Adjusted posterior t-value
sca.P.Value - Adjusted posterior p-value
sca.adj - sca.P.Value adjusted by BH method
}
\description{
This function is to generate DEqMS outputs in a data frame.
}

\examples{
library(ExperimentHub)
eh = ExperimentHub(localHub=TRUE)
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
    
DEqMS.results = outputResult(fit2, coef_col=3)
}
\author{
Yafeng Zhu
}
