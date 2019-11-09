\name{Residualplot}
\alias{Residualplot}
\title{plot the residuals against the number of quantified peptides/PSMs.}
\usage{
Residualplot(fit, xlab="log2(count)",
ylab="Variance(fitted - observed)", main="")
}
\arguments{
\item{fit}{an object returned from \code{spectraCounteBayes} function}
\item{xlab}{the title for x axis}
\item{ylab}{the title for y axis}
\item{main}{the title for the figure}
}

\description{
This function is to plot the residuals of fit model on the vertical axis and the
peptide or PSM count on the horizontal axis.
}

\value{
return a plot graphic
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
    
Residualplot(fit2,xlab="log2(PSM count)",main="TMT data PXD004163")
}
\author{
Yafeng Zhu
}
