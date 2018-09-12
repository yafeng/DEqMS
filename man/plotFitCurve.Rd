\name{plotFitCurve}
\alias{plotFitCurve}
\title{plot the fitted prior variance against the number of quantified peptides
or PSMs}
\usage{
plotFitCurve(fit,fit.method="loess",type="boxplot",xlab="",main="")
}
\arguments{
\item{fit}{an list object produced by spectraCounteBayes function}

\item{fit.method}{the method used to fit prior variance against the number of 
peptides. default loess.}

\item{type}{an character indicating the type of plot to be generated. 
options are boxplot and scatterplot. default is boxplot}

\item{xlab}{the title for x axis}

\item{main}{the title for the figure}
}

\description{
This function is to plot the fitted prior variance against the number of 
quantified peptides/PSMs.
}

\value{
return a plot graphic
}

\examples{
library(DEqMS)
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
    
plotFitCurve(fit4,type="boxplot",xlab="PSM count",
main="TMT data PXD004163")
}
\author{
Yafeng Zhu
}
