\name{medianSweeping}
\alias{medianSweeping}
\title{summarize peptide/PSM intensity into protein level relative abundance 
estimate by median sweeping method}
\usage{
medianSweeping(dat,group_col=2)
}
\arguments{
\item{dat}{an data frame with peptide/PSM intensities in log2 scale}

\item{group_col}{the column by which peptides/PSM intensity are grouped. 
Usually the gene/protein id column. Default is 2}
}
\value{
a data frame with protein relative abundance estimate in log2 scale
}
\description{
This function is to calculate proteins'relative abundance by 
median sweeping method}

\examples{
library(DEqMS)
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])

dat.gene.nm = medianSweeping(dat.psm.log,group_col = 2)
}
\author{
Yafeng Zhu
}
