\name{medpolishSummary}
\alias{medpolishSummary}
\title{summarize peptide/PSM intensity into protein level relative abundance 
estimate by Turkey median polish procedure}
\usage{
medpolishSummary(dat,group_col=2)
}
\arguments{
\item{dat}{an data frame containing peptide/psm intensities in log2 scale}

\item{group_col}{the column by which peptides/psm intensity are grouped. 
Usually the gene/protein column. Default is 2}
}
\value{
a data frame containing protein relative abundance estimate in log2 scale
}
\description{
This function is to calculate proteins'relative abundance by 
Turkey median polish}

\examples{
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])

dat.gene = medpolishSummary(dat.psm.log,group_col=2)
}
\author{
Yafeng Zhu
}
