\name{medianSummary}
\alias{medianSummary}
\title{summarize peptide/PSM intensity into protein level relative abundance 
estimate by taking the median}
\usage{
medianSummary(dat,group_col=2,ref_col)
}
\arguments{
\item{dat}{an data frame with peptide/psm intensities in log2 scale}

\item{group_col}{the column by which peptides/psm intensity are grouped. 
Usually the gene/protein id column. Default is 2}

\item{ref_col}{an integer vector indicating the column(s) used as denominator 
to calcualte relative petide ratio.}
}
\value{
a data frame containing protein relative abundance estimate in log2 scale
}
\description{
This function is to calculate proteins'relative abundance by median method}

\examples{
library(DEqMS)
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])
# use the 3 ctrl samples as reference channels to calculate log2 ratio
dat.gene = medianSummary(dat.psm.log,group_col = 2,ref_col =c(3,7,10))
}
\author{
Yafeng Zhu
}
