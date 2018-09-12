\name{farmsSummary}
\alias{farmsSummary}
\title{summarize peptide/PSM intensity into protein level relative abundance by 
factor analysis}
\usage{
farmsSummary(dat,group_col=2)
}
\arguments{
\item{dat}{an data frame with raw peptide or psm intensities}

\item{group_col}{the column by which peptides/psm intensity are grouped.
Usually it is the gene/protein id column. Default is 2}
}
\value{
a data frame containing protein relative abundance estimate in log2 scale
}
\description{
This function is to calculate proteins'relative abundance by factor analysis}

\examples{
library(DEqMS)
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]
# farms method does not tolerate missing values
dat.gene = farmsSummary(dat.psm,group_col=2) 
}
\author{
Yafeng Zhu
}
