\name{medianSweeping}
\alias{medianSweeping}
\title{summarize peptide/PSM intensity into protein level relative abundance estimate by median sweeping method}
\usage{
    medianSweeping(data,group_col=2)
}
\arguments{
\item{data}{an data frame with peptide/PSM intensities in log2 scale}

\item{group_col}{the column by which peptides/PSM intensity are grouped. Usually the gene/protein column. Default is 2}
}
\value{
    a data frame with protein relative abundance estimate in log2 scale
}
\description{
    This function is to calculate proteins'relative abundance by median sweeping method}
\details{

}
\examples{
    library(DEqMS)
    dat.gene = medianSweeping(dat.pep,group_col=2)

}
\author{
    Yafeng Zhu
}
