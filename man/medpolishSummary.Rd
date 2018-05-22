\name{medpolishSummary}
\alias{medpolishSummary}
\title{summarize peptide/PSM intensity into protein level relative abundance estimate by Turkey median polish procedure}
\usage{
    medpolishSummary(data,group_col=2)
}
\arguments{
\item{data}{an data frame containing peptide/psm intensities in log2 scale}

\item{group_col}{the column by which peptides/psm intensity are grouped. Usually the gene/protein column. Default is 2}
}
\value{
    a data frame containing protein relative abundance estimate in log2 scale
}
\description{
    This function is to calculate proteins'relative abundance by Turkey median polish}
\details{

}
\examples{
    dat.gene = medpolishSummary(dat.pep,group_col=2)

}
\author{
    Yafeng Zhu
}
