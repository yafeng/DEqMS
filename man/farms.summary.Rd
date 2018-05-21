\name{farms.summary}
\alias{farms.summary}
\title{summarize peptide/PSM intensity into protein level relative abundance estimate by factor analysis}
\usage{
farms.summary(data,group_col=2)
}
\arguments{
\item{data}{an data frame containing original peptide/psm intensities (not log transformed)}

\item{group_col}{the column by which peptides/psm intensity are grouped. Usually the gene/protein column. Default is 2}
}
\value{
a data frame containing protein relative abundance estimate in log2 scale
}
\description{
This function is to calculate proteins'relative abundance by factor analysis}
\details{

}
\examples{
Library(DEqMS)
dat.gene = farms.summary(dat.pep,group_col=2)

}
\author{
Yafeng Zhu
}
