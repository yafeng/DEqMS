\name{median.summary}
\alias{median summary}
\title{summarize peptide/PSM intensity into protein level relative abundance estimate by median method}
\usage{
median.summary(data,group_col=2,ref_col=3:5)
}
\arguments{
\item{data}{an data frame containing peptide/psm intensities in log2 scale}

\item{group_col}{the column by which peptides/psm intensity are grouped. Usually the gene/protein column. Default is 2}

\item{ref_col}{the column(s) used as denominator to calcualte relative petide ratio.}
}
\value{
a data frame containing protein relative abundance estimate in log2 scale
}
\description{
This function is to calculate proteins'relative abundance by median method}
\details{

}
\examples{
Library(DEqMS)
dat.gene = median.summary(dat.pep,group_col=2,ref_col=3:5)
dat.gene.nm = equal.median.normalization(dat.gene)

}
\author{
Yafeng Zhu
}
