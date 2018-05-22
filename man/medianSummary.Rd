\name{medianSummary}
\alias{medianSummary}
\title{summarize peptide/PSM intensity into protein level relative abundance estimate by taking the median}
\usage{
    medianSummary(data,group_col=2,ref_col)
}
\arguments{
\item{data}{an data frame with peptide/psm intensities in log2 scale}

\item{group_col}{the column by which peptides/psm intensity are grouped. Usually the gene/protein column. Default is 2}

\item{ref_col}{an integer vector indicating the column(s) used as denominator to calcualte relative petide ratio.}
}
\value{
    a data frame containing protein relative abundance estimate in log2 scale
}
\description{
    This function is to calculate proteins'relative abundance by median method}
\details{

}
\examples{
    dat.gene = medianSummary(dat.pep,group_col=2,ref_col=3:5)
    dat.gene.nm = equalMedianNormalization(dat.gene)

}
\author{
    Yafeng Zhu
}
