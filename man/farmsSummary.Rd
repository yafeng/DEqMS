\name{farmsSummary}
\alias{farmsSummary}
\title{summarize peptide/PSM intensity into protein level relative abundance by factor analysis}
\usage{
    farmsSummary(dat,group_col=2)
}
\arguments{
\item{dat}{an data frame with original peptide/psm intensities (not log transformed)}

\item{group_col}{the column by which peptides/psm intensity are grouped. Usually it is the gene/protein column. Default is 2}
}
\value{
    a data frame containing protein relative abundance estimate in log2 scale
}
\description{
    This function is to calculate proteins'relative abundance by factor analysis}
\details{

}
\examples{
    dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
    dat.pep[dat.pep==0] <- NA
    dat.gene = farmsSummary(na.omit(dat.pep),group_col=2) # farms method does not tolerate missing values

}
\author{
    Yafeng Zhu
}
