\name{medianSweeping}
\alias{medianSweeping}
\title{summarize peptide/PSM intensity into protein level relative abundance estimate by median sweeping method}
\usage{
    medianSweeping(dat,group_col=2)
}
\arguments{
\item{dat}{an data frame with peptide/PSM intensities in log2 scale}

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
    dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
    dat.pep[dat.pep==0] <- NA
    dat.pep.log = dat.pep
    dat.pep.log[,3:12] = log2(dat.pep[,3:12])
    dat.gene.nm = medianSweeping(dat.pep,group_col=2)

}
\author{
    Yafeng Zhu
}
