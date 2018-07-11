\name{equalMedianNormalization}
\alias{equalMedianNormalization}
\title{normalize to have equal medians in all samples}
\usage{
equalMedianNormalization(dat)
}
\arguments{
\item{dat}{an numeric data frame or matrix containing protein relative abundance in log2 scale}

}
\value{
a data frame or matrix with normalized protein relative abundance
}
\description{
This function is to normaliza out the differences of protein medians in different samples}

\examples{
dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
dat.pep[dat.pep==0] <- NA
dat.pep.log = dat.pep
dat.pep.log[,3:12] = log2(dat.pep[,3:12])
dat.gene = medianSummary(dat.pep.log,group_col = 2,ref_col =8:12)
dat.gene.nm = equalMedianNormalization(dat.gene)
}
\author{
Yafeng Zhu
}
