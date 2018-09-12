\name{equalMedianNormalization}
\alias{equalMedianNormalization}
\title{normalize to have equal medians in all samples}
\usage{
equalMedianNormalization(dat)
}
\arguments{
\item{dat}{an numeric data frame or matrix containing protein relative 
abundance in log2 scale}

}
\value{
a data frame or matrix with normalized protein relative abundance
}
\description{
This function is to normaliza out the differences of protein medians in 
different samples}

\examples{
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])
# use the 3 ctrl samples as reference channels to calculate log2 ratio
dat.gene = medianSummary(dat.psm.log,group_col = 2,ref_col =c(3,7,10))
dat.gene.nm = equalMedianNormalization(dat.gene)
}
\author{
Yafeng Zhu
}
