\name{equalMedianNormalization}
\alias{equalMedianNormalization}
\title{normalize to have equal medians in all samples}
\usage{
equalMedianNormalization(data)
}
\arguments{
\item{data}{an numeric data frame or matrix containing protein relative abundance in log2 scale}

}
\value{
    a data frame or matrix with normalized protein relative abundance
}
\description{
    This function is to normaliza out the differences of protein medians in different samples}
\details{

}
\examples{
    Library(DEqMS)
    dat.gene.nm = equalMedianNormalization(dat.gene)

}
\author{
    Yafeng Zhu
}
