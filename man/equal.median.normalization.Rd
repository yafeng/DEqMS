\name{equal.median.normalization}
\alias{equal.median.normalization}
\title{normalize to have equal medians in all samples}
\usage{
equal.median.normalization(data)
}
\arguments{
\item{data}{an data frame containing protein relative abundance in log2 scale}

}
\value{
a data frame containing normalized protein relative abundance
}
\description{
This function is to normaliza out the differences of protein medians in different samples due to loading}
\details{

}
\examples{
Library(DEqMS)
dat.gene.nm = equal.median.normalization(dat.gene)

}
\author{
Yafeng Zhu
}
