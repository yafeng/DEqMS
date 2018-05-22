\name{outputResult}
\alias{outputResult}
\title{output the DEqMS analysis results in a data frame}
\usage{
    outputResult(fit, coef_col=1, ...)
}
\arguments{
\item{fit}{an list object produced by spectraCounteBayes function}

\item{coef_col}{is an integer indicating the column of fit$coefficients for which corresponding t-statistics and p-values are extracted in the output}

}
\value{
    a data frame object with the last three columns being:
        sca.t	 - Peptide or Spectra Count Adjusted posterior t-value
        sca.P.Value	- Adjusted posterior p-value
        sca.adj - sca.P.Value adjusted by BH method
}
\description{
    This function is to generate DEqMS outputs in a data frame.
}
\details{

}
\examples{
    library(DEqMS)
    DEqMS.results = outputResult(fit, coef_col=1)

}
\author{
    Yafeng Zhu
}
