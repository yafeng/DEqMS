\name{output_results}
\alias{output_results}
\title{output the DEqMS analysis results in data frame}
\usage{
output_results(fit, coef_col=1, ...)
}
\arguments{
\item{fit}{an list object produced by spectra.count.eBayes function}

\item{coef_col}{is an integer indicating the column of fit$coefficients for which corresponding t-statistics and p-values to be extracted in the output}

}
\value{
a data frame object with the last three columns being:
  ##	sca.t	 - Peptide or Spectra Count Adjusted posterior t-value
  ##	sca.P.Value	- Adjusted posterior p-value
  ## sca.adj - sca.P.Value adjusted by BH method
}
\description{
This function is to generate DEqMS outputs in a data frame.
}
\details{

}
\examples{
Library(DEqMS)
#fit1 is an list object produced by Limma eBayes function.
fit1$count = pep.count.table[rownames(fit1$coefficients),"count"]
#Before running spectra.count.eBayes function, add gene peptide/spectra count for fit1
fit2=spectra.count.eBayes(fit1)	

DEqMS.results = output_results(fit2, coef_col=1)

}
\author{
Yafeng Zhu
}
