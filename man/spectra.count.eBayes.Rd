\name{spectra.count.eBayes}
\alias{spectra.count.eBayes}
\title{Peptide/Spectra Count Based Empirical Bayes Statistics for Differential Expression}
\usage{
spectra.count.eBayes(fit, fit.method="loess", ...)
}
\arguments{
\item{fit}{an list object produced by Limma eBayes function, it should at least have 
attributes named sigma, count, df.residual, coefficients, and stdev.unscaled.}

\item{fit.method}{the method used to fit prior variance against the number of peptides. Two available methods: loess and nls. default loess.}

\item{coef_col}{is an integer or vector indicating the column(s) of fit$coefficients for which the function is to be performed.}

}
\value{
a list object with the additions being:
  ##	sca.t	 - Spectra Count Adjusted posterior t-value
  ##	sca.p	 - Spectra Count Adjusted posterior p-value
  ##	sca.dfprior - Spectra Count Adjusted prior degrees of freedom
  ##	sca.priorvar- Spectra Count Adjusted estimated prior variance
  ##	sca.postvar - Spectra Count Adjusted posterior variance
  ##    loess.model - fitted loess model
}
\description{
This function is to calculate peptide/spectra count adjusted t-statistics, p-values.
}
\details{

}
\examples{
Library(DEqMS)
#fit1 is an list object produced by Limma eBayes function.
fit1$count = pep.count.table[rownames(fit1$coefficients),"count"]
#Before running spectra.count.eBayes function, add gene peptide/spectra count for fit1
fit2=spectra.count.eBayes(fit1)	

}
\author{
Yafeng Zhu
}
