\name{spectraCounteBayes}
\alias{spectraCounteBayes}
\title{Peptide/Spectra Count Based Empirical Bayes Statistics for Differential Expression}
\usage{
spectraCounteBayes(fit, fit.method="loess", ...)
}
\arguments{
\item{fit}{an list object produced by Limma eBayes function, it should have one
additional attribute named count, which stored the peptide or PSM count quantified for the gene in label-free or isobaric labelled data}

\item{fit.method}{the method used to fit prior variance against the number of peptides. Two available methods: loess and nls. default loess.}

\item{coef_col}{an integer vector indicating the column(s) of fit$coefficients for which the function is to be performed.}

}
\value{
a list object with the additions being:
  ##	sca.t	 - Spectra Count Adjusted posterior t-value
  ##	sca.p	 - Spectra Count Adjusted posterior p-value
  ##	sca.dfprior - Spectra Count Adjusted prior degrees of freedom
  ##	sca.priorvar- Spectra Count Adjusted estimated prior variance
  ##	sca.postvar - Spectra Count Adjusted posterior variance
  ##loess.model - fitted loess model
}
\description{
This function is to calculate peptide/spectra count adjusted t-statistics, p-values.
}
\details{

}
\examples{
Library(DEqMS)
#fit1 is an list object produced by Limma eBayes function.
#Before running spectraCounteBayes function, add gene's peptide/spectra count to fit1
fit1$count = pep.count.table[rownames(fit1$coefficients),"count"]
fit2= spectraCounteBayes(fit1)	

}
\author{
Yafeng Zhu
}
