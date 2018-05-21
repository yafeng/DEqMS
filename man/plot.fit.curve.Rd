\name{plot.fit.curve}
\alias{plot.fit.curve}
\title{plot the fitted prior variance to the peptides/spectra count}
\usage{
plot.fit.curve(fit,type="boxplot",xlab="peptide count")
}
\arguments{
\item{fit}{an list object produced by spectra.count.eBayes function}

\item{fit.method}{the method used to fit prior variance against the number of peptides. default loess.}

\item{type}{an character indicating the type of plot to be generated. options are boxplot and scatterplot. default is boxplot}

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
This function is to plot the fitted prior variance to the peptides/spectra count
}
\details{

}
\examples{
Library(DEqMS)
fit2=spectra.count.eBayes(fit1)	
plot.fit.curve(fit2,type="boxplot",xlab="peptide count")

}
\author{
Yafeng Zhu
}
