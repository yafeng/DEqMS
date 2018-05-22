\name{plotFitCurve}
\alias{plotFitCurve}
\title{plot the fitted prior variance against the number of quantified peptide/spectra}
\usage{
    plotFitCurve(fit,type="boxplot",fit.method="loess",xlab="peptide count",main="")
}
\arguments{
\item{fit}{an list object produced by spectraCounteBayes function}

\item{fit.method}{the method used to fit prior variance against the number of peptides. default loess.}

\item{type}{an character indicating the type of plot to be generated. options are boxplot and scatterplot. default is boxplot}

\item{xlab}{the title for x axis}

\item{main}{the title for the figure}
}

\description{
    This function is to plot the fitted prior variance against the number of quantified peptide/spectra.
}
\details{

}
\examples{
    library(DEqMS)
    fit2=spectraCounteBayes(fit1)	
    plotFitCurve(fit2,type="boxplot",xlab="peptide count")

}
\author{
    Yafeng Zhu
}
