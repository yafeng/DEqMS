\name{peptideProfilePlot}
\alias{peptideProfilePlot}
\title{plot log2 intensities of all peptides for one gene in different samples}
\usage{
peptideProfilePlot(data, col=2, gene)
}
\arguments{
\item{data}{a data frame with peptide/psm log2 intensities}

\item{col}{an integer indicates the column number where the gene protein name/id is. default is 2, asumming the gene/protein is is in the second column}

\item{gene}{an character indicates the gene name/id to be plotted}

}
\value{

}
\description{
This function is to plot log2 intensities of all peptides for one gene in different samples.
}
\details{

}
\examples{
Library(DEqMS)
peptideProfilePlot(dat.pep,col=2,gene="TGFBR2")

}
\author{
Yafeng Zhu
}
