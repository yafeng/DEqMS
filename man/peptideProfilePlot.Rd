\name{peptideProfilePlot}
\alias{peptideProfilePlot}
\title{plot log2 intensities of all peptides for one gene in different samples}
\usage{
peptideProfilePlot(dat, col=2, gene)
}
\arguments{
\item{dat}{a data frame with peptide/psm log2 intensities}

\item{col}{an integer indicates the column number where the gene protein 
id is. default is 2, asumming the gene/protein is in the second column}

\item{gene}{an character indicates the gene name/id to be plotted}
}

\description{
This function is to plot log2 intensities of all peptides for one gene 
in different samples.
}
\value{
return a ggplot2 object
}
\examples{
dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
dat.pep[dat.pep==0] <- NA
dat.pep.log = na.omit(dat.pep)
dat.pep.log[,3:12] = log2(dat.pep.log[,3:12])
    
peptideProfilePlot(dat.psm.log,col=2,gene="IPI00001159")

}
\author{
Yafeng Zhu
}
