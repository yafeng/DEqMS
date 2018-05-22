\name{peptideProfilePlot}
\alias{peptideProfilePlot}
\title{plot log2 intensities of all peptides for one gene in different samples}
\usage{
    peptideProfilePlot(dat, col=2, gene)
}
\arguments{
\item{dat}{a data frame with peptide/psm log2 intensities}

\item{col}{an integer indicates the column number where the gene protein name/id is. default is 2, asumming the gene/protein is is in the second column}

\item{gene}{an character indicates the gene name/id to be plotted}

}

\description{
    This function is to plot log2 intensities of all peptides for one gene in different samples.
}
\details{

}
\examples{
    dat.psm = readRDS(system.file("extdata","PXD004163.rds",package="DEqMS"))
    dat.psm[dat.psm == 0] <- NA # convert 0 to NA
    dat.psm = na.omit(dat.psm) # remove rows with NAs

    dat.psm.log = dat.psm
    dat.psm.log[,3:12] =  log2(dat.psm[,3:12])  # log2 transformation
    
    peptideProfilePlot(dat.psm.log,col=2,gene="TGFBR2")

}
\author{
    Yafeng Zhu
}
