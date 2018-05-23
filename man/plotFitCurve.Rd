\name{plotFitCurve}
\alias{plotFitCurve}
\title{plot the fitted prior variance against the number of quantified peptide/spectra}
\usage{
    plotFitCurve(fit,fit.method="loess",type="boxplot",xlab="peptide count",main="")
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
    dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
    dat.pep[dat.pep==0] <- NA
    
    dat.pep$cond1_na_count  = apply(dat.pep,1, function(x) sum(is.na(x[3:7])))
    dat.pep$cond2_na_count  = apply(dat.pep,1, function(x) sum(is.na(x[8:12])))

    #require minimum two observations in each condition
    dat.pep.filter =dat.pep[dat.pep$cond1_na_count<3 & dat.pep$cond2_na_count <3,1:12]
    
    dat.pep.log = dat.pep.filter
    dat.pep.log[,3:12] = log2(dat.pep.filter[,3:12])
    
    dat.gene.nm = medianSweeping(dat.pep.log,group_col=2)
    pep.count.table = as.data.frame(table(dat.gene.nm$Protein))
    
    cond = as.factor(c(rep("AF",5),rep("ANC",5)))
    design = model.matrix(~0+cond) # fitting without intercept
    colnames(design) = c("AF","ANC")
    
    protein.matrix = as.matrix(dat.gene.nm)
    
    fit1 = lmFit(protein.matrix,design = design)
    cont <- makeContrasts(AF-ANC, levels = design)
    fit2 = contrasts.fit(fit1,contrasts = cont)
    fit3 <- eBayes(fit2)

    fit3$count = pep.count.table[rownames(fit3$coefficients),2]
    fit4 = spectraCounteBayes(fit3)
    
    plotFitCurve(fit4,type="boxplot",xlab="peptide count",main="label-free data PXD007725")

}
\author{
    Yafeng Zhu
}
