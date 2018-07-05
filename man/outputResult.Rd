\name{outputResult}
\alias{outputResult}
\title{output the DEqMS analysis results in a data frame}
\usage{
    outputResult(fit, coef_col=1)
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
    dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
    dat.pep[dat.pep==0] <- NA
    
    dat.pep$cond1_na_count  = apply(dat.pep,1, function(x) sum(is.na(x[3:7])))
    dat.pep$cond2_na_count  = apply(dat.pep,1, function(x) sum(is.na(x[8:12])))

    #require minimum two observations in each condition
    dat.pep.filter =dat.pep[dat.pep$cond1_na_count<3 & dat.pep$cond2_na_count <3,1:12]
    
    dat.pep.log = dat.pep.filter
    dat.pep.log[,3:12] = log2(dat.pep.filter[,3:12])
    
    dat.gene.nm = medianSweeping(dat.pep.log,group_col=2)
    
    pep.count.table = as.data.frame(table(dat.pep.log$Leading.razor.protein))
    rownames(pep.count.table) = pep.count.table$Var1
    
    cond = as.factor(c(rep("AF",5),rep("ANC",5)))
    design = model.matrix(~0+cond) # without intercept
    colnames(design) = c("AF","ANC")
    
    protein.matrix = as.matrix(dat.gene.nm)
    
    fit1 = lmFit(protein.matrix,design = design)
    cont <- makeContrasts(AF-ANC, levels = design)
    fit2 = contrasts.fit(fit1,contrasts = cont)
    fit3 <- eBayes(fit2)

    fit3$count = pep.count.table[rownames(fit3$coefficients),2]
    fit4 = spectraCounteBayes(fit3)
    
    DEqMS.results = outputResult(fit4, coef_col=1)

}
\author{
    Yafeng Zhu
}
