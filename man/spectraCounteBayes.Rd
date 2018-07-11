\name{spectraCounteBayes}
\alias{spectraCounteBayes}
\title{Peptide/Spectra Count Based Empirical Bayes Statistics for Differential 
Expression}
\usage{
spectraCounteBayes(fit, fit.method="loess", coef_col)
}
\arguments{
\item{fit}{an list object produced by Limma eBayes function, it should have one
additional attribute named count, which stored the peptide or PSM count 
quantified for the gene in label-free or isobaric labelled data}

\item{fit.method}{the method used to fit prior variance against the number of 
peptides. Two available methods: "loess" and "nls". default "loess".}

\item{coef_col}{an integer vector indicating the column(s) of fit$coefficients 
for which the function is to be performed. if not specified, all contrasts are 
calculated}

}
\value{
a list object with the additional attributes being:
sca.t - Spectra Count Adjusted posterior t-value
sca.p - Spectra Count Adjusted posterior p-value
sca.dfprior - Spectra Count Adjusted prior degrees of freedom
sca.priorvar- Spectra Count Adjusted estimated prior variance
sca.postvar - Spectra Count Adjusted posterior variance
loess.model - fitted loess model
}
\description{
This function is to calculate peptide/PSM count adjusted t-statistics, p-values.
}
\details{
This function adjusts the T-statistics and p-values for quantitative MS 
proteomics experiment according to the number of peptides/PSMs used for 
quantification. The method is similar in nature to intensity-based Bayes 
method (Maureen A. Sartor et al BMC Bioinformatics 2006).
}
\examples{
dat.pep = readRDS(system.file("extdata","PXD007725.rds",package="DEqMS"))
dat.pep[dat.pep==0] <- NA
    
dat.pep$cond1_na_count  = apply(dat.pep,1, function(x) sum(is.na(x[3:7])))
dat.pep$cond2_na_count  = apply(dat.pep,1, function(x) sum(is.na(x[8:12])))

#require minimum two observations in each condition
dat.pep.filter =dat.pep[dat.pep$cond1_na_count<3 & 
dat.pep$cond2_na_count <3,1:12]
    
dat.pep.log = dat.pep.filter
dat.pep.log[,3:12] = log2(dat.pep.filter[,3:12])
    
dat.gene.nm = medianSweeping(dat.pep.log,group_col=2)
    
pep.count.table = as.data.frame(table(dat.pep.log$Leading.razor.protein))
rownames(pep.count.table) = pep.count.table$Var1
    
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
}
\author{
Yafeng Zhu
}
