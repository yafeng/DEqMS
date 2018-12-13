\name{spectraCounteBayes}
\alias{spectraCounteBayes}
\title{Peptide/Spectra Count Based Empirical Bayes Statistics for Differential 
Expression}
\usage{
spectraCounteBayes(fit, fit.method="loess", coef_col)
}
\arguments{
\item{fit}{an list object produced by Limma \code{eBayes} function, it should 
have one additional attribute \code{$count}, which stored the peptide or PSM count 
quantified for the gene in label-free or isobaric labelled data.}

\item{fit.method}{the method used to fit variance against the number of 
peptides/PSM count quantified. Two available methods: "loess","nls" and 
"spline". default "loess"."loess" uses \code{loess} and span = 0.75, "nls"" 
uses a explicit formula \code{y~a+b/x}. "spline" uses \code{smooth.spline} and 
"generalized cross-validation" for smoothing parameter computation. For "nls", 
independent variable x is peptide/PSM count, response y is pooled variance 
(\code{fit$sigma^2}). For "loess" and "spline" method, both x and y are log
transformed before applying the two methods. In most of time, "loess" is sufficient.
To quickly assess the fit model, use \code{VarianceScatterplot} and
\code{Residualplot} functions.}

\item{coef_col}{an integer vector indicating the column(s) of fit$coefficients 
for which the function is to be performed. if not specified, all coefficients 
are used.}

}
\value{
a list object with the additional attributes
\item{count} {Peptide or PSM count used for quantification}
\item{sca.t} {Spectra Count Adjusted posterior t-value}
\item{sca.p} {Spectra Count Adjusted posterior p-value}
\item{sca.dfprior} {Spectra Count Adjusted prior degrees of freedom}
\item{sca.priorvar} {Spectra Count Adjusted prior variance}
\item{sca.postvar} {Spectra Count Adjusted posterior variance}
\item{model} {fitted model}
\item{fit.method} {The method used to fit the model}
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
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "DEqMS")
dat.psm = eh[["EH1663"]]

dat.psm.log = dat.psm
dat.psm.log[,3:12] =  log2(dat.psm[,3:12])

dat.gene.nm = medianSweeping(dat.psm.log,group_col = 2)
    
psm.count.table = as.data.frame(table(dat.psm$gene)) # generate PSM count table
rownames(psm.count.table)=psm.count.table$Var1
    
cond = c("ctrl","miR191","miR372","miR519","ctrl",
"miR372","miR519","ctrl","miR191","miR372")

sampleTable <- data.frame(
row.names = colnames(dat.psm)[3:12],
cond = as.factor(cond)
)
    
gene.matrix = as.matrix(dat.gene.nm)
design = model.matrix(~cond,sampleTable)

fit1 <- eBayes(lmFit(gene.matrix,design))
# add PSM count for each gene
fit1$count <- psm.count.table[rownames(fit1$coefficients),2]  

fit2 = spectraCounteBayes(fit1)
}
\author{
Yafeng Zhu
}
