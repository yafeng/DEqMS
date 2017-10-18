spectra.count.eBayes<-function(mdata,testcol) {
  
  ##########################################################################
  #  Functions used in DEqMS package to calculate spectra count adjusted p-values
  #  Written by Yafeng Zhu.
  ##########################################################################
  ##
  ##  This function adjusts the T-statistics and p-values from proteomics experiment.  
  ##  The method contains elements similar in nature to intensity-based Bayes method 
  ##  (Maureen A. Sartor et al BMC Bioinformatics 2006).  
  ##  
  ##  Inputs:
  ##  2 objects: mdata and testcol
  ##  "mdata" should be a list object from the lmFit or eBayes fcn. in  
  ##       limma, or at least have attributes named sigma, count,  
  ##     df.residual, coefficients, and stdev.unscaled.
  ##  "testcol" is an integer or vector indicating the column(s) of
  ##       mdata$coefficients for which the function is to be performed.
  ##
  ##  Outputs:
  ##  object is augmented form of "mdata" (the input), with the additions being:
  ##	sca.t	 - Spectra Count Adjusted posterior t-value
  ##	sca.p	 - Spectra Count Adjusted posterior p-value
  ##	sca.dfprior - estimated prior degrees of freedom
  ##	sca.priorvar- estimated prior variance
  ##	sca.postvar - estimated posterior variance
  ##  nls.model - fitted non-linear model
  ##
  ##  Example Function Call:
  ##      sca.fit <- spectra.count.eBayes(eBayes.output,1:4)
  ## 
  ###########################################################################
  #  The function was adapted from:
  #  Function for IBMT (Intensity-based Moderated T-statistic)
  #  Written by: Maureen Sartor,sartorma@ucmail.uc.edu
  #  University of Cincinnati, 2006
  ###########################################################################
  
  
  library("stats")
  library("limma")
  
  logVAR<-log(mdata$sigma^2)
  df<-mdata$df.residual
  numgenes<-length(logVAR[df>0])	
  df[df==0]<-NA
  eg<-logVAR-digamma(df/2)+log(df/2)
  
  y=mdata$sigma^2
  x=mdata$count
  # start values have negelectable influence on the fitted model
  nls.model = nls(y ~ const+A/x,start = list(const=0.1,A=0.05))
  y.pred = predict(nls.model)
  
  egpred<-log(y.pred)-digamma(df/2)+log(df/2)
  
  myfct<- (eg-egpred)^2 - trigamma(df/2)
  print("non linear regression (Var ~ const+A/(count))")
  
  mean.myfct<-mean(myfct,na.rm=TRUE)
  priordf<-vector(); testd0<-vector()
  for (i in 1:(numgenes*10)) {
    testd0[i]<-i/10
    priordf[i]= abs(mean.myfct-trigamma(testd0[i]/2))
    if (i>2) {
      if (priordf[i-2]<priordf[i-1]) { break }
    }
  }
  d0<-testd0[match(min(priordf),priordf)]
  print("Prior degrees freedom found")
  
  s02<-exp(egpred + digamma(d0/2) - log(d0/2))
  
  post.var<- (d0*s02 + df*mdata$sigma^2)/(d0+df)
  post.df<-d0+df
  # sca.t and scc.p stands for spectra count adjusted t and p values.  
  sca.t<-mdata$coefficients[,testcol]/(mdata$stdev.unscaled[,testcol]*sqrt(post.var)) # divided by standard error
  sca.p<-2*(1-pt(abs(sca.t),post.df))
  print("P-values calculated")
  
  output<-mdata
  output$sca.t<-sca.t
  output$sca.p<-sca.p
  output$sca.postvar<-post.var
  output$sca.priorvar<-s02
  output$sca.dfprior<-d0
  output$nls.model = nls.model
  
  output
}

plot.nls.fit <- function (fit) {
  x = fit$count
  y = fit$sigma^2
  plot(log2(x),log2(y),xlab = "log2(count)",ylab="log2(pooled variance)")
  model = fit$nls.model
  y.pred <- predict(model)
  k = order(x)
  lines(log2(x[k]),log2(y.pred[k]),col='red',lwd=3)
}



make.profile.plot <- function(dat){
  
  library(reshape2)
  library(ggplot2)
  
  m = melt(dat) 
  m$PSM_id =  rep(1:nrow(dat),10)
  ggplot(m, aes(x=variable,y=value))+
    geom_point()+
    geom_line(aes(group=PSM_id,col=Sequence))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(dat[1,2])+
    xlab("samples")+
    ylab("log2 intensity")+
    theme(plot.title = element_text(hjust = 0.5))
}


median.summary <- function(dat,group_col,ref_col) {
  library(plyr)
  library(matrixStats)
  
  dat.ratio = dat
  dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - rowMeans(dat.ratio[,ref_col])
  dat.summary = ddply(dat.ratio,colnames(dat)[group_col],
                      function(x) colMedians(as.matrix(x[,3:ncol(dat)])))
  colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
  
  return (dat.summary)
}


median_polish <- function (df) {
  dat = medpolish(df,trace.iter=FALSE)$col
  return (dat)
}

mepolish.summary <- function(dat,group_col,ref_col) {
  library(plyr)
  library(matrixStats)
  
  dat.ratio = dat
  dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - rowMeans(dat.ratio[,ref_col])
  dat.summary = ddply(dat.ratio,colnames(dat)[group_col],
                      function(x) median_polish(as.matrix(x[,3:ncol(dat)])))
  colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
  return (dat.summary)
}

equal.median.normalization <- function(df) {
  sizefactor = colMedians(as.matrix(df[,2:ncol(df)]))
  dat.nm = sweep(df[,2:ncol(df)],2,sizefactor)
  return (dat.nm)
}

median.sweeping <- function(dat,group_col) {
  library(plyr)
  library(matrixStats)
  
  dat.ratio = dat
  dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - rowMedians(as.matrix(dat.ratio[,3:ncol(dat)]))
  dat.summary = ddply(dat.ratio,colnames(dat)[group_col],
                      function(x) colMedians(as.matrix(x[,3:ncol(dat)])))
  colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
  
  dat.summary.nm = equal.median.normalization(dat.summary)
  return (dat.summary.nm)
}
