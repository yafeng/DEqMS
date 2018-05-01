spectra.count.eBayes<-function(mdata,coef_col,fit.method="loess") {
  
  ##########################################################################
  #  Functions used in DEqMS package to calculate spectra count adjusted p-values
  #
  ##########################################################################
  ##
  ##  This function adjusts the T-statistics and p-values from proteomics experiment.  
  ##  The method is similar in nature to intensity-based Bayes method
  ##  (Maureen A. Sartor et al BMC Bioinformatics 2006).  
  ##  
  ##  Inputs:
  ##  2 objects: mdata and ceof_col
  ##  "mdata" should be a list object from the lmFit or eBayes fcn. in  
  ##       limma, or at least have attributes named sigma, count,  
  ##     df.residual, coefficients, and stdev.unscaled.
  ##  "coef_col" is an integer or vector indicating the column(s) of
  ##       mdata$coefficients for which the function is to be performed.
  ##
  ##  Outputs:
  ##  object is augmented form of "mdata" (the input), with the additions being:
  ##	sca.t	 - Spectra Count Adjusted posterior t-value
  ##	sca.p	 - Spectra Count Adjusted posterior p-value
  ##	sca.dfprior - estimated prior degrees of freedom
  ##	sca.priorvar- estimated prior variance
  ##	sca.postvar - estimated posterior variance
  ##    nls.model - fitted non-linear model
  ##
  ##  Example Function Call:
  ##      sca.fit <- spectra.count.eBayes(eBayes.output,1:4)
  ## 
  ###########################################################################
  #  The function was adapted from:
  #  Function for IBMT (Intensity-based Moderated T-statistic)
  #  Written by: Maureen Sartor
  #  University of Cincinnati, 2006
  ###########################################################################
  
  
  require("stats")
  require("limma")
  
  logVAR<-log(mdata$sigma^2)
  df<-mdata$df.residual
  numgenes<-length(logVAR[df>0])	
  df[df==0]<-NA
  eg<-logVAR-digamma(df/2)+log(df/2)
  names(mdata$count) = rownames(mdata$coefficients)
  
  output<-mdata
  
  if (fit.method == "loess"){
    x=log2(mdata$count)
    loess.model = loess(logVAR~x,span = 0.75)
    y.pred = predict(loess.model)
    output$loess.model = loess.model
    }else if (fit.method == "nls"){
      x=mdata$count
      y=mdata$sigma^2
      nls.model =  nls(y~a+b/x,start = (list(a=0.1,b=0.05)))
      y.pred = log(predict(nls.model))
      output$nls.model = nls.model
    }
  
  egpred<-y.pred-digamma(df/2)+log(df/2)
  
  myfct<- (eg-egpred)^2 - trigamma(df/2)
  
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
  #print("Prior degrees freedom found")
  
  s02<-exp(egpred + digamma(d0/2) - log(d0/2))
  
  post.var<- (d0*s02 + df*mdata$sigma^2)/(d0+df)
  post.df<-d0+df
  # sca.t and scc.p stands for spectra count adjusted t and p values.  
  sca.t<-as.matrix(mdata$coefficients[,coef_col]/(mdata$stdev.unscaled[,coef_col]*sqrt(post.var))) # divided by standard error
  sca.p<-as.matrix(2*(1-pt(abs(sca.t),post.df)))
  #print("P-values calculated")
  
  output$sca.t<-sca.t
  output$sca.p<-sca.p
  output$sca.postvar<-post.var
  output$sca.priorvar<-s02
  output$sca.dfprior<-d0
  
  return (output)
}

output_result <-function(sca.fit,coef_col){
  results.table = topTable(sca.fit,coef = coef_col,n= Inf)
  
  results.table$gene = rownames(results.table)
  results.table$PSMcount = sca.fit$count[results.table$gene]
  
  results.table$sca.t = sca.fit$sca.t[results.table$gene,coef_col]
  results.table$sca.P.Value = sca.fit$sca.p[results.table$gene,coef_col]
  results.table$sca.adj.pval = p.adjust(results.table$sca.P.Value,method = "BH")
  results.table = results.table[order(results.table$sca.P.Value), ]
}

plot.fit.curve <- function (fit,main="", fit.method="nls",xlab="feature count",type = "boxplot") {
  x = fit$count
  y = fit$sigma^2
  
  if (fit.method=="nls"){
    model = fit$nls.model
    
    if (type=="scatterplot"){
      plot(x,log(y),ylab="log(Variance)",main= main)
      
      y.pred <- log(predict(model))
      k = order(x)
      lines((x[k]),y.pred[k],col='red',lwd=3)
      
    }else if (type=="boxplot"){
      df.temp = data.frame(pep_count =x, variance = y )
      df.temp.filter = df.temp[df.temp$pep_count<21,]
      boxplot(log(variance)~pep_count,df.temp.filter, xlab=xlab,
              ylab = "log(Variance)", main=main)
      
      y.pred <- log(predict(model,data.frame(x=1:20)))
      lines(1:20,y.pred,col='red',lwd=3)
    }else{stop("only scatterplot and boxplot are supported")}
    
  }else if (fit.method=="loess"){
    model = fit$loess.model
  if (type=="scatterplot"){
    plot(x,log(y),ylab="log(Variance)",main= main)
    
    y.pred <- predict(model)
    k = order(x)
    lines((x[k]),y.pred[k],col='red',lwd=3)
    
  }else if (type=="boxplot"){
    df.temp = data.frame(pep_count =x, variance = y )
    df.temp.filter = df.temp[df.temp$pep_count<21,]
    boxplot(log(variance)~pep_count,df.temp.filter, xlab=xlab,
            ylab = "log(Variance)", main=main)
    
    y.pred <- predict(model,log2(1:20))
    lines(1:20,y.pred,col='red',lwd=3)
  }else{stop("only scatterplot and boxplot are supported")}
  }  
}

make.profile.plot <- function(dat){
  
  require(reshape2)
  require(ggplot2)
  
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


median.summary <- function(dat,group_col=2,ref_col) {
  require(plyr)
  require(matrixStats)
  
  dat.ratio = dat
  dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - rowMeans(dat.ratio[,ref_col],na.rm = T)
  dat.summary = ddply(dat.ratio,colnames(dat)[group_col],
                      function(x) colMedians(as.matrix(x[,3:ncol(dat)]),na.rm = T))
  colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
  
  dat.new = dat.summary[,-1]
  rownames(dat.new) = dat.summary[,1]
  return (dat.new)
}


median_polish <- function (m) {
  require(matrixStats)
  dat = medpolish(m,trace.iter=FALSE)$col
  return (dat)
}

medpolish.summary <- function(dat,group_col=2) {
  require(plyr)
  
  dat.summary = ddply(dat,colnames(dat)[group_col],
                      function(x) median_polish(as.matrix(x[,3:ncol(dat)])))
  colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
  
  dat.new = dat.summary[,-1]
  rownames(dat.new) = dat.summary[,1]
  return (dat.new)
}



equal.median.normalization <- function(dat) {
  require(matrixStats)
  sizefactor = colMedians(as.matrix(dat),na.rm = T)
  dat.nm = sweep(dat,2,sizefactor)
  return (dat.nm)
}

median.sweeping <- function(dat,group_col=2) {
  require(plyr)
  require(matrixStats)
  
  dat.ratio = dat
  dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - rowMedians(as.matrix(dat.ratio[,3:ncol(dat)]),na.rm = T)
  dat.summary = ddply(dat.ratio,colnames(dat)[group_col],
                      function(x) colMedians(as.matrix(x[,3:ncol(dat)]),na.rm = T))
  colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
  
  dat.new = dat.summary[,-1]
  rownames(dat.new) = dat.summary[,1]
  
  dat.nm = equal.median.normalization(dat.new)
  return (dat.nm)
}

farms.method <- function(df){
  library(farms)
  if (nrow(df)==1){
    dat = log2(as.matrix(df))
  }else {dat = generateExprVal.method.farms(as.matrix(na.omit(df)))$expr}
  return (dat)
}

farms.summary <- function(dat,group_col) {
  require(plyr)
  dat.log = ddply(dat,colnames(dat)[group_col],
                      function(x) farms.method(x[,3:ncol(dat)]))
  
  colnames(dat.log)[2:ncol(dat.log)]=colnames(dat)[3:ncol(dat)]
  
  dat.ratio = dat.log
  dat.ratio[,2:ncol(dat.ratio)]= dat.log[,2:ncol(dat.log)] - rowMedians(as.matrix(dat.log[,2:ncol(dat.log)]),na.rm = T)
  
  dat.summary = dat.ratio[,-1]
  rownames(dat.summary) = dat.ratio[,1]
  
  return (dat.summary)
}
