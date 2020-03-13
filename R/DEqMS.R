spectraCounteBayes<-function(fit,fit.method="loess",coef_col) {
    
    ################################################
    #  The function was adapted from:
    #  Function for IBMT (Intensity-based Moderated 
    #  T-statistic) Written by Maureen Sartor
    #  University of Cincinnati, 2006
    ################################################
    logVAR<-log(fit$sigma^2)
    df<-fit$df.residual
    numgenes<-length(logVAR[df>0])
    df[df==0]<-NA
    eg<-logVAR-digamma(df/2)+log(df/2)
    names(fit$count) <- rownames(fit$coefficients)
    output<-fit
    output$fit.method <- fit.method

    if (fit.method == "loess"){
        x<-log2(fit$count)
        loess.model <- loess(logVAR~x,span = 0.75)
        y.pred <- fitted(loess.model)
        output$model <- loess.model
    }else if (fit.method == "nls"){
        x<-fit$count
        y<-fit$sigma^2
        nls.model <-  nls(y~a+b/x,start = (list(a=0.1,b=0.05)))
        y.pred <- log(fitted(nls.model))
        output$model <- nls.model
    }else if (fit.method == "spline"){
        x<-log2(fit$count)
        spline.model <- smooth.spline(x,logVAR,cv=FALSE)
        y.pred <- fitted(spline.model)
        output$model <- spline.model
    }
    
    egpred<-y.pred-digamma(df/2)+log(df/2)
    
    myfct<- (eg-egpred)^2 - trigamma(df/2)
    
    mean.myfct<-mean(myfct,na.rm=TRUE)
    
    priordf<-vector()
    testd0<-vector()
    
    for (i in seq(1,numgenes*10)) {
        testd0[i]<-i/10
        priordf[i]= abs(mean.myfct-trigamma(testd0[i]/2))
        if (i>2) {
            if (priordf[i-2]<priordf[i-1]) { break }
        }
    }
    d0<-testd0[match(min(priordf),priordf)] # prior degree found
    
    s02<-exp(egpred + digamma(d0/2) - log(d0/2)) # calculate prior variance
    
    post.var<- (d0*s02 + df*fit$sigma^2)/(d0+df)
    post.df<-d0+df
    # sca.t and scc.p stands for spectra count adjusted t and p values.  
    sca.t<-as.matrix(fit$coefficients[,coef_col]/(fit$stdev.unscaled[,coef_col]
    *sqrt(post.var)))
    sca.p<-as.matrix(2*pt(abs(sca.t),post.df,lower.tail = FALSE))
    
    output$sca.t<-sca.t
    output$sca.p<-sca.p
    output$sca.postvar<-post.var
    output$sca.priorvar<-s02
    output$sca.dfprior<-d0
    
    return (output)
}

outputResult <-function(fit,coef_col=1){
    results.table = limma::topTable(fit,coef = coef_col,n= Inf)
    
    results.table$gene = rownames(results.table)
    results.table$count = fit$count[results.table$gene]
    
    results.table$sca.t = fit$sca.t[results.table$gene,coef_col]
    results.table$sca.P.Value = fit$sca.p[results.table$gene,coef_col]
    results.table$sca.adj.pval = p.adjust(results.table$sca.P.Value,
    method = "BH")
    results.table = results.table[order(results.table$sca.P.Value), ]
}


VarianceBoxplot <- function (fit, n=20, xlab="count",
                             ylab = "log(Variance)", main=""){
  x <- fit$count
  y <- fit$sigma^2
  
  df.temp <- data.frame(pep_count =x, variance = y )
  df.temp.filter <- df.temp[df.temp$pep_count<=n,]
  
  if (fit$fit.method=="nls"){
    y.pred <- log(predict(fit$model,data.frame(x=seq(1,n))))
  }else if (fit$fit.method=="loess"){
    y.pred <- predict(fit$model,data.frame(x=log2(seq(1,n))))}
  else if (fit$fit.method=="spline"){
    y.pred <- predict(fit$model,x=log2(seq(1,n)))$y
  }
  
  boxplot(log(variance)~pep_count,df.temp.filter, xlab=xlab,
          ylab = ylab, main=main)
  lines(seq(1,n),y.pred,col='red',lwd=3)
  
}

VarianceScatterplot <- function (fit, xlab="log2(count)",
                                 ylab = "log(Variance)",main=""){
  x <- fit$count
  y <- fit$sigma^2
  
  if (fit$fit.method=="nls"){
    y.pred <- log(fitted(fit$model))
  }else if (fit$fit.method=="loess" | fit$fit.method=="spline" ) {
    y.pred <- fitted(fit$model)}
  
  plot(log2(x),log(y),ylab=ylab,xlab=xlab,main= main)
  k = order(x)
  lines(log2(x[k]),y.pred[k],col='red',lwd=3)
  
}

Residualplot <- function (fit, xlab="log2(count)",
                          ylab="Variance(fitted - observed)", main=""){
  x <- fit$count
  
  if (fit$fit.method=="nls"){
    y <- log(fitted(fit$model)) - log(fit$sigma^2)
  }else if (fit$fit.method=="loess" | fit$fit.method=="spline") {
    y <- residuals(fit$model)}
  
  plot(log2(x),y, pch=20, cex=0.5,ylab=ylab,xlab=xlab,main= main)
}

peptideProfilePlot <- function(dat,col=2,gene){
    
    dat.sub = dat[dat[,col]==gene,]
    
    colnames(dat.sub)[1]="Peptide"
    sample_size = ncol(dat.sub)-2
    
    m = reshape2::melt(dat.sub) 
    m$PSM_id =  rep(seq(1,nrow(dat.sub)),sample_size)
    ggplot(m, aes(x=variable,y=value))+
        geom_point()+
        geom_line(aes(group=PSM_id,col=Peptide))+
        theme(axis.text.x = element_text(angle = 90,hjust = 1))+
        ggtitle(dat.sub[1,2])+
        xlab("samples")+
        ylab("log2 intensity")+
        theme(plot.title = element_text(hjust = 0.5))
}


medianSummary <- function(dat,group_col=2,ref_col) {
    dat.ratio = dat
    if (length(ref_col)>1){
        dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - 
            rowMeans(dat.ratio[,ref_col],na.rm = TRUE)  
    }else {
        dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - dat.ratio[,ref_col]
    }
    
    dat.summary = plyr::ddply(dat.ratio,colnames(dat)[group_col],
    function(x)matrixStats::colMedians(as.matrix(x[,3:ncol(dat)]),na.rm = TRUE))
    
    colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
    dat.new = dat.summary[,-1]
    rownames(dat.new) = dat.summary[,1]
    return (dat.new)
}


medianPolish <- function (m) {
    dat = medpolish(m,trace.iter=FALSE)$col
    return (dat)
}

medpolishSummary <- function(dat,group_col=2) {
    dat.summary = plyr::ddply(dat,colnames(dat)[group_col],
    function(x) medianPolish(as.matrix(x[,3:ncol(dat)])))
    
    colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
    dat.new = dat.summary[,-1]
    rownames(dat.new) = dat.summary[,1]
    return (dat.new)
}

equalMedianNormalization <- function(dat, optional_outfile) {
    sizefactor = matrixStats::colMedians(as.matrix(dat),na.rm = TRUE)
    dat.nm = sweep(dat,2,sizefactor)
    if (!missing(optional_outfile)) {
      write.table(data.frame(channels=colnames(x), medians=sizefactor), optional_outfile, sep='\t', header=T)
    }
    return (dat.nm)
}

medianSweeping <- function(dat,group_col=2,channelmedian_outfile) {
    dat.ratio = dat
    dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - 
        matrixStats::rowMedians(as.matrix(dat.ratio[,3:ncol(dat)]),na.rm = TRUE)
    dat.summary = plyr::ddply(dat.ratio,colnames(dat)[group_col],
    function(x)matrixStats::colMedians(as.matrix(x[,3:ncol(dat)]),na.rm = TRUE))
    
    colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
    dat.new = dat.summary[,-1]
    rownames(dat.new) = dat.summary[,1]
    
    dat.nm = equalMedianNormalization(dat.new, optional_outfile=channelmedian_outfile)
    return (dat.nm)
}

farmsMethod <- function(df){
    if (nrow(df)==1){
        dat = log2(as.matrix(df))
    }else {dat = farms::generateExprVal.method.farms(
        as.matrix(na.omit(df)))$expr}
    return (dat)
}

farmsSummary <- function(dat,group_col=2) {
    dat.log = plyr::ddply(dat,colnames(dat)[group_col],
    function(x) farmsMethod(x[,3:ncol(dat)]))
    
    colnames(dat.log)[2:ncol(dat.log)]=colnames(dat)[3:ncol(dat)]
    
    dat.ratio = dat.log
    dat.ratio[,2:ncol(dat.ratio)]= dat.log[,2:ncol(dat.log)] - 
    matrixStats::rowMedians(as.matrix(dat.log[,2:ncol(dat.log)]),na.rm = TRUE)
    
    dat.summary = dat.ratio[,-1]
    rownames(dat.summary) = dat.ratio[,1]
    
    return (dat.summary)
}
