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
    names(fit$count) = rownames(fit$coefficients)
    output<-fit

    if (fit.method == "loess"){
        x=log2(fit$count)
        loess.model = stats::loess(logVAR~x,span = 0.75)
        y.pred = stats::predict(loess.model)
        output$loess.model = loess.model
    }else if (fit.method == "nls"){
        x=fit$count
        y=fit$sigma^2
        nls.model =  stats::nls(y~a+b/x,start = (list(a=0.1,b=0.05)))
        y.pred = log(stats::predict(nls.model))
        output$nls.model = nls.model
    }
    
    egpred<-y.pred-digamma(df/2)+log(df/2)
    
    myfct<- (eg-egpred)^2 - trigamma(df/2)
    
    mean.myfct<-mean(myfct,na.rm=TRUE)
    priordf<-vector(); testd0<-vector()
    for (i in seq(1,numgenes*10)) {
        testd0[i]<-i/10
        priordf[i]= abs(mean.myfct-trigamma(testd0[i]/2))
        if (i>2) {
            if (priordf[i-2]<priordf[i-1]) { break }
        }
    }
    d0<-testd0[match(min(priordf),priordf)]
    #print("Prior degrees freedom found")
    
    s02<-exp(egpred + digamma(d0/2) - log(d0/2))
    
    post.var<- (d0*s02 + df*fit$sigma^2)/(d0+df)
    post.df<-d0+df
    # sca.t and scc.p stands for spectra count adjusted t and p values.  
    sca.t<-as.matrix(fit$coefficients[,coef_col]/(fit$stdev.unscaled[,coef_col]
    *sqrt(post.var)))
    sca.p<-as.matrix(2*(1-pt(abs(sca.t),post.df)))
    
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
    results.table$sca.adj.pval = stats::p.adjust(results.table$sca.P.Value,
    method = "BH")
    results.table = results.table[order(results.table$sca.P.Value), ]
}

plotFitCurve <- function (fit,fit.method="loess",type = "boxplot",
    xlab="peptide count",main="") {
    x = fit$count
    y = fit$sigma^2
    
    if (fit.method=="nls"){
        model = fit$nls.model
        
        if (type=="scatterplot"){
            plot(x,log(y),ylab="log(Variance)",main= main)
            
            y.pred <- log(stats::predict(model))
            k = order(x)
            lines((x[k]),y.pred[k],col='red',lwd=3)
            
        }else if (type=="boxplot"){
            df.temp = data.frame(pep_count =x, variance = y )
            df.temp.filter = df.temp[df.temp$pep_count<21,]
            boxplot(log(variance)~pep_count,df.temp.filter, xlab=xlab,
                    ylab = "log(Variance)", main=main)
            
            y.pred <- log(stats::predict(model,data.frame(x=seq(1,20))))
            lines(seq(1,20),y.pred,col='red',lwd=3)
        }else{stop("only scatterplot and boxplot are supported")}
        
    }else if (fit.method=="loess"){
        model = fit$loess.model
        if (type=="scatterplot"){
            plot(x,log(y),ylab="log(Variance)",main= main)
            
            y.pred <- stats::predict(model)
            k = order(x)
            lines((x[k]),y.pred[k],col='red',lwd=3)
            
        }else if (type=="boxplot"){
            df.temp = data.frame(pep_count =x, variance = y )
            df.temp.filter = df.temp[df.temp$pep_count<21,]
            boxplot(log(variance)~pep_count,df.temp.filter, xlab=xlab,
                    ylab = "log(Variance)", main=main)
            
            y.pred <- stats::predict(model,log2(seq(1,20)))
            lines(seq(1,20),y.pred,col='red',lwd=3)
        }else{stop("only scatterplot and boxplot are supported")}
    }  
}

peptideProfilePlot <- function(dat,col=2,gene){
    
    dat.sub = dat[dat[,col]==gene,]
    
    colnames(dat.sub)[1]="Peptide"
    sample_size = ncol(dat.sub)-2
    
    m = reshape2::melt(dat.sub) 
    m$PSM_id =  rep(seq(1,nrow(dat.sub)),sample_size)
    ggplot2::ggplot(m, ggplot2::aes(x=variable,y=value))+
        ggplot2::geom_point()+
        ggplot2::geom_line(ggplot2::aes(group=PSM_id,col=Peptide))+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
        hjust = 1))+ggplot2::ggtitle(dat[1,2])+
        ggplot2::xlab("samples")+
        ggplot2::ylab("log2 intensity")+
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}


medianSummary <- function(dat,group_col=2,ref_col) {
    dat.ratio = dat
    dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - 
        rowMeans(dat.ratio[,ref_col],na.rm = TRUE)
    dat.summary = plyr::ddply(dat.ratio,colnames(dat)[group_col],
    function(x)matrixStats::colMedians(as.matrix(x[,3:ncol(dat)]),na.rm = TRUE))
    
    colnames(dat.summary)[2:ncol(dat.summary)]=colnames(dat)[3:ncol(dat)]
    dat.new = dat.summary[,-1]
    rownames(dat.new) = dat.summary[,1]
    return (dat.new)
}


medianPolish <- function (m) {
    dat = stats::medpolish(m,trace.iter=FALSE)$col
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


equalMedianNormalization <- function(dat) {
    sizefactor = matrixStats::colMedians(as.matrix(dat),na.rm = TRUE)
    dat.nm = sweep(dat,2,sizefactor)
    return (dat.nm)
}

medianSweeping <- function(dat,group_col=2) {
    dat.ratio = dat
    dat.ratio[,3:ncol(dat)] = dat.ratio[,3:ncol(dat)] - 
        matrixStats::rowMedians(as.matrix(dat.ratio[,3:ncol(dat)]),na.rm = TRUE)
    dat.summary = plyr::ddply(dat.ratio,colnames(dat)[group_col],
    function(x)matrixStats::colMedians(as.matrix(x[,3:ncol(dat)]),na.rm = TRUE))
    
    colnames(dat.summary)[2:ncol(dat.summary)] = colnames(dat)[3:ncol(dat)]
    dat.new = dat.summary[,-1]
    rownames(dat.new) = dat.summary[,1]
    
    dat.nm = equalMedianNormalization(dat.new)
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
