### lifetime.mle
## failure.threshold is the percentage
################################################################################
addt.fit=function(formula, data, initial.val=100, proc="All", failure.threshold, time.rti=100000, method="Nelder-Mead", subset, na.action, starts=NULL,fail.thres.vec=c(70,80), semi.control= list(cor=F,...) ,...)
{
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  if (m[1]==0) stop("a formula argument is required")
  mf <- mf[c(1, m)]
  mf[[1]] <- as.name("model.frame")
  mdat <- eval(mf, parent.frame())
  n=nrow(mdat)
  mt <- attr(mdat, "terms")

  #y <- model.extract(mdat, "response")
  #X <- model.matrix(mt, mdat)
  #ll <- ncol(X)

  #browser()
  names(mdat)=c("Response", "Time", "Temp")
  dat0=as.data.frame(cbind(Temp=mdat[,"Temp"], Time=mdat[,"Time"], Response=mdat[,"Response"]))


  LS.obj=ML.obj=SemiPara.obj=NULL

  if(proc=="LS"|proc=="All"){
    # LSA
    LS.obj=lsa.fit(dat=dat0, initial.val=initial.val,failure.threshold=failure.threshold, time.rti=time.rti)
  }

  if(proc=="ML"|proc=="All"){
    # MLA
    if(is.null(starts))
    {
      #fail.thres=50
      #ok=0
      #while(ok<1){
      #  temp=try(lsa.fit(dat=dat0, initial.val=initial.val,failure.threshold=fail.thres, time.rti=time.rti), silent=T)
      #  ok=ifelse(class(temp)=="try-error", 0, 1)
      #  fail.thres=fail.thres+10
      #}
      #fail.thres1=fail.thres-10
      #fail.thres2=ifelse(fail.thres>100, 99, fail.thres)

      starts=addt.mle.initial.val(dat=dat0, time.rti=time.rti, initial.val=initial.val, failure.threshold=fail.thres.vec[1], failure.threshold1=fail.thres.vec[2])
    }
    ML.obj=lifetime.mle(dat=dat0, minusloglik=minus.loglik.kinetics, starts=starts, method = method , control=list(maxit=100000))
    ML.obj=mle.transform(obj=ML.obj)
    ML.obj=c(ML.obj, rti=true.ti.compute(pars=ML.obj$coef,failure.threshold=failure.threshold,initial.val=initial.val,time.rti=time.rti)$rti, beta = true.ti.compute(pars=ML.obj$coef,failure.threshold=failure.threshold,initial.val=initial.val,time.rti=time.rti)$coef)
  }

  if(proc=="SemiPara"|proc=="All"){
    #print("Succeeded")
    #dat=data.frame(TempC=c(50, mdat[,"Temp"]), TimeH=c(0, mdat[,"Time"]),Response=c(100, mdat[,"Response"]))
    dat=data.frame(TempC= mdat[,"Temp"], TimeH=mdat[,"Time"],Response=mdat[,"Response"])
    #print(mdat)
    #print(dat)

    ##############
    # initial value of beta
    # pineapple: This could cause problem
    #beta.ini(dat)
    #beta.ini(dat, smooth.spline=T)

    colnames(mdat) <- c("response", "time", "temp")


    if (missing(semi.control))
    {
      res=monotone.bspline.fit.nocor.knotselection.output(dat,cov.fun.type="cov.fun.REMLc")

      bestmodel=bspline.fit.addt.nocor(test.dat=dat, n.knots=NULL, knots=as.numeric(unlist(strsplit(as.character(res[,"Knots"]), ",")))
                                       , degree=res[,"Order"],  cov.fun.type="cov.fun.REMLc")


      TI.semi=try(TI.bspline.nocor(dat=dat, model.fit.obj=bestmodel,dd=time.rti, failure.threshold = failure.threshold/100 ), silent = T)

      #TI.values <- semi.para.TI(bestmodel, dat, threshold = failure.threshold, td = time.rti)

      #print(TI.semi)

      SemiPara.obj <-  c(bestmodel, TI= TI.semi)

      time.rti <- 1e+05
      initial.val <- NULL
      failure.threshold <- failure.threshold
      dat0 <- dat

    }
    else
    {
      if(semi.control$cor=="TRUE")
      {

        #print("corr succeed")
        # knot selection result
        res=monotone.bspline.fit.knotselection.output(dat,cov.fun.type="cov.fun.REMLc")

        bestmodel=bspline.fit.addt(test.dat=dat, n.knots=NULL, knots=as.numeric(unlist(strsplit(as.character(res[,"Knots"]), ",")))
                                   , degree=res[,"Order"],  cov.fun.type="cov.fun.REMLc")

        TI.semi=try(TI.bspline.nocor(dat=dat, model.fit.obj=bestmodel, dd=time.rti, failure.threshold = failure.threshold/100), silent = T)

        SemiPara.obj <-  c(bestmodel, TI= TI.semi)

        time.rti <- 1e+05
        initial.val <- NULL
        failure.threshold <- failure.threshold
        dat0 <- dat

      }
    }

  }
  addt.fit=list(LS.obj=LS.obj, ML.obj=ML.obj,SemiPara.obj=SemiPara.obj, dat=dat0, time.rti=time.rti, initial.val=initial.val,
                failure.threshold=failure.threshold)
  class(addt.fit)="addt.fit"

  return(addt.fit)
}

################################################################################
power.exponential.fun=function(tt, temp, alpha, beta0,beta1,gamma)
{
  x=11605/(temp+273.16)
  bb=exp(beta0+beta1*x)
  res=alpha*exp(-(tt/bb)^gamma)
  return(res)
}

################################################################################
true.ti.compute=function(pars, failure.threshold, initial.val, time.rti)
{
  #browser()
  pp=failure.threshold/100
  beta0=pars[2]/log(10)
  beta1=pars[3]/log(10)
  gamma=exp(pars[4])
  gamma.term=(1/gamma)*log((1-pp)/pp)/log(10)
  dK=273.16

  rti=beta1/(log10(time.rti)-beta0-gamma.term)-dK

  beta0.log10=beta0+gamma.term
  beta1.log10=beta1


  names(beta0.log10)="beta0"
  names(beta1.log10)="beta1"


  res=list(coef=c(beta0.log10, beta1.log10),rti=rti,time.rti=time.rti,failure.threshold=failure.threshold)

  return(res)
}

####################################################################################
print.addt.fit=function(x,...){
  obj=x
  if(!is.null(obj$LS.obj)){
    cat("Least Squares Approach: \n")
    print(round(x$LS.obj$coef0, 4))
    cat("est.TI:" , round(obj$LS.obj$rti), "\n")
  }

  if(!is.null(obj$ML.obj)){
    cat("\n", "Maximum Likelihood Approach:", "\n")

    alpha=exp(obj$ML.obj$coef[1])
    beta0=obj$ML.obj$coef[2]
    beta1=obj$ML.obj$coef[3]
    gamma=exp(obj$ML.obj$coef[4])

    sigma=sqrt(exp(obj$ML.obj$coef[5])^2+exp(obj$ML.obj$coef[6])^2)
    rho=exp(obj$ML.obj$coef[6])^2/(sigma^2)

    coef=c(alpha, beta0, beta1, gamma, sigma, rho)
    names(coef)=c("alpha", "beta0", "beta1", "gamma", "sigma", "rho")

    cat("Call:\n")
    print(obj$ML.obj$call)
    cat("\n")
    cat("Coefficients:\n")
    print(round(coef,4))
    cat("est.TI:", round(obj$ML.obj$rti), "\n")
    cat("\n")
    cat("Loglikelihod:\n")
    print(as.numeric(-obj$ML.obj$min))
  }

  if(!is.null(obj$SemiPara.obj)){
    cat("\n", "Semi-Parametric Approach:", "\n")
    betahat <- obj$SemiPara.obj$betahat
    print("beta estimates:")
    print(betahat)
    rho <- obj$SemiPara.obj$fit.monotone.bspline$coef[13]
    print("rho")
    print(rho)
    knots <- obj$SemiPara.obj$spline.inf$knots
    print("knots")
    print(knots)
    Loglike <- obj$SemiPara.obj$fit.monotone.bspline$loglik
    print("Log-likelihood")
    print(Loglike)
    print("TI")
    TI.semi <- obj$SemiPara.obj$TI
    print(TI.semi)
  }
}
####################################################################################
summary.addt.fit=function(object,...)
{
  obj=object
  if(!is.null(obj$ML.obj)){
    alpha=exp(obj$ML.obj$coef[1])
    beta0=obj$ML.obj$coef[2]
    beta1=obj$ML.obj$coef[3]
    gamma=exp(obj$ML.obj$coef[4])

    sigma=sqrt(exp(obj$ML.obj$coef[5])^2+exp(obj$ML.obj$coef[6])^2)
    rho=exp(obj$ML.obj$coef[6])^2/(sigma^2)

    sigma.dev=deriv(~sqrt(exp(x)^2+exp(y)^2),c("x","y"), function(x,y){})
    rho.dev=deriv(~exp(y)^2/(exp(x)^2+exp(y)^2),c("x","y"), function(x,y){})


    ll=length(obj$ML.obj$coef)

    coef.dev=diag(rep(1, ll))
    coef.dev[1,1]=alpha
    coef.dev[4,4]=gamma
    coef.dev[(ll-1),(ll-1):ll]=attr(sigma.dev(obj$ML.obj$coef[5], obj$ML.obj$coef[6]), "gradient")
    coef.dev[ll,(ll-1):ll]=attr(rho.dev(obj$ML.obj$coef[5], obj$ML.obj$coef[6]), "gradient")

    vcov=coef.dev%*%obj$ML.obj$vcov%*%t(coef.dev)

    coef=c(alpha, beta0, beta1, gamma, sigma, rho)
    names(coef)=c("alpha", "beta0", "beta1", "gamma", "sigma", "rho")

    mat=matrix(0, ll, 4)
    mat[,1]=coef
    mat[,2]=sqrt(diag(vcov))
    mat[,3]=mat[,1]*exp(-1.96*mat[,2]/mat[,1])
    mat[,4]=mat[,1]*exp(+1.96*mat[,2]/mat[,1])

    mat[c(2:3,ll),3]=mat[c(2:3,ll),1]-1.96*mat[c(2:3,ll),2]
    mat[c(2:3,ll),4]=mat[c(2:3,ll),1]+1.96*mat[c(2:3,ll),2]
    colnames(mat)=c("mean", "std", "95% Lower", "95% Upper")
    rownames(mat)=names(coef)

    ti.CI=addt.confint.ti.mle(obj=obj, conflevel=0.95)
    names(ti.CI)=c("est", "std", "95% Lower", "95% Upper")

    beta0 <- round(obj$ML.obj$beta.beta0 , 4)
    beta1 <- round(obj$ML.obj$beta.beta1 , 4)

    Temperature_time <- c(beta0, beta1)
    names(Temperature_time) <- c("beta0", "beta1")

    obj=c(obj, list(coef.mle.mat=mat), list(ti.CI=ti.CI), list(Temperature_time =Temperature_time))
  }

  if(!is.null(obj$SemiPara.obj)){
      betahat <- obj$SemiPara.obj$betahat
      rho <- as.numeric(obj$SemiPara.obj$fit.monotone.bspline$coef["rho"])
      knots <- as.numeric(obj$SemiPara.obj$spline.inf$knots)
      Boundary <- obj$SemiPara.obj$spline.inf$Boundary.knots
      Loglike <- obj$SemiPara.obj$fit.monotone.bspline$loglik
      aic <- obj$SemiPara.obj$aic
      aicc <- obj$SemiPara.obj$aicc
      TI.semi <- obj$SemiPara.obj$TI.TI
      beta0 <- obj$SemiPara.obj$TI.beta0
      beta1 <- obj$SemiPara.obj$TI.beta1

      TI.value <- c(TI.semi, beta0, beta1)
      names(TI.value) <- c("TI.semi", "beta0", "beta1")

      Semi.knots <- c(knots, Boundary)
      names(Semi.knots) <- c("knots", "Left Boundary", "Right Boundary")

      parameter <- c(betahat, rho)
      names(parameter) <- c("betahat", "rho")

      evaluation <- c(Loglike , aic, aicc)
      names(evaluation) <- c("Loglikelihood","AIC","AICC")

      obj = c(obj , betahat, knots, Loglike, aic, aicc, rho, TI.semi, beta0, beta1)
  }

  class(obj)="summary.addt.fit"
  return(obj)
}

####################################################################################
print.summary.addt.fit=function(x,...)
{
  if(!is.null(x$LS.obj)){
    cat("Least Squares Approach: \n")
    print(round(x$LS.obj$coef0, 4))
    cat("est.TI:" , round(x$LS.obj$rti), "\n")
    cat("Interpolation time: \n")
    print(x$LS.obj$interp.mat)
  }

  if(!is.null(x$ML.obj)){
    cat("\n", "Maximum Likelihood Approach:", "\n")
    cat("Call:\n")
    print(x$ML.obj$call)
    cat("\n")
    cat("Parameters:\n")
    print(round(x$coef.mle.mat, 4))
    cat("\n")
    cat("Temperature-Time Relationship: \n")
    print(x$Temperature_time)
    cat("\n")
    #print(c("beta0", "beta1"))
    #print(c(round(x$ML.obj$beta.beta0), round(x$ML.obj$beta.beta1)))
    cat("TI: \n")
    print(round(x$ti.CI, 4))
    cat("\n")
    cat("Loglikelihood:\n")
    print(as.numeric(-x$ML.obj$min))
  }

  if(!is.null(x$SemiPara.obj)){
    betahat <- round(x$SemiPara.obj$betahat,3)
    rho <- round(as.numeric(x$SemiPara.obj$fit.monotone.bspline$coef["rho"]) ,3)
    knots <- round(as.numeric(x$SemiPara.obj$spline.inf$knots) , 3)
    Boundary <-round(x$SemiPara.obj$spline.inf$Boundary.knots ,3)
    Loglike <- round(x$SemiPara.obj$fit.monotone.bspline$loglik,3)
    aic <- round(x$SemiPara.obj$aic ,3)
    aicc <- round(x$SemiPara.obj$aicc ,3)
    TI.semi <- round(x$SemiPara.obj$TI.TI ,3)
    beta0 <- round(x$SemiPara.obj$TI.beta0 , 3)
    beta1 <- round(x$SemiPara.obj$TI.beta1 , 3)
    TI.value <- c(TI.semi, beta0, beta1)
    names(TI.value) <- c("TI.semi", "beta0", "beta1")
    Semi.knots <- c(Boundary[1], knots, Boundary[2])
    names(Semi.knots) <- c("Left Boundary",rep("knots", length(knots)), "Right Boundary")

    if(!is.null(x$SemiPara.obj$fit.monotone.bspline$rho))
    {
      parameter <- c(betahat, rho)
      names(parameter) <- c("betahat", "rho")
    }
    if(is.null(x$SemiPara.obj$fit.monotone.bspline$rho))
    {
      parameter <- c(betahat)
      names(parameter) <- c("betahat")
    }
    evaluation <- c(Loglike , aicc)
    names(evaluation) <- c("Loglikelihood","AICC")
    cat("\n", "Semi-Parametric Approach:", "\n")
    cat("\n")
    cat("Parameters Estimates: \n")
    print(parameter)
    cat("\n")
    cat("TI estimates: \n")
    print(TI.value)
    cat("\n")
    cat("Model Evaluations: \n")
    print(evaluation)
    cat("\n")
    cat("B-spline: \n")
    print(Semi.knots)
  }

}

####################################################################################
plot.addt.fit=function(x, type,...){
  if(type=="data"){
    addt.data.plot(x$ML.obj$dat)
  }

  if(type=="ML"){
    addt.par.fitted.plot(x$ML.obj,mean.fun="kinetics",zero.correction=F,timeline.show=F,legend.pos=1)
  }

  if(type=="LS"){
    polynomial.interpolation(dat=x$dat,initial.val=x$initial.val,failure.threshold=x$failure.threshold)
  }
  if(type=="SEMI"){
    data0 <- x$dat
    colnames(data0) <- c("Temp", "Time", "Response")
    addt.data.plot(data0, xlab="Time in Weeks", ylab=expression(paste(Log, " of Strength (Newton)")), main="", legend.pos=1)

    colnames(data0) <- c("temp", "time", "response")
    monotone.bspline.data0.fitted=aggregate(x$SemiPara.obj$fit.monotone.bspline$yhat, by=list(data0$temp, data0$time), mean)
    names(monotone.bspline.data0.fitted)=c("temp", "time", "fitted")

    temp.uni=unique(data0$temp)

    #monotone.bspline.data0.fitted=monotone.bspline.fit.obj$fit.mat
    monotone.bspline.baseline=monotone.bspline.data0.fitted[monotone.bspline.data0.fitted$time==0, "fitted"]

    for(i in 1:length(temp.uni)){
      lines(c(0, monotone.bspline.data0.fitted[monotone.bspline.data0.fitted$temp==temp.uni[i], "time"]),
            c(monotone.bspline.baseline, monotone.bspline.data0.fitted[monotone.bspline.data0.fitted$temp==temp.uni[i], "fitted"]),
            lty=i, col=i)
    }

  }
}

################################################################################
addt.par.fitted.plot=function(obj,mean.fun,zero.correction=F,timeline.show=T,legend.pos=1,initial.val=100,...)
{
  dat=obj$dat
  addt.data.plot(dat=dat,zero.correction=zero.correction,timeline.show=timeline.show,legend.pos=legend.pos,...)

  mtime=max(dat[,"Time"])
  stime=ifelse(zero.correction,1,0)
  ww=seq(stime,mtime,,1000)

  coefs=obj$coef

  switch(mean.fun,
         "kinetics"={
           mfun=kinetics.fun
           alpha=exp(coefs[1])
           beta0=coefs[2]
           beta1=coefs[3]
           gamma=exp(coefs[4])
         },
         "kinetics0"={
           mfun=kinetics.fun
           alpha=initial.val
           beta0=coefs[1]
           beta1=coefs[2]
           gamma=exp(coefs[3])
         },
         "power.exponential"={
           mfun=power.exponential.fun
           alpha=exp(coefs[1])
           beta0=coefs[2]
           beta1=coefs[3]
           gamma=exp(coefs[4])
         },
         "power.exponential0"={
           mfun=power.exponential.fun
           alpha=initial.val
           beta0=coefs[1]
           beta1=coefs[2]
           gamma=exp(coefs[3])
         },

  )


  cc2=unique(dat[,"Temp"])

  for(i in 1:(length(cc2)))
  {
    yy=mfun(tt=ww, temp=cc2[i], alpha=alpha, beta0=beta0, beta1=beta1, gamma=gamma)
    lines(ww,yy,lwd=3,col=i)
  }

}
################################################################################
addt.data.plot=function(dat,zero.correction=F,timeline.show=T,legend.pos=1,xlab,ylab,...)
{
  cc=dat[,"Temp"]
  cc=as.factor(cc)
  cc1=as.numeric(cc)

  if(zero.correction)
  {
    dat[dat[,"Time"]==0,"Time"]=1
  }

  #browser()

  mm=max(dat[,"Response"])
  mmin=min(dat[,"Response"])
  mtime=max(dat[,"Time"])
  mintime=min(dat[,"Time"])

  if(!zero.correction)
  {
    mintime=0
  }

  plot(dat[,"Time"], dat[,"Response"],col=cc1,xlim=c(mintime,1.05*mtime),ylim=c(.8*mmin,1.2*mm),pch=cc1,xlab="Time", ylab="Responses Values",las=1,lwd=2,...)

  #browser()

  if(legend.pos==1)
  {
    legend(0.7*mtime,1.23*mm,paste(levels(cc),"C",sep=""),col=unique(cc1),pch=unique(cc1),lwd=2,lty=NA,bty="n")
  }

  if(legend.pos==2)
  {
    legend(mintime,.5*mm,paste(levels(cc),"C",sep=""),col=unique(cc1),pch=unique(cc1),lwd=2,lty=NA,bty="n")
  }


  if(timeline.show)
  {
    tt=unique(dat[,"Time"])
    abline(v=tt,col="grey",lty=2)
    text(tt,rep(.9*mmin,length(tt)),tt,cex=.9)
  }
}

################################################################################
lsa.fit=function(dat, initial.val, failure.threshold, time.rti)
{
  dK=273.16

  # polynominal fit to the data
  dat0=polynomial.interpolation(dat=dat,initial.val=initial.val,failure.threshold=failure.threshold,plot=F,plot.all=F)
  dat0=dat0[!is.na(dat0[,2]),]

  # least square fit to the data
  fit0=lm(I(log10(dat0[,2]))~I(1/(dat0[,1]+dK)))
  coef0=fit0$coef
  names(coef0)=c("beta0", "beta1")
  rti0=coef0[2]/(log10(time.rti)-coef0[1])-dK

  res=list(coef0=coef0, interp.mat=dat0, rti0=rti0)
  class(res)="addt.fit.lsa"
  return(res)

}
################################################################################
lifetime.mle=function(dat, minusloglik, starts, method = "BFGS",hessian = TRUE,...)
{
  call=match.call()
  f = function(p) {
    minusloglik(dat,p)
  }
  oout = optim(starts, f, method = method, hessian = hessian,...)#,control=list(trace=T))
  coef = oout$par
  #browser()
  if(hessian)
  {
    vcov =solve(oout$hessian)
  }else{
    vcov=NULL
  }
  min = oout$value
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat,minusloglik=minusloglik))
}

################################################################################
# sample mean at each combo of time and temp
addt.mean.summary=function(dat)
{
  aa=tapply(dat[,3],dat[,1:2],"mean")

  bb=rownames(aa)
  bb=as.numeric(bb)
  cc=colnames(aa)
  cc=as.numeric(cc)

  mm=dim(aa)[1]
  nn=dim(aa)[2]

  res=data.frame(Temp=rep(bb,nn), Time=rep(cc,rep(mm,nn)), Response=as.vector(aa))
  res=res[!is.na(res[,3]),]
  res=res[order(res[,1],res[,2]),]
  rownames(res)=NULL
  return(res)
}
################################################################################
addt.confint.ti.mle=function(obj,conflevel)
{
  if(is.null(obj$ML.obj)) {
    stop("this function can only use for maximum likelihood approach")
  }

  #browser()

  failure.threshold=obj$failure.threshold
  initial.val=obj$initial.val
  time.rti=obj$time.rti

  obj=obj$ML.obj
  Sigma.part=obj$vcov[c(2:4),c(2:4)]
  beta0.hat=obj$coef[2]
  beta1.hat=obj$coef[3]
  pp=failure.threshold/100
  gamma.hat=exp(obj$coef[4])
  gamma.term=(1/gamma.hat)*log((1-pp)/pp)
  ti.hat=true.ti.compute(pars=obj$coef,failure.threshold=failure.threshold,initial.val=initial.val,time.rti=time.rti)$rti

  beta0.log10=beta0.hat/log(10)
  beta1.log10=beta1.hat/log(10)
  gamma.log10=gamma.hat*log(10)
  gamma.term.log10=log((1-pp)/pp)/gamma.log10

  den=log10(time.rti)-beta0.log10-gamma.term.log10
  partial.b0=beta1.log10/(den^2)
  partial.b1=1/den
  partial.gamma=-(beta1.log10/(den^2))*(gamma.term.log10/gamma.log10)

  Sigma.part.trans=diag(c(1/log(10),1/log(10),log(10)*gamma.hat))%*%Sigma.part%*%diag(c(1/log(10),1/log(10),log(10)*gamma.hat))

  Sigma.ti=t(c(partial.b0,partial.b1,partial.gamma))%*%Sigma.part.trans%*%c(partial.b0,partial.b1,partial.gamma)
  CI=c(ti.hat-qnorm(0.5+conflevel/2)*sqrt(Sigma.ti),ti.hat+qnorm(0.5+conflevel/2)*sqrt(Sigma.ti))
  #cat(100*(1-conflevel), "% CI is ", "(", round(CI[1], 3), ",", round(CI[2], 3), ") \n", sep="")

  CI[1]=ifelse(CI[1]<0, 0, CI[1])

  res=c(ti.hat, sqrt(Sigma.ti), CI)
  names(res)=c("est.", "s.e.", "lower", "upper")
  return(res)
}

################################################################################
addt.predint.ybar.mle=function(obj,conflevel,num.fut.obs=5,temp,tt)
{
  if(is.null(obj$ML.obj)) {
    stop("this function can only use for maximum likelihood approach")
  }

  obj=obj$ML.obj
  alpha.hat=exp(obj$coef[1])
  beta0.hat=obj$coef[2]
  beta1.hat=obj$coef[3]
  gamma.hat=exp(obj$coef[4])

  sigma.sq.hat=exp(obj$coef[5])^2+exp(obj$coef[6])^2
  rho.hat=exp(obj$coef[6])^2/sigma.sq.hat

  Sigma.betaf=diag(c(alpha.hat,1,1,gamma.hat))%*%obj$vcov[c(1:4),c(1:4)]%*%diag(c(alpha.hat,1,1,gamma.hat))

  ybar.hat=kinetics.fun(tt=tt, temp=temp, alpha=alpha.hat, beta0=beta0.hat, beta1=beta1.hat, gamma=gamma.hat)
  partial.vec=first.dev.kinetics.fun(tt=tt, temp=temp, alpha=alpha.hat, beta0=beta0.hat, beta1=beta1.hat, gamma=gamma.hat)

  Sigma.ybar=sigma.sq.hat*(rho.hat+(1-rho.hat)/num.fut.obs)+
    t(partial.vec)%*%Sigma.betaf%*%partial.vec

  return(c(ybar.hat-qnorm(0.5+conflevel/2)*sqrt(Sigma.ybar),ybar.hat+qnorm(0.5+conflevel/2)*sqrt(Sigma.ybar)))
}


################################################################################
first.dev.kinetics.fun=function(tt, temp, alpha, beta0, beta1, gamma){
  x=1/(temp+273.16)

  partial.beta0 = alpha*gamma*(tt^gamma)*
    (1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-2)*exp(-(beta0+beta1*x)*gamma)
  partial.beta1 = alpha*x*gamma*(tt^gamma)*
    (1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-2)*exp(-(beta0+beta1*x)*gamma)
  partial.gamma = -alpha*(tt^gamma)*
    (1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-2)*exp(-(beta0+beta1*x)*gamma)*
    (log(tt)-(beta0+beta1*x))
  partial.alpha=(1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-1)

  return(c(partial.alpha, partial.beta0, partial.beta1, partial.gamma))
}


################################################################################
# initial.val is specified if no records available at time 0
addt.data.normalization=function(dat, initial.val)
{
  temps=unique(dat[,"Temp"])
  nn=length(temps)

  if(!any(dat[,"Time"]==0))
  {
    tmp1=cbind(Temp=temps, Time=0, Response=100)
    dat[,"Response"]=dat[,"Response"]/initial.val*100
    dat=rbind(tmp1,dat)
    dat=dat[order(dat[,"Temp"],dat[,"Time"]),]
    rownames(dat)=NULL
    res=dat
  }

  if(any(dat[,"Time"]==0))
  {
    initial.val=dat[dat[,"Temp"]==min(temps) & dat[,"Time"]==0,"Response"]

    res=NULL
    for(i in 1:nn)
    {
      xtmp=dat[dat[,"Temp"]==temps[i],]

      if(any(xtmp[,"Time"]==0))
      {
        xtmp[,"Response"]=xtmp[,"Response"]/xtmp[xtmp[,"Time"]==0,"Response"]*100
      }else{
        xtmp[,"Response"]=xtmp[,"Response"]/initial.val*100
        xtmp=rbind(xtmp,c(temps[i],0,100))
      }
      res=rbind(res,xtmp)
    }

    res=res[order(res[,"Temp"],res[,"Time"]),]
    rownames(res)=NULL
  }
  return(res)
}


################################################################################
polynomial.interpolation=function(dat,initial.val=100,failure.threshold=80,plot.all=T,plot=T)
{
  dat=addt.mean.summary(dat=dat)
  dat=addt.data.normalization(dat=dat, initial.val=initial.val)

  temps=unique(dat[,"Temp"])
  nn=length(temps)
  tres=rep(NA,nn)

  if(plot)
  {
    if(plot.all)
    {
      plot(0,0,type="n",xlim=c(0,max(dat[,"Time"])),ylim=c(0,150),ylab="Response (Relative %)",xlab="Time (hours)")
      abline(h=failure.threshold)
    }else{
      par(mfrow=c(ceiling(nn/round(sqrt(nn))),round(sqrt(nn))))
    }
  }

  for(i in 1:nn)
  {
    idx=(dat[,"Temp"]==temps[i])

    if(sum(idx)>=3)
    {
      yy=dat[idx,"Response"]
      xx=dat[idx,"Time"]
      #browser()
      if(sum(idx)==3 & min(yy)<failure.threshold)
      {
        afit=lm(yy~I(xx)+I(xx^2))
        coefs=afit$coef

        tt=seq(0,max(xx),,1000)
        yyhat=coefs[1]+coefs[2]*tt+coefs[3]*tt^2


        ctmp=polyroot(c(coefs[1]-failure.threshold,coefs[2:3]))

        #print(ctmp)
        ctmp1=Re(ctmp[abs(Im(ctmp))<1e-5])
        ctmp2=ctmp1[ctmp1>0 & ctmp1<=max(xx)]
        if(length(ctmp2)>0)
        {
          tres[i]=min(ctmp2)
        }else{
          ff=approx(yyhat, tt, failure.threshold)
          tres[i]=min(ff$y)
        }
        #browser()

        if(plot)
        {
          if(!plot.all)
          {
            plot(tt,yyhat,type="l",ylim=c(0,100),xlim=c(0,max(xx)),xlab="Time",ylab="Response (Relative %)")
            abline(h=failure.threshold)
            points(xx,yy)
          }else{
            lines(tt,yyhat,lwd=2,col=i)
            points(xx,yy,col=i,pch=i,lwd=2)
          }
          points(tres[i],failure.threshold,pch=13,col=2,lwd=2)
        }
      }

      if(sum(idx)>3 & min(yy)<failure.threshold)
      {
        afit=lm(yy~I(xx)+I(xx^2)+I(xx^3))
        coefs=afit$coef

        tt=seq(0,max(xx),,1000)
        yyhat=coefs[1]+coefs[2]*tt+coefs[3]*tt^2+coefs[4]*tt^3
        #ff=approx(yyhat, tt, failure.threshold)
        #tres[i]=min(ff$y)
        #browser()
        ctmp=polyroot(c(coefs[1]-failure.threshold,coefs[2:4]))
        #print(ctmp)
        ctmp1=Re(ctmp[abs(Im(ctmp))<1e-5])
        ctmp2=ctmp1[ctmp1>0 & ctmp1<=max(xx)]
        if(length(ctmp2)>0)
        {
          tres[i]=min(ctmp2)
        }else{
          ff=approx(yyhat, tt, failure.threshold)
          tres[i]=min(ff$y)
        }


        if(plot)
        {
          if(!plot.all)
          {
            plot(tt,yyhat,type="l",ylim=c(0,100),xlim=c(0,max(xx)),xlab="Time",ylab="Response (Relative %)")
            abline(h=failure.threshold)
            points(xx,yy)
          }else{
            lines(tt,yyhat,lwd=2,col=i)
            points(xx,yy,col=i,pch=i,lwd=2)
          }
          points(tres[i],failure.threshold,pch=13,col=2,lwd=2)
        }

      }
    }
  }

  res=cbind(temps,tres)
  colnames(res)=c("Temp","Time")

  if(plot.all)
  {
    legend(0.7*max(dat[,"Time"]),150, paste("T_",res[,"Temp"],"=",round(res[,"Time"]),sep="") ,pch=1:nn,col=1:nn,lty=1,lwd=2,bty="n")
  }
  return(res)

}

################################################################################
addt.mle.initial.val=function(dat=dat, time.rti=100000, initial.val=100, failure.threshold=70, failure.threshold1=50)
{
  #browser()
  base.factor=log10(exp(1))
  dK=273.16

  yy=dat[,"Response"]
  tt=dat[,"Time"]
  temp=dat[,"Temp"]

  #initial value for alpha
  alpha=mean(yy[tt==0])
  #initial value for beta0, beta1 and gamma
  tmp.sfit1=lsa.fit(dat=dat,time.rti=time.rti,initial.val=initial.val,failure.threshold=failure.threshold)
  tmp.sfit2=lsa.fit(dat=dat,time.rti=time.rti,initial.val=initial.val,failure.threshold=failure.threshold1)

  pp1=failure.threshold/100
  pp2=failure.threshold1/100

  coef1=tmp.sfit1$coef0
  coef2=tmp.sfit2$coef0

  #initial value for beta1
  beta1=coef1[2]
  sbeta0=log10(time.rti)-beta1/(tmp.sfit2$rti0+dK)

  cc1=log((1-pp1)/pp1)
  cc2=log((1-pp2)/pp2)
  kk1=coef1[1]/base.factor
  kk2=sbeta0/base.factor

  kk3=(kk1-kk2)/(cc1-cc2)
  gamma=1/kk3

  if(gamma<0)
  {
    gamma=1
  }

  beta0=kk1-kk3*cc1
  beta1=beta1/11605/base.factor
  beta0=beta0+beta1*11605/(min(temp)+dK)

  aa=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  des=yy-aa
  sigma=sqrt(mean(des^2))

  start.val0=c(log(alpha), beta0, beta1, log(gamma), log(sigma))
  names(start.val0)=NULL

  #browser()

  tmp.fit=lifetime.mle(dat=dat, minusloglik=minus.loglik.kinetics.no.cor, starts=start.val0, method = "Nelder-Mead", control=list(maxit=10000))


  #cat("-log likelihood without correlation:", tmp.fit$min, "\n")

  tcoef=tmp.fit$coef
  alpha=exp(tcoef[1])
  beta0=tcoef[2]
  beta1=tcoef[3]
  gamma=exp(tcoef[4])

  aa1=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  des1=yy-aa1
  xtmp=tapply(des1,paste(temp,tt),"mean")
  sigma1=sd(xtmp)

  res=c(tcoef,log(sigma1))

  return(res)
}
################################################################################
minus.loglik.kinetics.no.cor=function(dat,pars)
{
  #print(pars)
  yy=dat[,"Response"]
  tt=dat[,"Time"]
  temp=dat[,"Temp"]

  alpha=exp(pars[1])
  beta0=pars[2]
  beta1=pars[3]
  gamma=exp(pars[4])
  sigma=exp(pars[5])
  sigma1=0#exp(pars[6]) #batch level
  ####

  #yy=yy/100
  #A=A/100

  #######
  aa=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  time.rti=yy-aa

  #plot(tt,aa)

  #browser()

  temp.time=paste(dat[,"Temp"],"C",dat[,"Time"],"H",sep="")
  cc=unique(temp.time)
  nn=length(cc)

  ###
  res=0
  for(i in 1:nn)
  {
    #browser()
    idx=(temp.time==cc[i])
    pp=sum(idx)
    SS=diag(sigma^2,pp)+sigma1^2*matrix(1,pp,pp)*as.numeric(pp>1)
    SS.inv=solve(SS)
    tdd=time.rti[idx]
    tres=-0.5*pp*log(2*pi)-0.5*log(det(SS))-0.5*t(tdd)%*%SS.inv%*%tdd
    #print(tres)
    res=res+tres
    #print(round(SS.inv,3))
    #browser()
  }

  res=as.vector(res)
  res=(-1)*res
  #print(res)
  return(res)
}


################################################################################
minus.loglik.kinetics=function(dat,pars)
{
  #print(pars)
  yy=dat[,"Response"]
  tt=dat[,"Time"]
  temp=dat[,"Temp"]

  alpha=exp(pars[1])
  beta0=pars[2]
  beta1=pars[3]
  gamma=exp(pars[4])
  sigma=exp(pars[5])
  sigma1=exp(pars[6]) #batch level
  ####

  #yy=yy/100
  #A=A/100

  #######
  aa=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  time.rti=yy-aa

  #plot(tt,aa)

  #browser()

  temp.time=paste(dat[,"Temp"],"C",dat[,"Time"],"H",sep="")
  cc=unique(temp.time)
  nn=length(cc)

  ###
  res=0
  for(i in 1:nn)
  {
    #browser()
    idx=(temp.time==cc[i])
    pp=sum(idx)
    SS=diag(sigma^2,pp)+sigma1^2*matrix(1,pp,pp)*as.numeric(pp>1)
    SS.inv=solve(SS)
    tdd=time.rti[idx]
    tres=-0.5*pp*log(2*pi)-0.5*log(det(SS))-0.5*t(tdd)%*%SS.inv%*%tdd
    #print(tres)
    res=res+tres
    #print(round(SS.inv,3))
    #browser()
  }

  res=as.vector(res)
  res=(-1)*res
  #print(res)
  return(res)
}

################################################################################
mle.transform=function(obj)
{
  #browser()

  dat=obj$dat
  coefs=obj$coef
  vcovs=obj$vcov

  temp=dat[,"Temp"]
  x=11605/(temp+273.16)
  mx=min(x)

  mat=diag(length(coefs))
  mat[2,3]=-mx
  mat[3,3]=11605
  coefs=mat%*%coefs
  coefs=as.vector(coefs)
  vcovs=mat%*%vcovs%*%t(mat)

  coefs->obj$coef
  vcovs->obj$vcov

  return(obj)
}

################################################################################
kinetics.fun=function(tt, temp, alpha, beta0, beta1, gamma)
{

  #browser()
  x=1/(temp+273.16)

  mu=beta0+beta1*x
  nu=1/gamma
  zz=(log(tt)-mu)/nu
  res=alpha/(1+exp(zz))
  return(res)
}

################################################################################
#used for MLE
kinetics.fun1=function(tt, temp, alpha, beta0, beta1, gamma)
{

  #browser()
  x=11605/(temp+273.16)
  mx=min(x)
  mu=beta0+beta1*(x-mx)
  nu=1/gamma
  zz=(log(tt)-mu)/nu
  res=alpha/(1+exp(zz))
  return(res)
}



######################
# Semi-parametric model
######################



############################################################################
monotone.bspline.fit.knotselection.output=function(test.dat, n.knots.vec=1:5, degree=3, plot.loglike.fun=TRUE, use.aic=FALSE, use.aicc=TRUE, Boundary.knots=NULL, deg.vec=2:4, cov.fun.type){
  #deg.vec=2:4

  #browser()

  aicc.vec=lld.vec=edf.vec=NULL
  knot.list=NULL
  nknots.vec=NULL
  knotseq.vec=NULL

  for(m in 1:length(deg.vec)){

    tmp=monotone.bspline.fit.knotselection(test.dat=test.dat, degree=deg.vec[m], Boundary.knots=Boundary.knots, n.knots.vec=n.knots.vec, plot.loglike.fun=plot.loglike.fun, use.aic=use.aic, use.aicc=use.aicc, cov.fun.type=cov.fun.type)


    aicc.vec=c(aicc.vec, tmp$aicc)

    knot.list=c(knot.list, list(knots=tmp$knots))

  }

  for(k in 1:length(knot.list)){
    nknots.vec[k]=length(knot.list[[k]])
    knotseq.vec[k]=paste(round(knot.list[[k]], 2), collapse=", ")

    tmp=bspline.fit.addt(test.dat=test.dat, n.knots=NULL, knots=knot.list[[k]], degree=deg.vec[k], plot.loglike.fun=T, Boundary.knots=Boundary.knots,
                         cov.fun.type=cov.fun.type)

    lld.vec[k]=tmp$lld
    edf.vec[k]=tmp$fit.monotone.bspline$edf
  }

  res=data.frame(2:4, lld.vec, edf.vec, aicc.vec, nknots.vec, knotseq.vec)
  names(res)=c("Order", "Loglik", "edf", "AICC", "NumKnots", "Knots")

  #print(xtable(res), include.rownames=FALSE)

  return(invisible(res[which.min(res$AICC),]))
}

############################################################################
monotone.bspline.fit.nocor.knotselection.output=function(test.dat, n.knots.vec=1:5, degree=3, plot.loglike.fun=TRUE, use.aic=FALSE, use.aicc=TRUE, Boundary.knots=NULL, deg.vec=2:4, cov.fun.type=NULL){
  #deg.vec=2:4

  #browser()

  aicc.vec=lld.vec=edf.vec=NULL
  knot.list=NULL
  nknots.vec=NULL
  knotseq.vec=NULL

  for(m in 1:length(deg.vec)){


    tmp=monotone.bspline.fit.nocor.knotselection(test.dat=test.dat, degree=deg.vec[m], Boundary.knots=Boundary.knots, n.knots.vec=n.knots.vec, plot.loglike.fun=plot.loglike.fun, use.aic=use.aic, use.aicc=use.aicc, cov.fun.type=cov.fun.type)


    aicc.vec=c(aicc.vec, tmp$aicc)

    knot.list=c(knot.list, list(knots=tmp$knots))

  }

  for(k in 1:length(knot.list)){
    nknots.vec[k]=length(knot.list[[k]])
    knotseq.vec[k]=paste(round(knot.list[[k]], 2), collapse=", ")

    tmp=bspline.fit.addt.nocor(test.dat=test.dat, n.knots=NULL, knots=knot.list[[k]], degree=deg.vec[k], plot.loglike.fun=T, Boundary.knots=Boundary.knots)

    lld.vec[k]=tmp$lld
    edf.vec[k]=tmp$fit.monotone.bspline$edf
  }

  res=data.frame(2:4, lld.vec, edf.vec, aicc.vec, nknots.vec, knotseq.vec)
  names(res)=c("Order", "Loglik", "edf", "AICC", "NumKnots", "Knots")

  #print(xtable(res), include.rownames=FALSE)

  return(invisible(res[which.min(res$AICC),]))
}



###################################################################################
### fit monotone b-spline to the data
monotone.bspline.fit.knotselection=function(test.dat, n.knots.vec=1:5, degree=3, plot.loglike.fun=TRUE, use.aic=FALSE, use.aicc=TRUE, Boundary.knots=NULL, cov.fun.type){
  names(test.dat)=c("temp", "time", "Converted.Value")


  ## first step: find the number of knots that gives the minimum aic
  betahat.vec=aic=aicc=lld=rep(0, length(n.knots.vec))

  tmp.cat=function(syb, nsyb)
  {
    for(cat.count in 1:nsyb)
    {
      cat(syb)
    }
  }
  cat("|")
  tmp.cat(">", sum(n.knots.vec))
  cat("|\n|")

  for(nk in 1:length(n.knots.vec))
  {
    tmp.cat(">", nk)

    bspline.fit.addt.obj=bspline.fit.addt(test.dat=test.dat, n.knots=n.knots.vec[nk], knots=NULL, degree=degree, plot.loglike.fun=plot.loglike.fun, Boundary.knots=Boundary.knots, cov.fun.type=cov.fun.type)
    betahat.vec[nk]=bspline.fit.addt.obj$betahat
    lld[nk]=bspline.fit.addt.obj$lld
    aic[nk]=bspline.fit.addt.obj$aic
    aicc[nk]=bspline.fit.addt.obj$aicc

  }
  cat("|\n")

  #browser()

  #########################################################
  # use aic
  if(use.aic){
    nknotsbest=n.knots.vec[which.min(aic)]

    xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=betahat.vec[which.min(aic)], n.knots=nknotsbest, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)

    knots.nknotsbest=xt.mat$knots
    ## stepwise knot deletion
    aic.dele=rep(0, length(knots.nknotsbest))

    for(j in 1:length(knots.nknotsbest)){
      bspline.fit.addt.obj=bspline.fit.addt(test.dat=test.dat, n.knots=NULL, knots=knots.nknotsbest[-j], degree=degree, plot.loglike.fun=plot.loglike.fun, Boundary.knots=Boundary.knots, cov.fun.type=cov.fun.type)

      aic.dele[j]=bspline.fit.addt.obj$aic
    }

    # stepwise knot addition
  }


  #########################################################
  # use aicc
  if(use.aicc){
    nknotsbest=n.knots.vec[which.min(aicc)]

    xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=betahat.vec[which.min(aicc)], n.knots=nknotsbest, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)

    aicc.cur=min(aicc)

    aicc.diff=1

    knots.nknotsbest=xt.mat$knots

    # if number of knots equal 1, don't perform knot deletion
    if(nknotsbest>1){
      ## stepwise knot deletion
      while(aicc.diff>0){
        aicc.dele=rep(0, length(knots.nknotsbest))

        for(j in 1:length(knots.nknotsbest)){

          bspline.fit.addt.obj=try(bspline.fit.addt(test.dat=test.dat, n.knots=NULL, knots=knots.nknotsbest[-j], degree=degree, plot.loglike.fun=plot.loglike.fun, Boundary.knots=Boundary.knots, cov.fun.type=cov.fun.type), silent=T)

          if(!is.null(bspline.fit.addt.obj)){
            aicc.dele[j]=bspline.fit.addt.obj$aicc
          }




        }

        aicc.diff=aicc.cur-min(aicc.dele)

        if(aicc.diff>0){
          knots.nknotsbest=knots.nknotsbest[-which.min(aicc.dele)]
          aicc.cur=min(aicc.dele)
        }
      }

      # stepwise knot addition
    }


  }

  return(list(knots=knots.nknotsbest, aicc=aicc.cur))
}


###################################################################################
### fit monotone b-spline to the data, no correlation
monotone.bspline.fit.nocor.knotselection=function(test.dat, n.knots.vec=1:5, degree=3, plot.loglike.fun=TRUE, use.aic=FALSE, use.aicc=TRUE, Boundary.knots=NULL, cov.fun.type){
  names(test.dat)=c("temp", "time", "Converted.Value")


  ## first step: find the number of knots that gives the minimum aic
  betahat.vec=aic=aicc=lld=rep(0, length(n.knots.vec))

  tmp.cat=function(syb, nsyb)
  {
    for(cat.count in 1:nsyb)
    {
      cat(syb)
    }
  }
  cat("|")
  tmp.cat(">", sum(n.knots.vec))
  cat("|\n|")

  for(nk in 1:length(n.knots.vec))
  {
    tmp.cat(">", nk)

    bspline.fit.addt.obj=bspline.fit.addt.nocor(test.dat=test.dat, n.knots=n.knots.vec[nk], knots=NULL, degree=degree, plot.loglike.fun=plot.loglike.fun, Boundary.knots=Boundary.knots)
    betahat.vec[nk]=bspline.fit.addt.obj$betahat
    lld[nk]=bspline.fit.addt.obj$lld
    #aic[nk]=bspline.fit.addt.obj$aic
    aicc[nk]=bspline.fit.addt.obj$aicc

  }

  cat("|\n")
  #browser()


  #########################################################
  # use aicc
  if(use.aicc){
    nknotsbest=n.knots.vec[which.min(aicc)]

    xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=betahat.vec[which.min(aicc)], n.knots=nknotsbest, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)

    aicc.cur=min(aicc)

    aicc.diff=1

    knots.nknotsbest=xt.mat$knots

    # if number of knots equal 1, don't perform knot deletion
    if(nknotsbest>1){
      ## stepwise knot deletion
      while(aicc.diff>0){
        aicc.dele=rep(0, length(knots.nknotsbest))

        for(j in 1:length(knots.nknotsbest)){

          bspline.fit.addt.obj=bspline.fit.addt.nocor(test.dat=test.dat, n.knots=NULL, knots=knots.nknotsbest[-j], degree=degree, plot.loglike.fun=plot.loglike.fun, Boundary.knots=Boundary.knots)

          aicc.dele[j]=bspline.fit.addt.obj$aicc
        }

        aicc.diff=aicc.cur-min(aicc.dele)

        if(aicc.diff>0){
          knots.nknotsbest=knots.nknotsbest[-which.min(aicc.dele)]
          aicc.cur=min(aicc.dele)
        }
      }

      # stepwise knot addition
    }


  }

  return(list(knots=knots.nknotsbest, aicc=aicc.cur))
}


##############################################################################
bspline.fit.addt=function(test.dat, n.knots, knots, degree, plot.loglike.fun=T, Boundary.knots=NULL, rho.method="other", cov.fun.type){
  names(test.dat)=c("temp", "time", "Converted.Value")
  #browser()

  aa=seq(0.01, 1.93,len=50)
  #aa=seq(0.01, 3, len=100)
  #aa=seq(0.2, 1, len=60)
  #aa=seq(0.01, 0.5, len=50)
  raa=raa.lme=rep(NA, length(aa))

  if(is.null(n.knots)){ncoef=length(knots)+degree+1}else{ncoef=n.knots+degree+1}
  coef.mat.lme=matrix(NA, length(aa), ncoef+2)
  for(i in 1:length(aa))
  {
    #print(i)
    xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=aa[i], n.knots=n.knots, knots=knots, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)
    fit.dat.obj=xmat.obj.to.xmat(dat=test.dat,xmat.obj=xt.mat, intercept=FALSE)

    y=fit.dat.obj$dat[,"Y"]
    X=as.matrix(fit.dat.obj$dat[,3:(xt.mat$dfs+2)])
    ID=fit.dat.obj$dat[,"Batch"]

    #browser()
    # treat the results from B-spline as initial values
    #fit.gls <- try(gls(y~X-1, correlation=corCompSymm(form = ~ 1 | ID)), silent=T)

    # treat the results from B-spline as initial values
    fit.lme <- try( lme(y~X-1,  random = ~ 1|ID), silent=T)

    # treat the results from B-spline as initial values
    lme.flag=(attr(fit.lme,"class")=="try-error")

    pp=NULL

    if(!lme.flag){pp=as.vector(attr(fit.lme$apVar,"Pars"))}

    if(is.null(pp))
    {
      #res=list(conv=F)
      #return(res)
      ss=0.019
      rho=0.2
    }else{
      ts1=exp(pp[1])   #reStruct.Batch
      ts0=exp(pp[2])   #lSigma

      ss=sqrt(ts0^2+ts1^2)
      rho=ts1^2/(ts0^2+ts1^2)

    }

    if(!lme.flag){

      raa.lme[i]=loglik.compute(mm=length(unique(ID)), dat=fit.dat.obj$dat, ids=unique(ID), ss=ss, rho=rho,
                                yhat=as.numeric(X%*%fixef(fit.lme)))$loglik

      coef.mat.lme[i,]=c(fixef(fit.lme),ss,rho)


      #if(class(fit.gls)!="try-error"){
      #ss=fit.gls$sigma
      #rho=coef(summary(fit.gls)$modelStruct$corStruct, unconstrained=FALSE)

      #rho=ifelse(rho<0, 0.1, rho)

      #browser()

      #if(i==22 &cov.fun.type=="cov.fun.REMLc"){browser()}
      fit.monotone.bspline=try(bspline.lme.addt.cone(dat.obj=fit.dat.obj,beta.vec0=fixef(fit.lme),theta0=c(ss, rho), cov.fun.type=cov.fun.type, rho.method=rho.method), silent=T)


      #print(fit.monotone.bspline$coef["sigma"])
      #print(fit.monotone.bspline$coef["rho"])

      if(all(class(fit.monotone.bspline)!="try-error")){
        if(fit.monotone.bspline$conv==TRUE){
          if(fit.monotone.bspline$iter>=20){
            #print("greater than 20")
            fit.monotone.bspline=try(bspline.lme.addt.cone(dat.obj=fit.dat.obj,beta.vec0=fixef(fit.lme),theta0=c(ss, rho), cov.fun.type=cov.fun.type,rho.method=rho.method, adjust.stepsize=T), silent=T)
            #save(test.dat, file="test.dat")
          }
          raa[i]=fit.monotone.bspline$loglik
        }
      }

    }

    #print(i)
  }



  #max.id=which.max(raa)

  #print(raa)

  #browser()
  #pineapple
  #if(plot.loglike.fun){
  #  if(is.null(n.knots)){
  #    plot(aa,raa,type="l", xlab=expression(beta), ylab="log likelihood", main=paste("monotone b-spline, knots:",paste(round(knots, 2),  collapse = ",")))
  #    abline(v=aa[which.max(raa)], lty=2, col=2)
  #  } else{
  #    plot(aa,raa,type="l", xlab=expression(beta), ylab="log likelihood", main=paste("monotone b-spline, number of interior knots:", n.knots))
  #    abline(v=aa[which.max(raa)], lty=2, col=2)
  #  }
  #}



  # check whether it is a global maximum point for beta
  #if((all(diff(raa)[1:(max.id-1)]>0)) & (all(diff(raa)[max.id:(length(raa)-1)]<0))){
  #    betahat=aa[which.max(raa)]
  #} else{
  #  stop("cannot find the betahat")
  #}

  if(all(is.na(raa))){
    return(list(conv=F))
  }
  betahat=aa[which.max(raa)]

  coef.lme=c(aa[which.max(raa.lme)], coef.mat.lme[which.max(raa.lme),])

  xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=betahat, n.knots=n.knots, knots=knots, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)
  fit.dat.obj=xmat.obj.to.xmat(dat=test.dat,xmat.obj=xt.mat, intercept=FALSE)

  y=fit.dat.obj$dat[,"Y"]
  X=as.matrix(fit.dat.obj$dat[,3:(xt.mat$dfs+2)])
  ID=fit.dat.obj$dat[,"Batch"]

  fit.lme <- try( lme(y~X-1,  random = ~ 1|ID), silent=T)

  # treat the results from B-spline as initial values
  lme.flag=(attr(fit.lme,"class")=="try-error")

  if(!lme.flag){pp=as.vector(attr(fit.lme$apVar,"Pars"))}

  if(is.null(pp))
  {
    #res=list(conv=F)
    #return(res)
    ss=0.019
    rho=0.2
  }else{
    ts1=exp(pp[1])   #reStruct.Batch
    ts0=exp(pp[2])   #lSigma

    ss=sqrt(ts0^2+ts1^2)
    rho=ts1^2/(ts0^2+ts1^2)
  }
  # treat the results from B-spline as initial values
  #fit.gls <- try( gls(y~X-1,  correlation=corCompSymm(form = ~ 1 | ID)), silent=T)
  #if(class(fit.gls)!="try-error"){
  #  ss=fit.gls$sigma
  # rho=coef(summary(fit.gls)$modelStruct$corStruct, unconstrained=FALSE)

  if(!lme.flag){

    fit.monotone.bspline=try(bspline.lme.addt.cone(dat.obj=fit.dat.obj,beta.vec0=fixef(fit.lme),theta0=c(ss, rho), cov.fun.type=cov.fun.type, rho.method=rho.method), silent=T)


    if(all(class(fit.monotone.bspline)!="try-error")){
      if(fit.monotone.bspline$conv==TRUE){
        if(fit.monotone.bspline$iter>=20){
          #print("greater than 20")
          fit.monotone.bspline=try(bspline.lme.addt.cone(dat.obj=fit.dat.obj,beta.vec0=fixef(fit.lme),theta0=c(ss, rho), cov.fun.type=cov.fun.type,rho.method=rho.method, adjust.stepsize=T), silent=T)
          #save(test.dat, file="test.dat")
        }
      }
    }

    lld=fit.monotone.bspline$loglik

    n=dim(test.dat)[1]

    tmp=as.numeric(fit.monotone.bspline$dat$dat[,"Y"]-fit.monotone.bspline$yhat)

    aic=log(as.numeric(t(tmp)%*%solve(fit.monotone.bspline$Sigma)%*%tmp))+2*(2+fit.monotone.bspline$edf)/n

    aicc=-2*lld+2*(3+fit.monotone.bspline$edf)

    spline.inf=list(Boundary.knots=xt.mat$Boundary.knots, n.knots=n.knots, knots=knots, degree=degree)

    return(list(betahat=betahat, lld=lld, aic=aic, aicc=aicc, fit.monotone.bspline=fit.monotone.bspline, spline.inf=spline.inf, coef.lme=coef.lme, conv=T))


  }
}

##############################################################################
bspline.fit.addt.nocor=function(test.dat, n.knots, knots, degree, plot.loglike.fun=T, Boundary.knots=NULL,  cov.fun.type, aa=seq(0.01, 1.93,len=50)){
  names(test.dat)=c("temp", "time", "Converted.Value")
  #browser()

  aa=seq(0.01, 1.93,len=100)
  #aa=seq(0.01, 3, len=100)
  #aa=seq(0.2, 1, len=60)

  #aa=seq(0.01, 0.5, len=50)
  raa=raa.lm=rep(NA, length(aa))

  if(is.null(n.knots)){ncoef=length(knots)+degree+1}else{ncoef=n.knots+degree+1}
  coef.mat.lm=matrix(NA, length(aa), ncoef+1)
  for(i in 1:length(aa))
  {
    #print(i)
    xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=aa[i], n.knots=n.knots, knots=knots, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)
    fit.dat.obj=xmat.obj.to.xmat(dat=test.dat,xmat.obj=xt.mat, intercept=FALSE)

    y=fit.dat.obj$dat[,"Y"]
    X=as.matrix(fit.dat.obj$dat[,3:(xt.mat$dfs+2)])
    ID=fit.dat.obj$dat[,"Batch"]

    #browser()
    # treat the results from B-spline as initial values
    #fit.gls <- try(gls(y~X-1, correlation=corCompSymm(form = ~ 1 | ID)), silent=T)

    # treat the results from B-spline as initial values
    fit.lm <- try( lm(y~X-1), silent=T)

    # treat the results from B-spline as initial values
    lm.flag=(attr(fit.lm,"class")=="try-error")


    if(lm.flag)
    {
      #res=list(conv=F)
      #return(res)
    }else{

      raa.lm[i]=logLik(fit.lm)

      coef.mat.lm[i,]=c(coef(fit.lm), summary(fit.lm)$sigma)

      #if(i==22 &cov.fun.type=="cov.fun.REMLc"){browser()}
      fit.monotone.bspline=try(bspline.lme.addt.cone.nocor(dat.obj=fit.dat.obj,beta.vec0=coef(fit.lm),theta0=summary(fit.lm)$sigma), silent=T)

      #print(fit.monotone.bspline$coef["sigma"])
      #print(fit.monotone.bspline$coef["rho"])

      if(class(fit.monotone.bspline)!="try-error"){
        raa[i]=fit.monotone.bspline$loglik
      }



    }

    #print(i)
  }



  #max.id=which.max(raa)

  #print(raa)

  #browser()
  #pineapple
  #if(plot.loglike.fun){
  #  if(is.null(n.knots)){
  #    plot(aa,raa,type="l", xlab=expression(beta), ylab="log likelihood", main=paste("monotone b-spline, knots:",paste(round(knots, 2),  collapse = ",")))
  #    abline(v=aa[which.max(raa)], lty=2, col=2)
  #  } else{
  #    plot(aa,raa,type="l", xlab=expression(beta), ylab="log likelihood", main=paste("monotone b-spline, number of interior knots:", n.knots))
  #    abline(v=aa[which.max(raa)], lty=2, col=2)
  #  }
  #}



  # check whether it is a global maximum point for beta
  #if((all(diff(raa)[1:(max.id-1)]>0)) & (all(diff(raa)[max.id:(length(raa)-1)]<0))){
  #    betahat=aa[which.max(raa)]
  #} else{
  #  stop("cannot find the betahat")
  #}

  if(all(is.na(raa))){
    return(list(conv=F))
  }
  betahat=aa[which.max(raa)]

  coef.lm=c(aa[which.max(raa.lm)], coef.mat.lm[which.max(raa.lm),])

  xt.mat=addt.bsplines.xmat.obj.fun(dat=test.dat, beta=betahat, n.knots=n.knots, knots=knots, degree=degree, plot.splines=F, Boundary.knots=Boundary.knots)
  fit.dat.obj=xmat.obj.to.xmat(dat=test.dat,xmat.obj=xt.mat, intercept=FALSE)

  y=fit.dat.obj$dat[,"Y"]
  X=as.matrix(fit.dat.obj$dat[,3:(xt.mat$dfs+2)])
  ID=fit.dat.obj$dat[,"Batch"]

  # treat the results from B-spline as initial values
  fit.lm <- try( lm(y~X-1), silent=T)

  # treat the results from B-spline as initial values
  lm.flag=(attr(fit.lm,"class")=="try-error")


  if(lm.flag)
  {
    return(NULL)
  }else{
    #if(i==22 &cov.fun.type=="cov.fun.REMLc"){browser()}
    fit.monotone.bspline=try(bspline.lme.addt.cone.nocor(dat.obj=fit.dat.obj,beta.vec0=coef(fit.lm),theta0=summary(fit.lm)$sigma), silent=T)

    if(class(fit.monotone.bspline)!="try-error"){
      lld=fit.monotone.bspline$loglik

      tmp=as.numeric(fit.monotone.bspline$dat$dat[,"Y"]-fit.monotone.bspline$yhat)

      #aic=log(as.numeric(t(tmp)%*%solve(fit.monotone.bspline$Sigma)%*%tmp))+2*(2+fit.monotone.bspline$edf)/n

      aicc=-2*lld+2*(2+fit.monotone.bspline$edf)

      spline.inf=list(Boundary.knots=xt.mat$Boundary.knots, n.knots=n.knots, knots=knots, degree=degree)

      return(list(betahat=betahat, lld=lld,  aicc=aicc, fit.monotone.bspline=fit.monotone.bspline, spline.inf=spline.inf, coef.lm=coef.lm, conv=T))

    }

  }
}
##############################################################################

################################################################################
bspline.lme.addt.cone=function(dat.obj,theta0,beta.vec0, control.beta=0.001, control.theta=0.001, monotone=T,
                               cov.fun.type="lme", rho.method="other", adjust.stepsize=F)
{
  dat=dat.obj$dat
  lam=dat.obj$lam
  ncoef=length(lam)
  nn=dim(dat)[1]
  ids=unique(dat[,1])
  mm=length(ids)

  ss=theta0[1]
  rho=theta0[2]

  ss.vec0=c(ss,rho)
  beta.vec=beta.vec0

  #print(beta.vec0)

  ys=rep(0,nn)
  xmats=matrix(0,nrow=nn,ncol=ncoef)
  mat.inv=matrix(0,nrow=nn,ncol=nn)

  X=as.matrix(dat[, 3:(2+ncoef)])
  y=dat[, "Y"]

  convbeta=convtheta=1

  iter=1

  #compute log likelihood
  llk.obj0=loglik.compute(mm=mm, dat=dat, ids=ids, ss=ss, rho=rho, yhat=X%*%beta.vec0)

  while((convbeta>control.beta | convtheta>control.theta) & iter<20)
  {
    iter=iter+1

    ## construct covariance matrix
    Sigma=llk.obj0$Sigma

    #browser()
    Sigma.inv=solve(Sigma)

    #browser()

    # treat sigma matrix as known, define transformed X and y
    Ltilde=as.matrix(t(chol(Sigma.inv)))
    Ltildeinv=solve(Ltilde)
    ytilde=Ltildeinv%*%y
    Xtilde=Ltildeinv%*%X


    # use cone projection to solve the quadratic object function with constrains
    umat=chol(t(Xtilde)%*%Xtilde)
    uinv=solve(umat)

    A=cbind(diag(rep(1, ncoef-1)), rep(0, ncoef-1))+cbind(rep(0, ncoef-1), diag(rep(-1, ncoef-1)))
    bmata=A%*%uinv
    bmat=matrix(0, ncoef, ncoef)

    uvec=runif(ncoef)
    bmat[1,]=uvec-t(bmata)%*%solve(bmata%*%t(bmata))%*%bmata%*%uvec
    bmat[2:ncoef,]=bmata

    edges=t(solve(bmat))

    ysend=t(uinv)%*%t(Xtilde)%*%ytilde

    # use package coneproj
    coef=coneB(ysend,edges[2:ncoef,],matrix(t(edges[1,]),ncol=1))

    beta.vec=uinv%*%t(edges)%*%coef$coefs

    # find the effective degree of freedom
    sm=1e-8

    index=coef$coefs>sm
    index[1]=TRUE
    gmat=edges[index,]

    if(length(gmat)/ncoef==1){gmat=matrix(gmat,ncol=ncoef)}

    pcmat=Xtilde%*%uinv%*%t(gmat)%*%solve(gmat%*%t(gmat))%*%gmat%*%t(uinv)%*%t(Xtilde)

    edfc=sum(diag(pcmat))

    #
    yhat=X%*%beta.vec
    ww1=y-yhat

    #browser()


    # estimate covariance parameters
    if(cov.fun.type=="lme"){
      tdat=data.frame(Batch=dat[,c("Batch")],ww=ww1)
      tfit=try(lme(ww~-1,data=tdat,random=~1|Batch),silent=T)

      #browser()
      #tfit=try(gls(ww~1,data=tdat,correlation=corCompSymm(form = ~ 1 | Batch)),silent=T)

      lme.flag=(attr(tfit,"class")=="try-error")

      if(lme.flag)
      {
        res=list(conv=F)
        return(res)
      }


      pp=as.vector(attr(tfit$apVar,"Pars"))

      if(is.null(pp))
      {
        pp=c(-12, log(tfit$sigma))

        #res=list(conv=F)
        #return(res)
      }

      ts1=exp(pp[1])   #reStruct.Batch
      ts0=exp(pp[2])   #lSigma

      ss=sqrt(ts0^2+ts1^2)
      rho=ts1^2/(ts0^2+ts1^2)
      ss.vec=c(ss,rho)

      #fitted=yhat+tfit$coef$fixed+as.vector(fitted(tfit))
      fitted=yhat+as.vector(fitted(tfit))

      llk.obj=loglik.compute(mm=mm, dat=dat, ids=ids, ss=ss, rho=rho, yhat=yhat)
      #fitted=yhat
    }

    if(cov.fun.type=="cov.fun.gls"){
      tdat=data.frame(Batch=dat[,c("Batch")],ww=ww1)
      fit.gls=try(gls(ww~1,data=tdat,correlation=corCompSymm(form = ~ 1 | Batch)),silent=T)

      ss=fit.gls$sigma
      rho=coef(summary(fit.gls)$modelStruct$corStruct, unconstrained=FALSE)

      ss.vec=c(ss, rho)
      llk.obj=loglik.compute(mm=mm, dat=dat, ids=ids, ss=ss, rho=rho, yhat=yhat)
    }

    if(cov.fun.type=="cov.fun.REML"){

      #browser()

      ID=dat.obj$dat[,"Batch"]

      rho=glsnew(y~X-1, correlation=corCompSymm(form = ~ 1 | ID), coef.ini=as.numeric(round(beta.vec, 8)),
                 sigma.ini=ss)

      #print(rho)

      if(rho<0){
        #browser()
        rho=exp(try(nlminb(start=log(0.5), control=list(iter.max=20), objective=cov.fun.REML.rho, mm=mm, dat=dat, ids=ids, X_L=X_L, ww1=ww1, ss=ss)$par, silent=T))
      }

      for(j in 1:mm)
      {
        ww=sum((dat[,1]==ids[j])*1)
        RR=diag(ww)+rho*as.numeric(ww>1)-rho*diag(ww)*as.numeric(ww>1)

        tmp.ww1=as.numeric(ww1)[dat[,1]==ids[j]]

        if(j==1){qudsum=t(tmp.ww1)%*%solve(RR)%*%tmp.ww1} else{
          qudsum=t(tmp.ww1)%*%solve(RR)%*%tmp.ww1+qudsum
        }

      }

      ss=sqrt(as.numeric(qudsum)/(nrow(X)-ncol(X)))

      llk.obj=loglik.compute(mm=mm, dat=dat, ids=ids, ss=ss, rho=rho, yhat=yhat)

      ss.vec=c(ss, rho)
    }

    if(cov.fun.type=="cov.fun.REMLc"){
      #browser()
      X_L=NULL
      for(ll in unique(round(beta.vec, 8))){
        tmp=X[,round(beta.vec, 8)==ll]

        if(sum(round(beta.vec, 8)==ll)!=1){tmp=apply(tmp, 1, sum)}

        X_L=cbind(X_L, tmp)
      }

      #browser()
      colnames(X_L)=NULL

      ID=dat.obj$dat[,"Batch"]

      if(rho.method=="glsnew"){
        rho=glsnew(y~X_L-1, correlation=corCompSymm(form = ~ 1 | ID), coef.ini=as.numeric(unique(round(beta.vec, 8))),
                   sigma.ini=ss)

        #print(rho)

        if(rho<0){
          #browser()
          #if(iter==4){browser()}
          rho=exp(try(nlminb(start=qlogis(0.2), control=list(iter.max=20), objective=cov.fun.REML.rho, mm=mm, dat=dat, ids=ids, X_L=X_L, ww1=ww1, ss=ss)$par, silent=T))

          #print(rho)
        }
      }else{
        rho=plogis(try(nlminb(start=qlogis(0.5), control=list(iter.max=20), objective=cov.fun.REML.rho, mm=mm, dat=dat, ids=ids, X_L=X_L, ww1=ww1, ss=ss)$par, silent=T))
      }

      #browser()

      # compute sigma
      for(j in 1:mm)
      {
        ww=sum((dat[,1]==ids[j])*1)
        RR=diag(ww)+rho*as.numeric(ww>1)-rho*diag(ww)*as.numeric(ww>1)

        tmp.ww1=as.numeric(ww1)[dat[,1]==ids[j]]

        if(j==1){qudsum=t(tmp.ww1)%*%solve(RR)%*%tmp.ww1} else{
          qudsum=t(tmp.ww1)%*%solve(RR)%*%tmp.ww1+qudsum
        }

      }
      #browser()

      #print(rho)

      ss=sqrt(as.numeric(qudsum)/(nrow(X_L)-ncol(X_L)))

      #browser()
      llk.obj=loglik.compute(mm=mm, dat=dat, ids=ids, ss=ss, rho=rho, yhat=yhat)

      #print(llk.obj$loglik)

      if(adjust.stepsize==T){
        iter2=1
        while((llk.obj$loglik<llk.obj0$loglik) & iter2<50){
          iter2=iter2+1
          ss=(ss+ss.vec0[1])/2
          rho=(rho+ss.vec0[2])/2
          llk.obj=loglik.compute(mm=mm, dat=dat, ids=ids, ss=ss, rho=rho, yhat=yhat)

        }
      }



      ss.vec=c(ss, rho)
      #print(ss.vec)

      #ss.vec=REML(pars=c(ss.vec_L[1], ss.vec_L[2]),  mm=mm, dat=dat, ids=ids, X_L=X_L, ww1=ww1)$par


    }

    convbeta=max(abs(beta.vec0-beta.vec))
    convtheta=max(abs(ss.vec0-ss.vec))

    beta.vec0=beta.vec
    ss.vec0=ss.vec
    llk.obj0=llk.obj


    #print(round(beta.vec0, 4))
    #print(round(ss.vec0, 4))

    #print(paste("diff in beta", formatC(convbeta, digits=4, format="fg", flag="#")))
    #print(paste("diff in theta", formatC(convtheta, digits=4, format="fg", flag="#")))




  }

  if(cov.fun.type=="cov.fun.gls"|cov.fun.type=="cov.fun.REML"|cov.fun.type=="cov.fun.REMLc"){
    fitted=NULL
    error=NULL
    ts0=sqrt(ss.vec[1]^2*(1-ss.vec[2]))
    ts1=sqrt(ss.vec[1]^2*(ss.vec[2]))

    Sigma=llk.obj$Sigma
    loglik=llk.obj$loglik
  }else{
    error=dat[,"Y"]-fitted
    loglik=llk.obj$loglik
  }


  #std.error=error/ss

  coef=c(beta.vec,ss.vec, ts0, ts1)
  names(coef)=c(colnames(dat)[3:(ncoef+2)],"sigma","rho", "sigma_e", "sigma_u")

  #coef=c(beta.vec,ss.vec)
  #names(coef)=c(colnames(dat)[3:(ncoef+2)],"sigma","rho")


  res=list(dat=dat.obj,fitted=fitted,yhat=yhat,coef=coef,
           error=error,loglik=loglik,edf=edfc,Sigma=Sigma,conv=TRUE, iter=iter)
  return(res)
}


################################################################################
bspline.lme.addt.cone.nocor=function(dat.obj,theta0,beta.vec0, control.beta=0.001, control.theta=0.001, monotone=T)
{
  dat=dat.obj$dat
  lam=dat.obj$lam
  ncoef=length(lam)
  nn=dim(dat)[1]
  ids=unique(dat[,1])
  mm=length(ids)

  beta.vec=beta.vec0
  sigma=theta=theta0

  #print(beta.vec0)

  ys=rep(0,nn)
  xmats=matrix(0,nrow=nn,ncol=ncoef)
  mat.inv=matrix(0,nrow=nn,ncol=nn)

  X=as.matrix(dat[, 3:(2+ncoef)])
  y=dat[, "Y"]

  n=length(y)

  convbeta=convtheta=1

  iter=1

  #compute log likelihood

  while((convbeta>control.beta | convtheta>control.theta) & iter<20)
  {
    iter=iter+1

    ## construct covariance matrix
    Sigma=sigma^2*diag(rep(1, n))

    #browser()
    Sigma.inv=sigma^(-2)*diag(rep(1, n))

    #browser()

    # treat sigma matrix as known, define transformed X and y
    Ltilde=as.matrix(t(chol(Sigma.inv)))
    Ltildeinv=solve(Ltilde)
    ytilde=Ltildeinv%*%y
    Xtilde=Ltildeinv%*%X


    # use cone projection to solve the quadratic object function with constrains
    umat=chol(t(Xtilde)%*%Xtilde)
    uinv=solve(umat)

    A=cbind(diag(rep(1, ncoef-1)), rep(0, ncoef-1))+cbind(rep(0, ncoef-1), diag(rep(-1, ncoef-1)))
    bmata=A%*%uinv
    bmat=matrix(0, ncoef, ncoef)

    uvec=runif(ncoef)
    bmat[1,]=uvec-t(bmata)%*%solve(bmata%*%t(bmata))%*%bmata%*%uvec
    bmat[2:ncoef,]=bmata

    edges=t(solve(bmat))

    ysend=t(uinv)%*%t(Xtilde)%*%ytilde

    # use package coneproj
    coef=coneB(ysend,edges[2:ncoef,],matrix(t(edges[1,]),ncol=1))

    beta.vec=uinv%*%t(edges)%*%coef$coefs

    # find the effective degree of freedom
    sm=1e-8

    index=coef$coefs>sm
    index[1]=TRUE
    gmat=edges[index,]

    if(length(gmat)/ncoef==1){gmat=matrix(gmat,ncol=ncoef)}

    pcmat=Xtilde%*%uinv%*%t(gmat)%*%solve(gmat%*%t(gmat))%*%gmat%*%t(uinv)%*%t(Xtilde)

    edfc=sum(diag(pcmat))

    #
    yhat=X%*%beta.vec
    ww1=y-yhat

    #browser()

    theta=sigma=sqrt(sum(ww1^2)/(n-edfc))


    convbeta=max(abs(beta.vec0-beta.vec))
    convtheta=max(abs(theta-theta0))

    beta.vec0=beta.vec
    theta0=theta

    #print(paste("diff in beta", formatC(convbeta, digits=4, format="fg", flag="#")))
    #print(paste("diff in theta", formatC(convtheta, digits=4, format="fg", flag="#")))

  }


  loglik=sum(dnorm(ww1, mean=0, sd=sigma, log=TRUE))
  #std.error=error/ss

  coef=c(beta.vec, sigma)
  names(coef)=c(colnames(dat)[3:(ncoef+2)],"sigma")

  res=list(dat=dat.obj,fitted=fitted,yhat=yhat,coef=coef,
           loglik=loglik,edf=edfc, conv=TRUE, iter=iter)
  return(res)
}



################################################################################
addt.bsplines.xmat.obj.fun=function(dat,beta=0,n.knots=5,degree=3,plot.splines=F, knots=NULL, eq.alloc, Boundary.knots=NULL)
{
  time=dat[,"time"]
  temp=dat[,"temp"]
  max.temp=max(temp)
  dK=273.16

  temp.K.inv=11605/(temp+dK)-11605/(max.temp+dK)
  ss=exp(beta*temp.K.inv)
  time.ss=time/ss

  #print(time.ss)
  #browser()


  if(is.null(knots)){
    knots=quantile(time.ss, probs=(1:n.knots)/(n.knots+1))
  }
  #knots=(1:n.knots)*max(time.ss)/(n.knots+1)
  #knots=quantile(time.ss, probs=c(0.1, 1:(n.knots-1)/n.knots))


  if(is.null(Boundary.knots)){
    Boundary.knots=range(time.ss)
  }


  #browser()
  #spline.mat=bs(time.ss, knots=knots, degree = degree, Boundary.knots=Boundary.knots, intercept=TRUE)

  spline.mat=newbs.matfun(time.ss, knots=knots, degree = degree, Boundary.knots=Boundary.knots)





  #spline.mat=ns(time.ss, df=n.knots+degree)

  #spline.mat=(-1)*MIC.splines.basis(time.ss, df = n.knots+degree, knots = NULL,
  #                                  boundary.knots=NULL,type="Is",degree = degree,delta=1,eq.alloc=T)$mat

  #spline.mat=bs(time.ss, df=n.knots+degree, degree = degree)


  if(plot.splines)
  {
    #fig.paper(file="uv.splines.basis")
    #matplot(time.ss[order(time.ss)],spline.mat[order(time.ss),],col=1,type="l",xlab="",ylab="",las=1,cex.axis=1.5)
    matplot(time.ss[order(time.ss)],spline.mat[order(time.ss),],col=1,type="l",xlab="",ylab="",las=1)
    #dev.off()
  }

  a1=dim(spline.mat)[2]

  dfs=c(a1)

  res=list(dat=dat,dfs=dfs, spline.mat=spline.mat, time.ss=time.ss, knots=knots, degree=degree, Boundary.knots=Boundary.knots)

  return(invisible(res))
}

################################################################################
newbs.matfun=function(time.ss, knots, degree, Boundary.knots){
  spline.mat=matrix(0, length(time.ss), degree+length(knots)+1)

  for(mm in 1:length(time.ss)){
    spline.mat[mm,]=newbs(x=time.ss[mm], degree=degree, inner.knots=knots, Boundary.knots=Boundary.knots)
  }

  #pineapple


  return(spline.mat)
}

################################################################################
xmat.obj.to.xmat=function(dat,xmat.obj, intercept=TRUE)
{
  dfs=xmat.obj$dfs
  ids=paste(dat[,"temp"],dat[,"time"],sep="")

  if(intercept==TRUE){
    xmat=data.frame(ids, dat[,"Converted.Value"],1,xmat.obj$spline.mat)

    #browser()
    colnames(xmat)=c("Batch", "Y", paste("Time",0:(dfs[1]),sep=""))

    dat=as.data.frame(xmat)

    lam=c(0,rep(1,dfs[1]))
  } else {
    xmat=data.frame(ids, dat[,"Converted.Value"], xmat.obj$spline.mat)

    #browser()
    colnames(xmat)=c("Batch", "Y", paste("Time", 1:(dfs[1]),sep=""))

    dat=as.data.frame(xmat)

    lam=c(0,rep(1,dfs[1]-1))
  }


  res=list(dat=dat,lam=lam)

  return(res)
}




##############################################################################
newbs=function(x, degree, inner.knots, Boundary.knots) {
  Boundary.knots=sort(Boundary.knots);
  knots=c(rep(Boundary.knots[1], (degree+1)), sort(inner.knots),
          rep(Boundary.knots[2], (degree+1)));
  np=degree+length(inner.knots)+1
  s=rep(0, np)
  if(x==Boundary.knots[2]) {s[np]=1} else {for( i in 1: np)
    s[i]=basis(x, degree, i, knots)}
  return(s)}

##############################################################################
basis=function(x, degree, i, knots)
{ if(degree==0){ if((x<knots[i+1])&(x>=knots[i])) y=1 else
  y=0}else{
    if((knots[degree+i]-knots[i])==0) {temp1=0} else {temp1=
      (x-knots[i])/(knots[degree+i]-knots[i])};
    if((knots[i+degree+1]-knots[i+1])==0) {temp2=0} else {temp2=
      (knots[i+degree+1]-x)/(knots[i+degree+1]-knots[i+1])}
    y= temp1*basis(x, (degree-1), i, knots) +temp2*basis(x, (degree-1),
                                                         (i+1), knots)}
  return(y)}



################################################################################
#compute log likelihood
loglik.compute=function(mm, dat, ids, ss, rho, yhat){
  loglik=0
  for(j in 1:mm)
  {
    id.x=(dat[,1]==ids[j])
    tmp1=dat[id.x, "Y"]-yhat[id.x]
    ww=length(tmp1)
    #browser()
    SS=diag(ww)*ss^2+rho*ss^2*as.numeric(ww>1)-rho*ss^2*diag(ww)*as.numeric(ww>1)
    SS.inv=solve(SS)
    ll=-(ww/2)*log(2*pi)-.5*log(det(SS))-.5*as.vector(t(tmp1)%*%SS.inv%*%tmp1)

    #print(ll)

    if(j==1){ Sigma=SS } else{ Sigma=bdiag(Sigma, SS)}

    loglik=loglik+ll
  }

  return(list(Sigma=Sigma, loglik=loglik))
}



################################################################################
cov.fun.REML=function(pars, mm, dat, ids, X, ww1){
  ss=exp(pars[1])
  rho=pars[2]

  for(j in 1:mm)
  {
    ww=sum((dat[,1]==ids[j])*1)
    SS=diag(ww)*ss^2+rho*ss^2*as.numeric(ww>1)-rho*ss^2*diag(ww)*as.numeric(ww>1)

    if(j==1){ Sigma=SS } else{ Sigma=bdiag(Sigma, SS)}
    #cat("j=",j,"OK","\n")
  }

  Sigma.inv=solve(Sigma)

  det.sigma=determinant(Sigma)
  det.XinvX=determinant(t(X)%*%Sigma.inv%*%X)

  return(as.numeric(det.sigma$modulus*det.sigma$sign+det.XinvX$modulus*det.XinvX$sign+as.numeric(t(ww1)%*%Sigma.inv%*%ww1)))
}

################################################################################
cov.fun.REMLc=function(pars, X_L, mm, dat, ids,  ww1){
  ss=exp(pars[1])
  rho=pars[2]

  for(j in 1:mm)
  {
    ww=sum((dat[,1]==ids[j])*1)
    SS=diag(ww)*ss^2+rho*ss^2*as.numeric(ww>1)-rho*ss^2*diag(ww)*as.numeric(ww>1)

    if(j==1){ Sigma=SS } else{ Sigma=bdiag(Sigma, SS)}
    #cat("j=",j,"OK","\n")
  }

  Sigma.inv=solve(Sigma)

  det.sigma=determinant(Sigma)
  det.XinvX=determinant(t(X_L)%*%Sigma.inv%*%X_L)

  return(as.numeric(det.sigma$modulus*det.sigma$sign+det.XinvX$modulus*det.XinvX$sign+as.numeric(t(ww1)%*%Sigma.inv%*%ww1)))
}

################################################################################
cov.fun.REML.rho=function(pars, X_L, mm, dat, ids,  ww1, ss){
  rho=plogis(pars)

  for(j in 1:mm)
  {
    ww=sum((dat[,1]==ids[j])*1)
    SS=diag(ww)*ss^2+rho*ss^2*as.numeric(ww>1)-rho*ss^2*diag(ww)*as.numeric(ww>1)

    if(j==1){ Sigma=SS } else{ Sigma=bdiag(Sigma, SS)}
    #cat("j=",j,"OK","\n")
  }

  Sigma.inv=solve(Sigma)

  det.sigma=determinant(Sigma)
  det.XinvX=determinant(t(X_L)%*%Sigma.inv%*%X_L)

  return(as.numeric(det.sigma$modulus*det.sigma$sign+det.XinvX$modulus*det.XinvX$sign+as.numeric(t(ww1)%*%Sigma.inv%*%ww1)))
}

############################################################################
glsnew=function (model, data = sys.frame(sys.parent()), correlation = NULL,
                 weights = NULL, subset, method = c("REML", "ML"), na.action = na.fail,
                 control = list(), verbose = FALSE, coef.ini, sigma.ini)
{
  Call <- match.call()
  controlvals <- glsControl()
  if (!missing(control)) {
    controlvals[names(control)] <- control
  }
  if (!inherits(model, "formula") || length(model) != 3L) {
    stop("\nmodel must be a formula of the form \"resp ~ pred\"")
  }
  method <- match.arg(method)
  REML <- method == "REML"
  if (!is.null(correlation)) {
    groups <- getGroupsFormula(correlation)
  }
  else groups <- NULL
  glsSt <- glsStruct(corStruct = correlation, varStruct = varFunc(weights))
  model <- terms(model, data = data)
  mfArgs <- list(formula = asOneFormula(formula(glsSt), model,
                                        groups), data = data, na.action = na.action)
  if (!missing(subset)) {
    mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2L]]
  }
  mfArgs$drop.unused.levels <- TRUE
  dataMod <- do.call("model.frame", mfArgs)
  origOrder <- row.names(dataMod)
  if (!is.null(groups)) {
    groups <- eval(parse(text = paste("~1", deparse(groups[[2L]]),
                                      sep = "|")))
    grps <- getGroups(dataMod, groups, level = length(getGroupsFormula(groups,
                                                                       asList = TRUE)))
    ord <- order(grps)
    grps <- grps[ord]
    dataMod <- dataMod[ord, , drop = FALSE]
    revOrder <- match(origOrder, row.names(dataMod))
  }
  else grps <- NULL
  X <- model.frame(model, dataMod)
  contr <- lapply(X, function(el) if (inherits(el, "factor"))
    contrasts(el))
  contr <- contr[!unlist(lapply(contr, is.null))]
  X <- model.matrix(model, X)
  if (ncol(X) == 0L)
    stop("no coefficients to fit")
  y <- eval(model[[2L]], dataMod)
  N <- nrow(X)
  p <- ncol(X)
  parAssign <- attr(X, "assign")
  fTerms <- terms(as.formula(model), data = data)
  namTerms <- attr(fTerms, "term.labels")
  if (attr(fTerms, "intercept") > 0) {
    namTerms <- c("(Intercept)", namTerms)
  }
  namTerms <- factor(parAssign, labels = namTerms)
  parAssign <- split(order(parAssign), namTerms)
  attr(glsSt, "conLin") <- list(Xy = array(c(X, y), c(N, ncol(X) +
                                                        1L), list(row.names(dataMod), c(colnames(X), deparse(model[[2]])))),
                                dims = list(N = N, p = p, REML = as.integer(REML)), logLik = 0)
  glsEstControl <- controlvals["singular.ok"]
  #browser()
  glsSt <- Initialize(glsSt, dataMod, glsEstControl)

  attr(glsSt, "glsFit")[["beta"]]=coef.ini
  attr(glsSt, "glsFit")[["sigma"]]=sigma.ini

  parMap <- attr(glsSt, "pmap")

  numIter <- numIter0 <- 0L  # was commented out
  if (length(coef(glsSt))) {
    optRes <- if (controlvals$opt == "nlminb") {
      nlminb(c(coef(glsSt)), function(glsPars) -logLik(glsSt,
                                                       glsPars), control = list(trace = controlvals$msVerbose,
                                                                                iter.max = controlvals$msMaxIter))
    }
    else {
      cat("numIter=", numIter, "\n")
      optim(c(coef(glsSt)), function(glsPars) -logLik(glsSt, glsPars), method = controlvals$optimMethod, control = list(trace = controlvals$msVerbose, maxit = controlvals$msMaxIter, reltol = if(numIter ==0L) controlvals$msTol else 100 * .Machine$double.eps))
    }
    coef(glsSt) <- optRes$par
  }
  else {
    optRes <- list(convergence = 0)
  }
  return(coef(glsSt, unconstrained=FALSE))
}


############################################################################
AdhesiveBondB.data.read=function()
{
  tmp=read.csv(file="AdhesiveBondB.csv",header=T)
  tmp=tmp[tmp[,4]=="Exact",]
  tmp=tmp[,-4]
  tmp=tmp[!(tmp[,1]%in%c(60,70) & tmp[,2]==0),]
  tmp=tmp[order(tmp[,1],tmp[,2]),]
  dat=data.frame(TempC=tmp[,1],TimeH=tmp[,2]*24*7,Response=tmp[,3])
  return(dat)
}



############################################################################
############################################################################
# TI
TI.bspline.nocor=function(dat,  model.fit.obj, dd=100000, failure.threshold = 0.5){
  # browser()

  #browser()
  names(dat)[1:3]=c("temp", "time", "response")
  time=dat[,"time"]
  temp=dat[,"temp"]
  max.temp=max(temp)
  dK=273.16

  xmax=1/(max.temp+dK)

  beta=model.fit.obj$betahat

  Boundary.knots=model.fit.obj$spline.inf$Boundary.knots
  knots=model.fit.obj$spline.inf$knots
  degree=model.fit.obj$spline.inf$degree
  ncoef=length(knots)+degree+1

  time.ss.vec=seq(min(Boundary.knots), max(Boundary.knots), len=200)
  #time.ss.vec=seq(0, max(Boundary.knots), len=200)
  spline.mat=newbs.matfun(time.ss.vec, knots=knots, degree = degree, Boundary.knots=Boundary.knots)

  yests=as.numeric(spline.mat%*%model.fit.obj$fit.monotone.bspline$coef[1:ncoef])

  tt.max=approx(yests, time.ss.vec, xout=failure.threshold*yests[1])$y

  TI=11605/(log(dd/tt.max)/beta+11605/(max.temp+dK))-dK

  beta0=log10(tt.max)-beta*(11605)*xmax/log(10)
  beta1=11605*beta/log(10)
  #return(c("TI"=TI, "beta1" = beta1 , "beta0" = beta0))
  TI.values <- c(c("TI"=TI, "beta1" = beta1 , "beta0" = beta0))

  return("TI"=TI.values)
}




