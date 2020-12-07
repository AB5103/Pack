ParNPHCox=function(formula.scale, formula.shape,cluster,dist,data,expr,strata) 
{
  #######This function calculates negloglikelihood for Weibull or Gompertz survival model
  ptm <- proc.time()
  nr=nrow(data)
  obsdata <- NULL
  if (length(formula.scale[[2]]) == 3) {
    obsdata$trunc <- c(rep(0,nrow(data)))
    obsdata$time <- eval(formula.scale[[2]][[2]], envir = data)
    obsdata$event <- eval(formula.scale[[2]][[3]], envir = data)
  }   else if (length(formula.scale[[2]]) == 4) {
    obsdata$trunc <- eval(formula.scale[[2]][[2]], envir = data)
    obsdata$time <- eval(formula.scale[[2]][[3]], envir = data)
    obsdata$event <- eval(formula.scale[[2]][[4]], envir = data)
  }
  if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
    stop(paste("The status indicator 'event' in the Surv object", 
               "in the left-hand side of the formula object", "must be either 0 (no event) or 1 (event)."))
  }
  obsdata$x <- as.data.frame(model.matrix.lm(formula.scale, data = data,na.action='na.pass'))
  obsdata$xs <- as.data.frame(model.matrix.lm(formula.shape, data = data,na.action='na.pass')) #factors for shape
  names(obsdata$x)=paste(names(obsdata$x),"scale",sep=".")
  names(obsdata$xs)=paste(names(obsdata$xs),"shape",sep=".")
  ind.x=which(is.na(c(apply(obsdata$x,1,sum))))
  ind.xs=which(is.na(c(apply(obsdata$xs,1,sum))))

    if (is.null(cluster)) {
     obsdata$ncl <- 0
     obsdata$di <- sum(obsdata$event)
     obsdata$cluster <- c(rep(0,nrow(data)))
     ind.cl <- as.numeric({})
    }   else {
    if (!cluster %in% names(data)) {
      stop(paste0("object '", cluster, "' not found"))
          }
    obsdata$cluster <- eval(as.name(cluster), envir = data)
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    obsdata$di <- aggregate(obsdata$event, by = list(obsdata$cluster), 
                            FUN = sum)[, , drop = FALSE]
    cnames <- obsdata$di[, 1]
    obsdata$di <- as.vector(obsdata$di[, 2])
    names(obsdata$di) <- cnames
    ind.cl=which(is.na(obsdata$cl))
    }
  ind=sort(unique(c(ind.x,ind.xs,ind.cl))) 
  nx=rep(0,nr)
  nx[ind]=1
  if (is.factor(obsdata$cluster)) obsdata$cluster=as.character(obsdata$cluster)
  obs=data.frame(obsdata$xs[-1],obsdata$x[-1],obsdata$event,obsdata$trunc,obsdata$time,obsdata$cluster)[nx!=1,]
  nr=nrow(obs)
  namesk=names(obsdata$xs)[-1] 
  namesf=names(obsdata$x)[-1] 
  names(obs)=c(names(obsdata$xs)[-1],names(obsdata$x)[-1],names(obsdata$event),"event","trunc","time","cluster")
  ncl=obsdata$ncl
  nk=length(namesk)
  nf=length(namesf)
  D=obs
  par0=c(0,0)
  Result=ucminf(par=par0,LikGenNPH,gr=NULL,D=obs,nf=0,nk=0,ncl=0,dist=dist,hessian=1)
  if (ncl>0){
    par0=c(Result$par,rep(0,(1+nk+nf)))} else {
      par0=c(Result$par,rep(0,(nk+nf)))
    }  
  Result=ucminf(par=par0,LikGenNPH,gr=GrGenNPH,D=obs,nf=nf,nk=nk,ncl=ncl,dist=dist,hessian=1)
  par=Result$par
  if (any(!is.finite(as.matrix(Result$hessian))))
    stop("infinite or missing values in hessian. It is not possible to calculate the matrix of covariates. \n  Change the model and try again.")
if (any(suppressWarnings(diag(ginv(Result$hessian)))<0))
  stop("hessian cannot be correctly calculated. \n  Change the model and try again.")
  invHes=sqrt(diag(ginv(Result$hessian)))
  Lik=-Result$value
  HCT=data.frame(obs$time,negH,obs$event)
  colnames(HCT)=c("Time","H","Cens")
  CSE=survConcordance(Surv(Time,Cens) ~ H,HCT)
  Conc=c(as.numeric(CSE$concordance),as.numeric(CSE$std.err))
  Vnames1=c("a","b",namesk,namesf)
Vnames={}
    if ((nk+nf)>0)    Vnames=paste("exp(",c(namesk,namesf),")",sep="")
    if (ncl>0) Vnames1=c(Vnames1,"Sigma2")
  
  if (ncl==0 & dist=="Weibull"){
    Names=c("Sample size","Number of non-censored","a","b",Vnames,"Concordance (se)","Loglik","AIC")} else if (ncl>0 & dist=="Weibull"){
      Names=c("Sample size","Number of non-censored","a  ","b",Vnames,"Sigma2","Concordance (se)","Loglik","AIC")} else if (ncl==0 & dist=="Gompertz"){ 
        Names=c("Sample size","Number of non-censored","1000a ","100b",Vnames,"Concordance (se)","Loglik","AIC")} else if (ncl>0 & dist=="Gompertz"){ 
          Names=c("Sample size","Number of non-censored","1000a ","100b",Vnames,"Sigma2","Concordance (se)","Loglik","AIC")  
        }
set.seed(123)
vpar<<- mvrnorm(n = 1000000, Result$par, ginv(Result$hessian), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
p_h={}
p_hl={}
p_hu={}
for (lika in 1:length(par)){
hh0=ecdf(vpar[,lika])  
if (lika<3 | (lika==length(par) & ncl>0)) {
  p_h=c(p_h,round(hh0(-Inf),4))} else {
  p_h=c(p_h,round(2*min(hh0(0),1-hh0(0)),4))
  }
p_hl=c(p_hl,round(exp(quantile(hh0,probs=0.025)),3))
p_hu=c(p_hu,round(exp(quantile(hh0,probs=0.975)),3))
}
parval=c(as.character(nr),as.character(sum(obs$event)),round(exp(Result$par),3),paste0(round(Conc[1],3)," (",round(Conc[2],3),")"),round(-Result$value,2),round(2*(Result$value+length(par0)),2))
    parvalm=c("","",p_hl,"","","")
    parvalp=c("","",p_hu,"","","")
    CI=c("","",paste(parvalm[3:(2+length(par))],parvalp[3:(2+length(par))],sep=" - "),"","","")
  pval=c("","",p_h,"","","")
  pval0=as.numeric(pval[3:(2+length(par0))])
  Tab=data.frame(parval,CI,pval)
  colnames(Tab)=c("Estimates","CI","p-value")
  rownames(Tab)=Names
  capt=paste("Parameter estimates.",dist,"model.",sep=" ")
  print(xtable(Tab,caption=capt))
  pcon=NULL 
    if (!is.null(expr)){
    if(is.expression(expr)){
      D1<<-obs      
      nf1<<-nf
      nk1<<-nk
      ncl1<<-ncl
      dist1<<-dist
      pcon_={}
      for (iex in 1:length(expr)){
         for (ln in 1:length(par)){
       eval(parse(text=paste0("ru","=vpar[,ln]")))
       if (dist=="Weibull"){
       assign(Vnames1[ln],ru)
       if (ln<3 | (ln==length(par) & ncl>0))    assign(Vnames1[ln],exp(ru))
                }
           if (dist=="Gompertz"){
             assign(Vnames1[ln],ru)
             if (ln==length(par) & ncl>0)    assign(Vnames1[ln],exp(ru))
             if (ln==1)    assign(Vnames1[ln],1e-3*exp(ru))
             if (ln==2)    assign(Vnames1[ln],1e-2*exp(ru))
           }
         }
        hi=ecdf(eval(expr[iex]))
        expriex=eval(expr[iex])
        mea=mean(expriex)
        CIL=as.numeric(quantile(expriex,probs=0.025))
        CIU=as.numeric(quantile(expriex,probs=0.975))
        p_expr=2*min(hi(0),1-hi(0))
        pcon_=rbind(pcon_,c(round(mea,4),paste0(as.character(round(CIL,4)),"-",as.character(round(CIU,4))),as.character(round(p_expr,4)),deparse(expr[iex][[1]])))
      }
      pcon=data.frame(pcon_[,1:3])
      colnames(pcon)=c("contrast","CI","p-value")
      rownames(pcon)=pcon_[,4]
      capt=paste("Table of contrasts.",dist,"model.",sep=" ")
      print(xtable(pcon,caption=capt))
          } else {
    print(c("Variable 'expr' must have type 'expression'"))  
          }
  }
  vpar<<- mvrnorm(n = 10000, Result$par, ginv(Result$hessian), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  pstrata=NULL 
  if (!is.null(strata)){
      cn=colnames(D1)
      cs=names(strata)
      nf1<<-nf
      nk1<<-nk
      ncl1<<-ncl
      dist1<<-dist
      psurv_={}
      pcumhaz_={}
      phaz_={}
      D1<<-obs[,-nrow(obs)]
        if ((nk+nf)>0) {
          D1[,1:(nk+nf)]=0
        for (ik in 1:(nk+nf)){
          ix=which(cs==substr(cn[ik],1,(nchar(cn[ik])-6)))
          if (length(ix)>0) D1[,ik]=as.numeric(strata[ix])
        }
          for (ik in 1:(nk+nf)){
            ss=1
            iv=0
            for (ix in 1:length(cs)){
            if (regexpr(cs[ix],cn[ik])[1]>0) {
              ss=ss*as.numeric(strata[ix])
              iv=iv+1
            }
            }
            if (iv>0) D1[,ik]=ss
         }
      }
      D1<<-D1[order(D1$time),]
      Time=D1$time
        for (ln in 1:length(par)){
          eval(parse(text=paste0("ru","=vpar[,ln]")))
          if (dist=="Weibull"){
            assign(Vnames1[ln],ru)
            if (ln<3 | (ln==length(par) & ncl>0))    assign(Vnames1[ln],exp(ru))
          }
          if (dist=="Gompertz"){
            assign(Vnames1[ln],ru)
            if (ln==length(par) & ncl>0)    assign(Vnames1[ln],exp(ru))
            if (ln==1)    assign(Vnames1[ln],1e-3*exp(ru))
            if (ln==2)    assign(Vnames1[ln],1e-2*exp(ru))
          }
        }
        for (iy in 1:nrow(D1)){
        expriex=eval(expression(msurv(iy)))
        mea=mean(expriex)
        med=median(expriex)
        CIL=as.numeric(quantile(expriex,probs=0.025))
        CIU=as.numeric(quantile(expriex,probs=0.975))
        psurv_=rbind(psurv_,c(round(mea,4),round(med,4),round(CIL,4),round(CIU,4)))
        colnames(psurv_)=c("Mean survival","Median survival","Low survival","Upper survival")
        expriex=eval(expression(cumhazard(iy)))
        mea=mean(expriex)
        med=median(expriex)
        CIL=as.numeric(quantile(expriex,probs=0.025))
        CIU=as.numeric(quantile(expriex,probs=0.975))
        pcumhaz_=rbind(pcumhaz_,c(round(mea,4),round(med,4),round(CIL,4),round(CIU,4)))
        colnames(pcumhaz_)=c("Mean cumulative hazard","Median cumulative hazard","Low cumulative hazard","Upper cumulative hazard")
        expriex=eval(expression(hazard(iy)))
        mea=mean(expriex)
        med=median(expriex)
        CIL=as.numeric(quantile(expriex,probs=0.025))
        CIU=as.numeric(quantile(expriex,probs=0.975))
        phaz_=rbind(phaz_,c(round(mea,4),round(med,4),round(CIL,4),round(CIU,4)))
        colnames(phaz_)=c("Mean hazard","Median hazard","Low hazard","Upper hazard")
        }
      Time=sort(Time)
      pstrata=data.frame(Time,psurv_,pcumhaz_,phaz_)
      rownames(pstrata)<-NULL
  }
  Ti=print(c(paste0("Total calculation time=",as.character(round(as.numeric((proc.time() - ptm))[1]/60,4))," (min). ")))
  
  list(par=par,se=invHes,LogLik=Lik,Tab=Tab,Names=Vnames1,Conc=Conc,pval=pval0,p.contrast=pcon,pstrata=pstrata) 
}

#######################################
NamFact=function(data,formula.scale,formula.shape){
#This function returns the list of the factor names in the study after convering categorical variables in binary ones  
nsc=names(as.data.frame(model.matrix.lm(formula.scale, data = data,na.action='na.pass')))[-1]  
nsh=names(as.data.frame(model.matrix.lm(formula.shape, data = data,na.action='na.pass')))[-1]  
if ((length(nsc)+length(nsh))>0) {
return(sort(unique(c(nsc,nsh))))} else {
  return(NULL)}
}
#######################################
LikGenNPH=function(b,D,nf,nk,ncl,dist){
  #######This function calculates negloglikelihood for Weibull or Gompertz survival model
  if (dist=='Weibull'){
    lambda0=exp(b[1])  
    k0=exp(b[2])} else if (dist=='Gompertz') {
      lambda0=1e-3*exp(b[1])  
      k0=1e-2*exp(b[2])
    } 
  Coxk=c(rep(0,nrow(D)))
  Cox=c(rep(0,nrow(D)))
  if (nk>0){
    for (i in 1:nk){
      Coxk=Coxk+D[,i]*b[2+i]
    }
  }
  if (nf>0){
    for (i in 1:nf){
      Cox=Cox+D[,(nk+i)]*b[2+nk+i]
    }
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
  if (ncl>0) {
    G2=exp(b[3+nk+nf])  
    list=unique(D$cluster)  
    nl=length(list)    
  }
  k0=k0*Coxk
    Cens=D$event
    x1=D$time
    x0=D$trunc
    if (dist=='Weibull'){
    if (ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        nc1=sum(ic1)
        Hfull1= Cox[ind]*(x1[ind]/lambda0)^k0[ind] 
        Hfull0= Cox[ind]*(x0[ind]/lambda0)^k0[ind] 
        mufull1=Cox[ind]*k0[ind]*x1[ind]^(k0[ind]-1)/lambda0^k0[ind]#Weibull
        Lmufull1=LCox[ind]+((k0[ind]-1)*log(x1[ind]/lambda0)+log(k0[ind]/lambda0))
        if (nc1<=1) {
          ee1=0} else {
            ee1=sum(log(c(1:(nc1-1))*G2+1))
          }
        Lik=Lik+ee1+sum(Lmufull1*ic1)-(1/G2+nc1)*log(1+G2*sum(Hfull1))+(1/G2)*log(1+G2*sum(Hfull0))
      }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*(x1/lambda0)^k0 
      Hfull0= Cox*(x0/lambda0)^k0 
      Lmufull1=LCox+((k0-1)*log(x1/lambda0)+log(k0/lambda0))
      Lik=sum(Lmufull1*ic)-sum(Hfull1)+sum(Hfull0)
      
    }
      negH<<-c(Cox*((x0/lambda0)^k0-(x1/lambda0)^k0))
      
  }
  
  if (dist=='Gompertz'){
    if (ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        nc1=sum(ic1)
        Hfull1= Cox[ind]*(lambda0/k0[ind])*(exp(k0[ind]*x1[ind])-1)
        Hfull0= Cox[ind]*(lambda0/k0[ind])*(exp(k0[ind]*x0[ind])-1)
        mufull1=Cox[ind]*lambda0[ind]*exp(k0[ind]*x1[ind])
        Lmufull1=LCox[ind]+log(lambda0)+k0[ind]*x1[ind]
        if (nc1<=1) {
          ee1=0} else {
            ee1=sum(log(c(1:(nc1-1))*G2+1))
          }
        Lik=Lik+ee1+sum(Lmufull1*ic1)-(1/G2+nc1)*log(1+G2*sum(Hfull1))+(1/G2)*log(1+G2*sum(Hfull0))
      }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*(lambda0/k0)*(exp(k0*x1)-1)
      Hfull0= Cox*(lambda0/k0)*(exp(k0*x0)-1)
      Lmufull1=LCox+log(lambda0)+k0*x1
      Lik=sum(Lmufull1*ic)-sum(Hfull1)+sum(Hfull0)
      
    }
    negH<<-c(Cox*(lambda0/k0)*(exp(k0*x0)-exp(k0*x1)))
  }
  
  Lik=-Lik
  return(Lik)
}
#####################################
GrGenNPH=function(b,D,nf,nk,ncl,dist){
  #######This function calculates negloglikelihood for Weibull or Gompertz survival model
  if (dist=='Weibull'){
    lambda0=exp(b[1])  
    k0=exp(b[2])} else if (dist=='Gompertz') {
      lambda0=1e-3*exp(b[1])  
      k0=1e-2*exp(b[2])
    } 
  Coxk=c(rep(0,nrow(D)))
  Cox=c(rep(0,nrow(D)))
  if (nk>0){
    for (i in 1:nk){
      Coxk=Coxk+D[,i]*b[2+i]
    }
  }
  if (nf>0){
    for (i in 1:nf){
      Cox=Cox+D[,(nk+i)]*b[2+nk+i]
    }
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
    if (ncl>0) {
    G2=exp(b[3+nk+nf])  
    }
    list=unique(D$cluster)  
    nl=length(list)    

  k0=k0*Coxk
  Cens=D$event
  x1=D$time
  x0=D$trunc
  Gr=c(rep(0,length(b)))
  if (dist=='Weibull'){
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        nc1=sum(ic1)
        Hfull1= Cox[ind]*(x1[ind]/lambda0)^k0[ind] 
        Hfull0= Cox[ind]*(x0[ind]/lambda0)^k0[ind] 
        H_lambda_1=-Hfull1*k0[ind]
        lmu_lambda_1=-k0[ind]
        H_lambda_0=-Hfull0*k0[ind]
        H_k0_1=Hfull1*(log(x1[ind])-log(lambda0))*k0[ind]
        lmu_k0_1=k0[ind]*(log(x1[ind])-log(lambda0))+1
        H_k0_0=c(rep(0,nn))
        H_k0_0[x0[ind]>0]=Hfull0[x0[ind]>0]*(log(x0[ind][x0[ind]>0])-log(lambda0))*k0[ind][x0[ind]>0]
        H_beta_1={}
        if (nf>0)  H_beta_1=matrix(rep(Hfull1,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
        H_beta_0={}
        if (nf>0)  H_beta_0=matrix(rep(Hfull0,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
        H_betak_1={}
        if (nk>0)  {
          H_betak_1=matrix(rep(Hfull1*(log(x1[ind])-log(lambda0))*k0[ind],nk),nn,nk)*D[ind,1:nk]
        }
        H_betak_0={}
        if (nk>0)  {
          H_betak_0=matrix(0,nn,nk)
          if (sum(x0[ind]>0)>0) H_betak_0[x0[ind]>0,]=matrix(rep(Hfull0*(log(x1[ind])-log(lambda0))*k0[ind],nk),nn,nk)[x0[ind]>0,]*D[ind,1:nk][x0[ind]>0,]
        }
        lmu_beta_1=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
        lmu_beta_0=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
        lmu_betak_1={}
        if (nk>0)  {
        lmu_betak_1=matrix(rep(ic1*(k0[ind]*(log(x1[ind])-log(lambda0))+1),nk),nn,nk)*D[ind,1:nk]
        }
        if (ncl>0){
        Gr[1]=Gr[1]-(1+G2*nc1)*sum(H_lambda_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)/(1+G2*sum(Hfull0))        
        Gr[2]=Gr[2]-(1+G2*nc1)*sum(H_k0_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_k0_1)+sum(H_k0_0)/(1+G2*sum(Hfull0))
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-(1+G2*nc1)*c(apply(H_betak_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))/(1+G2*sum(Hfull0))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-(1+G2*nc1)*c(apply(H_beta_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(H_beta_0,2,sum))/(1+G2*sum(Hfull0))+c(apply(lmu_beta_1,2,sum))
        Gr[length(b)]=Gr[length(b)]+log(1+G2*sum(Hfull1))/G2-(1+G2*nc1)*sum(Hfull1)/(1+G2*sum(Hfull1))+(digamma(1/G2)-digamma(nc1+1/G2))/G2-log(1+G2*sum(Hfull0))/G2-sum(Hfull0)/(1+G2*sum(Hfull0))+nc1
        } else {        
          Gr[1]=Gr[1]-sum(H_lambda_1)+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)        
          Gr[2]=Gr[2]-sum(H_k0_1)+sum(ic1*lmu_k0_1)+sum(H_k0_0)
          if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-c(apply(H_betak_1,2,sum))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))
          if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-c(apply(H_beta_1,2,sum))+c(apply(lmu_beta_1,2,sum))+c(apply(H_beta_0,2,sum))
    }
  }
}
  
  if (dist=='Gompertz'){
    for (i in 1:nl){
      ID=list[i]
      ind=which(D$cluster==ID)
      nn=length(ind)
      ic1=1*(Cens[ind]==1)
      nc1=sum(ic1)
      Hfull1= Cox[ind]*(lambda0/k0[ind])*(exp(k0[ind]*x1[ind])-1)
      Hfull0= Cox[ind]*(lambda0/k0[ind])*(exp(k0[ind]*x0[ind])-1)
      mufull1=Cox[ind]*lambda0*exp(k0[ind]*x1[ind])
      mufull0=Cox[ind]*lambda0*exp(k0[ind]*x0[ind])
      Lmufull1=LCox[ind]+log(lambda0)+k0[ind]*x1[ind]
      H_lambda_1=Hfull1
      H_lambda_0=Hfull0
      lmu_lambda_1=1
      H_k0_1=-Hfull1+mufull1*x1[ind]
      H_k0_0=-Hfull0+mufull0*x0[ind]
      lmu_k0_1=k0[ind]*x1[ind]
      H_beta_1={}
      if (nf>0)  H_beta_1=matrix(rep(Hfull1,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
      H_beta_0={}
      if (nf>0)  H_beta_0=matrix(rep(Hfull0,nf),nn,nf)*D[ind,(1+nk):(nf+nk)]
      H_betak_1={}
      if (nk>0)  {
        H_betak_1=matrix(rep(-Hfull1+mufull1*x1[ind],nk),nn,nk)*D[ind,1:nk]
      }
      H_betak_0={}
      if (nk>0)  {
        H_betak_0=matrix(rep(-Hfull0+mufull0*x0[ind],nk),nn,nk)*D[ind,1:nk]
      }
      lmu_beta_1=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
      lmu_beta_0=matrix(rep(ic1,nf),nn,nf)*D[ind,(1+nk):(nk+nf)]
      lmu_betak_1={}
      if (nk>0)  {
        lmu_betak_1=matrix(rep(ic1*k0[ind]*x1[ind],nk),nn,nk)*D[ind,1:nk]
      }
      if (ncl>0){
        Gr[1]=Gr[1]-(1+G2*nc1)*sum(H_lambda_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)/(1+G2*sum(Hfull0))        
        Gr[2]=Gr[2]-(1+G2*nc1)*sum(H_k0_1)/(1+G2*sum(Hfull1))+sum(ic1*lmu_k0_1)+sum(H_k0_0)/(1+G2*sum(Hfull0))
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-(1+G2*nc1)*c(apply(H_betak_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))/(1+G2*sum(Hfull0))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-(1+G2*nc1)*c(apply(H_beta_1,2,sum))/(1+G2*sum(Hfull1))+c(apply(H_beta_0,2,sum))/(1+G2*sum(Hfull0))+c(apply(lmu_beta_1,2,sum))
        Gr[length(b)]=Gr[length(b)]+log(1+G2*sum(Hfull1))/G2-(1+G2*nc1)*sum(Hfull1)/(1+G2*sum(Hfull1))+(digamma(1/G2)-digamma(nc1+1/G2))/G2-log(1+G2*sum(Hfull0))/G2-sum(Hfull0)/(1+G2*sum(Hfull0))+nc1
      } else {        
        Gr[1]=Gr[1]-sum(H_lambda_1)+sum(ic1*lmu_lambda_1)+sum(H_lambda_0)        
        Gr[2]=Gr[2]-sum(H_k0_1)+sum(ic1*lmu_k0_1)+sum(H_k0_0)
        if (nk>0) Gr[3:(2+nk)]=Gr[3:(2+nk)]-c(apply(H_betak_1,2,sum))+c(apply(lmu_betak_1,2,sum))+c(apply(H_betak_0,2,sum))
        if (nf>0) Gr[(3+nk):(2+nk+nf)]=Gr[(3+nk):(2+nk+nf)]-c(apply(H_beta_1,2,sum))+c(apply(lmu_beta_1,2,sum))+c(apply(H_beta_0,2,sum))
      }
    }
  }
  return(-Gr)
}
#####################################
msurv=function(ID){
    #######This function calculates negloglikelihood for Weibull or Gompertz survival model
  if (dist1=='Weibull'){
    lambda0=exp(vpar[,1])  
    k0=exp(vpar[,2])} else if (dist1=='Gompertz') {
      lambda0=1e-3*exp(vpar[,1])  
      k0=1e-2*exp(vpar[,2])
    } 
  Coxk=0
  Cox=0
  if (nk1>0){
    for (i in 1:nk1){
      Coxk=Coxk+D1[ID,i]*vpar[,2+i]}  
  }
  if (nf1>0){
    for (i in 1:nf1){
      Cox=Cox+D1[ID,(nk1+i)]*vpar[,2+nk1+i]}
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
  if (ncl1>0) {
    G2=exp(vpar[,3+nk1+nf1])  
  }
  k0=k0*Coxk
  x1=D1$time[ID]
  x0=D1$trunc[ID]
  if (dist=='Weibull'){
    Hfull1= Cox*(x1/lambda0)^k0 
    Hfull0= Cox*(x0/lambda0)^k0 
    mufull1=Cox*k0*x1^(k0-1)/lambda0^k0
    if (ncl1>0){
      msurv=((1+G2*Hfull1)^(-1/G2))/((1+G2*Hfull0)^(-1/G2))
    } else {
     msurv=exp(Hfull0-Hfull1)
    }
  }
  
  if (dist1=='Gompertz'){
    Hfull1= (Cox*lambda0/k0)*(exp(k0*x1)-1) 
    Hfull0= (Cox*lambda0/k0)*(exp(k0*x0)-1) 
    mufull1= (Cox*lambda0)*exp(k0*x1) 
    
    if (ncl1>0){
      msurv=((1+G2*Hfull1)^(-1/G2))/((1+G2*Hfull0)^(-1/G2))
    } else {
      msurv=exp(Hfull0-Hfull1)
    }
  }
  return(msurv)
}
#####################################
cumhazard=function(ID){
  #######This function calculates negloglikelihood for Weibull or Gompertz survival model
  if (dist1=='Weibull'){
    lambda0=exp(vpar[,1])  
    k0=exp(vpar[,2])} else if (dist1=='Gompertz') {
      lambda0=1e-3*exp(vpar[,1])  
      k0=1e-2*exp(vpar[,2])
    } 
  Coxk=0
  Cox=0
  if (nk1>0){
    for (i in 1:nk1){
      Coxk=Coxk+D1[ID,i]*vpar[,2+i]}
  }
  if (nf1>0){
    for (i in 1:nf1){
      Cox=Cox+D1[ID,(nk1+i)]*vpar[,2+nk1+i]}
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
  if (ncl1>0) {
    G2=exp(vpar[,3+nk1+nf1])  
  }
  k0=k0*Coxk
  x1=D1$time[ID]
  x0=D1$trunc[ID]
  if (dist1=='Weibull'){
    Hfull1= Cox*(x1/lambda0)^k0 
    Hfull0= Cox*(x0/lambda0)^k0 
    mufull1=Cox*k0*x1^(k0-1)/lambda0^k0
    cumh=Hfull1-Hfull0
  }
  
  if (dist1=='Gompertz'){
    Hfull1= (Cox*lambda0/k0)*(exp(k0*x1)-1) 
    Hfull0= (Cox*lambda0/k0)*(exp(k0*x0)-1) 
    mufull1= (Cox*lambda0)*exp(k0*x1) 
    cumh=Hfull1-Hfull0
  }
  return(cumh)
}
#####################################
hazard=function(ID){
  #######This function calculates negloglikelihood for Weibull or Gompertz survival model
  if (dist1=='Weibull'){
    lambda0=exp(vpar[,1])  
    k0=exp(vpar[,2])} else if (dist1=='Gompertz') {
      lambda0=1e-3*exp(vpar[,1])  
      k0=1e-2*exp(vpar[,2])
    } 
  Coxk=0
  Cox=0
  if (nk1>0){
    for (i in 1:nk1){
      Coxk=Coxk+D1[ID,i]*vpar[,2+i]}
  }
  if (nf1>0){
    for (i in 1:nf1){
      Cox=Cox+D1[ID,(nk1+i)]*vpar[,2+nk1+i]}
  }
  Coxk=exp(Coxk)
  LCox=Cox
  Cox=exp(Cox)
  if (ncl1>0) {
    G2=exp(vpar[,3+nk1+nf1])  
  }
  k0=k0*Coxk
  x1=D1$time[ID]
  x0=D1$trunc[ID]
  if (dist1=='Weibull'){
    Hfull1= Cox*(x1/lambda0)^k0 
    Hfull0= Cox*(x0/lambda0)^k0 
    mufull1=Cox*k0*x1^(k0-1)/lambda0^k0
  }
  
  if (dist1=='Gompertz'){
    Hfull1= (Cox*lambda0/k0)*(exp(k0*x1)-1) 
    Hfull0= (Cox*lambda0/k0)*(exp(k0*x0)-1) 
    mufull1= (Cox*lambda0)*exp(k0*x1) 
    cumh=Hfull1-Hfull0
  }
  return(mufull1)
}
###################################
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}

#####################################
