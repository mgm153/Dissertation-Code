library(MASS)
library(gam)
library(mgcv)
library(data.table)
library(dplyr)

#specify the numbe of bootstraps and sample size of bootstrap here#
#this should be the same as the sample size you are estimating for the modeling procedure#
NUM.BOOTS = 200
#samp size and samp boot are the same#
SAMP.BOOT = 400
POP.SIZE = 20000
#this is the number of times you want to run this#
N.POPS <- 1000

#######################################
# STEP 1 - Generate the Population and Sample #
######################################

gen.pop.w.sample = function(seed, N, n, phi){
  set.seed(seed)
  #this function generates one sample according to the model specified in the paper
  
  e12i<-mvrnorm(N,mu=c(0,0),Sigma=matrix(c(1,phi,phi,1),2,2))
  
  x1i=rnorm(N)
  x2i=rnorm(N)     
  
  y1i_star = x1i + e12i[,1]
  y2i_star = x2i + e12i[,2]
  
  # need to put dataset in front of re-code#
  y1i <- as.numeric(y1i_star > 0)
  y2i= y2i_star*y1i
  
  z=rnorm(N,mean=1+y2i_star, sd=0.5)
  k=1/(1+exp(2.5-0.5*z))
  
  #poisson sampling
  Pi=n*k/sum(k) #inclusion probablities
  In=rbinom(N,1,Pi) 
  select=(In==1)
  pop <-  data.frame(x1i=x1i, x2i=x2i, y1i=y1i,y2i=y2i, d=1/Pi,logd=log((1/Pi)-1), Pi=Pi,
                     select=select, N=N, k=k)
  }

#################################
#Population Values #
#also creating a function to clean results into usable format#
################################

clean.results <- function(lm.fit, name.suffix="") {
  result <- as.data.frame(t(coef(lm.fit)))
  names(result) <- paste0(c("Int", "slope", "IMR"), name.suffix)
  return(result)
}

clean.results1=function(data){
  result <- as.data.frame(t(data))
  return(result)
}

get.pop.params <- function(Dat) {
  
  #step 1 - get pop IMR#
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Dat)  
  Bx10i = summary(probit)$coefficients[1,1]
  Bx11i = summary(probit)$coefficients[2,1]
  x1TB1=as.vector((Dat$x1i*Bx11i)+Bx10i)
  IMR=(dnorm(x1TB1, mean=0, sd=1)/pnorm(x1TB1, mean=0, sd=1))
  
  #step 2#
  Pop1=cbind(IMR,Dat)
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMR[y1i==1],  data=Pop1)
  
  return(clean.results(lm.fit, "Pop"))  
}

#######################################
# Creating Pseudo - Population - for use later#
#######################################
#varboots=function(Dat, weight){ #

pseudo.pop.samp=function(Dat,n)
{
  #using sample -getting U = Pseudo-Pop#
  
  weight=round(Dat$d,digits=0)
  Dat1=cbind(Dat,weight)
  U=data.frame(Dat1[rep(seq.int(1,nrow(Dat1)), weight),1:ncol(Dat1) ])
  
  #getting sample- poisson sampling
  U$Pi.PP=n*U$k/sum(U$k)
  N1=dim(U)[1]
  sI<-rbinom(N1,1, U$Pi.PP)
  cbind(sI, U)
  boot=subset(U, sI==1)
  return(boot)
}

################################
# STEP-2  - Creating Weight Functions#
#Using Popl and Psuedopop above to get all results#
################################

#######################
#unsmoothed Estimators#
#######################

# 1 - Design Weight #

DW=function(Dat,b){
  #This function returns Design-Weighted estimator
  Dat <- Dat[Dat$select, ]
  
  #Dat <- p.test[p.test$select, ]
  
  #step 1 - getting probit#
  
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Dat, weight=d)  
  Bx10iw = summary(probit)$coefficients[1,1]
  Bx11iw = summary(probit)$coefficients[2,1]
  x1TB1w=as.vector((Dat$x1i*Bx11iw)+Bx10iw)
  IMRw=(dnorm(x1TB1w, mean=0, sd=1)/pnorm(x1TB1w, mean=0, sd=1))
  
  #step 2 - Getting LM Results#
  Samp1=cbind(IMRw,Dat)
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRw[y1i==1],  data=Samp1, weights=d[y1i==1])
  lm.fit$df.null
  
   Set1=clean.results(lm.fit, "DW")
  
    #Step 3 - Getting Psuedo-Population Results
  
    bootresults=matrix(nrow=NUM.BOOTS,ncol=3)
    
    # Change to NUM.BOOTS
    for (b in 1:NUM.BOOTS){
    boot.samp=pseudo.pop.samp(Samp1, SAMP.BOOT)
    

    lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRw[y1i==1],  data=boot.samp, weights=d[y1i==1])
    bootresults[b,] <- coef(lm.fit)
    }
    
    colnames(bootresults) <- c("intercept","slope", "IMR")
    B0=mean(bootresults[,1])
    B1=mean(bootresults[,2])
    BIMR=mean(bootresults[,3])
    
    data.frame(bootresults)
    
    boot.results=data.frame(cbind(B0,B1,BIMR,bootresults))
    
    B0.sq=(boot.results$intercept-boot.results$B0)**2
    B1.sq=(boot.results$slope-boot.results$B1)**2
    BIMR.sq=(boot.results$IMR-boot.results$BIMR)**2
    
    Var.B0=mean(B0.sq)
    SD.B0=sqrt(Var.B0)
    UB.B0=B0+qt(0.975,lm.fit$df.null)*SD.B0
    LB.B0=B0-qt(0.975,lm.fit$df.null)*SD.B0
    
    Var.B1=mean(B1.sq)
    SD.B1=sqrt(Var.B1)
    UB.B1=B1+qt(0.975,lm.fit$df.null)*SD.B1
    LB.B1=B1-qt(0.975,lm.fit$df.null)*SD.B1
    
    Var.BIMR=mean(BIMR.sq)
    SD.BIMR=sqrt(Var.BIMR)
    UB.BIMR=BIMR+qt(0.975,lm.fit$df.null)*SD.BIMR
    LB.BIMR=BIMR-qt(0.975,lm.fit$df.null)*SD.BIMR
    
    #list1=list(IntDW=Set1$IntDW, slopeDW=Set1$SlopeDW, IMRDW=Set1$IMRDW, B0DW=B0, B1DW=B1, BIMRDW=BIMR,
              # SD.B0.DW=SD.B0, SD.B1.DW=SD.B1, SD.BIMR=SD.BIMR)
    
  return(list(IntDW=Set1$IntDW, slopeDW=Set1$slopeDW, IMRDW=Set1$IMRDW, B0DW=B0, B1DW=B1, BIMRDW=BIMR,
              UB.B0.DW=UB.B0, UB.B1.DW=UB.B1, UB.BIMR.DW=UB.BIMR, LB.B0.DW=LB.B0, LB.B1.DW=LB.B1, LB.BIMR.DW=LB.BIMR,
              SD.B0.DW=SD.B0, SD.B1.DW=SD.B1, SD.BIMR.DW=SD.BIMR))   
}

#DW(Samp)#

# 2 - PS Estimator #

PSGAM=function(Dat){
  Dat <- Dat[Dat$select, ]
  
  #Dat <- p.test[p.test$select, ]
  #This function returns PS-Weighted estimator
  fit.PS=gam(logd~s(x1i), data=Dat, fit=TRUE)
  
  mu.PS=fit.PS$fitted.values
  sigma.PS=sqrt(sum(fit.PS$residuals^2)/fit.PS$df.residual)
  PS=exp(mu.PS+sigma.PS^2/2)+1
  PS1=Dat$d/PS
  
  #step 2 - getting probit#
  Samp1=cbind(PS1,Dat)
  
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Samp1, weight=PS1)  
  Bx10iw = summary(probit)$coefficients[1,1]
  Bx11iw = summary(probit)$coefficients[2,1]
  x1TB1w=as.vector((Dat$x1i*Bx11iw)+Bx10iw)
  IMRPS=(dnorm(x1TB1w, mean=0, sd=1)/pnorm(x1TB1w, mean=0, sd=1))
  
  #step 3 - getting 2nd stage weight #
  Samp2=cbind(IMRPS,Samp1)
  
  fit.PS2=gam(logd[y1i==1]~s(x2i[y1i==1])+s(IMRPS[y1i==1]), data=Samp2, fit=TRUE)
  
  mu.PS2=fit.PS2$fitted.values
  sigma.PS2=sqrt(sum(fit.PS2$residuals^2)/fit.PS2$df.residual)
  mod2=exp(mu.PS2+sigma.PS2^2/2)+1
  Samp3=subset(Samp2, y1i==1)
  PS2=Samp3$d/mod2
  
  #step 4 - getting final estimates #
  Samp4=cbind(PS2 , Samp3)
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRPS[y1i==1],  data=Samp4, weights=PS2[y1i==1])
  
  Set1=clean.results(lm.fit, "PSGAM")
  
  #Step 5 - Getting Psuedo-Population Results
  
  bootresults=matrix(nrow=NUM.BOOTS,ncol=3)
  
  for (b in 1:NUM.BOOTS){
    boot.samp=pseudo.pop.samp(Samp4, SAMP.BOOT)
    
    
    lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRPS[y1i==1],  data=boot.samp, weights=PS2[y1i==1])
    bootresults[b,] <- coef(lm.fit)
  }
  
  colnames(bootresults) <- c("intercept","slope", "IMR")
  B0=mean(bootresults[,1])
  B1=mean(bootresults[,2])
  BIMR=mean(bootresults[,3])
  
  data.frame(bootresults)
  
  boot.results=data.frame(cbind(B0,B1,BIMR,bootresults))
  
  B0.sq=(boot.results$intercept-boot.results$B0)**2
  B1.sq=(boot.results$slope-boot.results$B1)**2
  BIMR.sq=(boot.results$IMR-boot.results$BIMR)**2
  
  Var.B0=mean(B0.sq)
  SD.B0=sqrt(Var.B0)
  UB.B0=B0+qt(0.975,lm.fit$df.null)*SD.B0
  LB.B0=B0-qt(0.975,lm.fit$df.null)*SD.B0
  
  Var.B1=mean(B1.sq)
  SD.B1=sqrt(Var.B1)
  UB.B1=B1+qt(0.975,lm.fit$df.null)*SD.B1
  LB.B1=B1-qt(0.975,lm.fit$df.null)*SD.B1
  
  Var.BIMR=mean(BIMR.sq)
  SD.BIMR=sqrt(Var.BIMR)
  UB.BIMR=BIMR+qt(0.975,lm.fit$df.null)*SD.BIMR
  LB.BIMR=BIMR-qt(0.975,lm.fit$df.null)*SD.BIMR
  
  return(list(IntPSGAM=Set1$IntPSGAM, slopePSGAM=Set1$slopePSGAM, IMRPSGAM=Set1$IMRPSGAM, B0PSGAM=B0, B1PSGAM=B1, BIMRPSGAM=BIMR,
              UB.B0.PSGAM=UB.B0, UB.B1.PSGAM=UB.B1, UB.BIMR.PSGAM=UB.BIMR, LB.B0.PSGAM=LB.B0, LB.B1.PSGAM=LB.B1, LB.BIMR.PSGAM=LB.BIMR,
              SD.B0.PSGAM=SD.B0, SD.B1.PSGAM=SD.B1, SD.BIMR.PSGAM=SD.BIMR))  
}

###################
# 3 - UOPT #
###################

#This function returns UnsmoothedWeighted estimator#

UOPTGAM=function(Dat){
  Dat <- Dat[Dat$select, ]
  
  #Dat <- p.test[p.test$select, ]
  # getting inital beta #
  fit.UOPT=gam(logd~s(x1i), data=Dat, fit=TRUE)
  
  #fitting GAM#
  mu.UOPT=fit.UOPT$fitted.values
  nUOPT=pnorm(mu.UOPT, mean=0, sd=1)
  p1=sum(Dat$y1i)/nrow(Dat)
  denominator=(Dat$d*(1-nUOPT)**2)*p1+Dat$d*((-nUOPT)**2)*(1-p1)
  
  sigma.UOPT=sqrt(sum(fit.UOPT$residuals^2)/fit.UOPT$df.residual)
  PS=exp(mu.UOPT+sigma.UOPT^2/2)+1
  UOPT1=(Dat$d*(nUOPT)*(1-nUOPT))/PS
  
  #step 2 - getting probit#
  Samp1=cbind(UOPT1,Dat)
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Samp1, weight=UOPT1)  
  Bx10iw = summary(probit)$coefficients[1,1]
  Bx11iw = summary(probit)$coefficients[2,1]
  x1TB1w=as.vector((Samp1$x1i*Bx11iw)+Bx10iw)
  IMRUOPT=(dnorm(x1TB1w, mean=0, sd=1)/pnorm(x1TB1w, mean=0, sd=1))
  
  #step 3 - getting 2nd stage weight #
  Samp2=cbind(IMRUOPT,Samp1)
  fit.UOPT2=gam(logd[y1i==1]~s(x2i[y1i==1])+s(IMRUOPT[y1i==1]), data=Samp2, fit=TRUE)
  mu.dt2=fit.UOPT2$fitted.values
  sigma.dt2=sqrt(sum(fit.UOPT2$residuals^2)/fit.UOPT2$df.residual)
  dt2=exp(mu.dt2+sigma.dt2^2/2)+1
  Samp3=subset(Samp2, y1i==1)
  UOPT2=Samp3$d/dt2

  #step 4 - getting final estimates #
  Samp4=cbind(UOPT2 , Samp3)
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRUOPT[y1i==1],  data=Samp3, weights=UOPT2[y1i==1])
  
  Set1=clean.results(lm.fit, "UOPTGAM")
  
  #Step 5 - Getting Psuedo-Population Results
  
  bootresults=matrix(nrow=NUM.BOOTS,ncol=3)
  
  for (b in 1:NUM.BOOTS){
    boot.samp=pseudo.pop.samp(Samp4, SAMP.BOOT)
    
    lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRUOPT[y1i==1],  data=boot.samp, weights=UOPT2[y1i==1])
    bootresults[b,] <- coef(lm.fit)
  }
  
  colnames(bootresults) <- c("intercept","slope", "IMR")
  B0=mean(bootresults[,1])
  B1=mean(bootresults[,2])
  BIMR=mean(bootresults[,3])
  
  data.frame(bootresults)
  
  boot.results=data.frame(cbind(B0,B1,BIMR,bootresults))
  
  B0.sq=(boot.results$intercept-boot.results$B0)**2
  B1.sq=(boot.results$slope-boot.results$B1)**2
  BIMR.sq=(boot.results$IMR-boot.results$BIMR)**2
  
  Var.B0=mean(B0.sq)
  SD.B0=sqrt(Var.B0)
  UB.B0=B0+qt(0.975,lm.fit$df.null)*SD.B0
  LB.B0=B0-qt(0.975,lm.fit$df.null)*SD.B0
  
  Var.B1=mean(B1.sq)
  SD.B1=sqrt(Var.B1)
  UB.B1=B1+qt(0.975,lm.fit$df.null)*SD.B1
  LB.B1=B1-qt(0.975,lm.fit$df.null)*SD.B1
  
  Var.BIMR=mean(BIMR.sq)
  SD.BIMR=sqrt(Var.BIMR)
  UB.BIMR=BIMR+qt(0.975,lm.fit$df.null)*SD.BIMR
  LB.BIMR=BIMR-qt(0.975,lm.fit$df.null)*SD.BIMR
  
  return(list(IntUOPTGAM=Set1$IntUOPTGAM, slopeUOPTGAM=Set1$slopeUOPTGAM, IMRUOPTGAM=Set1$IMRUOPTGAM, B0UOPTGAM=B0, B1UOPTGAM=B1, BIMRUOPTGAM=BIMR,
              UB.B0.UOPTGAM=UB.B0, UB.B1.UOPTGAM=UB.B1, UB.BIMR.UOPTGAM=UB.BIMR, LB.B0.UOPTGAM=LB.B0, LB.B1.UOPTGAM=LB.B1, LB.BIMR.UOPTGAM=LB.BIMR,
              SD.B0.UOPTGAM=SD.B0, SD.B1.UOPTGAM=SD.B1, SD.BIMR.UOPTGAM=SD.BIMR))  
}
  
#################################
#Smoothed Estimators#
################################


# 4 - Beaumont Estimator #

#Weight at Stage 1 and Stage 2 - nonparametric estimator#

SDWGAM=function(Dat){
  Dat <- Dat[Dat$select, ]
  
  #Dat <- p.test[p.test$select, ]
  
  #step 1 getting stage 1 weight #
  fit.SDW=gam(logd~s(x1i)+y1i, data=Dat)
  mu.SDW=fit.SDW$fitted.values
  sigma.SDW=sqrt(sum(fit.SDW$residuals^2)/fit.SDW$df.residual)
  SDW1=exp(mu.SDW+sigma.SDW^2/2)+1
  
  Samp1=cbind(SDW1,Dat)
  
  #step 2 - getting probit#
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Samp1, weight=SDW1)  
  Bx10iw = summary(probit)$coefficients[1,1]
  Bx11iw = summary(probit)$coefficients[2,1]
  x1TB1w=as.vector((Samp1$x1i*Bx11iw)+Bx10iw)
  IMRSDW=(dnorm(x1TB1w, mean=0, sd=1)/pnorm(x1TB1w, mean=0, sd=1))
  
  #Step 3 - compute d2 #
  
  Samp2=cbind(IMRSDW,Samp1)
  
  fit.SDW2=gam(logd~s(x2i)+s(y2i)+s(IMRSDW), data=Samp2)
  muSDW.dt2=fit.SDW2$fitted.values
  sigma.dt=sqrt(sum(fit.SDW2$residuals^2)/fit.SDW2$df.residual)
  SDW2=exp(muSDW.dt2+sigma.dt^2/2)+1
  Samp3=cbind(SDW2, Samp2)
  
  #step 4 #
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRSDW[y1i==1],  data=Samp3, weights=SDW2[y1i==1])
  
  #return(clean.results(lm.fit, "SDWGAM"))
  
  Set1=clean.results(lm.fit, "SDWGAM")
  
  #Step 5 - Getting Psuedo-Population Results
  
  bootresults=matrix(nrow=NUM.BOOTS,ncol=3)
  
  for (b in 1:NUM.BOOTS){
    boot.samp=pseudo.pop.samp(Samp3, SAMP.BOOT)
    
    
    lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRSDW[y1i==1],  data=boot.samp, weights=SDW2[y1i==1])
    bootresults[b,] <- coef(lm.fit)
  }
  
  colnames(bootresults) <- c("intercept","slope", "IMR")
  B0=mean(bootresults[,1])
  B1=mean(bootresults[,2])
  BIMR=mean(bootresults[,3])
  
  data.frame(bootresults)
  
  boot.results=data.frame(cbind(B0,B1,BIMR,bootresults))
  
  B0.sq=(boot.results$intercept-boot.results$B0)**2
  B1.sq=(boot.results$slope-boot.results$B1)**2
  BIMR.sq=(boot.results$IMR-boot.results$BIMR)**2
  
  Var.B0=mean(B0.sq)
  SD.B0=sqrt(Var.B0)
  UB.B0=B0+qt(0.975,lm.fit$df.null)*SD.B0
  LB.B0=B0-qt(0.975,lm.fit$df.null)*SD.B0
  
  Var.B1=mean(B1.sq)
  SD.B1=sqrt(Var.B1)
  UB.B1=B1+qt(0.975,lm.fit$df.null)*SD.B1
  LB.B1=B1-qt(0.975,lm.fit$df.null)*SD.B1
  
  Var.BIMR=mean(BIMR.sq)
  SD.BIMR=sqrt(Var.BIMR)
  UB.BIMR=BIMR+qt(0.975,lm.fit$df.null)*SD.BIMR
  LB.BIMR=BIMR-qt(0.975,lm.fit$df.null)*SD.BIMR
  
  return(list(IntSDWGAM=Set1$IntSDWGAM, slopeSDWGAM=Set1$slopeSDWGAM, IMRSDWGAM=Set1$IMRSDWGAM, B0SDWGAM=B0, B1SDWGAM=B1, BIMRSDWGAM=BIMR,
              UB.B0.SDWGAM=UB.B0, UB.B1.SDWGAM=UB.B1, UB.BIMR.SDWGAM=UB.BIMR, LB.B0.SDWGAM=LB.B0, LB.B1.SDWGAM=LB.B1, LB.BIMR.SDWGAM=LB.BIMR,
              SD.B0.SDWGAM=SD.B0, SD.B1.SDWGAM=SD.B1, SD.BIMR.SDWGAM=SD.BIMR)) 
  
}


# 5- SPS #

SPSGAM=function(Dat) {
  Dat <- Dat[Dat$select, ]
  
  #Dat <- p.test[p.test$select, ]
  
  #step 1 getting stage 1 weight #
  
  #numerator#
  
  fit.SDW=gam(logd~s(x1i)+y1i, data=Dat, fit=TRUE )
  mu.SDW=fit.SDW$fitted.values
  sigma.SDW=sqrt(sum(fit.SDW$residuals^2)/fit.SDW$df.residual)
  SDW1=exp(mu.SDW+sigma.SDW^2/2)+1
  
  #denominator#
  fit.PS=gam(logd~s(x1i), data=Dat, fit=TRUE)
  mu.PS=fit.PS$fitted.values
  sigma.PS=sqrt(sum(fit.PS$residuals^2)/fit.PS$df.residual)
  PS=exp(mu.PS+sigma.PS^2/2)+1
  SPS1=SDW1/PS
  
  Samp1=cbind(SPS1,Dat)
  
  #step 2 - getting probit#
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Samp1, weight=SPS1)#  
  Bx10iw = summary(probit)$coefficients[1,1]
  Bx11iw = summary(probit)$coefficients[2,1]
  x1TB1w=as.vector((Samp1$x1i*Bx11iw)+Bx10iw)
  IMRSPS=(dnorm(x1TB1w, mean=0, sd=1)/pnorm(x1TB1w, mean=0, sd=1))
  
  #Step 3 - compute d2 #
  
  Samp2=cbind(IMRSPS,Samp1)
  fit.dt=gam(logd~s(x2i)+s(y2i)+s(IMRSPS), data=Samp2)
  mu.dt=fit.dt$fitted.values
  sigma.dt=sqrt(sum(fit.dt$residuals^2)/fit.dt$df.residual)
  SPS2=exp(mu.dt+sigma.dt^2/2)+1
  
  Samp3=cbind(SPS2,Samp2)
  
  #step 4 #
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRSPS[y1i==1],  data=Samp3, weights=SPS2[y1i==1])
  
  Set1=clean.results(lm.fit, "SPSGAM")
  
  #Step 5 - Getting Psuedo-Population Results#
  
  bootresults=matrix(nrow=NUM.BOOTS,ncol=3)
  
  for (b in 1:NUM.BOOTS)  {
    boot.samp=pseudo.pop.samp(Samp3, SAMP.BOOT)
    
    
    lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRSPS[y1i==1],  data=boot.samp, weights=SPS2[y1i==1])
    bootresults[b,] <- coef(lm.fit)
  }
  
  colnames(bootresults) <- c("intercept","slope", "IMR")
  B0=mean(bootresults[,1])
  B1=mean(bootresults[,2])
  BIMR=mean(bootresults[,3])
  
  data.frame(bootresults)
  
  boot.results=data.frame(cbind(B0,B1,BIMR,bootresults))
  
  B0.sq=(boot.results$intercept-boot.results$B0)**2
  B1.sq=(boot.results$slope-boot.results$B1)**2
  BIMR.sq=(boot.results$IMR-boot.results$BIMR)**2
  
  Var.B0=mean(B0.sq)
  SD.B0=sqrt(Var.B0)
  UB.B0=B0+qt(0.975,lm.fit$df.null)*SD.B0
  LB.B0=B0-qt(0.975,lm.fit$df.null)*SD.B0
  
  Var.B1=mean(B1.sq)
  SD.B1=sqrt(Var.B1)
  UB.B1=B1+qt(0.975,lm.fit$df.null)*SD.B1
  LB.B1=B1-qt(0.975,lm.fit$df.null)*SD.B1
  
  Var.BIMR=mean(BIMR.sq)
  SD.BIMR=sqrt(Var.BIMR)
  UB.BIMR=BIMR+qt(0.975,lm.fit$df.null)*SD.BIMR
  LB.BIMR=BIMR-qt(0.975,lm.fit$df.null)*SD.BIMR
  
  return(list(IntSPSGAM=Set1$IntSPSGAM, slopeSPSGAM=Set1$slopeSPSGAM, IMRSPSGAM=Set1$IMRSPSGAM, B0SPSGAM=B0, B1SPSGAM=B1, BIMRSPSGAM=BIMR,
              UB.B0.SPSGAM=UB.B0, UB.B1.SPSGAM=UB.B1, UB.BIMR.SPSGAM=UB.BIMR, LB.B0.SPSGAM=LB.B0, LB.B1.SPSGAM=LB.B1, LB.BIMR.SPSGAM=LB.BIMR,
              SD.B0.SPSGAM=SD.B0, SD.B1.SPSGAM=SD.B1, SD.BIMR.SPSGAM=SD.BIMR)) 
  
}

# 6- SOPT #

SOPTGAM=function(Dat){
  Dat <- Dat[Dat$select, ]
  
  #Dat <- p.test[p.test$select, ]
  
  # Step 1a - getting inital beta #
  fit.SOPT=gam(logd~s(x1i), data=Dat, fit=TRUE)
  
  # step 1b - fitting GAM for denominator weight#
  mu.SOPT=fit.SOPT$fitted.values
  nSOPT=pnorm(mu.SOPT, mean=0, sd=1)
  p1=sum(Dat$y1i)/nrow(Dat)
  denominator=(Dat$d*(1-nSOPT)**2)*p1+Dat$d*((-nSOPT)**2)*(1-p1)
  
  sigma.SOPT=sqrt(sum(fit.SOPT$residuals^2)/fit.SOPT$df.residual)
  denom=exp(mu.SOPT+sigma.SOPT^2/2)+1
  
  #step 1c getting numerator weight #
  fit.SDW=gam(logd~s(x1i)+y1i, data=Dat)
  mu.SDW=fit.SDW$fitted.values
  sigma.SDW=sqrt(sum(fit.SDW$residuals^2)/fit.SDW$df.residual)
  Dat$SDW1=exp(mu.SDW+sigma.SDW^2/2)+1
  SOPT1=(Dat$SDW1*(nSOPT)*(1-nSOPT))/denom
  
  #step 2 - getting probit#
  Samp1=cbind(SOPT1,Dat)
  probit <- glm(y1i ~ x1i, family = quasibinomial(link = "probit"),data =Samp1, weight=SOPT1)  
  Bx10iw = summary(probit)$coefficients[1,1]
  Bx11iw = summary(probit)$coefficients[2,1]
  x1TB1w=as.vector((Samp1$x1i*Bx11iw)+Bx10iw)
  IMRSOPT=(dnorm(x1TB1w, mean=0, sd=1)/pnorm(x1TB1w, mean=0, sd=1))
  
  #step 3 - getting 2nd stage weight #
  Samp2=cbind(IMRSOPT,Samp1)
  
  fit.SOPT2=gam(logd[y1i==1]~s(x2i[y1i==1])+s(IMRSOPT[y1i==1]), data=Samp2, fit=TRUE)
  mu.dt2=fit.SOPT2$fitted.values
  sigma.dt2=sqrt(sum(fit.SOPT2$residuals^2)/fit.SOPT2$df.residual)
  denom2=exp(mu.dt2+sigma.dt2^2/2)+1
  Samp3=subset(Samp2, y1i==1)
  
  fit.SDW2=gam(logd[y1i==1]~s(x2i[y1i==1])+s(y2i[y1i==1])+s(IMRSOPT[y1i==1]), data=Samp2)
  muSDW.dt2=fit.SDW2$fitted.values
  sigma.dt=sqrt(sum(fit.SDW2$residuals^2)/fit.SDW2$df.residual)
  Samp3$SDW2=exp(muSDW.dt2+sigma.dt^2/2)+1
  
  SOPT2=Samp3$SDW2/denom2
  
  
  #step 4 - getting final estimates #
  Samp4=cbind(SOPT2 , Samp3)
  lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRSOPT[y1i==1],  data=Samp3, weights=SOPT2[y1i==1])
  
  #return(clean.results(lm.fit, "SOPTGAM"))
  
  Set1=clean.results(lm.fit, "SOPTGAM")
  
  #Step 5 - Getting Psuedo-Population Results
  
  bootresults=matrix(nrow=NUM.BOOTS,ncol=3)
  
  for (b in 1:NUM.BOOTS){
    boot.samp=pseudo.pop.samp(Samp4, SAMP.BOOT)
    
    
    lm.fit=glm(y2i[y1i==1] ~ x2i[y1i==1] + IMRSOPT[y1i==1],  data=boot.samp, weights=SOPT2[y1i==1])
    bootresults[b,] <- coef(lm.fit)
  }
  
  colnames(bootresults) <- c("intercept","slope", "IMR")
  B0=mean(bootresults[,1])
  B1=mean(bootresults[,2])
  BIMR=mean(bootresults[,3])
  
  data.frame(bootresults)
  
  boot.results=data.frame(cbind(B0,B1,BIMR,bootresults))
  
  B0.sq=(boot.results$intercept-boot.results$B0)**2
  B1.sq=(boot.results$slope-boot.results$B1)**2
  BIMR.sq=(boot.results$IMR-boot.results$BIMR)**2
  
  Var.B0=mean(B0.sq)
  SD.B0=sqrt(Var.B0)
  UB.B0=B0+qt(0.975,lm.fit$df.null)*SD.B0
  LB.B0=B0-qt(0.975,lm.fit$df.null)*SD.B0
  
  Var.B1=mean(B1.sq)
  SD.B1=sqrt(Var.B1)
  UB.B1=B1+qt(0.975,lm.fit$df.null)*SD.B1
  LB.B1=B1-qt(0.975,lm.fit$df.null)*SD.B1
  
  Var.BIMR=mean(BIMR.sq)
  SD.BIMR=sqrt(Var.BIMR)
  UB.BIMR=BIMR+qt(0.975,lm.fit$df.null)*SD.BIMR
  LB.BIMR=BIMR-qt(0.975,lm.fit$df.null)*SD.BIMR
  
  return(list(IntSOPTGAM=Set1$IntSOPTGAM, slopeSOPTGAM=Set1$slopeSOPTGAM, IMRSOPTGAM=Set1$IMRSOPTGAM, B0SOPTGAM=B0, B1SOPTGAM=B1, BIMRSOPTGAM=BIMR,
              UB.B0.SOPTGAM=UB.B0, UB.B1.SOPTGAM=UB.B1, UB.BIMR.SOPTGAM=UB.BIMR, LB.B0.SOPTGAM=LB.B0, LB.B1.SOPTGAM=LB.B1, LB.BIMR.SOPTGAM=LB.BIMR,
              SD.B0.SOPTGAM=SD.B0, SD.B1.SOPTGAM=SD.B1, SD.BIMR.SOPTGAM=SD.BIMR))  
}



#######################################################
# Step 3 - bootstrapping and getting Data for Analysis#
#######################################################


#library(boot)#
#resultsDW=boot(data=Samp, statistic=DW, R=500)


p.test <- gen.pop.w.sample(12345, POP.SIZE, SAMP.BOOT, 0.8)

#get.pop.params(p.test)
#DW(p.test)
#PSGAM(p.test)
#SDWGAM(p.test)
#SPSGAM(p.test)
#UOPTGAM(p.test)
#SOPTGAM(p.test)

custom.iter <- function(seed) {
  p.test <- gen.pop.w.sample(seed, POP.SIZE, SAMP.BOOT, 0.8)
  result <- cbind(data.frame(seed=seed),
                  get.pop.params(p.test),
                  clean.results1(DW(p.test)),
                  clean.results1(PSGAM(p.test)),
                  clean.results1(SDWGAM(p.test)),
                  clean.results1(SPSGAM(p.test)),
                  clean.results1(UOPTGAM(p.test)),
                  clean.results1(SOPTGAM(p.test)))
  
  # insert other methods#
  
  ### Fix for the columns showing up as lists (causing "non-numeric" errors later)
  result %<>% lapply(unlist) %>% as.data.frame
  
  return(result)
}

####################
#Getting Data - part 2#
####################

#Getting Data from all the weights and iterating number of desired times
# Generate a set of seeds = Desired number of Monte Carlo Sims#

set.seed(0)
seeds <- floor(runif(N.POPS, min=10000, max=100000))

results <- lapply(seeds, custom.iter)

########################################
# Stop Here #
#########################################



test <- rbindlist(results, use.names=T)

#data for table 1 - Bias for Different Weights#
myvars=c("IntPop","slopePop", "IMRPop", 
         "IntDW","slopeDW", "IMRDW", 
         "IntPSGAM","slopePSGAM", "IMRPSGAM",
         "IntUOPTGAM","slopeUOPTGAM", "IMRUOPTGAM", 
         "IntSDWGAM","slopeSDWGAM", "IMRSDWGAM",
         "IntSPSGAM","slopeSPSGAM", "IMRSPSGAM", 
         "IntSOPTGAM","slopeSOPTGAM", "IMRSOPTGAM")
test1=select(test, myvars)

#data for table 2 -Variace Values#
#myvars1=c("B0DW", "B1DW", "BIMRDW",
          #"B0PSGAM", "B1PSGAM", "BIMRPSGAM" ,
         # "B0UOPTGAM", "B1UOPTGAM", "BIMRUOPTGAM",
         # "B0SDW", "B1SDW", "BIMRSDW",
         # "B0SPSGAM", "B1SPSGAM", "BIMRSPSGAM",
         # "B0SOPTGAM", "B1SOPTGAM", "BIMRSOPTGAM"
         # )

test2=test[, c("B0DW", "B1DW", "BIMRDW",
              "B0PSGAM", "B1PSGAM", "BIMRPSGAM" ,
              "B0UOPTGAM", "B1UOPTGAM", "BIMRUOPTGAM",
              "B0SDWGAM", "B1SDWGAM", "BIMRSDWGAM",
              "B0SPSGAM", "B1SPSGAM", "BIMRSPSGAM",
              "B0SOPTGAM", "B1SOPTGAM", "BIMRSOPTGAM") := NULL]

test2=na.omit(test2)

##########################
# Step 4 - Getting Results#
###########################
cov.test=mean(ifelse(test2$IntPop >test2$LB.B0.DW & test2$IntPop >test2$LB.B0.DW,1,0))
cov.test1=mean(ifelse(test2$slopePop >test2$LB.B1.DW & test2$slopePop >test2$LB.B1.DW,1,0))
cov.test2=mean(ifelse(test2$IMRPop >test2$LB.BIMR.DW & test2$IMRPop >test2$LB.BIMR.DW,1,0))

cov.test3=mean(ifelse(test2$IntPop >test2$LB.B0.PSGAM & test2$IntPop >test2$LB.B0.PSGAM,1,0))
cov.test4=mean(ifelse(test2$slopePop >test2$LB.B1.PSGAM & test2$slopePop >test2$LB.B1.PSGAM,1,0))
cov.test5=mean(ifelse(test2$IMRPop >test2$LB.BIMR.PSGAM & test2$IMRPop >test2$LB.BIMR.PSGAM,1,0))

cov.test6=mean(ifelse(test2$IntPop >test2$LB.B0.UOPTGAM & test2$IntPop >test2$LB.B0.UOPTGAM,1,0))
cov.test7=mean(ifelse(test2$slopePop >test2$LB.B1.UOPTGAM & test2$slopePop >test2$LB.B1.UOPTGAM,1,0))
cov.test8=mean(ifelse(test2$IMRPop >test2$LB.BIMR.UOPTGAM & test2$IMRPop >test2$LB.BIMR.UOPTGAM,1,0))

cov.test9=mean(ifelse(test2$IntPop >test2$LB.B0.SDWGAM & test2$IntPop >test2$LB.B0.SDWGAM,1,0))
cov.test10=mean(ifelse(test2$slopePop >test2$LB.B1.SDWGAM & test2$slopePop >test2$LB.B1.SDWGAM,1,0))
cov.test11=mean(ifelse(test2$IMRPop >test2$LB.BIMR.SDWGAM & test2$IMRPop >test2$LB.BIMR.SDWGAM,1,0))

cov.test12=mean(ifelse(test2$IntPop >test2$LB.B0.SPSGAM & test2$IntPop >test2$LB.B0.SPSGAM,1,0))
cov.test13=mean(ifelse(test2$slopePop >test2$LB.B1.SPSGAM & test2$slopePop >test2$LB.B1.SPSGAM,1,0))
cov.test14=mean(ifelse(test2$IMRPop >test2$LB.BIMR.SPSGAM & test2$IMRPop >test2$LB.BIMR.SPSGAM,1,0))

cov.test15=mean(ifelse(test2$IntPop >test2$LB.B0.SOPTGAM & test2$IntPop >test2$LB.B0.SOPTGAM,1,0))
cov.test16=mean(ifelse(test2$slopePop >test2$LB.B1.SOPTGAM & test2$slopePop >test2$LB.B1.SOPTGAM,1,0))
cov.test17=mean(ifelse(test2$IMRPop >test2$LB.BIMR.SOPTGAM & test2$IMRPop >test2$LB.BIMR.SOPTGAM,1,0))

cov.overall=data.frame(rbind(cov.test, cov.test1,cov.test2,cov.test3, cov.test4,cov.test5,
                  cov.test6, cov.test7,cov.test8,cov.test9, cov.test10,cov.test11,
                  cov.test12, cov.test13,cov.test14,cov.test15, cov.test16,cov.test17))

AL=mean(test2$UB.B0.DW-test2$LB.B0.DW)
AL1=mean(test2$UB.B1.DW-test2$LB.B1.DW)
AL2=mean(test2$UB.BIMR.DW-test2$LB.BIMR.DW)

AL3=mean(test2$UB.B0.PSGAM-test2$LB.B0.PSGAM)
AL4=mean(test2$UB.B1.PSGAM-test2$LB.B1.PSGAM)
AL5=mean(test2$UB.BIMR.PSGAM-test2$LB.BIMR.PSGAM)

AL6=mean(test2$UB.B0.UOPTGAM-test2$LB.B0.UOPTGAM)
AL7=mean(test2$UB.B1.UOPTGAM-test2$LB.B1.UOPTGAM)
AL8=mean(test2$UB.BIMR.UOPTGAM-test2$LB.BIMR.UOPTGAM)

AL9=mean(test2$UB.B0.SDWGAM-test2$LB.B0.SDWGAM)
AL10=mean(test2$UB.B1.SDWGAM-test2$LB.B1.SDWGAM)
AL11=mean(test2$UB.BIMR.SDWGAM-test2$LB.BIMR.SDWGAM)

AL12=mean(test2$UB.B0.SPSGAM-test2$LB.B0.SPSGAM)
AL13=mean(test2$UB.B1.SPSGAM-test2$LB.B1.SPSGAM)
AL14=mean(test2$UB.BIMR.SPSGAM-test2$LB.BIMR.SPSGAM)

AL15=mean(test2$UB.B0.SOPTGAM-test2$LB.B0.SOPTGAM)
AL16=mean(test2$UB.B1.SOPTGAM-test2$LB.B1.SOPTGAM)
AL17=mean(test2$UB.BIMR.SOPTGAM-test2$LB.BIMR.SOPTGAM)

AL.overall=data.frame(rbind(AL, AL1,AL2,AL3, AL4,AL5, AL6, AL7,AL8,AL9,AL10, AL11, AL12, AL13,
                 AL14,AL15, AL16,AL17))

RSE=mean((test2$SD.B0.DW)-mean(((test2$IntDW-test2$IntPop)**2)**0.5))/mean(((test2$IntDW-test2$IntPop)**2)**0.5)*100
RSE1=mean((test2$SD.B1.DW)-mean(((test2$slopeDW-test2$slopePop)**2)**0.5))/mean(((test2$slopeDW-test2$slopePop)**2)**0.5)*100
RSE2=mean((test2$SD.BIMR.DW)-mean(((test2$IMRDW-test2$IMRPop)**2)**0.5))/mean(((test2$IMRDW-test2$IMRPop)**2)**0.5)*100


RSE3=mean((test2$SD.B0.PSGAM)-mean(((test2$IntPSGAM-test2$IntPop)**2)**0.5))/mean(((test2$IntPSGAM-test2$IntPop)**2)**0.5)*100
RSE4=mean((test2$SD.B1.PSGAM)-mean(((test2$slopePSGAM-test2$slopePop)**2)**0.5))/mean(((test2$slopePSGAM-test2$slopePop)**2)**0.5)*100
RSE5=mean((test2$SD.BIMR.PSGAM)-mean(((test2$IMRPSGAM-test2$IMRPop)**2)**0.5))/mean(((test2$IMRPSGAM-test2$IMRPop)**2)**0.5)*100

RSE6=mean((test2$SD.B0.UOPTGAM)-mean(((test2$IntUOPTGAM-test2$IntPop)**2)**0.5))/mean(((test2$IntUOPTGAM-test2$IntPop)**2)**0.5)*100
RSE7=mean((test2$SD.B1.UOPTGAM)-mean(((test2$slopeUOPTGAM-test2$slopePop)**2)**0.5))/mean(((test2$slopeUOPTGAM-test2$slopePop)**2)**0.5)*100
RSE8=mean((test2$SD.BIMR.UOPTGAM)-mean(((test2$IMRUOPTGAM-test2$IMRPop)**2)**0.5))/mean(((test2$IMRUOPTGAM-test2$IMRPop)**2)**0.5)*100

RSE9=mean((test2$SD.B0.SDWGAM)-mean(((test2$IntSDWGAM-test2$IntPop)**2)**0.5))/mean(((test2$IntSDWGAM-test2$IntPop)**2)**0.5)*100
RSE10=mean((test2$SD.B1.SDWGAM)-mean(((test2$slopeSDWGAM-test2$slopePop)**2)**0.5))/mean(((test2$slopeSDWGAM-test2$slopePop)**2)**0.5)*100
RSE11=mean((test2$SD.BIMR.SDWGAM)-mean(((test2$IMRSDWGAM-test2$IMRPop)**2)**0.5))/mean(((test2$IMRSDWGAM-test2$IMRPop)**2)**0.5)*100


RSE12=mean((test2$SD.B0.SPSGAM)-mean(((test2$IntSPSGAM-test2$IntPop)**2)**0.5))/mean(((test2$IntSPSGAM-test2$IntPop)**2)**0.5)*100
RSE13=mean((test2$SD.B1.SPSGAM)-mean(((test2$slopeSPSGAM-test2$slopePop)**2)**0.5))/mean(((test2$slopeSPSGAM-test2$slopePop)**2)**0.5)*100
RSE14=mean((test2$SD.BIMR.SPSGAM)-mean(((test2$IMRSPSGAM-test2$IMRPop)**2)**0.5))/mean(((test2$IMRSPSGAM-test2$IMRPop)**2)**0.5)*100

RSE15=mean((test2$SD.B0.SOPTGAM)-mean(((test2$IntSOPTGAM-test2$IntPop)**2)**0.5))/mean(((test2$IntSOPTGAM-test2$IntPop)**2)**0.5)*100
RSE16=mean((test2$SD.B1.SOPTGAM)-mean(((test2$slopeSOPTGAM-test2$slopePop)**2)**0.5))/mean(((test2$slopeSOPTGAM-test2$slopePop)**2)**0.5)*100
RSE17=mean((test2$SD.BIMR.SOPTGAM)-mean(((test2$IMRSOPTGAM-test2$IMRPop)**2)**0.5))/mean(((test2$IMRSOPTGAM-test2$IMRPop)**2)**0.5)*100


RSE.overall=data.frame(rbind(RSE, RSE1,RSE2,RSE3, RSE4,RSE5,RSE6, RSE7,RSE8,RSE9,RSE10, RSE11,
                  RSE12, RSE13,RSE14,RSE15, RSE16,RSE17))


Table1=cbind(RSE.overall, cov.overall,AL.overall)
rownames(Table1)=c( "IntDW","slopeDW", "IMRDW", 
                   "IntPSGAM","slopePSGAM", "IMRPSGAM",
                   "IntUOPTGAM","slopeUOPTGAM", "IMRUOPTGAM", 
                   "IntSDWGAM","slopeSDWGAM", "IMRSDWGAM",
                   "IntSPSGAM","slopeSPSGAM", "IMRSPSGAM", 
                   "IntSOPTGAM","slopeSOPTGAM", "IMRSOPTGAM")
colnames(Table1)=c( "RSE(%)","Coverage", "Avg Length")






















#############################################
# Getting Bias Avg-Relative Bias SE and RMSE#
#############################################

fun2 <- function(v1, v2, data) {
   #[[]] means to extract the contents (when specified as a string)
  bias <- mean(data[[v2]] - data[[v1]])
  arb <- mean((data[[v2]] - data[[v1]]) / data[[v1]])
  rmse <- sqrt(mean((data[[v2]]-data[[v1]])^2)+(sd(data[[v2]]^2)))
  se <- sd(data[[v2]])
  return(data.frame(Bias=bias, SE=se,  RMSE=rmse,  ARelBias=arb))
}

fun3 <- function(wttype, data) {
  params <- c("Int", "slope", "IMR")
  var.names.pop <- paste0(params, "Pop")
  var.names.samp <- paste0(params, wttype)
  for (i in 1:3) {
    print(var.names.pop[i])
    print(var.names.samp[i])
    temp <- cbind(data.frame(weight=wttype, estimate=params[i]),
                  fun2(var.names.pop[i], var.names.samp[i], data))
    
    if (i == 1) {
      result <- temp
    } else {
      result <- rbind(result, temp)
    }
  }
  return(result)
}

rbind(fun3("DW", test1),
      fun3("PSGAM", test1),
      fun3("UOPTGAM", test1),
      fun3("SDWGAM", test1),
      fun3("SPSGAM", test1),
      fun3("SOPTGAM", test1))


# Getting the SE-Rel Bias, Coverage, Average Length#

#fun4 <- function(v1, v2,v3, data) {
#  # [[]] means to extract the contents (when specified as a string)
#  SERelBias <- mean((data[[v1]])**(0.5) - mean((data[[v2]]-data[[v3]]**2)**0.5))
#  Averagelength <- mean((data[[v2]] - data[[v1]]))
#  Coverage <- mean(ifelse(data[[v3]] >data[[v1]] & data[[v3]] >data[[v2]],1,0))
#  return(data.frame(SERelBias=SERelBias, Coverage=Coverage,  Averagelength=Averagelength))
#}

check.est <- function(dat, par.name, wt.type) {
  # Construct dynamic variable names based on parameter of interest and weight method.
  # See naming conventions above.
  par.name2 <- ifelse(par.name == "Int", "B0",
                      ifelse(par.name == "slope", "B1", "BIMR"))
  lb.name <- paste0("LB.", par.name2, ".", wt.type)
  ub.name <- paste0("UB.", par.name2, ".", wt.type)
  pop.name <- paste0(par.name, "Pop")
  sd.name <- paste0("SD.", par.name2, ".", wt.type)
  est.name <- paste0(par.name, wt.type)
  
  # Then, extract results of interest.
  # Store in shorted variable names for ease of readability
  p <- dat[[pop.name]]
  ub <- dat[[ub.name]]
  lb <- dat[[lb.name]]
  sd.est <- dat[[sd.name]]
  est <- dat[[est.name]]
  
  # Set up results list.
  result <- list()
  
  # Compute diagnostics.
  result$coverage <- mean(p >= lb & p <= ub)
  result$averagelength <- mean(ub - lb)
  bs <- mean(sqrt(((est - p)^2)))
  result$serelbias <- (mean(sd.est) - bs) / bs
  result$bias <- mean(est - p)
  result$SE <- sd(est)
  result$rmse <- sqrt(mean((est - p)^2))
  result$arb <- mean(est - p) / mean(p)
  
  # Done!
  return(result)
}

# Set up data frame for each combination of parameter and weight method.
est.df <- expand.grid(parameter=c("Int", "slope", "IMR"),
                      method=c("DW", "PSGAM", "UOPTGAM", "SDWGAM", "SPSGAM", "SOPTGAM"))
est.df$bias <-est.df$SE<- est.df$rmse <- est.df$arb<-
  est.df$cov <- est.df$averagelength <- est.df$serelbias <- NA

# Iterate over parameters and weight methods.
for (r in 1:nrow(est.df)) {
  # Get diagnostics for each combination.
  temp <- check.est(test2, par.name=est.df$parameter[r],
                    wt.type=est.df$method[r])
  
  # Store in data frame in appropriate columns.
  est.df$bias[r]<- temp$bias
  est.df$SE[r]<- temp$SE
  est.df$rmse[r]<- temp$rmse
  est.df$arb[r]<- temp$arb
  est.df$cov[r] <- temp$coverage
  est.df$averagelength[r] <- temp$averagelength
  est.df$serelbias[r] <- temp$serelbias
}

print(est.df)

###############################################
# Data Checking #
###############################################


################################################
# check results#
#################################################


