## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(EpiEstim, tidyverse, incidence, plotrix, distr,
       surveillance,mgcv,bbmle,BayesianTools,grDevices)

## Get data together and estimate R0
covid.data <- read_csv("./COVID_Surveillance_Impact_20200810-20201204.csv") %>%
    select(CRU_COMPLETED_DATE,Test_Type, Description, Positive) 
dates_onset <- rep(covid.data$CRU_COMPLETED_DATE, covid.data$Positive)
location.temp <- ifelse(covid.data$Description == "student","local","imported")
location <- rep(location.temp, covid.data$Positive)
location[1] <- "imported"
imported.prop <- sum(location == "imported")/length(location)
change.indices <- sample((1:length(location))[location == "local"],
                         imported.prop*sum(location == "local"))
location[change.indices] <- "imported"
incid <- incidence(dates_onset, groups = location)
res_with_imports <- estimate_R(incid, method = "parametric_si",
                   config = make_config(list(
                       mean_si = 5.11, std_si = 2.68)))##Zhang


### Plan for analyzing wastewater data
## 1. Generate an incidence curve from the Rt estimate, imports, etc.,
## 2. use a parametric distribution (gamma) of individual shedding, to estimate shedding over time
## 3. Calibrate the gamma parameters

## 1. set up params
inc.shape <- 1.88
inc.scale <- 7.97
delay.mean <- 1
delay.rate <- 1/delay.mean

## deconvolve using  surveillance package
dmax <- 50 # truncation of incubation plus pmf 
dates <- res_with_imports$dates
cases <- c(rep(0,dmax),res_with_imports$I)## case detection time series
names(cases) <- c(seq(min(dates)-dmax,
                      min(dates)-1,length.out=dmax),dates) #seq(dmax,by=1,length.out=length(cases))
infections <- c(rep(0,dmax-delay.mean),res_with_imports$I,rep(NA,delay.mean))
names(infections) <- names(cases)
dates <- as.Date(names(cases))
sts <- new("sts", epoch=1:length(cases),observed=matrix(cases,ncol=1))
## functions to convolve
Inc <- Weibull(shape=inc.shape, scale=inc.scale)
Delay <- Pois(delay.mean) ## Try a shift, and a poisson
## Delay <- Gammad(shape=delay.mean/1e6,scale=1e6)
p.inf2test <- p(convpow(Inc+Delay, 1))
d.inf2test <- d(convpow(Inc+Delay, 1))
## convert pdf to pmf
inc.pmf <- c(0,(p.inf2test(1:dmax) - p.inf2test(0:(dmax-1)))/p.inf2test(dmax))

#Call non-parametric back-projection function with hook function but
#without bootstrapped confidence intervals
plotIt <- function(cur.sts) {
  plot(cur.sts,xaxis.labelFormat=NULL, legend=NULL,ylim=c(0,350))
}
##bpnp.control <- list(k=2,eps=rep(0.005,2),iter.max=rep(250,2),B=-1)
bpnp.control <- list(k=10,eps=rep(1e-5,2),iter.max=rep(250,2),B=-1,
                     ## hookFun=plotIt,verbose=TRUE)
                     verbose=TRUE)

#Fast C version (use argument: eq3a.method="C")! 
sts.bp <- backprojNP(sts, incu.pmf=inc.pmf,
    control=modifyList(bpnp.control,list(eq3a.method="C")))

plot(sts.bp,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),main="")
#Do the convolution for the expectation
mu <- matrix(0,ncol=ncol(sts.bp),nrow=nrow(sts.bp))
#Loop over all series
for (j in 1:ncol(sts.bp)) { 
  #Loop over all time points
  for (t in 1:nrow(sts.bp)) {
    #Convolution, note support of inc.pmf starts at zero (move idx by 1)
    i <- seq_len(t)
    mu[t,j] <- sum(inc.pmf[t-i+1] * upperbound(sts.bp)[i,j],na.rm=TRUE)
  }
}
#Show the fit
lines(1:nrow(sts.bp)-0.5,mu[,1],col="green",lwd=3,lty="dashed")

## get data
ww.data <- read_csv("ND WW Data_11-30-2020.csv") %>%
    rename(date=X1) %>%
    mutate(date=as.Date(date,"%d-%b-%Y")) %>%
    filter(trust)
N1.med <- apply(ww.data[,2:4],1,median)
N1.min <- apply(ww.data[,2:4],1,min)
N1.max <- apply(ww.data[,2:4],1,max)

## 2. 3. ML fitting with optim, with profiled likelihood on the fitting parameter
n.exp <- upperbound(sts.bp)
date.exp <- dates
NLL.dispersion <- function(par) {
    shed.shape <- (exp(par[1])/exp(par[2]))^2
    shed.rate <- shed.shape/exp(par[1])
    shedding <- rep(0,length(n.exp))
    for (jj in 1:(length(n.exp)-1)) {
        shedding <- shedding + n.exp[jj]*c(rep(0,jj),dgamma(1:(length(n.exp)-jj),
                                                            shed.shape,shed.rate))
    }
    shedding.scaled <- shedding*exp(par[3])
    nbinom.size <- exp(par[4])
    if (sum(shedding.scaled)) {
        return(-sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,2]),size=nbinom.size,
                            mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
               sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,3]),size=nbinom.size,
                           mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
               sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,4]),size=nbinom.size,
                           mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)))
    } else {
        return(Inf)
    }
}
LL.dispersion <- function(par) {-NLL.dispersion(par)}
NLL.mle2.dispersion <- function(logmean=1,logsd=1,logshedscale=1,lognbinomsize=1) {
    return(NLL(c(logmean,logsd,logshedscale,lognbinomsize)))
}

optim.out <- optim(c(log(1),log(1),log(1)),NLL.dispersion)
mle2.out <- mle2(NLL.mle2.dispersion)
mle2.profile <- profile(mle2.out)
save(mle2.profile,mle2.out,file=paste0("mle2_delay_pois",floor(1/delay.rate+0.5),".RData"))
load(paste0("mle2_delay_pois",floor(1/delay.rate+0.5),".RData"),verbose=TRUE)

## MCMC - with dispersion
lower <- c(logmean=log(1e0),logsd=log(1e0),logshedscale=log(1e-1), lognbinomsize=log(1e-3))#confint(mle2.profile)[,1]*0.1
upper <- c(logmean=log(1e5),logsd=log(1e7),logshedscale=log(1e3), lognbinomsize=log(1e1))#confint(mle2.profile)[,2]*1.9
bayesianSetup = createBayesianSetup(LL.dispersion,lower=lower,upper=upper)
mcmc.out = runMCMC(bayesianSetup)
plot(mcmc.out,start = 500)
samples = getSample(mcmc.out,start=500)
save(mcmc.out,samples,file=paste0("mcmc_delay_pois",floor(1/delay.rate+0.5),".RData"))
load(paste0("mcmc_delay_pois",floor(1/delay.rate+0.5),".RData"))

DIC(mcmc.out)$DIC

marginalPlot(mcmc.out, prior = TRUE, start=500)

gelmanDiagnostics(mcmc.out, plot = T)

## 4. plot it
pars <- exp(samples)
shed.mean <- pars[,1]
shed.sd <- pars[,2]
shed.shape <- (shed.mean/shed.sd)^2
shed.rate <- shed.shape/shed.mean
shed.scale <- pars[,3]
maxplot.day <- 30
nbinom.size <- pars[,4]
shedding <- array(0,dim=c(length(shed.shape),length(n.exp)))
shed.dist <- array(NA,dim=c(length(shed.shape),maxplot.day+1))
for (jj in 1:length(shed.shape)) {
    for (ii in 1:(length(n.exp)-1)) {
        shedding[jj,] <- shedding[jj,] + shed.scale[jj]*n.exp[ii]*c(rep(0,ii),
                                                                    dgamma(1:(length(n.exp)-ii),
                                                                           shed.shape[jj],
                                                                           shed.rate[jj]))
    }
    shed.dist[jj,] <- dgamma(0:maxplot.day,shed.shape[jj],shed.rate[jj])*shed.scale[jj]
}

## Prediction interval
n.samples <- 1e6
sampled.rows <- sample(nrow(shedding),n.samples,replace=TRUE)
sampled.mu <- shedding[sampled.rows,]
pred <- matrix(rnbinom(n.samples*ncol(sampled.mu),
                       size=nbinom.size[sampled.rows],
                       mu=t(sampled.mu)),
               nrow=nrow(sampled.mu),
               ncol=ncol(sampled.mu),byrow=TRUE)

## plot it
pdf(paste0("shedding_distribution_pois_",floor(1/delay.rate+0.5),"_mean_delay.pdf"))
par(mfrow=c(2,1))
plot(0:30, colMeans(shed.dist),
     type='l', ylim=c(0,max(apply(shed.dist[,-1],2,function(x)quantile(x,0.975)))),
     ylab="Shedding intensity (N1 GC / 100ml)",xlab="Time since infection",yaxs="i",las=1,bty="n")
polygon(c(1:30,30:1),c(apply(shed.dist[,-1],2,function(x)quantile(x,0.025)),
                       rev(apply(shed.dist[,-1],2,function(x)quantile(x,0.975)))),
        col=adjustcolor("gray",0.5),border=F)
plot(date.exp,colMeans(pred),type='l',xlim=c(min(ww.data$date),max(ww.data$date)),
     ylim=c(0,max(ww.data[,2:4], apply(pred,2,function(x)quantile(x,0.975)))),
     ylab="RNA (N1 GC / 100ml)",xlab="Date",las=1,bty="n")
plotCI(ww.data$date,N1.med,li=N1.min,ui=N1.max,add=TRUE,gap=TRUE,sfrac=0.003)
polygon(c(date.exp,rev(date.exp)),c(apply(pred,2,function(x)quantile(x,0.025)),
                                    rev(apply(pred,2,function(x)quantile(x,0.975)))),
        col=adjustcolor("gray",0.5),border=F)
dev.off()

pdf(paste0("shedding_distribution_",floor(1/delay.rate+0.5),"_mean_delay_spaghettid.pdf"))
par(mfrow=c(2,1))
plot(0:30, colMeans(shed.dist),
     type='l', ylim=c(0,0.7),
     ylab="Shedding intensity (N1 GC / 100ml)",xlab="Time since infection",yaxs="i",las=1,bty="n")
polygon(c(1:30,30:1),c(apply(shed.dist[,-1],2,function(x)quantile(x,0.025)),
                       rev(apply(shed.dist[,-1],2,function(x)quantile(x,0.975)))),
        col=adjustcolor("gray",0.5),border=F)
plot(-1,-1,type='l',xlim=c(min(ww.data$date),max(ww.data$date)),
     ylim=c(0,max(ww.data[,2:4])),
     ylab="RNA (N1 GC / 100ml)",xlab="Date",las=1,bty="n")
for (ii in sample(nrow(pred),100)) {
    lines(date.exp,pred[ii,],col=adjustcolor("gray",alpha.f=0.5))
}
lines(date.exp,colMeans(pred),col="red",lwd=2)
plotCI(ww.data$date,N1.med,li=N1.min,ui=N1.max,add=TRUE,gap=TRUE,sfrac=0.003)
dev.off()


## ## Boneyard
## ## Method B for 2., variable shedding period, with constant shedding during that period
## time.exp <- rep(1:length(cases),cases)
## inc.shape <- 1.88
## inc.scale <- 7.97
## shed.mean <- 10
## shed.sd <- 5
## shed.shape <- (shed.mean/shed.sd)^2
## shed.rate <- shed.shape/shed.mean
## time.start <- time.exp ##+ rweibull(length(time.exp), shape=inc.shape, scale=inc.scale)
## time.end <- time.start + rgamma(length(time.start), shed.shape, shed.rate)

## shedding <- rep(0,length(cases))
## for (ii in 1:length(shedding)) {
##     shedding[ii] <- sum((time.start <= ii) & (time.end > ii))
## }

## shedding.scaled <- shedding*(max(N1.med) - min(N1.med))/(max(shedding) - min(shedding)) +
##     (max(N1.med)*min(shedding) - min(N1.med)*max(shedding))/(max(shedding) - min(shedding))
## par(mfrow=c(2,1))
## plot(dgamma(0:30,shed.shape,shed.rate),type='l',
##      ylab="Shedding length",xlab="Time since infection",yaxs="i",las=1,bty="n")
## plot(dates,shedding.scaled,type='l',xlim=c(min(ww.data$date),max(ww.data$date)),
##      ylim=c(0,max(ww.data[,2:4])),
##      ylab="Expected measurement",xlab="Date",las=1,bty="n")
## points(ww.data$date,N1.med)

## ## functions to convolve
## Inc <- Weibull(shape=inc.shape, scale=inc.scale)
## Delay <- Exp(delay.rate)
## d.inf2test <- d(convpow(Inc+Delay, 1))
## r.inf2test <- r(convpow(Inc+Delay, 1))
## cases <- res_with_imports$I
## dates <- res_with_imports$dates
## earliest.date <- 30
## date.exp <- seq(min(dates) - earliest.date,max(dates),"1 day")
## n.exp <- rep(0,length(date.exp))

## for (ii in 1:length(dates)) {
##     for (jj in 0:earliest.date) {
##         n.exp[ii+jj] = n.exp[ii+jj]  + cases[ii]*d.inf2test(earliest.date-jj)
##     }
## }

## plot(date.exp,n.exp,ylim=c(0,max(cases)))
## lines(dates,cases)

## ## 2. and 3. Calculate likelihood for a range of values
## ## Use upperbound(sts.bp) to get deconvoluted timeseries
## nbinom.size <- 10^seq(-4,4,1)
## shed.means <- 10^seq(0,1,0.1)
## shed.sds <- 10^seq(0,1,0.1)
## shed.scaling <- seq(0.1,1,0.1)
## par <- expand.grid(shed.mean=shed.means,shed.sd=shed.sds,nbinom.size=nbinom.size,
##                    shed.scaling=shed.scaling)
## LL <- rep(0,nrow(par))
## n.exp <- upperbound(sts.bp)
## date.exp <- dates
## for (ii in 1:nrow(par)) {
##     ## print(ii)
##     shed.shape <- (par[ii,1]/par[ii,2])^2
##     shed.rate <- shed.shape/par[ii,1]
##     shedding <- rep(0,length(n.exp))
##     for (jj in 1:(length(n.exp)-1)) {
##         shedding <- shedding + n.exp[jj]*c(rep(0,jj),dgamma(1:(length(n.exp)-jj),
##                                                             shed.shape,shed.rate))
##     }
##     shedding.scaled <- shedding*max(N1.med)*par[ii,4]/max(shedding)
##     LL[ii] <- sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,2]),size=par[ii,3],
##                           mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE))
##     LL[ii] <- LL[ii] + sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,3]),size=par[ii,3],
##                                    mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE))
##     LL[ii] <- LL[ii] + sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,4]),size=par[ii,3],
##                                    mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE))
## }

## ## 4. plot it
## shed.mean <- par[which.max(LL),1]
## shed.sd <- par[which.max(LL),2]
## shed.shape <- (shed.mean/shed.sd)^2
## shed.rate <- shed.shape/shed.mean
## shedding <- rep(0,length(n.exp))
## for (ii in 1:(length(n.exp)-1)) {
##     shedding <- shedding + n.exp[ii]*c(rep(0,ii),dgamma(1:(length(n.exp)-ii),shed.shape,shed.rate))
## }
## shedding.scaled <- shedding*max(N1.med)*par[which.max(LL),4]/max(shedding)
## pdf(paste0("shedding_distribution_",floor(1/delay.rate+0.5),"_mean_delay.pdf"))
## par(mfrow=c(2,1))
## plot(0:30, dgamma(0:30,shed.shape,shed.rate)*max(N1.med)*par[which.max(LL),4]/max(shedding),
##      type='l',
##      ylab="Shedding intensity",xlab="Time since infection",yaxs="i",las=1,bty="n")
## plot(date.exp,shedding.scaled,type='l',xlim=c(min(ww.data$date),max(ww.data$date)),
##      ylim=c(0,max(ww.data[,2:4])),
##      ylab="Expected measurement",xlab="Date",las=1,bty="n")
## plotCI(ww.data$date,N1.med,li=N1.min,ui=N1.max,add=TRUE,gap=TRUE,sfrac=0.003)
## dev.off()



## ### Boneyard
## covid.data <- read_csv("./COVID_Surveillance_Impact_20200810-20201204.csv") %>%
##     select(CRU_COMPLETED_DATE, Description, Positive)

## dates_onset <- rep(covid.data$CRU_COMPLETED_DATE, covid.data$Positive)
## location.temp <- ifelse(covid.data$Description == "student","local","imported")
## location <- rep(location.temp, covid.data$Positive)
## location[1] <- "imported"
## imported.prop <- sum(location == "imported")/length(location)
## change.indices <- sample((1:length(location))[location == "local"],
##                          imported.prop*sum(location == "local"))
## location[change.indices] <- "imported"

## incid <- incidence(dates_onset, groups = location)
## plot(incid)

## res_with_imports <- estimate_R(incid, method = "parametric_si",
##                    config = make_config(list(
##                        mean_si = 5.11, std_si = 2.68)))##Zhang

## ## pdf("./Rt_plot.pdf",7,7)
## plot(res_with_imports, add_imported_cases=TRUE)
## ## dev.off()

## ## ## epiestim vignette
## ## data(Flu2009)

## ## dates_onset <- Flu2009$incidence$dates[unlist(lapply(1:nrow(Flu2009$incidence), function(i) 
## ##     rep(i, Flu2009$incidence$I[i])))]

## ## location <- sample(c("local","imported"), length(dates_onset), replace=TRUE)
## ## location[1] <- "imported" # forcing the first case to be imported

## ## ## get incidence per group (location)
## ## incid <- incidence(dates_onset, groups = location)

## ## plot(incid)

## ## ## Estimate R with assumptions on serial interval
## ## res_with_imports <- estimate_R(incid, method = "parametric_si",
## ##                    config = make_config(list(
## ##                        mean_si = 2.6, std_si = 1.5)))

## ## plot(res_with_imports, add_imported_cases=TRUE)

## ## Poisson likelihood
## NLL.pois <- function(par) {
##     shed.shape <- (exp(par[1])/exp(par[2]))^2
##     shed.rate <- shed.shape/exp(par[1])
##     shedding <- rep(0,length(n.exp))
##     for (jj in 1:(length(n.exp)-1)) {
##         shedding <- shedding + n.exp[jj]*c(rep(0,jj),dgamma(1:(length(n.exp)-jj),
##                                                             shed.shape,shed.rate))
##     }
##     shedding.scaled <- shedding*exp(par[3])
##     LL <- sum(dpois(unlist(ww.data[ww.data$date %in% date.exp,2]),
##                     lambda=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) +
##         sum(dpois(unlist(ww.data[ww.data$date %in% date.exp,3]),
##                   lambda=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) +
##         sum(dpois(unlist(ww.data[ww.data$date %in% date.exp,4]),
##                   lambda=shedding.scaled[date.exp %in% ww.data$date],log=TRUE))
##     return(-LL)
## }

## ## Scatter plot
## plot(unlist(ww.data[ww.data$date %in% date.exp,2]),
##      colMeans(shedding.scaled)[date.exp %in% ww.data$date])
## points(unlist(ww.data[ww.data$date %in% date.exp,3]),
##      colMeans(shedding.scaled)[date.exp %in% ww.data$date])
## points(unlist(ww.data[ww.data$date %in% date.exp,4]),
##        colMeans(shedding.scaled)[date.exp %in% ww.data$date])

## ## Correlation
## cor.vec <- vector(mode="numeric",length=nrow(shedding.scaled))
## cor.data <- c(unlist(ww.data[ww.data$date %in% date.exp,2]),
##               unlist(ww.data[ww.data$date %in% date.exp,3]),
##               unlist(ww.data[ww.data$date %in% date.exp,4]))
## cor.output <- cbind(shedding.scaled[,date.exp %in% ww.data$date],
##                     shedding.scaled[,date.exp %in% ww.data$date],
##                     shedding.scaled[,date.exp %in% ww.data$date])
## correlation <- cor(t(cor.output),cor.data)[,1]
## summary(correlation)


## ### Try fiting an exponential instead
## NLL.exp <- function(par) {
##     shed.rate <- 1/exp(par[1])
##     shedding <- rep(0,length(n.exp))
##     for (jj in 1:(length(n.exp)-1)) {
##         shedding <- shedding + n.exp[jj]*c(rep(0,jj),dexp(1:(length(n.exp)-jj),
##                                                           shed.rate))
##     }
##     shedding.scaled <- shedding*exp(par[2])
##     optim.out <- optim(1,function(x) {
##         return(-sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,2]),size=x,
##                             mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
##                sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,3]),size=x,
##                            mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
##                sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,4]),size=x,
##                            mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)))},
##         method="Brent",lower=0,upper=1e6)
##     return(optim.out$value)
## }
## LL.exp <- function(par) {-NLL.exp(par)}
## NLL.exp.mle2 <- function(logmean=1,logshedscale=1) {return(NLL.exp(c(logmean,logshedscale)))}

## ## optim.out <- optim(c(log(1),log(1),log(1)),NLL)
## mle2.exp.out <- mle2(NLL.exp.mle2,trace=TRUE)
## mle2.exp.profile <- profile(mle2.exp.out)
## save(mle2.exp.profile,mle2.exp.out,file=paste0("mle2_exp_delay",floor(1/delay.rate+0.5),".RData"))
## load(paste0("mle2_exp_delay",floor(1/delay.rate+0.5),".RData"),verbose=TRUE)

## ## MCMC
## lower <- confint(mle2.exp.profile)[,1]*0.9
## upper <- confint(mle2.exp.profile)[,2]*1.1
## bayesianSetup = createBayesianSetup(LL.exp,lower=lower,upper=upper)
## mcmc.exp.out = runMCMC(bayesianSetup)
## samples.exp = getSample(mcmc.exp.out)
## plot(mcmc.exp.out)
## DIC(mcmc.exp.out)$DIC
## save(mcmc.exp.out,samples.exp,file=paste0("mcmc_exp_delay",floor(1/delay.rate+0.5),".RData"))
## load(paste0("mcmc_exp_delay",floor(1/delay.rate+0.5),".RData"))

## ## 4. plot it
## pars <- exp(samples.exp)
## shed.mean <- pars[,1]
## shed.rate <- 1/shed.mean
## shed.scale <- pars[,2]
## maxplot.day <- 30
## shedding <- array(0,dim=c(length(shed.shape),length(n.exp)))
## shed.dist <- array(NA,dim=c(length(shed.shape),maxplot.day+1))
## for (jj in 1:length(shed.shape)) {
##     for (ii in 1:(length(n.exp)-1)) {
##         shedding[jj,] <- shedding[jj,] + n.exp[ii]*c(rep(0,ii),
##                                                      dexp(1:(length(n.exp)-ii),shed.rate[jj]))
##     }
##     shed.dist[jj,] <- dexp(0:maxplot.day,shed.rate[jj])*shed.scale[jj]
## }
## shedding.scaled <- shedding*shed.scale
## pdf(paste0("shedding_distribution_exp_",floor(1/delay.rate+0.5),"_mean_delay.pdf"))
## par(mfrow=c(2,1))
## plot(0:30, colMeans(shed.dist),
##      type='l', ylim=c(0,0.7),
##      ylab="Shedding intensity (N1 GC / 100ml)",xlab="Time since infection",yaxs="i",las=1,bty="n")
## polygon(c(1:30,30:1),c(apply(shed.dist[,-1],2,function(x)quantile(x,0.025)),
##                        rev(apply(shed.dist[,-1],2,function(x)quantile(x,0.975)))),
##         col=adjustcolor("gray",0.5),border=F)
## plot(date.exp,colMeans(shedding.scaled),type='l',xlim=c(min(ww.data$date),max(ww.data$date)),
##      ylim=c(0,max(ww.data[,2:4], shedding.scaled)),
##      ylab="RNA (N1 GC / 100ml)",xlab="Date",las=1,bty="n")
## plotCI(ww.data$date,N1.med,li=N1.min,ui=N1.max,add=TRUE,gap=TRUE,sfrac=0.003)
## polygon(c(date.exp,rev(date.exp)),c(apply(shedding.scaled,2,function(x)quantile(x,0.025)),
##                        rev(apply(shedding.scaled,2,function(x)quantile(x,0.975)))),
##         col=adjustcolor("gray",0.5),border=F)
## dev.off()

## NLL <- function(par) {
##     shed.shape <- (exp(par[1])/exp(par[2]))^2
##     shed.rate <- shed.shape/exp(par[1])
##     shedding <- rep(0,length(n.exp))
##     for (jj in 1:(length(n.exp)-1)) {
##         shedding <- shedding + n.exp[jj]*c(rep(0,jj),dgamma(1:(length(n.exp)-jj),
##                                                             shed.shape,shed.rate))
##     }
##     shedding.scaled <- shedding*exp(par[3])
##     if (sum(shedding.scaled)) {
##         optim.out <- optim(1,function(x) {
##             return(-sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,2]),size=x,
##                                 mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
##                    sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,3]),size=x,
##                                mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
##                    sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,4]),size=x,
##                                mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)))},
##             method="Brent",lower=0,upper=1e6)
##         return(optim.out$value)
##     } else {
##         return(Inf)
##     }
## }
## LL <- function(par) {-NLL(par)}
## NLL.mle2 <- function(logmean=1,logsd=1,logshedscale=1) {return(NLL(c(logmean,logsd,logshedscale)))}
