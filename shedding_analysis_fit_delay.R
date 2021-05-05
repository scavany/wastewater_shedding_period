## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(EpiEstim, tidyverse, incidence, plotrix, distr,
       surveillance,mgcv,bbmle,BayesianTools,grDevices)

## Get data together and estimate R0
covid.data <- read_csv("./COVID_Surveillance_Impact_20200810-20201204.csv") %>%
    select(CRU_COMPLETED_DATE, Description, Positive)
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

## deconvolve using  surveillance package
dmax <- 50 # truncation of incubation plus pmf 
dates <- res_with_imports$dates
cases <- c(rep(0,dmax),res_with_imports$I)## case detection time series
names(cases) <- c(seq(min(dates)-dmax,
                      min(dates)-1,length.out=dmax),dates) #seq(dmax,by=1,length.out=length(cases))
dates <- as.Date(names(cases))
sts <- new("sts", epoch=1:length(cases),observed=matrix(cases,ncol=1))
## functions to convolve
Inc <- Weibull(shape=inc.shape, scale=inc.scale)

#Call non-parametric back-projection function with hook function but
#without bootstrapped confidence intervals
bpnp.control <- list(k=2,eps=rep(0.005,2),iter.max=rep(250,2),B=-1)

## get data
ww.data <- read_csv("ND WW Data_11-30-2020.csv") %>%
    rename(date=X1) %>%
    mutate(date=as.Date(date,"%d-%b-%Y")) %>%
    filter(trust)
N1.med <- apply(ww.data[,2:4],1,median)
N1.min <- apply(ww.data[,2:4],1,min)
N1.max <- apply(ww.data[,2:4],1,max)

## 2. 3. ML fitting with optim, with profiled likelihood on the fitting parameter
date.exp <- dates
NLL <- function(par) {
    delay.rate <- 1/exp(par[4])
    Delay <- Exp(delay.rate)
    p.inf2test <- p(convpow(Inc+Delay, 1))
    d.inf2test <- d(convpow(Inc+Delay, 1))
    inc.pmf <- c(0,(p.inf2test(1:dmax) - p.inf2test(0:(dmax-1)))/p.inf2test(dmax))
    sts.bp <- backprojNP(sts, incu.pmf=inc.pmf,
                         control=modifyList(bpnp.control,list(eq3a.method="C")))
    n.exp <- upperbound(sts.bp)
    shed.shape <- (exp(par[1])/exp(par[2]))^2
    shed.rate <- shed.shape/exp(par[1])
    shedding <- rep(0,length(n.exp))
    for (jj in 1:(length(n.exp)-1)) {
        shedding <- shedding + n.exp[jj]*c(rep(0,jj),dgamma(1:(length(n.exp)-jj),
                                                            shed.shape,shed.rate))
    }
    shedding.scaled <- shedding*exp(par[3])
    optim.out <- optim(1,function(x) {
        return(-sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,2]),size=x,
                            mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
               sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,3]),size=x,
                           mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)) -
               sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,4]),size=x,
                           mu=shedding.scaled[date.exp %in% ww.data$date],log=TRUE)))},
        method="Brent",lower=0,upper=1e6)
    return(optim.out$value)
}
LL <- function(par) {-NLL(par)}
NLL.mle2 <- function(logmean=1,logsd=1,logshedscale=1) {return(NLL(c(logmean,logsd,logshedscale)))}

## MCMC
## lower <- c(logmean=log(1),logsd=log(1),logshedscale=log(1),logdelayrate=log(1/24))
## upper <- c(logmean=log(365),logsd=log(365),logshedscale=log(100),logdelaymean=log(2*7))
## bayesianSetup = createBayesianSetup(LL,lower=lower,upper=upper)
## mcmc.out = runMCMC(bayesianSetup)
## samples = getSample(mcmc.out)
## save(mcmc.out,samples,file="mcmc_fitted_delay.RData")
load("mcmc_fitted_delay.RData")
plot(mcmc.out)

DIC(mcmc.out)$DIC

marginalPlot(mcmc.out, prior = TRUE)

gelmanDiagnostics(mcmc.out, plot = T)

## 4. plot it
pars <- exp(samples)
delay.rate <- 1/exp(pars[,4])
shed.mean <- pars[,1]
shed.sd <- pars[,2]
shed.shape <- (shed.mean/shed.sd)^2
shed.rate <- shed.shape/shed.mean
shed.scale <- pars[,3]
maxplot.day <- 30
n.exp.arr <- array(NA,dim=c(length(shed.shape),nrow(sts)))
shedding <- array(0,dim=c(length(shed.shape),nrow(sts)))
shed.dist <- array(NA,dim=c(length(shed.shape),maxplot.day+1))
nbinom.size <- vector(mode="numeric",length=length(shed.shape))
for (jj in 1:length(shed.shape)) {
    if (jj%%100 == 0) print(jj)
    Delay <- Exp(delay.rate[jj])
    p.inf2test <- p(convpow(Inc+Delay,1))
    d.inf2test <- d(convpow(Inc+Delay,1))
    inc.pmf <- c(0,(p.inf2test(1:dmax) - p.inf2test(0:(dmax-1)))/p.inf2test(dmax))
    sts.bp <- backprojNP(sts, incu.pmf=inc.pmf,
                         control=modifyList(bpnp.control,list(eq3a.method="C")))
    n.exp.arr[jj,] <- upperbound(sts.bp)
    for (ii in 1:(ncol(n.exp.arr)-1)) {
        shedding[jj,] <- shedding[jj,] + shed.scale[jj]*n.exp.arr[jj,ii]*c(rep(0,ii),
                                                                           dgamma(1:(ncol(n.exp.arr)-ii),
                                                                                  shed.shape[jj],
                                                                                  shed.rate[jj]))
    }
    shed.dist[jj,] <- dgamma(0:maxplot.day,shed.shape[jj],shed.rate[jj])*shed.scale[jj]
    nbinom.size[jj] <- (optim(1,function(x) {
        return(-sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,2]),size=x,
                            mu=shedding[jj,date.exp %in% ww.data$date],log=TRUE)) -
               sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,3]),size=x,
                           mu=shedding[jj,date.exp %in% ww.data$date],log=TRUE)) -
               sum(dnbinom(unlist(ww.data[ww.data$date %in% date.exp,4]),size=x,
                           mu=shedding[jj,date.exp %in% ww.data$date],log=TRUE)))},
        method="Brent",lower=0,upper=1e6))$par
}
save(shedding,shed.dist,nbinom.size,n.exp.arr,file="prediction_interval_objects.RData")
load("prediction_interval_objects.RData",verbose=TRUE)
## shedding.scaled <- shedding*shed.scale
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
pdf("shedding_distribution_fitted_delay.pdf")
par(mfrow=c(2,1))
plot(0:30, colMeans(shed.dist),
     type='l', ylim=c(0,0.7),
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

pdf("deconvoluted_cases.pdf")
barplot(cases)
lines(colMeans(n.exp.arr),col="red",lwd=2)
polygon(c(1:length(cases),length(cases):1),
        c(apply(n.exp.arr,2,function(x)quantile(x,0.025)),
          rev(apply(n.exp.arr,2,function(x)quantile(x,0.975)))),
        col=adjustcolor("red",0.3),border=FALSE)
dev.off()

## Scatter plot
plot(unlist(ww.data[ww.data$date %in% date.exp,2]),
     colMeans(shedding.scaled)[date.exp %in% ww.data$date])
points(unlist(ww.data[ww.data$date %in% date.exp,3]),
     colMeans(shedding.scaled)[date.exp %in% ww.data$date])
points(unlist(ww.data[ww.data$date %in% date.exp,4]),
       colMeans(shedding.scaled)[date.exp %in% ww.data$date])

## Correlation
cor.vec <- vector(mode="numeric",length=nrow(shedding.scaled))
cor.data <- c(unlist(ww.data[ww.data$date %in% date.exp,2]),
              unlist(ww.data[ww.data$date %in% date.exp,3]),
              unlist(ww.data[ww.data$date %in% date.exp,4]))
cor.output <- cbind(shedding.scaled[,date.exp %in% ww.data$date],
                    shedding.scaled[,date.exp %in% ww.data$date],
                    shedding.scaled[,date.exp %in% ww.data$date])
correlation <- cor(t(cor.output),cor.data)[,1]
summary(correlation)
