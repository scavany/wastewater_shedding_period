## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(EpiEstim, tidyverse, incidence, plotrix, distr, data.table, deSolve,
       surveillance,mgcv,bbmle,BayesianTools,grDevices,viridis,TeachingDemos)

## Get data together and estimate R0
covid.data <- fread("./COVID_Surveillance_Impact_20200810-20201204.csv")[Description=="student"]
dates.onset <- rep(covid.data$CRU_COMPLETED_DATE, covid.data$Positive)
test.type <- rep(covid.data$Test_Type, covid.data$Positive)
test.type[test.type != "Symptomatic"] = "Asymptomatic"
incid <- incidence(as.Date(dates.onset), groups = test.type)
dates <- incid$dates
## cases.ct <- incid$counts[,"Contact_Tracing"]
## cases.tarsurv <- incid$counts[,"ST_Targeted"]
## cases.ransurv <- incid$counts[,"ST_Random"]
cases.symp <- incid$counts[,"Symptomatic"]
cases <- rowSums(incid$counts)
cases.asymp <- incid$counts[,"Asymptomatic"]

### Plan for analyzing wastewater data
## 1. Generate an incidence curve from the Rt estimate, imports, etc.,
## 2. use a parametric distribution (gamma) of individual shedding, to estimate shedding over time
## 3. Calibrate the gamma parameters

## 1. set up params
inc.meanlog <- 1.621
inc.sdlog <- 0.418
inc.reduction <- 1 + log(0.5) / (inc.meanlog + inc.sdlog^2/2)
## deconvolve using  surveillance package
dmax <- 50 # truncation of incubation plus pmf
cases <- c(rep(0,dmax),cases)
cases.symp <- c(rep(0,dmax),cases.symp)
cases.asymp <- c(rep(0,dmax),cases.asymp)
names(cases) <- c(seq(min(dates)-dmax,
                      min(dates)-1,length.out=dmax),dates) #seq(dmax,by=1,length.out=length(cases))
names(cases.symp) <- names(cases)
names(cases.asymp) <- names(cases)
dates <- as.Date(names(cases))
sts.symp <- new("sts", epoch=1:length(cases.symp),observed=matrix(cases.symp,ncol=1))
sts.asymp <- new("sts", epoch=1:length(cases.asymp),observed=matrix(cases.asymp,ncol=1))
## functions to convolve
Inc <- Lnorm(inc.meanlog,inc.sdlog)
Inc.asymp <- Lnorm(inc.meanlog*inc.reduction,inc.sdlog*sqrt(inc.reduction))
pdf(paste0("shedding_examples.pdf"),width=10,height=10)
layout(matrix(c(1,1:5), nrow = 3,byrow=TRUE))
## plot it
shed.shape <- c(1,2,3,5,9,7.5,0.5)
shed.rate <- 1/c(2,2,2,1,0.5,1,1)
colors <- viridis(length(shed.shape))
plot((0:200)/10,dgamma((0:200)/10,shed.shape[1],shed.rate[1]),col=colors[1],type='l',
     las=1,xaxs='i',yaxs='i',bty='n',xlab="Days since infection",
     ylab="Relative shedding intensity",ylim=c(0,0.5),lwd=2)
for (ii in 2:length(shed.shape)) {
    lines((0:200)/10,dgamma((0:200)/10,shed.shape[ii],shed.rate[ii]),col=colors[ii],lwd=2)
}
legend("topright",
       legend=sapply(1:length(shed.shape),
                     ## function(i) as.expression(substitute(paste(mu==A,", ",sigma^2==B),
                     ##                                      list(A = as.name(shed.shape[i]/shed.rate[i]),
                     ##                                           B = as.name(shed.shape[i]/shed.rate[i]^2))))),
                     function(i) as.expression(substitute(paste(alpha==A,", ",beta==B),
                                                          list(A = as.name(shed.shape[i]),
                                                               B = as.name(shed.rate[i]))))),
       col=colors,lty=1,lwd=2,bty="n")
## Loop over different delays
for (delay.mean in c(0,1,2,5)) {
    delay.rate <- 1/delay.mean
    if (delay.mean) {
        Delay <- Pois(delay.mean) ## Try a shift, and a poisson
        p.inf2test <- p(Inc+Delay)
        d.inf2test <- d(Inc+Delay)
    } else {
        p.inf2test <- p(Inc)
        d.inf2test <- d(Inc)
    }
    p.inf2test.asymp <- p(Inc.asymp)
    d.inf2test.asymp <- d(Inc.asymp)
    ## convert pdf to pmf
    inc.symp.pmf <- c(0,(p.inf2test(1:dmax) - p.inf2test(0:(dmax-1)))/p.inf2test(dmax))
    inc.asymp.pmf <- c(0,(p.inf2test.asymp(1:dmax)-p.inf2test.asymp(0:(dmax-1)))/p.inf2test.asymp(dmax))
    ##bpnp.control <- list(k=2,eps=rep(0.005,2),iter.max=rep(250,2),B=-1)
    bpnp.control <- list(k=10,eps=rep(1e-5,2),iter.max=rep(250,2),B=-1)
    #Fast C version (use argument: eq3a.method="C")! 
    sts.symp.bp <- backprojNP(sts.symp, incu.pmf=inc.symp.pmf,
                              control=modifyList(bpnp.control,list(eq3a.method="C")))
    sts.asymp.bp <- backprojNP(sts.asymp, incu.pmf=inc.asymp.pmf,
                               control=modifyList(bpnp.control,list(eq3a.method="C")))
    ## 2. 3. ML fitting with optim, with profiled likelihood on the fitting parameter
    n.exp <- upperbound(sts.symp.bp)+upperbound(sts.asymp.bp)
    date.exp <- dates
    ## 4. plot it
    shedding <- array(0,dim=c(length(shed.shape),length(n.exp)))
    for (jj in 1:length(shed.shape)) {
        for (ii in 1:(length(n.exp)-1)) {
            cum.dens.0 <- pgamma(0:(length(n.exp)-ii),shed.shape[jj],shed.rate[jj])
            cum.dens.1 <- pgamma(1:(length(n.exp)-ii+1),shed.shape[jj],shed.rate[jj])
            shedding[jj,] <- shedding[jj,] + n.exp[ii]*c(rep(0,ii-1),cum.dens.1-cum.dens.0)
        }
    }
    start.date <- which(date.exp == as.Date("2020-07-25"))
    plot(date.exp[start.date:length(date.exp)],shedding[1,start.date:length(date.exp)],
         col=colors[1],type='l',
         las=1,xaxs='i',yaxs='i',bty='n',xlab="Date",ylab="Relative RNA concentration in WW",
         ylim=c(0,125),lwd=2)
    title(main=paste0("Mean delay = ",delay.mean))
    for (ii in 2:length(shed.shape)) {
        lines(date.exp,shedding[ii,],col=colors[ii],lwd=2)
    }
}
dev.off()

## Plot with inset plot
pdf(paste0("shedding_examples_inset.pdf"),width=10,height=10)
layout(matrix(c(5,5,1:4), nrow = 3,byrow=TRUE))
## plot it
shed.shape <- c(1,2,3,5,9,7.5,0.5)
shed.rate <- 1/c(2,2,2,1,0.5,1,1)
colors <- viridis(length(shed.shape))
##Loop over different delays
for (delay.mean in c(0,1,2,5)) {
    delay.rate <- 1/delay.mean
    if (delay.mean) {
        Delay <- Pois(delay.mean) ## Try a shift, and a poisson
        p.inf2test <- p(Inc+Delay)
        d.inf2test <- d(Inc+Delay)
    } else {
        p.inf2test <- p(Inc)
        d.inf2test <- d(Inc)
    }
    p.inf2test.asymp <- p(Inc.asymp)
    d.inf2test.asymp <- d(Inc.asymp)
    ## convert pdf to pmf
    inc.symp.pmf <- c(0,(p.inf2test(1:dmax) - p.inf2test(0:(dmax-1)))/p.inf2test(dmax))
    inc.asymp.pmf <- c(0,(p.inf2test.asymp(1:dmax)-p.inf2test.asymp(0:(dmax-1)))/p.inf2test.asymp(dmax))
    ##bpnp.control <- list(k=2,eps=rep(0.005,2),iter.max=rep(250,2),B=-1)
    bpnp.control <- list(k=10,eps=rep(1e-5,2),iter.max=rep(250,2),B=-1)
    #Fast C version (use argument: eq3a.method="C")! 
    sts.symp.bp <- backprojNP(sts.symp, incu.pmf=inc.symp.pmf,
                              control=modifyList(bpnp.control,list(eq3a.method="C")))
    sts.asymp.bp <- backprojNP(sts.asymp, incu.pmf=inc.asymp.pmf,
                               control=modifyList(bpnp.control,list(eq3a.method="C")))
    ## 2. 3. ML fitting with optim, with profiled likelihood on the fitting parameter
    n.exp <- upperbound(sts.symp.bp)+upperbound(sts.asymp.bp)
    date.exp <- dates
    ## 4. plot it
    shedding <- array(0,dim=c(length(shed.shape),length(n.exp)))
    for (jj in 1:length(shed.shape)) {
        for (ii in 1:(length(n.exp)-1)) {
            cum.dens.0 <- pgamma(0:(length(n.exp)-ii),shed.shape[jj],shed.rate[jj])
            cum.dens.1 <- pgamma(1:(length(n.exp)-ii+1),shed.shape[jj],shed.rate[jj])
            shedding[jj,] <- shedding[jj,] + n.exp[ii]*c(rep(0,ii-1),cum.dens.1-cum.dens.0)
        }
    }
    start.date <- which(date.exp == as.Date("2020-07-25"))
    plot(date.exp[start.date:length(date.exp)],shedding[1,start.date:length(date.exp)],
         col=colors[1],type='l',
         las=1,xaxs='i',yaxs='i',bty='n',xlab="Date",ylab="Relative RNA concentration in WW",
         ylim=c(0,125),lwd=2)
    title(main=paste0("Mean delay = ",delay.mean))
    for (ii in 2:length(shed.shape)) {
        lines(date.exp,shedding[ii,],col=colors[ii],lwd=2)
    }
}
plot((0:200)/10,dgamma((0:200)/10,shed.shape[1],shed.rate[1]),col=colors[1],type='l',
     las=1,xaxs='i',yaxs='i',bty='n',xlab="Days since infection",
     ylab="Relative shedding intensity",ylim=c(0,0.5),lwd=2)
for (ii in 2:length(shed.shape)) {
    lines((0:200)/10,dgamma((0:200)/10,shed.shape[ii],shed.rate[ii]),col=colors[ii],lwd=2)
}
subplot(plot(shed.shape,shed.rate,col=colors,las=1,xaxs="i",yaxs="i",bty="n",cex=1.9,
             ylim=c(0,2.1),xlim=c(0,10),pch=19,xlab=expression(alpha),ylab=expression(beta)),
        c(14.5,17.4),c(0.225,0.55))
dev.off()

### Example incidence curves with shedding dists
delay.mean <- 2
delay.rate <- 1/delay.mean
load(paste0("mcmc_quarantine_recovery_lod_delay_pois",floor(1/delay.rate+0.5),".RData"),verbose=TRUE)
load(paste0("nexp_delay_pois",floor(1/delay.rate+0.5),".RData"),verbose=TRUE)
inc.meanlog <- 1.621
inc.sdlog <- 0.418

n.exp.ex1 <- c(n.exp[1:54],n.exp[54:39],rep(0,length(n.exp)-54-54+38))
n.exp.ex2 <- c(n.exp[1:which(dates == as.Date("2020-10-01"))],
               n.exp[44:(43 + length(n.exp) - which(dates == as.Date("2020-10-01")))])
n.exp.ex3 <- c(n.exp[1:60],rep(n.exp[60],45)*(1+0.05*sin((45/5+1:45)*2*pi/(45/2.5))),
               n.exp[48:(48+length(n.exp)-105-1)])
n.exp.ex4 <- c(n.exp[1:57],
               rep(n.exp[57],
                   length(n.exp)-57)*(1+0.05*sin((45/5+1:(length(n.exp)-57))*2*pi/(45/2.5))))
plot(dates,n.exp.ex4)

n.exp.list <- list(n.exp.ex1,n.exp.ex2,n.exp.ex3,n.exp.ex4)

## ## generate the n.exps
## ## 1. fixed R0
## R0 <- 3
## N <- 1000
## D <- 7
## beta.const <- function(R0,N,D,t,S) {
##     rep(R0/N/D,length(t))
## }
## beta.sin <- function(R0,N,D,t,S,period=120,amp=1-1/R0) {
##     R0/N/D*(1 + amp*cos(t*2*pi/period))
## }
## beta.step <- function(R0,N,D,t,S,step=30,Rt=1) {
##     ifelse(t < step,R0,Rt)/N/D
## }
## beta.step2 <- function(R0,N,D,t,S,step=30,Rt=1) {
##     ifelse(t < step,R0,Rt)/S/D
## }
## inc.period <- exp(inc.meanlog + inc.sdlog^2/2)
## I0 <- 1
## times <- seq(0,150,0.1)
## state <- c(S=N-I0,E=0,I=I0,R=0,E.cum=0)
## parms <- c(R0=R0,inc.period=inc.period,D=D,N=N)
## ode.eqns <- function(t,state,parameters,beta){
##     with(as.list(c(state, parameters)),{
##         dS = -I*S*beta(R0,N,D,t,S)
##         dE = I*S*beta(R0,N,D,t,S) - E/inc.period
##         dI = E/inc.period - I/D
##         dR = I/D
##         dE.cum = I*S*beta(R0,N,D,t,S)
##         list(c(dS,dE,dI,dR,dE.cum))
##     })
## }

## dat = ode(y = state, times = times,
##           func = ode.eqns, parms = parms,beta=beta.const)
## plot(dat)

## whole.days <- which(times %in% seq(min(times),max(times),1))
## n.exp <- diff(dat[whole.days,"E.cum"])
## plot(times[whole.days][-1]-0.5,n.exp)
pdf("shedding_dist_example_incidences.pdf")
par(mfrow=c(2,2))
for (iii in 1:length(n.exp.list)) {
    print(iii)
    n.exp <- n.exp.list[[iii]]
    pars <- samples
    shed.shape <- mean(pars[,1])
    shed.rate <- mean(pars[,2])
    shed.scale <- mean(pars[,3])
    maxplot.day <- 30
    shedding <- array(0,dim=c(1,length(n.exp)))
    for (ii in 1:(length(n.exp)-1)) {
        cum.dens.0 <- pgamma(0:(length(n.exp)-ii),shed.shape,shed.rate)
        cum.dens.1 <- pgamma(1:(length(n.exp)-ii+1),shed.shape,shed.rate)
        shedding[1,] <- shedding[1,] + shed.scale*n.exp[ii]*c(rep(0,ii-1),
                                                              cum.dens.1-cum.dens.0)
    }
    par(mar=0.1+c(5,4,4,5))
    plot(dates[which(dates==as.Date("2020-07-20")):length(dates)],
         n.exp[which(dates==as.Date("2020-07-20")):length(dates)],
         type='l',col=viridis(3)[1],lwd=2,
         bty="n",xaxs="i",yaxs="i",las=1,ylim=c(0,max(n.exp)),
         ylab="Infections",xlab="Days")
    if(iii==1) {
        legend("topright",legend=c("Infections","Viral RNA" ),lty=1,col=viridis(3)[1:2],
               bty="n",lwd=2)
    }
    par(new=TRUE)
    plot(dates[which(dates==as.Date("2020-07-20")):length(dates)],
         colMeans(shedding)[which(dates==as.Date("2020-07-20")):length(dates)]/10^5,
         type='l',col=viridis(3)[2],
         lwd=2,bty="n",xaxs="i",yaxs="i",las=1,yaxt="n",xlab="",ylab="",
         xaxt="n",ylim=c(0,1e-5*max(colMeans(shedding))))
    axis(4,las=1)
    mtext(bquote("RNA (10"^.(5)
                 ~ "N1 GC/l)"),4,3,cex=0.75)
}
dev.off()
