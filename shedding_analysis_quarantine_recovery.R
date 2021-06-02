## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(tidyverse,
       incidence,
       plotrix,
       distr,
       data.table,
       viridis,
       surveillance,
       BayesianTools,
       grDevices,
       plotfunctions)

delay.means <- c(0,1,2,5)
for (delay.mean in delay.means) {
    rm(list=ls()[-which(ls() %in% c("delay.mean","delay.means"))])
    ## Get data together
    covid.data <- fread("./COVID_Surveillance_Impact_20200810-20201204.csv")[Description=="student"]
    dates.onset <- rep(covid.data$CRU_COMPLETED_DATE, covid.data$Positive)
    test.type <- rep(covid.data$Test_Type, covid.data$Positive)
    test.type[test.type != "Symptomatic"] = "Asymptomatic"
    incid <- incidence(as.Date(dates.onset), groups = test.type)
    ## Open figure 1
    pdf(paste0("Fig1_recovery_",delay.mean,".pdf"),height=6.92,width=6.92,pointsize=10)
    layout(matrix(c(1,2,3,2),nrow=2))
    barplot.out <- barplot(t(incid$counts)[2:1,],las=1,col=c("black","dark gray"),
                           border=FALSE)
    axis(1,label=FALSE, tcl=0)
    tick.locs <- which(incid$dates %in% seq(as.Date("2020-09-01"),as.Date("2021-01-01"),"1 month"))
    axis(1,at=barplot.out[tick.locs],labels=FALSE)
    label.locs <- which(incid$dates %in% seq(as.Date("2020-08-16"),as.Date("2021-11-16"),"1 month"))
    axis(1,at=barplot.out[label.locs],
         labels=c("Aug","Sep","Oct","Nov"),tick=FALSE)
    mtext("Date",side=1,line=3,cex=3/4)
    mtext("Cases",side=2,line=3,cex=3/4)
    legend("center",bty="n",legend=c("Symptomatic","Non-symptomatic"),
           fill=c("black","gray"),border=FALSE)
    mtext("A",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    abline(v=0,lty="dotted",lwd=2)
    text(-3,75,"Students returning",srt=90)
    cooling.period <- which(incid$dates %in% c(as.Date("2020-08-19"),as.Date("2020-09-02")))
    polygon(c(barplot.out[cooling.period],rev(barplot.out[cooling.period])),
            c(0,0,200,200),border=FALSE,col=adjustcolor("lightblue",0.3))
    text(mean(barplot.out[cooling.period]),150,"Online\ninstruction",
         adj=0.5,srt=90)
    max(rowSums(incid$counts))
    sum(rowSums(incid$counts))
    incid$dates[which.max(rowSums(incid$counts))]
    ## cases.ct <- incid$counts[,"Contact_Tracing"]
    ## cases.tarsurv <- incid$counts[,"ST_Targeted"]
    ## cases.ransurv <- incid$counts[,"ST_Random"]
    dates <- incid$dates
    cases.symp <- incid$counts[,"Symptomatic"]
    cases <- rowSums(incid$counts)
    cases.asymp <- incid$counts[,"Asymptomatic"]
### Plan for analyzing wastewater data
    ## 1. Generate an incidence curve from the Rt estimate, imports, etc.,
    ## 2. use a parametric distribution (gamma) of individual shedding, to estimate shedding over time
    ## 3. Calibrate the gamma parameters
    ## 1. set up params (inc from Lauer, lnorm)
    ## inc.shape <- 5.807
    ## inc.scale <- 0.948
    inc.meanlog <- 1.621
    inc.sdlog <- 0.418
    inc.mean <- exp(inc.meanlog+inc.sdlog^2/2)
    inc.reduction <- 1 + log(0.5) / (inc.meanlog + inc.sdlog^2/2) 
    delay.rate <- 1/delay.mean
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
    Delay <- Pois(delay.mean) ## Try a shift, and a poisson
    ## Delay <- Gammad(shape=delay.mean/1e6,scale=1e6)
    if (delay.mean) {
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
                                        #Call non-parametric back-projection function with hook function but
                                        #without bootstrapped confidence intervals
    plotIt <- function(cur.sts) {
        plot(cur.sts,xaxis.labelFormat=NULL, legend=NULL,ylim=c(0,100))
    }
    ##bpnp.control <- list(k=2,eps=rep(0.005,2),iter.max=rep(250,2),B=-1)
    bpnp.control <- list(k=10,eps=rep(1e-5,2),iter.max=rep(250,2),B=-1)
                                        #Fast C version (use argument: eq3a.method="C")! 
    sts.symp.bp <- backprojNP(sts.symp, incu.pmf=inc.symp.pmf,
                              control=modifyList(bpnp.control,list(eq3a.method="C")))
    sts.asymp.bp <- backprojNP(sts.asymp, incu.pmf=inc.asymp.pmf,
                               control=modifyList(bpnp.control,list(eq3a.method="C")))
    n.exp.symp <- upperbound(sts.symp.bp)
    n.exp.asymp <- upperbound(sts.asymp.bp)
    n.exp <- n.exp.symp+n.exp.asymp
    save(dates,n.exp,n.exp.symp,n.exp.asymp,
         file=paste0("nexp_delay_pois",floor(1/delay.rate+0.5),".RData"))
    ## pdf("infections.pdf")
    ## par(mfrow=c(2,1))
    plot(dates,n.exp.symp,bty="n",yaxs="i",xaxs="i",las=1,type="l",col="blue",lwd=1,
         xlim=c(dates[30],dates[nrow(sts.symp.bp)]),xaxt="n",xlab="Date",
         ylim=c(0,max(incid$counts[,"Symptomatic"])),ylab="No. Infected")
    ## plot(sts.symp.bp,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),main="",
    ##      bty="n",yaxs="i",xaxs="i",las=1,xlim=c(30,nrow(sts.symp.bp)),xaxt="n",xlab="Date")
    tick.locs <- dates[dates %in% seq(as.Date("2020-07-01"),as.Date("2021-01-01"),"1 month")]
    axis(1,at=tick.locs,label=FALSE)
    label.locs <- dates[dates %in% seq(as.Date("2020-07-16"),as.Date("2021-01-16"),"1 month")]
    axis(1,at=label.locs,label=c("Jul","Aug","Sep","Oct","Nov"),tick=FALSE)
    ## add_bars(incid$dates,rowSums(incid$counts),col="white")
    add_bars(incid$dates,incid$counts[,"Symptomatic"],col="white")
    lines(dates,n.exp.symp,col=viridis(3)[1],lwd=3)
##Do the convolution for the expectation
    mu.symp <- matrix(0,ncol=ncol(sts.symp.bp),nrow=nrow(sts.symp.bp))
                                        #Loop over all series
    for (j in 1:ncol(sts.symp.bp)) { 
                                        #Loop over all time points
        for (t in 1:nrow(sts.symp.bp)) {
                                        #Convolution, note support of inc.pmf starts at zero (move idx by 1)
            i <- seq_len(t)
            mu.symp[t,j] <- sum(inc.symp.pmf[t-i+1] * upperbound(sts.symp.bp)[i,j],na.rm=TRUE)
        }
    }
    mu.asymp <- matrix(0,ncol=ncol(sts.asymp.bp),nrow=nrow(sts.asymp.bp))
                                        #Loop over all series
    for (j in 1:ncol(sts.asymp.bp)) { 
                                        #Loop over all time points
        for (t in 1:nrow(sts.asymp.bp)) {
                                        #Convolution, note support of inc.pmf starts at zero (move idx by 1)
            i <- seq_len(t)
            mu.asymp[t,j] <- sum(inc.asymp.pmf[t-i+1] * upperbound(sts.asymp.bp)[i,j],na.rm=TRUE)
        }
    }
                                        #Show the fit
    ## lines(dates,mu.symp[,1]+mu.asymp[,1],col=viridis(3)[2],lwd=3,lty="dashed")
    lines(dates,mu.symp[,1],col=viridis(3)[2],lwd=3,lty="dashed")
    cooling.period <- dates[dates %in% c(as.Date("2020-08-19"),as.Date("2020-09-02"))]
    polygon(c(cooling.period,rev(cooling.period)),
            c(0,0,200,200),border=FALSE,col=adjustcolor("lightblue",0.4))
    abline(v=as.Date("2020-08-09"),lty="dotted",lwd=2)
    legend("center", lty=c(NA,"solid","dashed"),col=c(NA,viridis(3)[1:2]),fill=c("white",NA,NA),
           legend=c("Symptomatic cases","Infections","Predicted Cases"),bty="n",border=c("black",NA,NA),lwd=c(NA,3,3))
    mtext("B",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    ## plot(sts.asymp.bp,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),main="",
    ##      xaxs="i",bty="n",yaxs="i",las=1,xlim=c(30,nrow(sts.symp.bp)),xaxt="n",xlab="Date")
    ## tick.locs <- which(dates %in% seq(as.Date("2020-07-01"),as.Date("2021-01-01"),"1 month"))
    ## axis(1,at=tick.locs,label=FALSE)
    ## label.locs <- which(dates %in% seq(as.Date("2020-07-16"),as.Date("2021-01-16"),"1 month"))
    ## axis(1,at=label.locs,label=c("Jul","Aug","Sep","Oct","Nov"),tick=FALSE)
    ##                                     #Do the convolution for the expectation
                                        #Show the fit
    ## lines(1:nrow(sts.asymp.bp)-0.5,mu[,1],col="green",lwd=3,lty="dashed")
    ## cooling.period <- which(dates %in% c(as.Date("2020-08-19"),as.Date("2020-09-02")))
    ## abline(v=which(dates == as.Date("2020-08-09")),lty="dotted",lwd=2)
    ## polygon(c(cooling.period,rev(cooling.period)),
    ##         c(0,0,200,200),border=FALSE,col=adjustcolor("lightblue",0.4))
    ## mtext("C",side=3,line=0, 
    ##       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
    ##       cex=1.2)
    ##dev.off()
    ## get data
    LoD95 <- 220 # N1 GC/L
    ww.data <- read_csv("NDCampus_ddPCRdata_cleanedrawdata_0s.csv") %>%
        mutate(Date=as.Date(Date,"%m/%d/%Y")) %>%
        select(-WW_flow) %>%
        ## select(Date,N1R1RC,N1R2RC,N1R3RC) %>%
        drop_na()
    N1.med <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,median)
    N1.mean <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,mean)
    N1.min <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,min)
    N1.max <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,max)
    N1R1RC.D <- ww.data$N1R1RC[ww.data$N1R1RC>0.1]
    N1R2RC.D <- ww.data$N1R2RC[ww.data$N1R2RC>0.1]
    N1R3RC.D <- ww.data$N1R3RC[ww.data$N1R3RC>0.1]
    LoD95RC <- LoD95 / ww.data$recovery_efficiency
    ww.dates.R1D <- ww.data$Date[ww.data$N1R1RC>0.1]
    ww.dates.R2D <- ww.data$Date[ww.data$N1R2RC>0.1]
    ww.dates.R3D <- ww.data$Date[ww.data$N1R3RC>0.1]
    ww.dates.R1ND <- ww.data$Date[ww.data$N1R1RC<0.1]
    ww.dates.R2ND <- ww.data$Date[ww.data$N1R2RC<0.1]
    ww.dates.R3ND <- ww.data$Date[ww.data$N1R3RC<0.1]    
    LoD95RC.R1ND <- LoD95RC[ww.data$N1R1RC<0.1]
    LoD95RC.R2ND <- LoD95RC[ww.data$N1R2RC<0.1]
    LoD95RC.R3ND <- LoD95RC[ww.data$N1R3RC<0.1]
    N1.med.RC <- ifelse(N1.med<LoD95RC,(N1.med+LoD95RC)/2,N1.med)
    N1.min.RC <- ifelse(N1.min<LoD95RC,(N1.min+LoD95RC)/2,N1.min)
    N1.max.RC <- ifelse(N1.max<LoD95RC,(N1.max+LoD95RC)/2,N1.max)
    N1.mean.RC <- (N1.med.RC + N1.min.RC + N1.max.RC)/3

    ## pdf("wastewater_data.pdf")
    plotCI(ww.data$Date,N1.med/1e6,li=N1.min/1e6,ui=N1.max/1e6,gap=TRUE,sfrac=0.003,
           xlab="Date",ylab="",
           xaxs="i",bty="n",las=1,xaxt="n")
    mtext(bquote("RNA (10"^6 ~ "N1 GC/l)"),2,3,cex=3/4)
    axis(1,label=FALSE, tcl=0)
    tick.locs <- which(ww.data$Date %in% seq(as.Date("2020-07-01"),as.Date("2021-01-01"),"1 month"))
    axis(1,at=ww.data$Date[tick.locs],label=FALSE)
    label.locs <- which(ww.data$Date %in% seq(as.Date("2020-07-16"),as.Date("2021-01-16"),"1 month"))
    axis(1,at=ww.data$Date[label.locs],label=c("Aug","Sep","Oct","Nov"),tick=FALSE)
    abline(v=as.Date("2020-08-09"),lty="dotted",lwd=2)
    polygon(c(cooling.period,rev(cooling.period)),
            c(0,0,1.5e6,1.5e6),border=FALSE,col=adjustcolor("lightblue",0.3))
    ## The next bit is the logged version
    ## plotCI(ww.data$Date,pmax(N1.med,LoD95),li=pmax(N1.min,LoD95),ui=pmax(LoD95,N1.max),
    ##        gap=TRUE,sfrac=0.003,
    ##        xlab="Date",ylab="RNA (N1 GC / l)",xaxs="i",bty="n",las=1,
    ##        xaxt="n",log="y")
    ## abline(h=LoD95,lty="dashed")
    ## axis(1,label=FALSE, tcl=0)
    ## tick.locs <- which(ww.data$Date %in% seq(as.Date("2020-07-01"),as.Date("2021-01-01"),"1 month"))
    ## axis(1,at=ww.data$Date[tick.locs],label=FALSE)
    ## label.locs <- which(ww.data$Date %in% seq(as.Date("2020-07-16"),as.Date("2021-01-16"),"1 month"))
    ## axis(1,at=ww.data$Date[label.locs],label=c("Aug","Sep","Oct","Nov"),tick=FALSE)
    ## abline(v=as.Date("2020-08-09"),lty="dotted",lwd=2)
    ## polygon(c(cooling.period,rev(cooling.period)),
    ##         c(LoD95,LoD95,1.5e6,1.5e6),border=FALSE,col=adjustcolor("lightblue",0.3))
    mtext("C",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    dev.off()

    ## 2. MCMC fitting with BayesianTools
    date.exp <- dates
    quarantine.length <- 10
    ## Set up proportion quarantining
    load("./prop_quarantine.RData")
    data.dates <- as.Date(names(prop.q))
    p.q.vec <- c(rep(0,min(data.dates)-min(date.exp)),
                 as.numeric(prop.q),
                 rep(0,max(date.exp)-max(data.dates)))
    NLL.dispersion <- function(par) {
        shed.shape <- par[1]
        shed.rate <- par[2]
        p.d <- par[5]
        p.scale.asymp <- 0#par[6]
        shedding <- rep(0,length(n.exp))
        for (jj in 1:(length(n.exp)-1)) {
            cum.dens.0 <- pgamma(0:(length(n.exp)-jj),shed.shape,shed.rate)
            cum.dens.1 <- pgamma(1:(length(n.exp)-jj+1),shed.shape,shed.rate)
            individual.shedding <- cum.dens.1-cum.dens.0
            p.enter <- p.inf2test(0:(length(n.exp)-jj))
            p.enter.asymp <- p.inf2test.asymp(0:(length(n.exp)-jj))
            if (length(n.exp) - jj >= quarantine.length) {
                p.exit <- c(rep(0,quarantine.length),
                            p.inf2test(0:(length(n.exp)-jj-quarantine.length)))
                p.exit.asymp <- c(rep(0,quarantine.length),
                                  p.inf2test.asymp(0:(length(n.exp)-jj-quarantine.length)))
            } else {
                p.exit <- rep(0,length(n.exp)-jj+1)
                p.exit.asymp <- rep(0,length(n.exp)-jj+1)
            }
            ## Make sure interpret the p.q correctly in writing! Is it more like an odds than a proportion?
            p.q <- p.q.vec[jj+1:length(individual.shedding)-1]
            individual.shedding.symp <- individual.shedding * (1 + p.d * p.q * (p.exit - p.enter))
            individual.shedding.asymp <- individual.shedding * p.scale.asymp * (1 + p.d * p.q * (p.exit.asymp - p.enter.asymp))
            shedding <- shedding + n.exp.asymp[jj]*c(rep(0,jj-1),individual.shedding.asymp)
            shedding <- shedding + n.exp.symp[jj]*c(rep(0,jj-1),individual.shedding.symp)
        }
        shedding.scaled <- shedding*par[3]/p.d
        nbinom.size <- par[4]
        ## tryCatch(-sum(dnbinom(N1R1RC.D,size=nbinom.size,
        ##                         mu=shedding.scaled[date.exp %in% ww.dates.R1D],log=TRUE)) -
        ##            sum(dnbinom(N1R2RC.D,size=nbinom.size,
        ##                        mu=shedding.scaled[date.exp %in% ww.dates.R2D],log=TRUE)) -
        ##            sum(dnbinom(N1R3RC.D,size=nbinom.size,
        ##                        mu=shedding.scaled[date.exp %in% ww.dates.R3D],log=TRUE)) -
        ##            sum(pnbinom(LoD95RC.R1ND,size=nbinom.size,
        ##                        mu=shedding.scaled[date.exp %in% ww.dates.R1ND],log.p=TRUE)) -
        ##            sum(pnbinom(LoD95RC.R2ND,size=nbinom.size,
        ##                        mu=shedding.scaled[date.exp %in% ww.dates.R2ND],log.p=TRUE)) -
        ##            sum(pnbinom(LoD95RC.R3ND,size=nbinom.size,
        ##                        mu=shedding.scaled[date.exp %in% ww.dates.R3ND],log.p=TRUE)),
        ##          warning = function(w) {
        ##              message(c("Warning at parms",par))
        ##          })
        if (sum(shedding.scaled)) {
            return(-sum(dnbinom(N1R1RC.D,size=nbinom.size,
                                mu=shedding.scaled[date.exp %in% ww.dates.R1D],log=TRUE)) -
                   sum(dnbinom(N1R2RC.D,size=nbinom.size,
                               mu=shedding.scaled[date.exp %in% ww.dates.R2D],log=TRUE)) -
                   sum(dnbinom(N1R3RC.D,size=nbinom.size,
                               mu=shedding.scaled[date.exp %in% ww.dates.R3D],log=TRUE)) -
                   sum(pnbinom(LoD95RC.R1ND,size=nbinom.size,
                               mu=shedding.scaled[date.exp %in% ww.dates.R1ND],log.p=TRUE)) -
                   sum(pnbinom(LoD95RC.R2ND,size=nbinom.size,
                               mu=shedding.scaled[date.exp %in% ww.dates.R2ND],log.p=TRUE)) -
                   sum(pnbinom(LoD95RC.R3ND,size=nbinom.size,
                               mu=shedding.scaled[date.exp %in% ww.dates.R3ND],log.p=TRUE)))
        } else {
            print("model all zero")
            return(Inf)
        }
    }
    LL.dispersion <- function(par) {-NLL.dispersion(par)}
    ## optim.out <- optim(c(log(1),log(1),log(1),log(1)),NLL.dispersion)
    ## mle2.out <- mle2(NLL.mle2.dispersion)
    ## mle2.profile <- profile(mle2.out)
    ## save(mle2.profile,mle2.out,file=paste0("mle2_delay_pois",floor(1/delay.rate+0.5),".RData"))
    ## load(paste0("mle2_delay_pois",floor(1/delay.rate+0.5),".RData"),verbose=TRUE)

    ## MCMC - with dispersion
    lower <- c(shape=1,rate=0,shedscale=0,
               nbinomsize=0,p.d=0)#,p.inf.asymp=0)#confint(mle2.profile)[,1]*0.1
    upper <- c(shape=1e3,rate=1e3,shedscale=1e5,
               nbinomsize=1e5,p.d=1)#,p.inf.asymp=1)#confint(mle2.profile)[,2]*1.9
    bayesianSetup = createBayesianSetup(LL.dispersion,lower=lower,upper=upper)
    ## Do 3e5 for 2,5. Or 9e5 for 1 or 0 (converge more slowly)
    mcmc.out = runMCMC(bayesianSetup,
                       settings=list(iterations=9e5))
    pdf(paste0("mcmc_out_quarantine_recovery_",delay.mean,".pdf"))
    plot(mcmc.out,start = 1.5e5) # 3e4 for 2,5; 1.5e5 for 0,1
    dev.off()
    samples = getSample(mcmc.out,start=1.5e5)
    save(mcmc.out,samples,file=paste0("mcmc_quarantine_recovery_lod_delay_pois",floor(1/delay.rate+0.5),".RData"))
    load(paste0("mcmc_quarantine_recovery_lod_delay_pois",floor(1/delay.rate+0.5),".RData"))

    DIC(mcmc.out)$DIC

    ## marginalPlot(mcmc.out, prior = TRUE, start=30000)
    pdf(paste0("mcmc_correlation_quarantine_recovery_",delay.mean,".pdf"))
    correlationPlot(mcmc.out,start=1.5e5)
    dev.off()
    ## gelmanDiagnostics(mcmc.out, plot = T)

    ## 4. plot it
    shed.shape <- samples[,1]
    shed.rate <- samples[,2]
    shed.scale <- samples[,3]
    maxplot.day <- 30
    nbinom.size <- samples[,4]
    p.d <- samples[,5]
    p.inf.asymp <- rep(0,nrow(samples))#samples[,6]
    shedding <- array(0,dim=c(length(shed.shape),length(n.exp)))
    shed.dist <- array(NA,dim=c(length(shed.shape),maxplot.day+1))

    ## n.samples <- 1e6
    ## sampled.rows <- sample(nrow(samples),n.samples,replace=FALSE)
    for (jj in 1:length(shed.shape)) {
        for (ii in 1:(length(n.exp)-1)) {
            cum.dens.0 <- pgamma(0:(length(n.exp)-ii),shed.shape[jj],shed.rate[jj])
            cum.dens.1 <- pgamma(1:(length(n.exp)-ii+1),shed.shape[jj],shed.rate[jj])
            individual.shedding <- cum.dens.1-cum.dens.0
            p.enter <- p.inf2test(0:(length(n.exp)-ii))
            p.enter.asymp <- p.inf2test.asymp(0:(length(n.exp)-ii))
            if (length(n.exp) - ii >= quarantine.length) {
                p.exit <- c(rep(0,quarantine.length),
                            p.inf2test(0:(length(n.exp)-ii-quarantine.length)))
                p.exit.asymp <- c(rep(0,quarantine.length),
                                  p.inf2test.asymp(0:(length(n.exp)-ii-quarantine.length)))
            } else {
                p.exit <- rep(0,length(n.exp)-ii+1)
                p.exit.asymp <- rep(0,length(n.exp)-ii+1)
            }
            p.q <- p.q.vec[ii+1:length(individual.shedding)-1]
            individual.shedding.symp <- individual.shedding * (1 + p.d[jj] * p.q * (p.exit - p.enter))
            individual.shedding.asymp <- individual.shedding * p.inf.asymp[jj] * (1 + p.d[jj] * p.q * (p.exit.asymp - p.enter.asymp))
            shedding[jj,] <- shedding[jj,] + shed.scale[jj]*n.exp.asymp[ii]*c(rep(0,ii-1),
                                                                              individual.shedding.asymp) / p.d[jj]
            shedding[jj,] <- shedding[jj,] + shed.scale[jj]*n.exp.symp[ii]*c(rep(0,ii-1),
                                                                             individual.shedding.symp) / p.d[jj]
        }
        shed.dist[jj,] <- dgamma(0:maxplot.day,shed.shape[jj],shed.rate[jj])*shed.scale[jj]
    }
    save(shed.dist,shedding,file=paste0("shedding_dists_pois",delay.mean,".RData"))

    load(paste0("shedding_dists_pois",delay.mean,".RData"))

    ## Prediction interval
    n.samples <- 1e5 
    sampled.rows <- sample(nrow(shedding),n.samples,replace=TRUE)
    sampled.mu <- shedding[sampled.rows,]
    pred <- matrix(rnbinom(n.samples*ncol(sampled.mu),
                           size=nbinom.size[sampled.rows],
                           mu=t(sampled.mu)),
                   nrow=nrow(sampled.mu),
                   ncol=ncol(sampled.mu),byrow=TRUE)

    ## plot it
    ## pdf(paste0("shedding_distribution_pois_",floor(1/delay.rate+0.5),"_mean_delay.pdf"))
    pdf(paste0("Fig2_quarantine_recovery_lod_",delay.mean,".pdf"),width=6.92,height=6.92,pointsize=10)
    layout(matrix(c(1,2,3,3),byrow=TRUE,nrow=2))
    p.enter <- p.inf2test(0:(ncol(shed.dist)-1))
    p.exit <- c(rep(0,quarantine.length),
                p.inf2test(0:(ncol(shed.dist)-quarantine.length-1)))
    shed.dist.mean <- colMeans(shed.dist)
    plot(0:maxplot.day, shed.dist.mean * (1 - p.enter + p.exit),
         type='l', ylim=c(0,max(apply(shed.dist[,2:(maxplot.day+1)],2,function(x)quantile(x,0.975)))),
         ylab="Shedding intensity (N1 GC / l)",xlab="Time since infection",yaxs="i",las=1,bty="n",
         col="red",lwd=2)
    ## plot(0:maxplot.day, shed.dist.mean * (1 - p.enter + p.exit),
    ##      type='l', ylim=c(0.9*LoD95,
    ##                       max(apply(shed.dist[,2:(maxplot.day+1)],2,function(x)quantile(x,0.975)))),
    ##      ylab="Shedding intensity (N1 GC / l)",xlab="Time since infection",yaxs="i",las=1,bty="n",
    ##      col="red",lwd=2,log="y")
    ## abline(h=LoD95,lty="dashed")
    lines(0:maxplot.day,shed.dist.mean,lwd=2)
    ## lines(0:maxplot.day,apply(shed.dist,2,median),lwd=2)
    mtext("A",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    polygon(c(0:30,30:0),c(apply(shed.dist[,],2,function(x)quantile(x,0.025)),
                           rev(apply(shed.dist[,],2,function(x)quantile(x,0.975)))),
            col=adjustcolor("gray",0.5),border=F)
    abline(v=which.max(colMeans(shed.dist[,2:(maxplot.day+1)])),lty="dashed")
    shed.conc.max <- max(colMeans(shedding))
    ## barplot(dnbinom(0:10,size=mean(nbinom.size),mu=shed.conc.max),ylab="Proportion of samples",
    ##         names.arg=0:10,
    ##         xlab="Measured viral concentration at peak (N1 GC / l)")
    shed.conc.max <- max(colMeans(shedding))
    increment <- 10^floor(log(shed.conc.max,10))
    mean.group <- floor(shed.conc.max/increment) + 1
    barplot(c(pnbinom(increment-1,size=mean(nbinom.size),mu=shed.conc.max),
              diff(pnbinom(1:9*increment-1,size=mean(nbinom.size),mu=shed.conc.max)),
              1 - pnbinom(9*increment-1,size=mean(nbinom.size),mu=shed.conc.max)),
            ylab="Proportion of samples",
            names.arg=c(paste0("[",seq(0,8,1),", ",seq(1,9,1),"]"),
                        expression(paste("[9, ", infinity,"]"))), 
            col=c(rep("grey",mean.group - 1),"red",rep("grey",10 - mean.group)),
            xlab=bquote("Predicted recovery-corrected viral concentration (10"^.(log(increment,10))
                        ~ "N1 GC/l)"))

    ## dnbinom(0,size=mean(nbinom.size),mu=shed.conc.max)
    pnbinom(LoD95RC[ww.data$Date == dates[which.max(colMeans(shedding))]],
            size=mean(nbinom.size),mu=shed.conc.max)
    mtext("B",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    ##dev.off()
    1+(0:30)[which.max(shed.dist.mean[2:(maxplot.day+1)])]
    quantile((0:30)[apply(shed.dist,1,which.max)],c(0.025,0.25,0.75,0.975))
    quantile(nbinom.size,c(0.025,0.975))
    mean(nbinom.size)

    ##pdf(paste0("shedding_timeseries_pois_",floor(1/delay.rate+0.5),"_mean_delay.pdf"))
    lower95 <- apply(pred,2,function(x)quantile(x,0.025))
    upper95 <- apply(pred,2,function(x)quantile(x,0.975))
    lower50 <- apply(pred,2,function(x)quantile(x,0.25))
    upper50 <- apply(pred,2,function(x)quantile(x,0.75))
    plot(date.exp,colMeans(pred),type='l',xlim=c(min(ww.data$Date),max(ww.data$Date)),
         ylim=c(0,max(N1.max)),
         ylab="RNA (N1 GC / l)",xlab="Date",las=1,bty="n",col='red',lwd=2)
    ## plot(date.exp,colMeans(pred),type='l',xlim=c(min(ww.data$Date),max(ww.data$Date)),
    ##      ylim=c(0.9*LoD95,max(N1.max)),log="y",
    ##      ylab="RNA (N1 GC / 100ml)",xlab="Date",las=1,bty="n",col='red',lwd=2)
    ## abline(h=LoD95,lty="dashed")
    mtext("C",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    plotCI(ww.data$Date,N1.med.RC,li=N1.min.RC,ui=N1.max.RC,add=TRUE,gap=TRUE,sfrac=0.003)    
    polygon(c(date.exp,rev(date.exp)),c(lower95,
                                        rev(upper95)),
            col=adjustcolor("gray",0.25),border=F)
    polygon(c(date.exp,rev(date.exp)),c(lower50,
                                        rev(upper50)),
            col=adjustcolor("gray",0.5),border=F)
    ## for (ii in sample(nrow(pred),1,TRUE)) {
    ##     other.measurements <- sample(nrow(pred),2,TRUE)
    ##     points(date.exp,(pred[ii,]+colSums(pred[other.measurements,]))/3,
    ##            col=adjustcolor("red",0.5))
    ## }
    ## for (ii in sample(nrow(pred),10)) {
    ##     lines(date.exp,pred[ii,],col=adjustcolor("red",alpha.f=0.5))
    ## }
    dev.off()

    ww.data$Date[which.max(N1.max)]
    date.exp[which.max(colMeans(pred))]
    
    indices <- (dates %in% ww.data$Date)
    ## Proportion in 95% interval... but this complicated by non-detects
    (sum((N1.med.RC <= upper95[indices]) & (N1.med.RC >= lower95[indices])) + sum((N1.min.RC <= upper95[indices]) & (N1.min.RC >= lower95[indices])) + sum((N1.max.RC <= upper95[indices]) & (N1.max.RC >= lower95[indices])))/3/sum(indices)
    (sum((N1.med.RC <= upper50[indices]) & (N1.med.RC >= lower50[indices])) + sum((N1.min.RC <= upper50[indices]) & (N1.min.RC >= lower50[indices])) + sum((N1.max.RC <= upper50[indices]) & (N1.max.RC >= lower50[indices])))/3/sum(indices)

    ## Instead do proportion less than 97.5% and 75%, including non-detects
    (sum((N1.med <= upper95[indices])) + sum((N1.min <= upper95[indices])) + sum((N1.max <= upper95[indices])))/3/sum(indices)
    (sum((N1.med <= upper50[indices])) + sum((N1.min <= upper50[indices])) + sum((N1.max <= upper50[indices])))/3/sum(indices)

    ## correlations
    cor(colMeans(pred)[indices],N1.mean,method="spearman")
    quantile(apply(pred[,indices],1,function(x)cor(x,N1.mean)),c(0.025,0.975))


    pdf(paste0("Fig2_quarantine_recovery_lod_",delay.mean,"_log.pdf"),width=6.92,height=6.92,
        pointsize=10)
    layout(matrix(c(1,2,3,3),byrow=TRUE,nrow=2))
    plot(0:maxplot.day, shed.dist.mean * (1 - p.enter + p.exit),
         type='l', ylim=c(min(shed.dist.mean[2:maxplot.day] * (1 - p.enter + p.exit)),
                          max(apply(shed.dist[,2:(maxplot.day+1)],2,function(x)quantile(x,0.975)))),
         ylab="Shedding intensity (N1 GC / l)",xlab="Time since infection",yaxs="i",las=1,bty="n",
         col="red",lwd=2,log="y")
    ##lines(0:maxplot.day,LoD95RC,lty="dashed")
    lines(0:maxplot.day,shed.dist.mean,lwd=2)
    mtext("A",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    polygon(c(0:30,30:0),c(pmax(apply(shed.dist[,],2,function(x)quantile(x,0.025)),1e-5),
                           rev(apply(shed.dist[,],2,function(x)quantile(x,0.975)))),
            col=adjustcolor("gray",0.5),border=F)
    abline(v=which.max(colMeans(shed.dist[,2:(maxplot.day+1)])),lty="dashed")
    barplot(c(pnbinom(increment-1,size=mean(nbinom.size),mu=shed.conc.max),
              diff(pnbinom(1:9*increment-1,size=mean(nbinom.size),mu=shed.conc.max)),
              1 - pnbinom(9*increment-1,size=mean(nbinom.size),mu=shed.conc.max)),
            ylab="Proportion of samples",
            names.arg=c(paste0("[",seq(0,8,1),", ",seq(1,9,1),"]"),
                        expression(paste("[9, ", infinity,"]"))), 
            col=c(rep("grey",mean.group - 1),"red",rep("grey",10 - mean.group)),
            xlab=bquote("Predicted recovery-corrected viral concentration (10"^.(log(increment,10))
                        ~ "N1 GC/l)"))
    ## barplot(c(pnbinom(99999,size=mean(nbinom.size),mu=shed.conc.max),
    ##           diff(pnbinom(1:10*1e5-1,size=mean(nbinom.size),mu=shed.conc.max))),
    ##         ylab="Proportion of samples",
    ##         names.arg=paste0("[",seq(0,9,1),", ",seq(1,10,1),"]"),
    ##         col=c(rep("grey",2),"red",rep("grey",7)),
    ##         xlab=expression(paste("Measured recovery-corrected viral concentration at peak (10"
    ##                               ^5, " N1 GC / l)")))
    mtext("B",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    plot(date.exp,colMeans(pred),type='l',xlim=c(min(ww.data$Date),max(ww.data$Date)),
         ylim=c(min(N1.min.RC),max(N1.max)),log="y",
         ylab="RNA (N1 GC / l)",xlab="Date",las=1,bty="n",col='red',lwd=2)
    lines(ww.data$Date,LoD95RC,lty="dashed")
    mtext("C",side=3,line=0, 
          at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
          cex=1)
    plotCI(ww.data$Date,N1.med.RC,
           li=N1.min.RC,ui=N1.max.RC,
           add=TRUE,gap=TRUE,sfrac=0.003)
    lower95 <- ifelse(lower95<LoD95,LoD95/2,lower95)
    upper95 <- ifelse(upper95<LoD95,LoD95/2,upper95)
    lower50 <- ifelse(lower50<LoD95,LoD95/2,lower50)
    upper50 <- ifelse(upper95<LoD95,LoD95/2,upper50)
    polygon(c(date.exp,rev(date.exp)),
            c(lower95,rev(upper95)),
            col=adjustcolor("gray",0.25),border=F)
    polygon(c(date.exp,rev(date.exp)),
            c(lower50,rev(upper50)),
            col=adjustcolor("gray",0.5),border=F)
    ## for (ii in sample(nrow(pred),1,TRUE)) {
    ##     other.measurements <- sample(nrow(pred),2,TRUE)
    ##     points(date.exp,(pred[ii,]+colSums(pred[other.measurements,]))/3,
    ##            col=adjustcolor("red",0.5))
    ## }
    dev.off()    

    ## Calculate the probability that the incubation period is completed before the peak in shedding
    mean(plnorm((shed.shape-1)/shed.rate,inc.meanlog,inc.sdlog))
    quantile(plnorm((shed.shape-1)/shed.rate,inc.meanlog,inc.sdlog),c(0.025,0.975))

    ## Calculate the proportion of shedding which is after the incubation period
    1 - mean(pgamma(rlnorm(1e6,inc.meanlog,inc.sdlog),mean(shed.shape),mean(shed.rate)))
    1 - quantile(pgamma(rlnorm(1e6,inc.meanlog,inc.sdlog),mean(shed.shape),mean(shed.rate)),c(0.025,0.975))
}
