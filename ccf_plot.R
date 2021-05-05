## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(EpiEstim, tidyverse, incidence, plotrix, distr, data.table,viridis,
       surveillance,mgcv,bbmle,BayesianTools,grDevices,grid,scales)

## covid.data
covid.data <- fread("./COVID_Surveillance_Impact_20200810-20201204.csv")[Description=="student"]
covid.data.total <- covid.data[,.(totalPos=sum(Positive)),by=CRU_COMPLETED_DATE]
    

## ww data
ww.data <- read_csv("ND WW Data_11-30-2020.csv") %>%
    rename(date=X1) %>%
    mutate(date=as.Date(date,"%d-%b-%Y")) %>%
    filter(trust)
N1.med <- apply(ww.data[,2:4],1,median)
N1.mean <- apply(ww.data[,2:4],1,mean)
N1.min <- apply(ww.data[,2:4],1,min)
N1.max <- apply(ww.data[,2:4],1,max)

ww.in.inc <- (ww.data$date %in% covid.data.total$CRU_COMPLETED_DATE)
inc.in.ww <- (covid.data.total$CRU_COMPLETED_DATE %in% ww.data$date)

## Plot it
pdf("cross-correlation_plot.pdf")
ccf(N1.mean[ww.in.inc],covid.data.total$totalPos[inc.in.ww],
    na.action=na.pass,main="",bty="n")
polygon(c(-0.5,0.5,0.5,-0.5),c(-1,-1,1,1),
        col=adjustcolor("red",alpha.f=0.2),border=FALSE)
dev.off()
