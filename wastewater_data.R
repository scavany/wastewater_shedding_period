library(data.table)
library(zoo)
library(grDevices)

## load data
elkhart.epi <- fread("Elkhart Epi Data through 10-25-2020.csv")
elkhart.ww <- fread("Elkhart WW Data through 10-25-2020.csv")
nd.epi <- fread("ND Epi Data_11-30-2020.csv")
nd.ww <- fread("ND WW Data_11-30-2020.csv")

## process data
colnames(elkhart.epi) <- c("date",
                           "incidence","incidence.rollmean",
                           "positivity","positivity.rollmean")
colnames(elkhart.ww) <- c("date",
                          "N1.1","N1.2","N1.3",
                          "PMMoV.1","PMMoV.2","PMMoV.3",
                          "recovery.1","recovery.2","recovery.3",
                          "ratioN1PMMoV.1","ratioN1PMMoV.2","ratioN1PMMoV.3",
                          "N1Load.1","N1Load.2","N1Load.3")
colnames(nd.epi) <- c("date",
                      "incidence","positivity",
                      "incidence.rollmean","positivity.rollmean")
colnames(nd.ww) <- colnames(elkhart.ww)[1:ncol(nd.ww)]
elkhart.epi[,date:=as.Date(date,format="%d-%b-%Y")]
nd.epi[,date:=as.Date(date,format="%d-%b-%Y")]
elkhart.ww[,date:=as.Date(date,format="%d-%b-%Y")]
nd.ww[,date:=as.Date(date,format="%d-%b-%Y")]

## ND analysis
par(mar=c(5.1,4.1,4.1,4.1))
plot(as.Date(nd.epi$date),nd.epi$incidence.rollmean,xlab="Date",type="l",
     ylab="Incidence (7-day rolling)",las=1,bty="n")
par(new=T)
plot(nd.ww$date,nd.ww$N1.1,col="red",xlab="",ylab="",xaxt="n",yaxt="n",pch=1,
     ylim=c(0,1.1*max(nd.ww$N1.1,nd.ww$N1.2,nd.ww$N1.3)),las=1,bty="n")
points(nd.ww$date,nd.ww$N1.2,col="red",pch=2)
points(nd.ww$date,nd.ww$N1.2,col="red",pch=3)
axis(4)
mtext("N1 GC per 100 ml",4,line=3)

plot(nd.ww$date,rowSums(nd.ww[,2:4])/3,ylim=c(0,600))

## fiddle with the positivity a bit
nd.epi[,total.tests:=100*incidence/positivity]
nd.epi$total.tests[7:nrow(nd.epi)] <- na.spline(nd.epi$total.tests[7:nrow(nd.epi)])
nd.epi[,positivity.rollmean.new:=incidence.rollmean/frollmean(total.tests,7)]
plot(nd.epi$positivity.rollmean/100, nd.epi$positivity.rollmean.new,
     ylim=c(0,max(c(nd.epi$positivity.rollmean/100, nd.epi$positivity.rollmean.new),na.rm=TRUE)))
cor(nd.epi$positivity.rollmean/100, nd.epi$positivity.rollmean.new,use="na")
abline(a=0,b=1)

## Cross correlations
nd.ww.N1 <- nd.ww[,1:4]
nd.ww.N1$N1.mean <- apply(nd.ww[,2:4],1,mean)
nd.ww.N1$N1.sd <- apply(nd.ww[,2:4],1,sd)
nd.ww.N1 <- nd.ww.N1[,-(2:4)]
nd.epi <- merge(nd.epi,nd.ww.N1,by="date",all=TRUE)

ccf(nd.epi$N1.mean,nd.epi$incidence.rollmean,
    na.action=na.pass,main="ND cross-correlation, incidence on N1")
polygon(c(-0.5,0.5,0.5,-0.5),c(-1,-1,1,1),
        col=adjustcolor("red",alpha.f=0.2),border=FALSE)

ccf(nd.epi$N1.mean,nd.epi$positivity.rollmean,
    na.action=na.pass,main="ND cross-correlation, incidence on N1")
polygon(c(-0.5,0.5,0.5,-0.5),c(-1,-1,1,1),
        col=adjustcolor("red",alpha.f=0.2),border=FALSE)

ccf(nd.epi$N1.mean,nd.epi$positivity.rollmean.new,
    na.action=na.pass,main="ND cross-correlation, incidence on N1")
polygon(c(-0.5,0.5,0.5,-0.5),c(-1,-1,1,1),
        col=adjustcolor("red",alpha.f=0.2),border=FALSE)

## What about comparing just to the surveillance tests?? Because it's less messy
## How can I download this data?
