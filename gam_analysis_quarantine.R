library(data.table)
library(mgcv)

data <- fread("./Sewage_Sampling_Fall_Data_20210507.csv")

## Create data
date.temp <- rep(data$CRU_COMPLETED_DATE,data$Total_Pos_Daily)
iso.off.binary.temp <- c()
for (ii  in 1:length(data$Total_Pos_Daily)){
    iso.off.binary.temp <- c(iso.off.binary.temp,
                             rep(c(1,0),
                                 c(data$Iso_Off_Campus[ii],data$Total_Pos_Daily[ii]-data$Iso_Off_Campus[ii])))
}
data.binary <- data.frame(list(Date=date.temp,
                               Iso_Off_Campus=as.logical(iso.off.binary.temp)))
setDT(data.binary)
data.binary[,time:=Date-min(Date)]
rm(date.temp, iso.off.binary.temp)

g.log <- gam(Iso_Off_Campus~s(time),data=data.binary,
             family=binomial,
             method = "REML")
summary(g.log)

gam.check(g.log)

pdf("./gam_output.pdf",width=10,height=5)
plot(g.log,trans=plogis,shift=coef(g.log)[1],
     seWithMean=TRUE,shade=TRUE,rug=FALSE,select=1,
     xlab="Date",ylim=c(0,1),las=1,yaxs="i",bty="n",xaxt="n",
     ylab="Proportion isolating off-campus")
axis(1,at=seq(0,100,20),labels=seq(min(data.binary$Date),by="20 days",length.out=6))
for (t in unique(data.binary$time)) {
    points(t,sum(data.binary[time==t]$Iso_Off_Campus) /
                length(data.binary[time==t]$Iso_Off_Campus))
}
dev.off()

prop.q <- plogis(predict(g.log, newdata=list(time=0:max(data.binary$time)), type="link"))
names(prop.q) <- seq(min(data.binary$Date),max(data.binary$Date),by="1 day")
save(g.log,prop.q,file="./prop_quarantine.RData")

## Create data
data[,Iso_Total:=Iso_Off_Campus+Iso_On_Campus]
date.temp <- rep(data$CRU_COMPLETED_DATE,data$Iso_Total)
iso.off.binary.temp <- c()
for (ii  in 1:nrow(data)){
    iso.off.binary.temp <- c(iso.off.binary.temp,
                             rep(c(1,0),
                                 c(data$Iso_Off_Campus[ii],data$Iso_On_Campus[ii])))
}
data.binary <- data.frame(list(Date=date.temp,
                               Iso_Off_Campus=as.logical(iso.off.binary.temp)))
setDT(data.binary)
data.binary[,time:=Date-min(Date)]
rm(date.temp, iso.off.binary.temp)

g.log <- gam(Iso_Off_Campus~s(time),data=data.binary,
             family=binomial,
             method = "REML")
summary(g.log)

gam.check(g.log)

pdf("./gam_output_iso_only.pdf",width=10,height=5)
plot(g.log,trans=plogis,shift=coef(g.log)[1],
     seWithMean=TRUE,shade=TRUE,rug=FALSE,select=1,
     xlab="Date",ylim=c(0,1),las=1,yaxs="i",bty="n",xaxt="n",
     ylab="Proportion isolating")
axis(1,at=seq(0,100,20),labels=seq(min(data.binary$Date),by="20 days",length.out=6))
for (t in unique(data.binary$time)) {
    points(t,sum(data.binary[time==t]$Iso_Off_Campus) /
                length(data.binary[time==t]$Iso_Off_Campus))
}
dev.off()

prop.q <- plogis(predict(g.log, newdata=list(time=0:max(data.binary$time)), type="link"))
names(prop.q) <- seq(min(data.binary$Date),max(data.binary$Date),by="1 day")
save(g.log,prop.q,file="./prop_quarantine_iso_only.RData")
