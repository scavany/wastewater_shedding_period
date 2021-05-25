## Load packages
if(!require(pacman)){install.packages('pacman'); library(pacman)}
p_load(tidyverse,data.table,viridis,grDevices)

## covid.data
covid.data <- fread("./COVID_Surveillance_Impact_20200810-20201204.csv")[Description=="student"]
covid.data.total <- covid.data[Test_Type=="Symptomatic",
                               .(totalPos=sum(Positive)),
                               by=CRU_COMPLETED_DATE]
    
## ww data
ww.data <- read_csv("NDCampus_ddPCRdata_cleanedrawdata_0s.csv") %>%
    mutate(Date=as.Date(Date,"%m/%d/%Y")) %>%
    select(-WW_flow) %>%
    ## select(Date,N1R1RC,N1R2RC,N1R3RC) %>%
    drop_na()
N1.med <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,median)
N1.mean <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,mean)
N1.min <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,min)
N1.max <- apply(ww.data[,c("N1R1RC","N1R2RC","N1R3RC")],1,max)

all.dates <- seq(min(c(ww.data$Date,covid.data.total$CRU_COMPLETED_DATE)),
                 max(c(ww.data$Date,covid.data.total$CRU_COMPLETED_DATE)),
                 "1 day")

cases.filled <- rep(NA,length(all.dates))
wwrna.filled <- rep(NA,length(all.dates))
cases.filled[all.dates %in% covid.data.total$CRU_COMPLETED_DATE] = covid.data.total$totalPos
wwrna.filled[all.dates %in% ww.data$Date] = N1.mean
## ww.in.inc <- (ww.data$Date %in% covid.data.total$CRU_COMPLETED_DATE)
## inc.in.ww <- (covid.data.total$CRU_COMPLETED_DATE %in% ww.data$Date)



## Plot it
pdf("cross-correlation_plot.pdf")
ccf(wwrna.filled,cases.filled,
    na.action=na.pass,
    main="",bty="n",
    ylab="Cross-correlation function")
## polygon(c(7.5,10.5,10.5,7.5),c(-1,-1,1,1),
##         col=adjustcolor("red",alpha.f=0.2),border=FALSE)
dev.off()
