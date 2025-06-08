
setwd("./")

library(survsim)
library(dplyr)

paper_simu=1

n_start=2000 

fup=30

r01=0.01
r02=0.01
r03=0.01
r04=0.01

r12=0.02
r13=0.02
r14=0.02

r23_0_1=0.03
r24_0_1=0.03

r23_0=r23_0_1
r24_0=r23_0_1

r34_0=0.04
r34_01=r34_0
r34_02=r34_0
r34_012=r34_0


for (ps in 1:paper_simu){
  
  df0 <- crisk.sim(n=n_start, foltime=fup, nsit=4, 
                   ##competing events
                   dist.ev=c("weibull","weibull","weibull","weibull"),
                   anc.ev=c(1, 1, 1, 1),  ##anc.ev is the shape
                   beta0.ev=c(-log(r01),-log(r02),-log(r03),-log(r04)), 
                   
                   ##censor
                   dist.cens="weibull", 
                   anc.cens=1,
                   beta0.cens=600  ##censor rate is exp(-600) very small
                   )


df0 = subset(df0, select = c(nid,cause,time) )
  colnames(df0)
  
  df0<-df0 %>% 
    dplyr::rename(ID = nid,
                  status1=cause,
                  tte1=time)
  
  df0$Tstart1=0
  df0$Tstop1=df0$tte1
  
  df0[which(is.na(df0$status1)==T),]$status1<-0
  
  table(df0$status1)
  summary(df0$tte1)
  
  df0[1,]
  
  ##simulate the next state
  
  df1_before<-df0[which(df0$status1==1),]
  
  df1_n=nrow(df1_before)

    df1_after <- crisk.sim(n=df1_n, foltime=fup, 
                         nsit=3, dist.ev=c("weibull","weibull","weibull"),
                         anc.ev=c(1,1,1),  ##anc.ev is the shape
                         beta0.ev=c(-log(r12), -log(r13), -log(r14)), ##beta0.ev is the scale, 1/rate 
                         dist.cens="weibull", anc.cens=1,beta0.cens=600)
                         

  df1_after = subset(df1_after, select = c(nid,cause,time) )
  
  df1_after<-df1_after %>% dplyr:: rename(ID = nid,
                                          tte2=time)
  
  df1_after$Tstart2=0
  df1_after$Tstop2=df1_after$tte2
  
  df1_after$status2=df1_after$cause+1
  
  df1_after<-df1_after%>%
    mutate(
      status2 = case_when(
        is.na(status2)==T ~1 ,
        TRUE ~ status2
      ),
    )
  
  df1_after = subset(df1_after, select = -c(cause, ID) )
  
  table(df1_after$status2)
  
  colnames(df1_before)
  colnames(df1_after)
  
  ##randomly combine past state and current state
  df1 <- cbind(df1_before, df1_after)      
  
  df1[1,]
  
  #library(tidyverse)
  #mesa_all=df_list %>% reduce(full_join, by='ID')
  
  ##simulate y3
  
  df2_before1=df0[which(df0$status1==2),]
  
  df2_before2=df1[which(df1$status2==2),]
  
  colnames(df2_before1)
  
  colnames(df2_before2)
  
  df2_n0=nrow(df2_before1)
  
  df2_n01=nrow(df2_before2)
  
  df2_after_0 <- crisk.sim(n=df2_n0, foltime=fup, 
                           nsit=2, dist.ev=c("weibull","weibull"),
                           anc.ev=c(1,1),  ##anc.ev is the shape
                           beta0.ev=c(-log(r23_0), -log(r24_0)), ##beta0.ev is the scale, 1/rate 
                           dist.cens="weibull", anc.cens=1,beta0.cens=600)
  
  colnames(df2_after_0) =c("nid","cause","time","status","start","stop")
  
  df2_after_0
  
  df2_after_01 <- crisk.sim(n=df2_n01, foltime=fup, 
                            nsit=2, dist.ev=c("weibull","weibull"),
                            anc.ev=c(1,1),  ##anc.ev is the shape
                            beta0.ev=c(-log(r23_0_1), -log(r24_0_1)), ##beta0.ev is the scale, 1/rate 
                            dist.cens="weibull", anc.cens=1,beta0.cens=600)
  
  colnames(df2_after_01) =c("nid","cause","time","status","start","stop")
  
  df2_after_01
  
  df2_after_0 = subset(df2_after_0, select = c(cause,time) )
  
  df2_after_0<-df2_after_0 %>% dplyr:: rename(tte3=time)
  
  df2_after_0$Tstart3=0
  df2_after_0$Tstop3=df2_after_0$tte3
  
  df2_after_0$status3=df2_after_0$cause+2
  
  df2_after_0<-df2_after_0%>%
    mutate(
      status3 = case_when(
        is.na(status3)==T ~2 ,
        TRUE ~ status3
      ),
    )
  
  df2_after_01 = subset(df2_after_01, select = c(cause,time) )
  
  df2_after_01<-df2_after_01 %>% dplyr:: rename(tte3=time)
  
  df2_after_01$Tstart3=0
  df2_after_01$Tstop3=df2_after_01$tte3
  
  df2_after_01$status3=df2_after_01$cause+2
  
  df2_after_01<-df2_after_01%>%
    mutate(
      status3 = case_when(
        is.na(status3)==T ~2 ,
        TRUE ~ status3
      ),
    )
  
  df2_after_0 = subset(df2_after_0, select = -c(cause) )
  
  df2_after_01 = subset(df2_after_01, select = -c(cause) )
  
  df2_0 = cbind(df2_before1, df2_after_0 )
  colnames(df2_0)
  
  df2_01 = cbind(df2_before2, df2_after_01 )
  colnames(df2_01)
  
  df2 <- data.frame(bind_rows(df2_0, df2_01))
  
  summary(df2)
  
  df2[1:2,]
  
  nrow(df2)
  nrow(df2_before2)+nrow(df2_before1)
  
  summary(df2[which(df2$status1==2),]$tte1)
  summary(df2[which(df2$status1==2),]$tte2)
  
  summary(df2[which(df2$status2==2),]$tte1)
  summary(df2[which(df2$status2==2),]$tte2)
  
  ##simulate y4
  
  df3_before_0=df0[which(df0$status1==3),]
  ##status1=3
  
  df3_before_01=df1[which(df1$status2==3),]  
  ##status1=1, status2=3
  
  df3_before_02=df2[which(df2$status1==2 & is.na(df2$status2)==T & df2$status3==3),]  
  ##status1=2, status2=NA, status3=3
  
  df3_before_012=df2[which(df2$status1==1 & df2$status2==2 & df2$status3==3),]  
  ##status1=1, status2=2, status3=3  
  
  colnames(df3_before_0)
  
  colnames(df3_before_01)
  
  colnames(df3_before_02)
  
  colnames(df3_before_012)
  
  
  df3_after_0 <- crisk.sim(n=nrow(df3_before_0), foltime=fup, 
                           nsit=1, dist.ev=c("weibull"),
                           anc.ev=c(1),  ##anc.ev is the shape
                           beta0.ev=c(-log(r34_0)), ##beta0.ev is the scale, 1/rate 
                           dist.cens="weibull", anc.cens=1,beta0.cens=600)
  
  colnames(df3_after_0) =c("nid","cause","time","status","start","stop")
  
  
  df3_after_01 <- crisk.sim(n=nrow(df3_before_01), foltime=fup, 
                            nsit=1, dist.ev=c("weibull"),
                            anc.ev=c(1),  ##anc.ev is the shape
                            beta0.ev=c(-log(r34_01)), ##beta0.ev is the scale, 1/rate 
                            dist.cens="weibull", anc.cens=1,beta0.cens=600)

  
  colnames(df3_after_01) =c("nid","cause","time","status","start","stop")
  
  
  df3_after_02 <- crisk.sim(n=nrow(df3_before_02), foltime=fup, 
                            nsit=1, dist.ev=c("weibull"),
                            anc.ev=c(1),  ##anc.ev is the shape
                            beta0.ev=c(-log(r34_02)), ##beta0.ev is the scale, 1/rate 
                            dist.cens="weibull", anc.cens=1,beta0.cens=600)

  
  colnames(df3_after_02) =c("nid","cause","time","status","start","stop")
  
  
  df3_after_012 <- crisk.sim(n=nrow(df3_before_012), foltime=fup, 
                             nsit=1, dist.ev=c("weibull"),
                             anc.ev=c(1),  ##anc.ev is the shape
                             beta0.ev=c(-log(r34_012)), ##beta0.ev is the scale, 1/rate 
                             dist.cens="weibull", anc.cens=1,beta0.cens=600)
  
  
  colnames(df3_after_012) =c("nid","cause","time","status","start","stop")
  
  
  df3_after_0 = subset(df3_after_0, select = c(nid,cause,time) )
  
  df3_after_0<-df3_after_0 %>% dplyr:: rename(ID = nid,
                                              tte4=time)
  
  df3_after_0$Tstart4=0
  df3_after_0$Tstop4=df3_after_0$tte4
  
  df3_after_0$status4=df3_after_0$cause+3
  
  df3_after_0<-df3_after_0%>%
    mutate(
      status4 = case_when(
        is.na(status4)==T ~3 ,
        TRUE ~ status4
      ),
    )
  
  table(df3_after_0$status4)
  
  
  df3_after_0 = subset(df3_after_0, select = -c(cause, ID) )
  
  df3_after_01 = subset(df3_after_01, select = c(nid,cause,time) )
  
  df3_after_01<-df3_after_01 %>% dplyr:: rename(ID = nid,
                                                tte4=time)
  
  df3_after_01$Tstart4=0
  df3_after_01$Tstop4=df3_after_01$tte4
  
  df3_after_01$status4=df3_after_01$cause+3
  
  df3_after_01<-df3_after_01%>%
    mutate(
      status4 = case_when(
        is.na(status4)==T ~3 ,
        TRUE ~ status4
      ),
    )
  
  table(df3_after_01$status4)
  
  
  df3_after_01 = subset(df3_after_01, select = -c(cause, ID) )
  
  df3_after_02 = subset(df3_after_02, select = c(nid,cause,time) )
  
  df3_after_02<-df3_after_02 %>% dplyr:: rename(ID = nid,
                                                tte4=time)
  
  df3_after_02$Tstart4=0
  df3_after_02$Tstop4=df3_after_02$tte4
  
  df3_after_02$status4=df3_after_02$cause+3
  
  df3_after_02<-df3_after_02%>%
    mutate(
      status4 = case_when(
        is.na(status4)==T ~3 ,
        TRUE ~ status4
      ),
    )
  
  table(df3_after_02$status4)
  
  
  df3_after_02 = subset(df3_after_02, select = -c(cause, ID) )
  
  df3_after_012 = subset(df3_after_012, select = c(nid,cause,time) )
  
  df3_after_012<-df3_after_012 %>% dplyr:: rename(ID = nid,
                                                  tte4=time)
  
  df3_after_012$Tstart4=0
  df3_after_012$Tstop4=df3_after_012$tte4
  
  df3_after_012$status4=df3_after_012$cause+3
  
  df3_after_012<-df3_after_012%>%
    mutate(
      status4 = case_when(
        is.na(status4)==T ~3 ,
        TRUE ~ status4
      ),
    )
  
  table(df3_after_012$status4)
  
  
  df3_after_012 = subset(df3_after_012, select = -c(cause, ID) )
  
  df3_0 = cbind(df3_before_0, df3_after_0 )
  df3_01 = cbind(df3_before_01, df3_after_01 )
  df3_02 = cbind(df3_before_02, df3_after_02 )
  df3_012 = cbind(df3_before_012, df3_after_012 )
  
  df3 <- data.frame(bind_rows(df3_0, df3_01,df3_02,df3_012))
  colnames(df3)
  colnames(df3_0)
  colnames(df3_01)
  colnames(df3_02)
  colnames(df3_012)
  
  nrow(df3_0)
  nrow(df3_01)
  nrow(df3_02)
  nrow(df3_012)
  
  nrow(df3_0)+nrow(df3_01)
  nrow(df3_0)+nrow(df3_02)
  
  df3[1,]

  ##combine the data
  
  colnames(df0)
  colnames(df1)
  colnames(df2)
  colnames(df3)
  
  df0_for_c<-data.frame(df0)
  df1_for_c<-subset(data.frame(df1), select=-c(status1,tte1, Tstart1,Tstop1))
  df2_for_c<-subset(data.frame(df2), select=-c(status1,tte1, Tstart1,Tstop1,
                                              status2,tte2, Tstart2,Tstop2))
  df3_for_c<-subset(data.frame(df3), select=-c(status1,tte1, Tstart1,Tstop1,
                                               status2,tte2, Tstart2,Tstop2,
                                               status3,tte3, Tstart3,Tstop3))
  
  
  df0_for_c$marker1=1
  df1_for_c$marker2=1
  df2_for_c$marker3=1
  df3_for_c$marker4=1
  
  colnames(df0_for_c)
  colnames(df1_for_c)
  colnames(df2_for_c)
  colnames(df3_for_c)
  
  library(tidyverse)
  df_list <- list(df0_for_c, df1_for_c,df2_for_c,df3_for_c)
  
  df_all=Reduce(function(x, y) merge(x, y,by="ID",all.x=T), df_list)
  
  
  write.table(df_all, "./multi_state_data.txt")
  
}








