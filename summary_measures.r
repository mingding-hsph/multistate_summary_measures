
library(riskRegression, verbose = FALSE, quietly = TRUE)
library(survival)
library(Hmisc)
library(rms)
library(splines)
library(ggplot2)
#install.packages("splitstackshape")
library(splitstackshape) ##expandrows function

setwd("./")

simu_data_all<-read.table("./multi_state_data.txt")

############Extract baseline data

simu_data_base <-select(simu_data_all, ID)

###make multiple copies of each ID by person-time (such as age)
###the purpose is for prediction use to estimate the age-specific transiton rates 

break_interval=2
fup=30
analysis_interval=fup/break_interval
tm=analysis_interval

simu_data_long <- data.frame(expandRows(simu_data_base, count=tm, count.is.col=F))

for (m in 0:(nrow(simu_data_base)-1)){
  for (i in (tm*m+1):(tm*m+tm)){
    if (i==1+tm*m) {
      simu_data_long$interval[i]=1
    } else if (i>1+tm*m & i<=tm+tm*m)  {simu_data_long$interval[i]=1+i-tm*m-1
    } 
  }
}


###Estimate age-specific TR and estimate MALY and disease path 

TR_var <- list()  
SOP_var <- list()  
SSLY<- list() 
LE<- list() 
LE_total<- list() 
MALY<- list() 
MALY_total<- list() 
path_total<- list() 
path_year<- list() 

I=diag(nrow=9)

var_loop=20     ##bootstrap loops

##start of bootstrap

for (simu in 1:var_loop){
  
simu_data_all_boot<-simu_data_all[sample(1:nrow(simu_data_all), size=nrow(simu_data_all), replace=TRUE),]

##subset datasets by disease states

##starting from state 0 

data0<-simu_data_all_boot[which(simu_data_all_boot$marker1==1),]
data0$tte<-data0$tte1
data0$S0_y<-data0$status1

data0<-data0%>%
  mutate(
    S0_ylong = case_when(
      S0_y>=1 ~ 1,
      S0_y==0 ~ 0,
      TRUE ~ NA
    )  
  )

table(data0$S0_ylong)  
table(data0$S0_y)

summary(data0$S0_ylong)  
summary(data0$S0_y)

##starting from state 1 

data1<-simu_data_all_boot[which(simu_data_all_boot$marker2==1),]
data1$tte<-data1$tte2
summary(data1$tte2)

data1$S1_y<-data1$status2

data1<-data1%>%
  mutate(
    S1_ylong = case_when(
      S1_y>=2 ~ 1,
      S1_y==1 ~ 0,
      TRUE ~ NA
    )
  )

table(data1$S1_ylong)  
table(data1$S1_y)


##starting from state 2

data2<-simu_data_all_boot[which(simu_data_all_boot$marker3==1),]

data2 <- data2 %>%
  mutate(tte = tte3,
         
         S2_y =status3,
         
         S2_ylong = case_when(
           S2_y>=3 ~ 1,
           S2_y==2 ~ 0,
           TRUE ~ NA
         )
  )

summary(data2$tte)

table(data2$S2_ylong)  
table(data2$S2_y)

##starting from state 3

data3<-simu_data_all_boot[which(simu_data_all_boot$marker4==1),]

data3 <- data3 %>%
  mutate(tte = tte4,
         
         S3_y =status4,
         
         S3_ylong = case_when(
           S3_y>=4 ~ 1,
           S3_y==3 ~ 0,
           TRUE ~ NA
         )
  )

##change data structure from wide to long

library(survival)

data0_t=round(max(data0$tte), digits = 0)+1
data1_t=round(max(data1$tte), digits = 0)+1
data2_t=round(max(data2$tte), digits = 0)+1
data3_t=round(max(data3$tte), digits = 0)+1

data0_long <- survSplit(Surv(tte, S0_ylong)~.,
                        data = data0,
                        cut = c(0:(data0_t/break_interval))*break_interval,   
                        episode="interval",
                        start="tstart",
                        end = "tstop")
data0_long=data.frame(data0_long)
data0_long$interval=data0_long$interval-1

table(data0_long$interval)

data1_long <- survSplit(Surv(tte, S1_ylong)~.,
                        data = data1,
                        cut = c(0:(data1_t/break_interval))*break_interval,   
                        episode="interval",
                        start="tstart",
                        end = "tstop")

data1_long=data.frame(data1_long)
data1_long$interval=data1_long$interval-1

data2_long <- survSplit(Surv(tte, S2_ylong)~.,
                        data = data2,
                        cut = c(0:(data2_t/break_interval))*break_interval,   
                        episode="interval",
                        start="tstart",
                        end = "tstop")

data2_long=data.frame(data2_long)
data2_long$interval=data2_long$interval-1

data3_long <- survSplit(Surv(tte, S3_ylong)~.,
                        data = data3,
                        cut = c(0:(data3_t/break_interval))*break_interval,   
                        episode="interval",
                        start="tstart",
                        end = "tstop")

data3_long=data.frame(data3_long)
data3_long$interval=data3_long$interval-1


table(data0_long$interval)
table(data1_long$interval)
table(data2_long$interval)
table(data3_long$interval)


##define competing outcomes in each data subset

data0_long$S0_mylong<-0
data0_long[which(data0_long$S0_y==1 & data0_long$S0_ylong==1),]$S0_mylong<-1
data0_long[which(data0_long$S0_y==2 & data0_long$S0_ylong==1),]$S0_mylong<-2
data0_long[which(data0_long$S0_y==3 & data0_long$S0_ylong==1),]$S0_mylong<-3
data0_long[which(data0_long$S0_y==4 & data0_long$S0_ylong==1),]$S0_mylong<-4


data1_long$S1_mylong<-0
data1_long[which(data1_long$S1_y==2 & data1_long$S1_ylong==1),]$S1_mylong<-2
data1_long[which(data1_long$S1_y==3 & data1_long$S1_ylong==1),]$S1_mylong<-3
data1_long[which(data1_long$S1_y==4 & data1_long$S1_ylong==1),]$S1_mylong<-4


data2_long$S2_mylong<-0
data2_long[which(data2_long$S2_y==3 & data2_long$S2_ylong==1),]$S2_mylong<-3
data2_long[which(data2_long$S2_y==4 & data2_long$S2_ylong==1),]$S2_mylong<-4

data3_long$S3_mylong<-0
data3_long[which(data3_long$S3_y==4 & data3_long$S3_ylong==1),]$S3_mylong<-4

data0_long = subset(data0_long, select = c(ID, S0_mylong, tstart, tstop, interval, status1, tte1) )
data1_long = subset(data1_long, select = c(ID, S1_mylong, tstart, tstop, interval, status1, tte1, status2, tte2) )
data2_long = subset(data2_long, select = c(ID, S2_mylong, tstart, tstop, interval, status1, tte1, status2, tte2, status3, tte3) )
data3_long = subset(data3_long, select = c(ID, S3_mylong, tstart, tstop, interval, status1, tte1, status2, tte2, status3, tte3, status4, tte4) )

library(dplyr)


data0_long <- data0_long %>%
  mutate(
    
    tte_cmr=tstop-tstart,
    
    S0_y_1 = case_when(
      S0_mylong==1 ~ 1,
      TRUE ~ 0
    ),
    S0_y_2 = case_when(
      S0_mylong==2 ~ 1,
      TRUE ~ 0
    ),
    S0_y_3 = case_when(
      S0_mylong==3 ~ 1,
      TRUE ~ 0
    ),
    S0_y_4 = case_when(
      S0_mylong==4 ~ 1,
      TRUE ~ 0
    )
  )


data1_long<-data1_long %>%
  mutate(
    
    tte_cmr=tstop-tstart,
    
    S1_y_2 = case_when(
      S1_mylong==2 ~ 1,
      TRUE ~ 0
    ),
    
    S1_y_3 = case_when(
      S1_mylong==3 ~ 1,
      TRUE ~ 0
    ),
    
    
    S1_y_4 = case_when(
      S1_mylong==4 ~ 1,
      TRUE ~ 0
    ) 
  )

data2_long<-data2_long %>%
  mutate(
    
    tte_cmr=tstop-tstart,
    
    S2_y_3 = case_when(
      S2_mylong==3 ~ 1,
      TRUE ~ 0
    ),
    
    S2_y_4 = case_when(
      S2_mylong==4 ~ 1,
      TRUE ~ 0
    ),
    
    sub2=case_when (
      status1==1 & status2==2 ~ "htn-chd",
      status1==2 & is.na(status2)==T ~ "well-chd"),
  )


data3_long<-data3_long %>%
  mutate(
    
    tte_cmr=tstop-tstart,
    
    S3_y_4=case_when(
      S3_mylong==4 ~ 1,
      TRUE~0
    ),
    
    sub3=case_when (
      status1==1 & status2==2 & status3==3 ~ "htn-chd-chf",
      status1==1 & status2==3 & is.na(status3)==T ~ "htn-chf",
      status1==2 & is.na(status2)==T & status3==3 ~ "chd-chf",
      status1==3 ~ "chf"
    )
    
  )

  ##fit cause-specific Cox models within each data subset
cfit0 <- CSC(formula = list(Hist(tte_cmr,S0_mylong) ~ rcs(interval,3)), data = data0_long)

cfit1 <- CSC(formula = list(Hist(tte_cmr,S1_mylong) ~rcs(interval,3)), data = data1_long)

cfit2 <- CSC(formula = list(Hist(tte_cmr,S2_mylong) ~ rcs(interval,3)+as.factor(sub2)), data = data2_long)

cfit3 <- coxph(Surv(tte_cmr,S3_y_4) ~ rcs(interval,3)+as.factor(sub3), data = data3_long, x = TRUE)

##for models cfit2 and cfit3, test the markov assumption in the original population before running bootstrap, 
## by including interaction terms of time and past states into the model, 
## apply likelihood ratio test, with a p<0.05 indicating that the markov assumption is violated

#estimate time(age)-specific transition rates
  
##data subset 0
pfit01 <- predict(cfit0, newdata = simu_data_long, cause = 1, times =1, se = F, keep.newdata = FALSE)
pfit02 <- predict(cfit0, newdata = simu_data_long, cause = 2, times =1, se = F, keep.newdata = FALSE)
pfit03 <- predict(cfit0, newdata = simu_data_long, cause = 3, times =1, se = F, keep.newdata = FALSE)
pfit04 <- predict(cfit0, newdata = simu_data_long, cause = 4, times =1, se = F, keep.newdata = FALSE)

predict0<-data.frame(pfit01$absRisk,  
pfit02$absRisk, 
pfit03$absRisk, 
pfit04$absRisk)

colnames(predict0)=c("r01", "r02", "r03", "r04")

##data subset 1

pfit12 <- predict(cfit1, newdata = simu_data_long, cause = 2, times =1, se = F, keep.newdata = FALSE)
pfit13 <- predict(cfit1, newdata = simu_data_long, cause = 3, times =1, se = F, keep.newdata = FALSE)
pfit14 <- predict(cfit1, newdata = simu_data_long, cause = 4, times =1, se = F, keep.newdata = FALSE)

predict1<-data.frame(
                pfit12$absRisk,  
                pfit13$absRisk,  
                pfit14$absRisk)

colnames(predict1)=c("r12","r13","r14")

##data subset 2

new_data2_1<-simu_data_long
new_data2_1$sub2="htn-chd"

new_data2_2<-simu_data_long
new_data2_2$sub2="well-chd"

new_data2=rbind(new_data2_1,new_data2_2)

pfit23 <- predict(cfit2, newdata = new_data2, cause = 3, times =1, se = F, keep.newdata = FALSE)
pfit24 <- predict(cfit2, newdata = new_data2, cause = 4, times =1, se = F, keep.newdata = FALSE)

predict2<-data.frame(new_data2$sub2,
                pfit23$absRisk, 
                pfit24$absRisk)

colnames(predict2)=c("sub2","r23", "r24")

is.numeric(predict2$r23)

predict2_0=predict2[which(predict2$sub2=="well-chd"),]

predict2_1=predict2[which(predict2$sub2=="htn-chd"),]

##data subset 3

new_data3_1<-simu_data_long
new_data3_1$sub3="chd-chf"

new_data3_2<-simu_data_long
new_data3_2$sub3="chf"

new_data3_3<-simu_data_long
new_data3_3$sub3="htn-chd-chf"

new_data3_4<-simu_data_long
new_data3_4$sub3="htn-chf"

new_data3=rbind(new_data3_1,new_data3_2, new_data3_3, new_data3_4)


##in case the data from this dataset is empty, use the following code 

if (nrow(as.matrix(table(new_data3$sub3)))==4){

pfit34 <- predictCox(cfit3, newdata=new_data3, times=1, se = F)

} else {

  pfit34$survival=matrix(NA, nrow=nrow(new_data3), ncol=1)

}


predict3<-data.frame(new_data3$sub3,
                pfit34$survival)

colnames(predict3)=c("sub3","r34_survival")


predict3$r34=1-predict3$r34_survival


predict3_0=predict3[which(predict3$sub3=="chf"),]

predict3_01=predict3[which(predict3$sub3=="htn-chf"),]

predict3_02=predict3[which(predict3$sub3=="chd-chf"),]

predict3_012=predict3[which(predict3$sub3=="htn-chd-chf"),]


##estimate MALY and disease path for each individual 

  for (i in 1:nrow(simu_data_base)) {
    
##construct the transition rate matrix at each time 
    
    for (t in 1:tm) {    

      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=matrix(0, nrow=9, ncol=9)
      
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,2]=predict0$r01[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,3]= predict0$r02[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,5]= predict0$r03[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,9]= predict0$r04[t+(i-1)*tm]
      
      
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][2,4]=predict1$r12[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][2,6]=predict1$r13[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][2,9]=predict1$r14[t+(i-1)*tm]
      
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][3,7]=predict2_0$r23[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][3,9]=predict2_0$r24[t+(i-1)*tm]
      
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][4,8]=predict2_1$r23[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][4,9]=predict2_1$r24[t+(i-1)*tm]
      
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][5,9]=predict3_0$r34[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][6,9]=predict3_01$r34[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][7,9]=predict3_02$r34[t+(i-1)*tm]
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][8,9]=predict3_012$r34[t+(i-1)*tm]
      
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,1]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,2:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][2,2]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][2,3:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][3,3]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][3,4:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][4,4]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][4,5:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][5,5]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][5,6:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][6,6]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][6,7:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][7,7]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][7,8:9])
      TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][8,8]=-sum(TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][8,9:9])
      
    }
    
    #state occupational probability at the start of follow up
    
    SOP_var_0=c(1, 0, 0, 0, 0, 0, 0, 0, 0)
  
    ##state occupational probability at the end of the first year t=1
    
    SOP_var[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=SOP_var_0%*%(I+TR_var[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]])
    
    for (t in 2:tm) {
      SOP_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=SOP_var[[t-1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]%*%(I+TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]])
    }
    
    
    #as to death (S4) before the end of follow up: 
    ###if death (S4) probability>0.95, we assume death and SOP and TR are coded as 0 after death, which will not affect the estimation of MALY and PATH
    
    for (t in 1:tm) {
      if (SOP_var[[t +(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][9]>=0.95){
        for (t_in_t in 1:(tm-t)) {
          SOP_var[[t+t_in_t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=c(0, 0, 0, 0, 0, 0, 0, 0, 0)
          TR_var[[t+t_in_t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=matrix(0, nrow=9, ncol=9)
        }
      }
    }
    
    ##path at the end of year 1

    path_year[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=matrix(0, nrow=1, ncol=9)
    
    #path equals to the probability of being transitioned to death from each state at start of year 1 to end of year 1
    
    for (k_s in 1:9){
      path_year[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,k_s]=SOP_var_0[k_s]%*%TR_var[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][k_s,9]
    }
    
    
    for (t in 2:tm) {
      path_year[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=matrix(0, nrow=1, ncol=9)
      
      for (k_s in 1:9){
        path_year[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,k_s]=SOP_var[[t-1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][k_s]%*%TR_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][k_s,9]
      }
      
    }
    
##SSLY and MALY
    
    SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=SOP_var[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]
    
##assign weight to each disease, here we used the weight of heart disease as an example 
    
DW1=1                    #healthy
DW2=1-0.049              #at metabolic risk
DW3=1-0.072              #CHD without metabolic risk factors
DW4=(1-0.072)*(1-0.049)  #multimorbidity: CHD with metabolic risk factors
DW5=1-0.074              #Heart failure without metabolic risk factors or CHD
DW6=(1-0.074)*(1-0.049)  #multimorbidity: Heart failure with metabolic risk factors
DW7=(1-0.074)*(1-0.072)  #multimorbidity: Heart failure with CHD
DW8=(1-0.074)*(1-0.072)*(1-0.049)  #multimorbidity: Heart failure with metabolic risk factors and CHD

    MALY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,1]*DW1+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,2]*DW2+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,3]*DW3+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,4]*DW4+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,5]*DW5+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,6]*DW6+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,7]*DW7+
      SSLY[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,8]*DW8
    
    path_total[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]= path_year[[1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]
    
    for (t in 2:tm) {
      
      SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=SSLY[[t-1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]+SOP_var[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]
      
      MALY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,1]*DW1+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,2]*DW2+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,3]*DW3+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,4]*DW4+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,5]*DW5+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,6]*DW6+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,7]*DW7+
        SSLY[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,8]*DW8
      
      path_total[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]=path_total[[t-1+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]+ path_year[[t+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]
      
    }     
    
    
  }
  
}


###End of bootstrap


##organize the results

SSLY_1=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_2=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_3=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_4=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_5=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_6=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_7=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_8=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
SSLY_9=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))

MALY_1=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))

##path 1:9 calculates probabilities of death transited from 
##S0, S0→S1, S0→S2, S0→S1→S2, S0→S3, S0→S1→S3, S0→S2→S3, S0→S1→S2→S3, and S4, respectively. 

path_1=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_2=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_3=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_4=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_5=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_6=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_7=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_8=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))
path_9=matrix(NA, nrow=var_loop,ncol=nrow(simu_data_base))


SSLY_1_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_2_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_3_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_4_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_5_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_6_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_7_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_8_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
SSLY_9_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)

MALY_1_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)

path_1_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_2_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_3_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_4_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_5_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_6_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_7_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_8_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)
path_9_pop=matrix(NA, nrow=var_loop*nrow(simu_data_base),ncol=1)


for (simu in 1:var_loop) {
  
  for (i in 1:nrow(simu_data_base)) {
    
    ##individual level
    
    SSLY_1[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,1]
    SSLY_2[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,2]
    SSLY_3[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,3]
    SSLY_4[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,4]
    SSLY_5[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,5]
    SSLY_6[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,6]
    SSLY_7[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,7]
    SSLY_8[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,8]
    SSLY_9[simu,i]=SSLY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,9]
    
    MALY_1[simu,i]=MALY[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]]

    path_1[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,1]
    path_2[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,2]
    path_3[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,3]
    path_4[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,4]
    path_5[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,5]
    path_6[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,6]
    path_7[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,7]
    path_8[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,8]
    path_9[simu,i]=path_total[[tm+(i-1)*tm+(simu-1)*tm*nrow(simu_data_base)]][1,9]
    
  }
  
}


SSLY_1_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_2_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_3_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_4_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_5_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_6_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_7_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_8_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
SSLY_9_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))

MALY_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))

path_1_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_2_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_3_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_4_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_5_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_6_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_7_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_8_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))
path_9_final=matrix(NA, ncol=3,nrow=nrow(simu_data_base))


for (h in 1:nrow(simu_data_base)) {
  
  ##individual level
  
  SSLY_1_final[h,]=quantile(SSLY_1[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_2_final[h,]=quantile(SSLY_2[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_3_final[h,]=quantile(SSLY_3[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_4_final[h,]=quantile(SSLY_4[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_5_final[h,]=quantile(SSLY_5[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_6_final[h,]=quantile(SSLY_6[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_7_final[h,]=quantile(SSLY_7[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_8_final[h,]=quantile(SSLY_8[,h], probs = c(0.025, 0.5, 0.975))
  SSLY_9_final[h,]=quantile(SSLY_9[,h], probs = c(0.025, 0.5, 0.975))
  
  
  MALY_final[h,]=quantile(MALY_1[,h], probs = c(0.025, 0.5, 0.975))
  
  path_1_final[h,]=quantile(path_1[,h], probs = c(0.025, 0.5, 0.975))
  path_2_final[h,]=quantile(path_2[,h], probs = c(0.025, 0.5, 0.975))
  path_3_final[h,]=quantile(path_3[,h], probs = c(0.025, 0.5, 0.975))
  path_4_final[h,]=quantile(path_4[,h], probs = c(0.025, 0.5, 0.975))
  path_5_final[h,]=quantile(path_5[,h], probs = c(0.025, 0.5, 0.975))
  path_6_final[h,]=quantile(path_6[,h], probs = c(0.025, 0.5, 0.975))
  path_7_final[h,]=quantile(path_7[,h], probs = c(0.025, 0.5, 0.975))
  path_8_final[h,]=quantile(path_8[,h], probs = c(0.025, 0.5, 0.975))
  path_9_final[h,]=quantile(path_9[,h], probs = c(0.025, 0.5, 0.975))
  
}


##population level

SSLY_1_pop=data.frame(y=unlist(data.frame(SSLY_1)))  ##into one column
SSLY_2_pop=data.frame(y=unlist(data.frame(SSLY_2))) 
SSLY_3_pop=data.frame(y=unlist(data.frame(SSLY_3))) 
SSLY_4_pop=data.frame(y=unlist(data.frame(SSLY_4))) 
SSLY_5_pop=data.frame(y=unlist(data.frame(SSLY_5))) 
SSLY_6_pop=data.frame(y=unlist(data.frame(SSLY_6))) 
SSLY_7_pop=data.frame(y=unlist(data.frame(SSLY_7))) 
SSLY_8_pop=data.frame(y=unlist(data.frame(SSLY_8))) 
SSLY_9_pop=data.frame(y=unlist(data.frame(SSLY_9))) 

MALY_1_pop=data.frame(y=unlist(data.frame(MALY_1)))


##population level by baseline status

simu_data_long_all=cbind(t(path_1), t(path_2),t(path_3),
               t(path_4), t(path_5),t(path_6),
               t(path_7), t(path_8), t(path_9),
               simu_data_base)


path_1_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,(1:var_loop)])))  
path_2_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((1*var_loop+1):(2*var_loop))])))  
path_3_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((2*var_loop+1):(3*var_loop))])))   
path_4_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((3*var_loop+1):(4*var_loop))])))  
path_5_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((4*var_loop+1):(5*var_loop))])))   
path_6_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((5*var_loop+1):(6*var_loop))])))  
path_7_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((6*var_loop+1):(7*var_loop))])))  
path_8_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((7*var_loop+1):(8*var_loop))])))  
path_9_pop=data.frame(y=unlist(data.frame(simu_data_long_all[,((8*var_loop+1):(9*var_loop))])))  


##population level, summarize

##NOTE: exclude 0 values, as 0 indicate participants never stay in the state

SSLY_1_pop[,1][which(SSLY_1_pop[,1]==0)]<-NA
SSLY_2_pop[,1][which(SSLY_2_pop[,1]==0)]<-NA
SSLY_3_pop[,1][which(SSLY_3_pop[,1]==0)]<-NA
SSLY_4_pop[,1][which(SSLY_4_pop[,1]==0)]<-NA
SSLY_5_pop[,1][which(SSLY_5_pop[,1]==0)]<-NA
SSLY_6_pop[,1][which(SSLY_6_pop[,1]==0)]<-NA
SSLY_7_pop[,1][which(SSLY_7_pop[,1]==0)]<-NA
SSLY_8_pop[,1][which(SSLY_8_pop[,1]==0)]<-NA
SSLY_9_pop[,1][which(SSLY_9_pop[,1]==0)]<-NA


SSLY_1_final_pop=quantile(SSLY_1_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_2_final_pop=quantile(SSLY_2_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_3_final_pop=quantile(SSLY_3_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_4_final_pop=quantile(SSLY_4_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_5_final_pop=quantile(SSLY_5_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_6_final_pop=quantile(SSLY_6_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_7_final_pop=quantile(SSLY_7_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_8_final_pop=quantile(SSLY_8_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)
SSLY_9_final_pop=quantile(SSLY_9_pop[,1], probs = c(0.025, 0.5, 0.975), na.rm=T)

MALY_final_pop=quantile(MALY_1_pop[,1], probs = c(0.025, 0.5, 0.975))

path_1_final_pop=quantile(path_1_pop[,1], probs = c(0.025, 0.5, 0.975))
path_2_final_pop=quantile(path_2_pop[,1], probs = c(0.025, 0.5, 0.975))
path_3_final_pop=quantile(path_3_pop[,1], probs = c(0.025, 0.5, 0.975))
path_4_final_pop=quantile(path_4_pop[,1], probs = c(0.025, 0.5, 0.975))
path_5_final_pop=quantile(path_5_pop[,1], probs = c(0.025, 0.5, 0.975))
path_6_final_pop=quantile(path_6_pop[,1], probs = c(0.025, 0.5, 0.975))
path_7_final_pop=quantile(path_7_pop[,1], probs = c(0.025, 0.5, 0.975))
path_8_final_pop=quantile(path_8_pop[,1], probs = c(0.025, 0.5, 0.975))
path_9_final_pop=quantile(path_9_pop[,1], probs = c(0.025, 0.5, 0.975))



##final tables for each individual 

SSLY_table=cbind(SSLY_1_final,
                 SSLY_2_final,
                 SSLY_3_final,
                 SSLY_4_final,
                 SSLY_5_final,
                 SSLY_6_final,
                 SSLY_7_final,
                 SSLY_8_final,
                 SSLY_9_final)

colnames(SSLY_table)=c("S0_ll","S0_mean","S0_ul",
                       "S1_ll", "S1_mean", "S1_ul", 
                       "S2_0_ll", "S2_0_mean", "S2_0_ul",
                       "S2_1_ll", "S2_1_mean", "S2_1_ul",
                       "S3_0_ll", "S3_0_mean", "S3_0_ul",
                       "S3_1_ll", "S3_1_mean", "S3_1_ul",
                       "S3_2_ll", "S3_2_mean", "S3_2_ul",
                       "S3_12_ll",  "S3_12_mean", "S3_12_ul",
                       "S4_ll", "S4_mean", "S4_ul")

MALY_table=MALY_final
colnames(MALY_table)=c("year_ll","year_mean","year_ul")

path_table=cbind(path_1_final,
                 path_2_final,
                 path_3_final,
                 path_4_final,
                 path_5_final,
                 path_6_final,
                 path_7_final,
                 path_8_final,
                 path_9_final)

colnames(path_table)=c(
  "path1_ll", "path1_mean", "path1_ul", 
  "path2_ll", "path2_mean", "path2_ul",
  "path3_ll", "path3_mean", "path3_ul",
  "path4_ll", "path4_mean", "path4_ul",
  "path5_ll", "path5_mean", "path5_ul",
  "path6_ll", "path6_mean", "path6_ul",
  "path7_ll",  "path7_mean", "path7_ul",
  "path8_ll",  "path8_mean", "path8_ul",
  "path9_ll",  "path9_mean", "path9_ul")

##final tables in the total population

SSLY_table_pop=rbind(SSLY_1_final_pop,
                     SSLY_2_final_pop,
                     SSLY_3_final_pop,
                     SSLY_4_final_pop,
                     SSLY_5_final_pop,
                     SSLY_6_final_pop,
                     SSLY_7_final_pop,
                     SSLY_8_final_pop,
                     SSLY_9_final_pop)


MALY_table_pop=rbind(MALY_final_pop)


path_table_pop=rbind(path_1_final_pop,
                             path_2_final_pop,
                             path_3_final_pop,
                             path_4_final_pop,
                             path_5_final_pop,
                             path_6_final_pop,
                             path_7_final_pop,
                             path_8_final_pop,
                             path_9_final_pop)

rownames(path_table_pop)=c("S0→S4","S0→S1→S4","S0→S2→S4","S0→S1→S2→S4","S0→S3→S4",
                                   "S0→S1→S3→S4","S0→S2→S3→S4","S0→S1→S2→S3→S4", "S4")

##estimate the probability of each path among those who were dead (S4) by the end of follow up

path_table_pop_rescale=path_table_pop[,2]/sum(path_table_pop[,2])


##output

SSLY_table_pop  ##SSLY

MALY_table_pop   ##MALY

t(t(path_table_pop_rescale))   ##PATH

