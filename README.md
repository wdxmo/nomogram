Nomogram
========
This is about how to establish a nomogram model to predict patients' prognosis<br>
### Our outcome: overall survival (OS)<br>
### Methods: Cox proportional hazards (PH) regression model

#### Preparation before start  
locate your work file  
`setwd("C:\\Users\\29416\\OneDrive\\document\\Rstudio\\gc")`  
load corresponding packages needed  
`library(survival)`  
`library(rms)`  
`library(foreign)`  

read the data  
`seer<-read.csv("simple_training.csv")`  
change the categorical variables into factors<br> 
`seer$Age<-factor(seer$Age,labels=c("<70",">=70"))`<br>
`seer$Sex<-factor(seer$Sex,labels=c("Male","Female"))`<br>
`seer$Race<-factor(seer$Race,labels=c("White","Black","Other"))`<br>
`seer$Location<-factor(seer$Location,labels=c("Cardia","Non-cardia"))`<br>
`seer$tumor_size<-factor(seer$tumor_size,labels=c("<2cm","2-10cm",">10cm/diffuse"))`<br>
`seer$Grade<-factor(seer$Grade,labels=c("Well","Moderately","Poorly"))`<br>
`seer$stage_T<-factor(seer$stage_T,labels=c("T1","T2","T3","T4"))`<br>
`seer$stage_N<-factor(seer$stage_N,labels=c("N0","N1","N2","N3"))`<br>
`seer$stage_M<-factor(seer$stage_M,labels=c("M0","M1"))`<br>
`seer$Chemotherapy<-factor(seer$Chemotherapy,labels=c("No/Unknown","Yes"))`<br>
package the data<br>
`ddist <- datadist(seer)`<br>
`options(datadist='ddist')`<br>

### Let's start!
<br>

## To build Cox model
`fmla1 <- as.formula(Surv(Survival_month,Status) ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy)`<br>
`cox2 <- coxph(fmla1,data=seer)`<br>
`summary(cox2)`<br>
The C-index of our model are in the summary<br>

### AJCC 8th TNM C-index in training set
    fmla2 <- as.formula(Surv(Survival_month,Status) ~ stage_T + stage_N + stage_M)
    cox3 <- coxph(fmla2,data=seer)
    summary(cox3)

compare the two C-indexes<br>  

    coxpe <- predict(cox2)  
    coxpe1 <- predict(cox3)
    
    install.packages("compareC")
    library(compareC)
    compareC(seer$Survival_month,seer$Status,coxpe,coxpe1)

### C-index in validation set  
##### *change categorical variables as factors in the training set*
    validation <- read.csv("simple_validate.csv")
    validation$Age<-factor(validation$Age,labels=c("<70",">=70"))
    validation$Sex<-factor(validation$Sex,labels=c("Male","Female"))
    validation$Race<-factor(validation$Race,labels=c("White","Black","Other"))
    validation$Location<-factor(validation$Location,labels=c("Cardia","Non-cardia"))
    validation$tumor_size<-factor(validation$tumor_size,labels=c("<2cm","2-10cm",">10cm/diffuse"))
    validation$Grade<-factor(validation$Grade,labels=c("Well","Moderately","Poorly"))
    validation$stage_T<-factor(validation$stage_T,labels=c("T1","T2","T3","T4"))
    validation$stage_N<-factor(validation$stage_N,labels=c("N0","N1","N2","N3"))
    validation$stage_M<-factor(validation$stage_M,labels=c("M0","M1"))
    validation$Chemotherapy<-factor(validation$Chemotherapy,labels=c("No/Unknown","Yes"))

package the data<br>
`ddist <- datadist(validation)`<br>
`options(datadist='ddist')`<br>

#### Method 1: rcorr.cens
`library(Hmisc)`<br>

`fit <- cph(Surv(Survival_month,Status)~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data=seer)`<br>
`fp <- predict(fit,validation)`<br>
`fox <- rcorr.cens(fp,Surv(validation$Survival_month,validation$Status))`<br>
`cindex.orig=1-fox`<br>
`cindex.orig`<br>
The C-index in validatin set is presented.

### AJCC 8th TNM C-index in validation set
Just the same method to make it.

    fit1 <- cph(Surv(Survival_month,Status)~ stage_T + stage_N + stage_M, data=seer)
    fp1 <- predict(fit1,validation)
    dog <- rcorr.cens(fp1,Surv(validation$Survival_month,validation$Status))
    dog
    cindex.orig1=1-rcorr.cens(fp1,Surv(validation$Survival_month,validation$Status))
    cindex.orig1
compare the two C-indexes in validation set to acquire the statistical significance  
`compareC(validation$Survival_month,validation$Status,fp,fp1)`<br>


#### Method 2: survConcordance
`library(survival)`<br>
`fit <- cph(Surv(Survival_month,Status)~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data=seer)`<br>
`c_index <- survConcordance(Surv(validation$Survival_month,validation$Status)~predict(fit,validation))$concordance`<br>
`c_index`<br>

`fit1 <- cph(Surv(Survival_month,Status)~ stage_T + stage_N + stage_M, data=seer)`<br>
`c_index1 <- survConcordance(Surv(validation$Survival_month,validation$Status)~predict(fit1,validation))$concordance`<br>
`c_index1`<br>
<br>

## To make the Nomogram

`cox <- cph(Surv(Survival_month,Status) ~Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy,surv=T,x=T, y=T,data=seer) `<br>
`surv <- Survival(cox)`<br>
### 3-year survival
`sur_3_year<-function(x)surv(1*12*3,lp=x)`<br>
### 5-year survival
`sur_5_year<-function(x)surv(1*12*5,lp=x)`<br>
`nom_sur <- nomogram(cox,fun=list(sur_3_year,sur_5_year),lp= F,funlabel=c('3-year survival','5-year survival'),maxscale=100,fun.at=c('0.95','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1'))`
### raw nomogram
`pdf("nom.pdf")`<br>
`par(mar=c(2,4,1,1)+.1)`<br>
`plot(nom_sur,xfrac=0.4)`<br>
`dev.off()`<br>
<br>

## ROC curve (training set)
read the data<br>
`seer<-read.csv("simple_training.csv")`<br>
### add the risk score<br>
`cox_m <- coxph(Surv(Survival_month,Status) ~Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data = seer)`<br>
`cox_m1<-step(cox_m,direction = "both")`<br>
`risk_score<-predict(cox_m1,type="risk",newdata=seer)`<br>
`risk_level<-as.vector(ifelse(risk_score>median(risk_score),"High","Low"))`<br>
generate a new text file
`write.table(cbind(id=rownames(cbind(seer[,1:2],risk_score,risk_level)),cbind(seer[,1:2],risk_score,risk_level)),"risk_score.txt",sep="\t",quote=F,row.names=F)`<br>
<br>

### Plot the ROC curve<br>
`library(survivalROC)`<br>
read the text generated in last step<br>
`seer<-read.table("risk_score.txt",header=T,sep="\t")`<br>

#### *3 year ROC*
    predict_time<-12*3
    myroc<-survivalROC(Stime=seer$Survival_month, status=seer$Status, marker=seer$risk_score, predict.time=predict_time,method="KM")
    a<-round(myroc$AUC,3)
    a`
Plot the ROC<br>

    pdf("ROC_3year.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(myroc$FP,myroc$TP,type="l",xlim=c(0,1),ylim=c(0,1),xaxs = "i", yaxs ="i",
         cex.axis=1.5,cex.lab=1.5, col="blue",
         xlab="1-Specificity",ylab="Sensitivity",main=paste("3-year Survival","AUC=",round(myroc$AUC,3)))
    abline(0,1)
    dev.off()

#### *5 year ROC*
    predict_time<-12*5 
    myroc<-survivalROC(Stime=seer$Survival_month, status=seer$Status, marker=seer$risk_score, predict.time=predict_time,method="KM")
    a<-round(myroc$AUC,3)
    a
Plot the ROC<br>

    pdf("ROC_5year.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(myroc$FP,myroc$TP,type="l",xlim=c(0,1),ylim=c(0,1),xaxs = "i", yaxs ="i",
         col="blue",cex.axis=1.5,cex.lab=1.5,
         xlab="1-Specificity",ylab="Sensitivity",main=paste("5-year Survival","AUC=",round(myroc$AUC,3)))
    abline(0,1)
    dev.off()
<br>


## ROC curve (validation set)
read the two sets of data<br>
`seer<-read.csv("simple_training.csv")`<br>
`seer1<-read.csv("simple_validate.csv")`<br>
### add the risk score<br>
`cox_m <- coxph(Surv(Survival_month,Status) ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data = seer)`<br>
`cox_m1<-step(cox_m,direction = "both")`<br>
`risk_score<-predict(cox_m1,type="risk",newdata=seer1)`<br>
`risk_level<-as.vector(ifelse(risk_score>median(risk_score),"High","Low"))`<br>
generate a new text file<br>
`write.table(cbind(id=rownames(cbind(seer1[,1:2],risk_score,risk_level)),cbind(seer1[,1:2],risk_score,risk_level)),"risk_score_val.txt",sep="\t",quote=F,row.names=F)`<br>


### Plot the ROC curve
read the text generated in last step<br>
`seer<-read.table("risk_score_val.txt",header=T,sep="\t")`<br>

#### *3 year ROC*
    predict_time<-12*3
    myroc<-survivalROC(Stime=seer$Survival_month, status=seer$Status, marker=seer$risk_score, predict.time=predict_time,method="KM")
    a<-round(myroc$AUC,3)
    a
Plot the ROC<br>

    pdf("ROC_val_3year.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(myroc$FP,myroc$TP,type="l",xlim=c(0,1),ylim=c(0,1),xaxs = "i", yaxs ="i",
         col="blue",cex.axis=1.5,cex.lab=1.5,
         xlab="1-Specificity",ylab="Sensitivity",main=paste("3-year Survival","AUC=",round(myroc$AUC,3)))
    abline(0,1)
    dev.off()

#### *5 year ROC*
    predict_time<-12*5
    myroc<-survivalROC(Stime=seer$Survival_month, status=seer$Status, marker=seer$risk_score, predict.time=predict_time,method="KM")
    a<-round(myroc$AUC,3)
    a
Plot the ROC<br>

    pdf("ROC_val_5year.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(myroc$FP,myroc$TP,type="l",xlim=c(0,1),ylim=c(0,1),xaxs = "i", yaxs ="i",
         col="blue",cex.axis=1.5,cex.lab=1.5,
         xlab="1-Specificity",ylab="Sensitivity",main=paste("5-year Survival","AUC=",round(myroc$AUC,3)))
    abline(0,1)
    dev.off()
<br>

## Decision curve analysis (DCA)

`attach(seer)`<br>
`str(seer)`<br>

First we need to download a "stdca.R" file<br>

download from [HERE](https://www.mskcc.org/sites/default/files/node/4509/documents/downloadrcode.zip)<br>

`source("stdca.R")`<br>
`library(survival)`<br>

`Srv = Surv(seer$Survival_month, seer$Status)`<br>

### *3 year in training set*
TNM model & Our nomogram

    coxmod1 <- coxph(Srv ~ stage_T + stage_N + stage_M, data=seer)
    coxmod2 <- coxph(Srv ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data=seer)

    seer$TNM_stage <- c(1 - (summary(survfit(coxmod1,newdata=seer), times=3)$surv))
    seer$Nomogram <- c(1 - (summary(survfit(coxmod2,newdata=seer), times=3)$surv))
Plot the DCA curve<br>

    pdf("dca3_tra.pdf")
    stdca(data=seer, outcome="Status", ttoutcome="Survival_month", timepoint=3,predictors=c("TNM_stage","Nomogram"), xstop=0.6, smooth=TRUE)
    dev.off()

### *5 year in training set*
    coxmod1 <- coxph(Srv ~ stage_T + stage_N + stage_M, data=seer)
    coxmod2 <- coxph(Srv ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data=seer)

    seer$TNM_stage <- c(1 - (summary(survfit(coxmod1,newdata=seer), times=5)$surv))
    seer$Nomogram <- c(1 - (summary(survfit(coxmod2,newdata=seer), times=5)$surv))
Plot the DCA curve<br>

    pdf("dca5_tra2.pdf")
    stdca(data=seer, outcome="Status", ttoutcome="Survival_month", timepoint=5,predictors=c("TNM_stage","Nomogram"), xstop=0.6, smooth=TRUE)
    dev.off()


`str(validation)`<br>

### *3 year in validation set*
TNM model & Our nomogram

    coxmod1 <- coxph(Srv ~ stage_T + stage_N + stage_M, data=seer)
    coxmod2 <- coxph(Srv ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data=seer)

    validation$TNM_stage <- c(1 - (summary(survfit(coxmod1,newdata=validation), times=3)$surv))
    validation$Nomogram <- c(1 - (summary(survfit(coxmod2,newdata=validation), times=3)$surv))
Plot the DCA curve

    pdf("dca3_val.pdf")
    stdca(data=validation, outcome="Status", ttoutcome="Survival_month", timepoint=3,predictors=c("TNM_stage","Nomogram"), xstop=0.6, smooth=TRUE)
    dev.off()

### *5 year in validation set*
    coxmod1 <- coxph(Srv ~ stage_T + stage_N + stage_M, data=seer)
    coxmod2 <- coxph(Srv ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy, data=seer)

    validation$TNM_stage <- c(1 - (summary(survfit(coxmod1,newdata=validation), times=5)$surv))
    validation$Nomogram <- c(1 - (summary(survfit(coxmod2,newdata=validation), times=5)$surv))
Plot the DCA curve

    pdf("dca5_val.pdf")
    stdca(data=validation, outcome="Status", ttoutcome="Survival_month", timepoint=5,predictors=c("TNM_stage","Nomogram"), xstop=0.6, smooth=TRUE)
    dev.off()
<br>

## Calibration curve

`library(rms)`<br>
`library(foreign)`<br>
`library(survival)`<br>

### *3-year in training set*
    cox_3 <- cph(Surv(Survival_month,Status) ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy,surv=T,x=T, y=T,time.inc = 1*12*3,data=seer) 
    cal <- calibrate(cox_3, cmethod="KM", method="boot", u=1*12*3, m=1500, B=1000)
Plot the curve<br>

    pdf("calibrate3.pdf",10,8)
    par(mar = c(8,6,3,3),lwd=3)
    plot(cal,lwd=2,lty=1,xlim = c(0,1),ylim = c(0,1),xaxs = "i", yaxs ="i",
     xlab ="Nomogram-Predicted 3-Year Survival",ylab="Actual 3-Year Survival (proportion)",
     cex.lab=1.5,col="blue",axes=F)
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.5)
    axis(1,at=c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.5)
    dev.off()

### *5-year in training set*
    cox_5 <- cph(Surv(Survival_month,Status) ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + Chemotherapy,surv=T,x=T, y=T,time.inc = 1*12*5,data=seer) 
    cal <- calibrate(cox_5, cmethod="KM", method="boot", u=1*12*5, m=1500, B=1000)
Plot the curve<br>

    pdf("calibrate5.pdf",10,8)
    par(mar = c(8,6,3,3),lwd=3)
    plot(cal,lwd=2,lty=1,xlim = c(0,1),ylim = c(0,1),xaxs = "i", yaxs ="i",
     xlab ="Nomogram-Predicted 5-Year Survival",ylab="Actual 5-Year Survival (proportion)",
     cex.lab=1.5,col="blue",axes=F)
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.5)
    axis(1,at=c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.5)
    dev.off()


### *3-year in external validation set*

    fev3 <- cph(Surv(Survival_month,Status) ~ fp, x=T, y=T, surv=T, data=validation, time.inc=36)
    calev3 <- calibrate(fev3, cmethod="KM", method="boot", u=36, m=600, B=1000)
Plot the curve<br>

    pdf("calibrate3_vali.pdf",10,8)
    par(mar = c(8,6,3,3),lwd=3)
    plot(calev3,lwd=2,lty=1,xlim = c(0,1),ylim = c(0,1),xaxs = "i", yaxs ="i",
     xlab ="Nomogram-Predicted 3-Year Survival",ylab="Actual 3-Year Survival (proportion)",
     cex.lab=1.5,col="blue",axes=F)
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.5)
    axis(1,at=c(0.0,0.2,0.4,0.6,0.8,1.0),cex.axis=1.5)
    dev.off()
<br>


## KM survival curve

`library(survival)`<br>
read the data<br>
`inputdata1<- read.csv("simple_training.csv")`<br>

### by Age
    kms<-survfit(Surv(Survival_month,Status)~Age,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~Age,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-age.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","goldenrod1"),lwd = 2,mark.time = T,
     xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
     main=paste("Age (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("<70",">=70"),lty=2,col=c("dodgerblue3","goldenrod1"))
    dev.off()

### by Sex
    kms<-survfit(Surv(Survival_month,Status)~Sex,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~Sex,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-sex.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","goldenrod1"),lwd = 2,mark.time = T,
     xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
     main=paste("Sex (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("Male","Female"),lty=2,col=c("dodgerblue3","goldenrod1"))
    dev.off()

### by Race
    kms<-survfit(Surv(Survival_month,Status)~Race,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~Race,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-race.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("goldenrod1","firebrick","dodgerblue3"),lwd = 2,mark.time = T,
     xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
     main=paste("Race (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("Other","White","Black"),lty=2,col=c("dodgerblue3","goldenrod1","firebrick"))
    dev.off()

### by Location
    kms<-survfit(Surv(Survival_month,Status)~Location,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~Location,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-location.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("goldenrod1","dodgerblue3"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("Location (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("Non-cardia","Cardia"),lty=2,col=c("dodgerblue3","goldenrod1"))
    dev.off()

### by tumor size
    kms<-survfit(Surv(Survival_month,Status)~tumor_size,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~tumor_size,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-tumor_size.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","gold","firebrick"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("Tumor size (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("<2cm","2-10cm",">10cm/diffuse"),lty=2,col=c("dodgerblue3","gold","firebrick"))
    dev.off()

### by Grade
    kms<-survfit(Surv(Survival_month,Status)~Grade,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~Grade,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-grade.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","gold","firebrick"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("Grade (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("Well","Moderately","Poorly"),lty=2,col=c("dodgerblue3","gold","firebrick"))
    dev.off()

### by AJCC 8th Stage
    kms<-survfit(Surv(Survival_month,Status)~AJCC_8th,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~AJCC_8th,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-AJCC_8th_Stage.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","gold","firebrick","black"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("AJCC 8th (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("Stage I","Stage II","Stage III","Stage IV"),lty=2,col=c("dodgerblue3","gold","firebrick","black"))
    dev.off()

### by Stage T
    kms<-survfit(Surv(Survival_month,Status)~stage_T,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~stage_T,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-stage_T.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","gold","firebrick","black"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("Stage T (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("T1","T2","T3","T4"),lty=2,col=c("dodgerblue3","gold","firebrick","black"))
    dev.off()

### by Stage N
    kms<-survfit(Surv(Survival_month,Status)~stage_N,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~stage_N,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-stage_N.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","gold","firebrick","black"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("Stage N (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("N0","N1","N2","N3"),lty=2,col=c("dodgerblue3","gold","firebrick","black"))
    dev.off()

### by Stage M
    kms<-survfit(Surv(Survival_month,Status)~stage_M,data=inputdata1)
    kmdffexp=survdiff(Surv(Survival_month,Status)~stage_M,data=inputdata1)
    pValue=round(1-pchisq(kmdffexp$chisq,df=1),3)
Plot the curve<br>

    pdf("survival-stage_M.pdf")
    par(mar=c(5,5,4,2)+.1)
    plot(kms,lty=2,col=c("dodgerblue3","goldenrod1"),lwd = 2,mark.time = T,
         xlab="Survival time (years)",ylab="Overall Survival",cex.lab=1.5,cex.main=1.5,axes=F,
         main=paste("Stage M (P=", pValue ,")",sep=""))
    axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),las=1,pos=0,cex.axis=1.5)
    axis(1,at=c(0,24,48,72,96,120,144),labels=c(0,2,4,6,8,10,12),pos=0,cex.axis=1.5)
    abline(h = 0.5,lty = 3)
    abline(v = 36,lty = 3)
    abline(v = 60,lty = 3)
    legend("topright",c("M0","M1"),lty=2,col=c("dodgerblue3","goldenrod1"))
    dev.off()

### by chemotherapy (adjusted before and after)

`cox2<-coxph(Surv(Survival_month,Status) ~ Age + Race + Location + tumor_size + Grade + stage_T + stage_N + stage_M + strata(Chemotherapy), data=seer)`<br>
`cox.zph(cox2,transform=rank)`<br>

`library("ggplot2")`<br>
`library("survival")`<br>
`library("survminer")`<br>

#### plot curve after adjusted

    pdf("survival-chemotherapy_corr.pdf")
    ggadjustedcurves(cox2,pval=TRUE,data=seer,method = "marginal", variable = "Chemotherapy",legend.title="",
                     xlab="Survival time (months)",ylab="Overall Survival",
                     main=paste("Chemotherapy Adjusted"))
    dev.off()
#### plot curve before adjusted

    pdf("survival-chemotherapy.pdf")
    ggadjustedcurves(cox2,pval=TRUE,data=seer,method = "average", variable = "Chemotherapy",legend.title="",
                     xlab="Survival time (months)",ylab="Overall Survival",
                     main=paste("Chemotherapy"))
    dev.off()


## make table1

`install.packages("table1")`<br>
`library(table1)`<br>
`library(lubridate)`<br>

`seer<-read.csv("table1_simple.csv")`<br>
`names(seer)`<br>
Change categorical variables into factors<br>

    seer$Age_group<-factor(seer$Age_group,labels=c("<70",">=70"))
    seer$Sex<-factor(seer$Sex,labels=c("Male","Female"))
    seer$Race<-factor(seer$Race,labels=c("White","Black","Other"))
    seer$Set<-factor(seer$Set,labels=c("training set","validation set"))
    seer$Location<-factor(seer$Location,labels=c("Cardia","Non-cardia"))
    seer$tumor_size<-factor(seer$tumor_size,labels=c("<2cm","2-10cm",">10cm/diffuse"))
    seer$Grade<-factor(seer$Grade,labels=c("Well","Moderately","Poorly"))
    seer$AJCC_8th<-factor(seer$AJCC_8th,labels=c("I","II","III","IV"))
    seer$stage_T<-factor(seer$stage_T,labels=c("T1","T2","T3","T4"))
    seer$stage_N<-factor(seer$stage_N,labels=c("N0","N1","N2","N3"))
    seer$stage_M<-factor(seer$stage_M,labels=c("M0","M1"))
    seer$Radiation<-factor(seer$Radiation,labels=c("No radiation or surgery","Yes"))
    seer$Chemotherapy<-factor(seer$Chemotherapy,labels=c("No/Unknown","Yes"))

### make table 1 by group (Set)
`units(seer$Age) <- "years"`<br>
`table1(~ Age + Age_group + Sex + Race + Location + Grade + AJCC_8th + stage_T + stage_N + stage_M + tumor_size + nodes_examined + nodes_positive + Radiation + Chemotherapy | Set, data=seer, overall = "Total")`<br>

### add more statistics to table 1
    units(seer$nodes_examined) <- "No."
    units(seer$nodes_positive) <- "No."
    table1(~ Age + Age_group + Sex + Race + Location + Grade + AJCC_8th + stage_T + stage_N + stage_M + tumor_size + nodes_examined + nodes_positive + Radiation + Chemotherapy | Set, data=seer,
           render.continuous=c(.="Mean (SD)",.="Median [Q1, Q3]"))

### add P value
`library(MatchIt)`<br>
`seer$Set <- factor(seer$Set,levels=c(0, 1, 2),labels=c("traning set", "validation set","P-value"))`<br>

### creat P values
    data <- seer
    outcome <- seer$Set
    
    rndr <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- data[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
          p <- t.test(y ~ outcome)$p.value
        } else {
          p <- chisq.test(table(y, droplevels(outcome)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
      } else {
        render.default(x=x, name=name, ...)
      }
    }
    
    rndr.strat <- function(label, n, ...) {
      ifelse(n==0, label, render.strat.default(label, n, ...))
    }

### make table 1 new

    table1(~ Age + Age_group + Sex + Race + Location + Grade + AJCC_8th + stage_T + stage_N + stage_M + tumor_size + nodes_examined + nodes_positive + Radiation + Chemotherapy | Set, data=seer, 
       render.continuous=c(.="Mean (SD)",.="Median [Q1, Q3]"),droplevels=F, render=rndr, render.strat=rndr.strat, overall="Total")
       
