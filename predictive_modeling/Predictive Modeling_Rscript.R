
### set working directory manually
# Session > Set Working Directory > To Source File Location


### loading R libraries
library(MASS)
library(survival)
library(penalized)
library(glmnet)
library(survAUC)
library(timeROC)
library(Hmisc)
library(ClassComparison)


### load TCGA BRCA R dataset
load("brca.Rda")
### load clustering group information
load("brca.clustering.group.Rda")




######################################################################
### 
### build a new dataset
###
######################################################################

### filter out genes with 0 count for all samples
zerorow <- apply(logBrcaData,1,function(x) sum(x==log2(3)))
logBrcaData1 <- logBrcaData[zerorow<dim(logBrcaData)[2],]

### process data
logBrcaData1T <- t(logBrcaData1)
temp <- sweep(logBrcaData1T, 2, apply(logBrcaData1T, 2, mean), "-")
colsd <- apply(logBrcaData1T,2,sd)
geneExpr <- sweep(temp, 2, colsd, "/")
Idmatching <- match(rownames(logBrcaData1T),rownames(brcaClin))

### build a new dataset
data <- list()
data$name  <- "TCGA breast cancer dataset"
data$geneExpr <- geneExpr
data$geneId <- colnames(logBrcaData1T)
data$sampleId <- rownames(logBrcaData1T)
data$age <- brcaClin$Diagnosis.Age[Idmatching]
data$race <- brcaClin$Race.Category[Idmatching]
data$stage <- brcaClin$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code[Idmatching]
data$DFS.time <- brcaClin$Disease.Free..Months.[Idmatching]
data$DFS.status <- brcaClin$Disease.Free.Status[Idmatching]
data$OS.time <- brcaClin$Overall.Survival..Months.[Idmatching]
data$OS.status <- brcaClin$Overall.Survival.Status[Idmatching]
data$OS.status <- 1*(data$OS.status=="DECEASED")
data$DFS.status <- 1*(data$DFS.status=="Recurred/Progressed")





######################################################################
### 
### Survival curves by clustering groups
###
######################################################################

### remove NA and 0 in overall survival time
patient.sel <- which(data$OS.time>0) 
time <- data$OS.time[patient.sel]
status <- data$OS.status[patient.sel]
S <- Surv(time,status)
group <- group[patient.sel]
logBrcaData1 <- logBrcaData[,patient.sel]

### Kaplan-Meier curve
fit.surv <- survfit(S ~ group)
plot(fit.surv, lwd = 2, lty = 1, col = c("red","blue","green","purple"), xlab = 'Time', ylab = 'Estimated Survival Function')
legend("bottomleft", legend=c("G1","G2","G3","G4"), lty = c(1,1), col = c("red","blue","green","purple"), lwd = 2)
title("Kaplan-Meier Curve")


#### gene-by-gene t-test
group.24 <- group %in% c("G2","G4")
group.test <- factor(group[group.24])
mtt <- MultiTtest(logBrcaData1[,group.24], group.test)
# model to estimate FDR
bum <- Bum(mtt@p.values)
hist(bum)
countSignificant(bum, alpha=1e-6, by="FDR")
sum(filt & selectSignificant(bum, alpha=1e-6, by="FDR"))

# note: filt was defined in the 07-clustering.R script; it filters based
# on a combination of mean expression and variance
# for later use in gene set enrichment analysis, we save the filtered
# gene that are significantly differentially expressed
temp <- as.data.frame(mtt)
temp <- temp[,-3]
sigGenes <- temp[filt & selectSignificant(bum, alpha=1e-6, by="FDR"),]
genename <- rownames(sigGenes)
genelist <- as.character(genename[!is.na(match(genename,data$geneId))])





######################################################################
### 
### predictive modeling with the obtained genelist
###
######################################################################

set.seed(8)
AUC.Uno <- matrix(0,nrow=1,ncol=11)
AUC.Blanche <- matrix(0,nrow=1,ncol=11)
C.Harrell <- matrix(0,nrow=1,ncol=11)
C.Uno <- matrix(0,nrow=1,ncol=11)



### split data into 2/3 training set and 1/3 test set
# subset expression matrix with the gene list
X <- data$geneExpr[patient.sel,genelist]
dim(X) #884 patients and 339 genes
K <- 3
fold <- split(sample(1:nrow(X)), rep(1:K, length=nrow(X)))
id.train <- c(fold[[1]],fold[[2]])
id.test <- fold[[K]]
X.train <- X[id.train,]
S.train <- S[id.train,]
time.train <- time[id.train]
status.train <- status[id.train]
X.test <- X[id.test,]
S.test <- S[id.test,]
time.test <- time[id.test]
status.test <- status[id.test]



### fitting univariate Cox model 
cox.uni.p <- NULL
m <- dim(X.train)[2]
cox.uni.p <- sapply(1:m, function(i) summary(coxph(S.train ~ X.train[,i], method="efron"))$logtest['pvalue'])

### histogram of univariate p-values
cox.uni.p.bum <- Bum(cox.uni.p)
hist(cox.uni.p.bum)

### select top genes for building predictive modeling
G <- 50
sel.order <- order(cox.uni.p)
top.id <- genelist[sel.order[1:G]]
filter1.id <- genelist[cox.uni.p<0.2]
filter2.id <- genelist


### stepwise (backward) selection with AIC
cox1 <- coxph(S.train ~ ., data=as.data.frame(X.train[,top.id]))
s1 <- step(cox1,direction="backward")

#step AIC plot
steps <- 1:length(s1$anova$Step)
plot(steps,s1$anova$AIC,xaxt="n",ylab="AIC",xlab="step",type='o')
axis(1,at=steps,labels=FALSE)
text(x=steps,par("usr")[3]-5, labels=s1$anova$Step, srt=60,pos=1, xpd=TRUE)

#genes selected with optimal model fitting
sel.step <- names(coef(s1)) 
#fit model with optimal genes
cox.step <- coxph(S.train ~ ., data=as.data.frame(X.train[,sel.step]),method="breslow")
#risk score
rs.step.train <- predict(cox.step,type="lp")
rs.step.test <- predict(cox.step,newdata=as.data.frame(X.test[,sel.step]),type="lp")



### stepwise selection with cross-validation
cv.l <- NULL
for (i in 1:length(top.id))
   {
	cv.l[i] <- try(cvl(S.train, X.train[, sel.order[1:i]], fold = 5)$cvl)
   }

#optimal tuning parameter
lambda.cv <- which.max(cv.l)  #number of genes selected to reach the maximum loglikelihood; 
#cross-validation likelihood plot
plot(1:length(top.id),cv.l,xlab="number of genes in order",ylab="CV likelihood",type='o')
abline(v=c(lambda.cv),col="red",lwd=2)

#genes selected with optimal model fitting
sel.cv <- genelist[sel.order[1:lambda.cv]]
#fit model with optimal lambda.uni
cox.cv <- coxph(S.train ~ ., data=as.data.frame(X.train[,sel.cv]),method="breslow")
#risk score
rs.cv.train <- predict(cox.cv,type="lp")
rs.cv.test <- predict(cox.cv,newdata=as.data.frame(X.test[,sel.cv]),type="lp")



### penalized likelihood method using glmnet
pen.lik <- function(geneid,alpha)
{
  cv.pen <- try(cv.glmnet(X.train[,geneid], S.train, family = "cox", alpha=alpha))
  #genes selected with optimal model fitting
  coef.min <-  coef(cv.pen, s = "lambda.min")
  active.min <- which(coef.min != 0)
  index.min <- coef.min[active.min]
  sel.pen <- rownames(coef.min)[active.min]
  #risk score
  rs.train <- predict(cv.pen,newx=X.train[,geneid],s="lambda.min")
  rs.test <- predict(cv.pen,newx=X.test[,geneid],s="lambda.min")
  #lasso fit
  coxfit <- glmnet(X.train[,geneid], S.train, family = "cox", alpha=alpha)
  return(list(cv.pen=cv.pen,coxfit=coxfit,rs.train=rs.train,rs.test=rs.test,sel.pen=sel.pen))
}

plot.pen <- function(pen.method)
{
  par(mfrow=c(2,1))
  plot(pen.method$cv.pen)
  plot(pen.method$coxfit,xvar="lambda",label=TRUE)
  abline(v=log(pen.method$cv.pen$lambda.min),col="black",lty=2)
  par(mfrow=c(1,1))
}
                     
f <- "pen.Rda"       
if (file.exists(f)) {
  load(f)
} else {
  #lasso method with top 50 genes
  pen.lasso <- pen.lik(top.id,1)
  #elastic net method with top 50 genes
  pen.elanet <- pen.lik(top.id,0.5)
  #ridge method with top 50 genes
  pen.ridge <- pen.lik(top.id,0)
  
  #lasso method with genes p<0.2
  pen.lasso1 <- pen.lik(filter1.id,1)
  #elastic net method with genes p<0.2
  pen.elanet1 <- pen.lik(filter1.id,0.5)
  #ridge method with genes p<0.2
  pen.ridge1 <- pen.lik(filter1.id,0)
  
  #lasso method with all genes
  pen.lasso2 <- pen.lik(filter2.id,1)
  #elastic net method with all genes
  pen.elanet2 <- pen.lik(filter2.id,0.5)
  #ridge method with all genes
  pen.ridge2 <- pen.lik(filter2.id,0)
  save(pen.lasso,pen.elanet,pen.ridge,pen.lasso1,pen.elanet1,pen.ridge1,pen.lasso2,pen.elanet2,pen.ridge2,file=f)
}
rm(f)

#lasso method with top 50 genes
plot.pen(pen.lasso)
#elastic net method with top 50 genes
plot.pen(pen.elanet)
#ridge method with top 50 genes
plot.pen(pen.ridge)

#lasso method with genes p<0.2
plot.pen(pen.lasso1)
#elastic net method with genes p<0.2
plot.pen(pen.elanet1)
#ridge method with genes p<0.2
plot.pen(pen.ridge1)

#lasso method with all genes
plot.pen(pen.lasso2)
#elastic net method with all genes
plot.pen(pen.elanet2)
#ridge method with all genes
plot.pen(pen.ridge2)






######################################################################
### 
### evaluate survival prediction
###
######################################################################

### plot risk prediction scores by 4 groups
plot(group[id.test],rs.step.test,main="Step-AIC",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],rs.cv.test,main="Step-CV",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.lasso$rs.test,main="Lasso",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.elanet$rs.test,main="Elastic Net",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.ridge$rs.test,main="Ridge",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.lasso1$rs.test,main="Lasso p<0.2",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.elanet1$rs.test,main="Elastic Net p<0.2",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.ridge1$rs.test,main="Ridge p<0.2",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.lasso2$rs.test,main="Lasso All",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.elanet2$rs.test,main="Elastic Net All",ylab="risk score",col=c("red","blue","green","purple"))
plot(group[id.test],pen.ridge2$rs.test,main="Ridge All",ylab="risk score",col=c("red","blue","green","purple"))
                                                  

pred.measure <- function(type,rs.train,rs.test,sel,times=times)
{
  if (type=="AUC.Blanche") 
  {
    ROC <- timeROC(T=time.test, delta=status.test, marker=rs.test, cause=1, 
                   times=times,iid=TRUE)
    w <- ROC$survProb/sum(ROC$survProb)
    value <- w[!is.na(ROC$AUC)]%*%ROC$AUC[!is.na(ROC$AUC)]
  }
  if (type=="AUC.Uno") value <- AUC.uno(S.train,S.test,rs.test,times)$iauc
  if (type=='C.Harrell') value <- rcorr.cens(-rs.test,S.test)["C Index"]
  if (type=='C.Uno') value <- UnoC(S.train,S.test,rs.test)
  return(value)
}
times <- quantile(time,probs=seq(0.1,0.9,0.05))
par(cex.lab=1,cex.axis=1)
methods <- c("Step-AIC top50","Step-CV top50","Lasso top50","ElasticNet top50","Ridge top50","Lasso p<0.2","ElasticNet p<0.2","Ridge p<0.2","Lasso all","ElasticNet all","Ridge all")

### Uno et al.'s AUC
AUC.Uno[1] <- pred.measure('AUC.Uno',rs.step.train,rs.step.test,sel.step,times)
AUC.Uno[2] <- pred.measure('AUC.Uno',rs.cv.train,rs.cv.test,sel.cv,times)
AUC.Uno[3] <- pred.measure('AUC.Uno',pen.lasso$rs.train,pen.lasso$rs.test,pen.lasso$sel.pen,times)
AUC.Uno[4] <- pred.measure('AUC.Uno',pen.elanet$rs.train,pen.elanet$rs.test,pen.elanet$sel.pen,times)
AUC.Uno[5] <- pred.measure('AUC.Uno',pen.ridge$rs.train,pen.ridge$rs.test,pen.ridge$sel.pen,times)
AUC.Uno[6] <- pred.measure('AUC.Uno',pen.lasso1$rs.train,pen.lasso1$rs.test,pen.lasso1$sel.pen,times)
AUC.Uno[7] <- pred.measure('AUC.Uno',pen.elanet1$rs.train,pen.elanet1$rs.test,pen.elanet1$sel.pen,times)
AUC.Uno[8] <- pred.measure('AUC.Uno',pen.ridge1$rs.train,pen.ridge1$rs.test,pen.ridge1$sel.pen,times)
AUC.Uno[9] <- pred.measure('AUC.Uno',pen.lasso2$rs.train,pen.lasso2$rs.test,pen.lasso2$sel.pen,times)
AUC.Uno[10] <- pred.measure('AUC.Uno',pen.elanet2$rs.train,pen.elanet2$rs.test,pen.elanet2$sel.pen,times)
AUC.Uno[11] <- pred.measure('AUC.Uno',pen.ridge2$rs.train,pen.ridge2$rs.test,pen.ridge2$sel.pen,times)
boxplot(AUC.Uno,col=rep(c("red","green","blue"),c(5,3,3)),xaxt="n",at=c(1:5,7:9,11:13),ylim=c(0.5,0.85),ylab="AUC.Uno")
axis(1,at=c(1:5,7:9,11:13),labels=FALSE)
text(x=c(1:5,7:9,11:13),par("usr")[3]-0.03,labels=methods, srt=20,pos=1, xpd=TRUE,cex=0.8)


### Blanche et al.'s AUC
AUC.Blanche[1] <- pred.measure('AUC.Blanche',rs.step.train,rs.step.test,sel.step,times)
AUC.Blanche[2] <- pred.measure('AUC.Blanche',rs.cv.train,rs.cv.test,sel.cv,times)
AUC.Blanche[3] <- pred.measure('AUC.Blanche',pen.lasso$rs.train,pen.lasso$rs.test,pen.lasso$sel.pen,times)
AUC.Blanche[4] <- pred.measure('AUC.Blanche',pen.elanet$rs.train,pen.elanet$rs.test,pen.elanet$sel.pen,times)
AUC.Blanche[5] <- pred.measure('AUC.Blanche',pen.ridge$rs.train,pen.ridge$rs.test,pen.ridge$sel.pen,times)
AUC.Blanche[6] <- pred.measure('AUC.Blanche',pen.lasso1$rs.train,pen.lasso1$rs.test,pen.lasso1$sel.pen,times)
AUC.Blanche[7] <- pred.measure('AUC.Blanche',pen.elanet1$rs.train,pen.elanet1$rs.test,pen.elanet1$sel.pen,times)
AUC.Blanche[8] <- pred.measure('AUC.Blanche',pen.ridge1$rs.train,pen.ridge1$rs.test,pen.ridge1$sel.pen,times)
AUC.Blanche[9] <- pred.measure('AUC.Blanche',pen.lasso2$rs.train,pen.lasso2$rs.test,pen.lasso2$sel.pen,times)
AUC.Blanche[10] <- pred.measure('AUC.Blanche',pen.elanet2$rs.train,pen.elanet2$rs.test,pen.elanet2$sel.pen,times)
AUC.Blanche[11] <- pred.measure('AUC.Blanche',pen.ridge2$rs.train,pen.ridge2$rs.test,pen.ridge2$sel.pen,times)
boxplot(AUC.Blanche,col=rep(c("red","green","blue"),c(5,3,3)),xaxt="n",at=c(1:5,7:9,11:13),ylim=c(0.5,0.85),ylab="AUC.Blanche")
axis(1,at=c(1:5,7:9,11:13),labels=FALSE)
text(x=c(1:5,7:9,11:13),par("usr")[3]-0.03,labels=methods, srt=20,pos=1, xpd=TRUE,cex=0.8)


### Harrell et al.'s C-Index
C.Harrell[1] <- pred.measure('C.Harrell',rs.step.train,rs.step.test,sel.step,times)
C.Harrell[2] <- pred.measure('C.Harrell',rs.cv.train,rs.cv.test,sel.cv,times)
C.Harrell[3] <- pred.measure('C.Harrell',pen.lasso$rs.train,pen.lasso$rs.test,pen.lasso$sel.pen,times)
C.Harrell[4] <- pred.measure('C.Harrell',pen.elanet$rs.train,pen.elanet$rs.test,pen.elanet$sel.pen,times)
C.Harrell[5] <- pred.measure('C.Harrell',pen.ridge$rs.train,pen.ridge$rs.test,pen.ridge$sel.pen,times)
C.Harrell[6] <- pred.measure('C.Harrell',pen.lasso1$rs.train,pen.lasso1$rs.test,pen.lasso1$sel.pen,times)
C.Harrell[7] <- pred.measure('C.Harrell',pen.elanet1$rs.train,pen.elanet1$rs.test,pen.elanet1$sel.pen,times)
C.Harrell[8] <- pred.measure('C.Harrell',pen.ridge1$rs.train,pen.ridge1$rs.test,pen.ridge1$sel.pen,times)
C.Harrell[9] <- pred.measure('C.Harrell',pen.lasso2$rs.train,pen.lasso2$rs.test,pen.lasso2$sel.pen,times)
C.Harrell[10] <- pred.measure('C.Harrell',pen.elanet2$rs.train,pen.elanet2$rs.test,pen.elanet2$sel.pen,times)
C.Harrell[11] <- pred.measure('C.Harrell',pen.ridge2$rs.train,pen.ridge2$rs.test,pen.ridge2$sel.pen,times)
boxplot(C.Harrell,col=rep(c("red","green","blue"),c(5,3,3)),xaxt="n",at=c(1:5,7:9,11:13),ylim=c(0.5,0.85),ylab="C.Harrell")
axis(1,at=c(1:5,7:9,11:13),labels=FALSE)
text(x=c(1:5,7:9,11:13),par("usr")[3]-0.03,labels=methods, srt=20,pos=1, xpd=TRUE,cex=0.8)


### Uno et al.'s C-statistics
C.Uno[1] <- pred.measure('C.Uno',rs.step.train,rs.step.test,sel.step,times)
C.Uno[2] <- pred.measure('C.Uno',rs.cv.train,rs.cv.test,sel.cv,times)
C.Uno[3] <- pred.measure('C.Uno',pen.lasso$rs.train,pen.lasso$rs.test,pen.lasso$sel.pen,times)
C.Uno[4] <- pred.measure('C.Uno',pen.elanet$rs.train,pen.elanet$rs.test,pen.elanet$sel.pen,times)
C.Uno[5] <- pred.measure('C.Uno',pen.ridge$rs.train,pen.ridge$rs.test,pen.ridge$sel.pen,times)
C.Uno[6] <- pred.measure('C.Uno',pen.lasso1$rs.train,pen.lasso1$rs.test,pen.lasso1$sel.pen,times)
C.Uno[7] <- pred.measure('C.Uno',pen.elanet1$rs.train,pen.elanet1$rs.test,pen.elanet1$sel.pen,times)
C.Uno[8] <- pred.measure('C.Uno',pen.ridge1$rs.train,pen.ridge1$rs.test,pen.ridge1$sel.pen,times)
C.Uno[9] <- pred.measure('C.Uno',pen.lasso2$rs.train,pen.lasso2$rs.test,pen.lasso2$sel.pen,times)
C.Uno[10] <- pred.measure('C.Uno',pen.elanet2$rs.train,pen.elanet2$rs.test,pen.elanet2$sel.pen,times)
C.Uno[11] <- pred.measure('C.Uno',pen.ridge2$rs.train,pen.ridge2$rs.test,pen.ridge2$sel.pen,times)
boxplot(C.Uno,col=rep(c("red","green","blue"),c(5,3,3)),xaxt="n",at=c(1:5,7:9,11:13),ylim=c(0.5,0.85),ylab="C.Uno")
axis(1,at=c(1:5,7:9,11:13),labels=FALSE)
text(x=c(1:5,7:9,11:13),par("usr")[3]-0.03,labels=methods, srt=20,pos=1, xpd=TRUE,cex=0.8)



