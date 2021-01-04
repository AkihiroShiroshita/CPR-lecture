#Clinical prediction model
library(tidyverse)
library(DescTools) #Winsorize(x)
library(gt)
library(mice)
library(rms)
library(mfp)
library(mgcv)
library(VIM)
library(naniar)
library(ROCR)
library(tableone)
library(ICC)
library(lmerTest)
library(glmnet)
library(Hmisc)
library(ROCR)
library(optimx)
pkgFile <- "norm2_2.0.3.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
library(norm2)
#Data preparation
dat <- read_csv("C:/Users/akihi/Downloads/Clinical_prediction_model_lecture/CPR-lecture/Sample_data.csv",
                locale = locale(encoding = "SHIFT-JIS"),
                col_types = cols(
                  id = col_double(),
                  age = col_double(),
                  gender = col_logical(),
                  wbc = col_double(),
                  eosi_p = col_double(),
                  sbp = col_double(),
                  dbp = col_double(),
                  bun = col_double(),
                  ams = col_logical(),
                  hr = col_double(),
                  death = col_logical(),
                  adl = col_logical(),
                  hospital = col_factor(),
                  hospitalterm = col_double(),
                  rr = col_double()
                ),
                guess_max = 1500, 
                na = "NA") 
dat %>% glimpse()
par(mfrow=c(1,1))
na.patterns <- naclus(dat)
plot(na.patterns, ylab="Fraction of NAs in common")
par(mfrow=c(2,2))
naplot(na.patterns)
na.pattern(dat)
aggr(dat)
ggplot(dat,
       aes(x = hospital,
           y = sbp)) +
  geom_miss_point()
summary(dat)
dat_e <- dat %>% 
  mutate(age_cat = ifelse((gender=="TRUE" & age>70) | (gender=="FALSE" & age >= 75), 1, 0)) %>% 
  mutate(bun_cat = ifelse(bun >= 21, 1, 0)) %>% 
  mutate(rr_cat = ifelse(rr >= 30, 1, 0)) %>% 
  mutate(bp_cat = ifelse(sbp <= 90 | dbp <= 50, 1, 0)) %>% 
  select(death, age_cat, bun_cat, rr_cat, bp_cat, ams) %>% 
  mutate_all(.funs = ~ as.logical(.))
dat_e %>% glimpse()
dat_c <- dat %>% 
  na.omit()
dat_ec <- dat_e %>% 
  na.omit()
summary(dat)
#Complete case analysis
##Development
full <- lrm(death ~ age+gender+hr+rr, data=dat_c, x = TRUE, y = TRUE)
full #Check c-statistics
full_int <- lrm(death ~ age*(gender+hr+rr), data=dat_c, x = TRUE, y = TRUE)
full_up <- lrm(death ~ age+gender+hr+rr+sbp, data=dat_c, x = TRUE, y = TRUE)
#update(full, .~. + age*gender + age*hr + age*rr)
#Breier score
B <- mean((full$y) * (1-plogis(full$linear.predictors))^2 + 
                (1-full$y) * plogis(full$linear.predictors)^2)
B
Bmax <- mean(full$y) * (1-mean(full$y))
Bmax
Bscaled <- 1-B/Bmax
Bscaled
#Pearson R2
cor(x=plogis(full$linear.predictors), y=full$y)^2
#Calibration
val.prob(logit=full$linear.predictor, y=full$y)
#Set up Hosmer-Lemeshow test
#p:predicted probability
#Y:outcome variable
#g:number of groups to calculate H-L (10 is default)
hl.ext2<-function(p,y,g=10, df=g-1)
{
  matres	<-matrix(NA,nrow=g,ncol=5)
  sor	<-order(p)
  p	<-p[sor]
  y	<-y[sor]
  groep	<-cut2(p,g=g)									#g more or less equal sized groups
  
  len		<-tapply(y,groep,length)					#n per group
  sump	<-tapply(p,groep,sum)						#expected per group
  sumy	<-tapply(y,groep,sum)						#observed per group
  meanp	<-tapply(p,groep,mean)						#mean probability per group
  meany	<-tapply(y,groep,mean)						#mean observed per group
  matres	<-cbind(len,meanp,meany, sump, sumy)
  contr<-((sumy-sump)^2)/(len*meanp*(1-meanp))		#contribution per group to chi square
  chisqr<-sum(contr)									#chi square total
  pval<-1-pchisq(chisqr,df)							#p-value corresponding to chi square with df degrees of freedom 
  cat("\nChi-square",chisqr," p-value", pval,"\n")
  dimnames(matres)	<-list(c(1:g),Cs(n,avg(p),avg(y), Nexp, Nobs))
  result  <- list(table(groep), matres,chisqr,pval) 
}
hl.ext2(p=plogis(full$linear.predictor),y=full$y, g=10, df=9) #Df: external -> g-1, internal -> g-2
#Disrimination
#Discrimination slope
mean(plogis(lp[full$y==1])) - mean(plogis(lp[full$y==0]))
#ROC curve
#Lorenz curve
pred.full <- prediction(full$linear.predictor, full$y)
perf <- performance(pred.full, "tpr", "fpr")
plot(perf)
auc.tmp <- performance(pred.full,"auc")
auc <- as.numeric(auc.tmp@y.values)
auc
perf <- performance(pred.full, "sens", "spec")
plot(perf)
perf.full <- performance(pred.full,"fnr","rnp") # Lorenz curve data+plot
plot(perf.full, xlab="Fraction not resected", ylab="Tumors missed")
abline(a=0,b=1)
#Box plot
boxplot(plogis(full$linear.predictors)~full$y,
        ylab="Predicted risk", xlab="Death",ylim=c(0,1))
#Scatter plot
plot(x= plogis(full$linear.predictors), y= plogis(full_up$linear.predictors),
     xlab="Predicted risk without SBP", ylab="Predicted risk with SBP",
     las=1,pty='s',xlim=c(0,1),ylim=c(0,1),
     pch=(full$y+.5)*2) # pch 1 and 3
abline(a=0,b=1)
#CI around c-statistics
cstatNo <- rcorr.cens(full$linear.predictors, full$y) 
cat(cstatNo[1], "[", cstatNo[1]-1.96/2*cstatNo[3], " - ", cstatNo[1]+1.96/2*cstatNo[3],"]")
cstatYes <- rcorr.cens(full_up$linear.predictors, full$y) 
cat(cstatYes[1], "[", cstatYes[1]-1.96/2*cstatYes[3], " - ", cstatYes[1]+1.96/2*cstatYes[3],"]")
cstatYes - cstatNo
cstatvAL <- rcorr.cens(lp.val, val$death) 
cat(cstatvAL[1], "[", cstatvAL[1]-1.96/2*cstatvAL[3], " - ", cstatvAL[1]+1.96/2*cstatvAL[3],"]")
par(mfrow = c(1,3), pty='s', mar=c(4.2,4,4,1),cex=0.95, font =1, col=1)
val.prob.ci(logit=full$linear.predictor,y=full$y, pl=T,smooth=T,logistic.cal=F, g=10,
            xlab="Predicted risk",
            ylab="Observed death",riskdist='predicted',
            d1lab="Death", d0lab="Alive", dist.label=-0.95, cutoff=.2)
#NRI&IDI
slopeNo <- mean(plogis(full$linear.predictors[full$y==1])) - mean(plogis(full$linear.predictors[full$y==0]))
slopeYes <- mean(plogis(full_up$linear.predictors[full_up$y==1])) - mean(plogis(full_up$linear.predictors[full_up$y==0]))
slopeYes - slopeNo #IDI
NoSBP <- cut(plogis(full$linear.predictors),breaks = c(0,.2,1), labels = c("No SBP, <=20%","No SBP, >20%"))
YesSBP <- cut(plogis(full_up$linear.predictors),breaks = c(0,.2,1), labels = c("With SBP, <=20%","With SBP, >20%"))
ValReclass <- cut(plogis(lp.val),breaks = c(0,.2,1), labels = c("Val, <=20%","Val, >20%"))
table(NoSBP)
table(YesSBP)
tabReclas <- table(NoSBP, YesSBP)
tabReclas
tabReclas[2]/(tabReclas[2]+tabReclas[4]) # reclassified among high risk
tabReclas[3]/(tabReclas[1]+tabReclas[3]) # reclassified among low risk
table(NoSBP, YesSBP, full$y) #By outcome
tab <- table(NoSBP, YesSBP, full$y )
tab[5:8]/(tab[1:4]+tab[5:8])
table(NoSBP, full$y )
tabNo <- table(NoSBP, full$y )
tabNo[3:4]+tabNo[1:2]
tabNo[3:4]/(tabNo[3:4]+tabNo[1:2])
table(YesSBP, full$y )
tabYes <- table(YesSBP, full$y )
tabYes[3:4]+tabYes[1:2]
tabYes[3:4]/(tabYes[3:4]+tabYes[1:2])
improveProb(x1=as.numeric(NoSBP)-1, x2=as.numeric(YesSBP)-1, y=full$y) #NRI and IDI calculations with Harrell's functions
##Internal validation
validate(full, B=200)
B <- 200
for (i in 1:B) {
  dat_c_b <- dat_c[sample(nrow(dat_c), replace=T),]
  full.gam.gcv.B <- gam(death ~ s(age) + gender+s(hr)+s(rr), data = dat_c, family = binomial)
  val.prob(y=dat_c_b$death, logit=predict(full.gam.gcv.B))
}
#Bootstrap discrimination slope
nrowB	<- nrow(dat_c)
B <- 200            # 200 bootstraps
matB <- matrix(NA,nrow=B,ncol=3) # Matrix for results
dimnames(matB) <- list(c(1:B), Cs(Slopeapp, Slopetest, optimism ))
for (i in 1:B) {
  if (i%%10==0) cat("Start Bootstrap sample nr", i, "\n")
  Brows <- sample(nrowB,replace=T)
  # Bsample is bootstrap sample from development set
  Bsample	<- dat_c[Brows,]
  devfull <- lrm(death ~ age+gender+rr+hr, data=Bsample,linear.predictors=T, x=T, y=T)
  matB[i,1] <- mean(plogis(devfull$linear.predictors[devfull$y==1])) - 
    mean(plogis(devfull$linear.predictors[devfull$y==0]))
  lp  <- full$x %*% devfull$coef[2:length(devfull$coef)] + devfull$coef[1] # lp with coefs from bootstrap
  matB[i,2] <- mean(plogis(lp[full$y==1])) - mean(plogis(lp[full$y==0]))  # Testing on original sample
  
}
cat("\n\nEinde bereikt!\n\n")
matB[,3] <- matB[,1] - matB[,2] # optimism
# matB bekijken mean en SD
apply(matB,2,mean)
apply(matB,2,sd)
##External validation
dat_ec <-dat_ec %>% 
  mutate(score = age_cat+bun_cat+rr_cat+bp_cat+ams) %>% 
  mutate(score = as.factor(score))
table1 <- CreateTableOne(vars = "death",
                         strata = "score",
                         data = dat_ec)
summary(table1)
ext.full <- lrm(death ~ age_cat+bun_cat+rr_cat+bp_cat+ams,
                data=dat_ec, x=TRUE, y=TRUE)
###Penalization
#Restricted cubic splines
lrm(death ~ rcs(age,4) + gender + rcs(hr,4)+rcs(rr,4), data = dat_c)
#Fractional polynomials
mfp(death ~ fp(age) + gender + fp(hr)+fp(rr), alpha = 1, data = dat_c, family = binomial)
#Generalized additive model
gam(death ~ s(age) + gender+s(hr)+s(rr), data = dat_c, family = binomial)
#LASSO
glmmod <- glmnet(full$x, y=full$y, alpha=1, family="binomial")
##Multilevel model
#Random intercept model
ICCest(as.factor(hospital),death, data=dat_c, alpha=0.05, CI.type="Smith")
dat_c$age.m <- ave(dat_c$age, dat_c$hospital)
dat_c$age.cwc <- dat_c$age - dat_c$age.m
dat_c$gender.m <- ave(dat_c$gender, dat_c$hospital)
dat_c$gender.cwc <- dat_c$gender - dat_c$gender.m
dat_c$hr.m <- ave(dat_c$hr, dat_c$hospital)
dat_c$hr.cwc <- dat_c$hr - dat_c$hr.m
dat_c$rr.m <- ave(dat_c$rr, dat_c$hospital)
dat_c$rr.cwc <- dat_c$rr - dat_c$rr.m
remodel <- glmer(death ~ age.cwc+gender.cwc+hr.cwc+rr.cwc+(1|hospital),
                 control=glmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
                 family = "binomial",
                 data=dat_c)
summary(remodel)
##Multiple imputation
dat_mi <- dat %>% 
  select(-id)
dat_mi0 <- mice(dat_mi, maxit = 0)
dat_mi0$method
dat_mi0$predictorMatrix
dat_mi10 <- mice(dat, m = 10, maxit = 20, printFlag = FALSE, seed = 1234) 
plot(dat_mi10)
densityplot(dat_mi0)
mira <- with(dat_mi0, glm(death ~ age + gender + rr + hr,  family = binomial))
result <- summary(pool(mira))
exp(result[2,2])
exp(result[2,3])
#Another way
est <- se <- vector(length = dat_mi10$m, mode = "list")
for (i in 1:dat_mi10$m) {
  com <- complete(dat_mi10, i)
  fit2 <- glm(death ~ age+ gender + rr + hr, com, family = binomial(link = "logit"))
  s <- summary(fit2)
  est[[i]] <- s$coefficients[2:5, 1]
  se[[i]] <- s$coefficients[2:5, 2]
}
miinf <- miInference(est, se)
print(miinf)
exp(miinf$est)
exp(miinf$std.err)
##Multilevel model with MI dataset
#Random intercept model
est <- se <- vector(length = dat_mi10$m, mode = "list")
for (i in 1:dat_mi10$m) {
  com <- complete(dat_mi10, i)
  com$age.m <- ave(com$age, com$hospital)
  com$age.cwc <- com$age - com$age.m
  com$gender.m <- ave(com$gender, com$hospital)
  com$gender.cwc <- com$gender - com$gender.m
  com$hr.m <- ave(com$hr, com$hospital)
  com$hr.cwc <- com$hr - com$hr.m
  com$rr.m <- ave(com$rr, com$hospital)
  com$rr.cwc <- com$rr - com$rr.m
  remodel <- glmer(death ~ age.cwc+gender.cwc+hr.cwc+rr.cwc+(1|hospital),
                   control=glmerControl(optimizer="optimx",optCtrl=list(method="nlminb")),
                   family = "binomial",
                   data=com)
  s <- summary(fit2)
  est[[i]] <- s$coefficients[2:5, 1]
  se[[i]] <- s$coefficients[2:5, 2]
}
miinf <- miInference(est, se)
print(miinf)
exp(miinf$est)
exp(miinf$std.err)
##Try to creat MI data sets and apply bootstrapping to each of them