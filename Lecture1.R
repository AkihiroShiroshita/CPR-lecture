#Clinical prediction model
library(tidyverse)
library(gt)
library(mice)
library(rms)
library(mfp)
library(mgcv)
library(foreign)
library(ROCR)
library(tableone)
library(ICC)
library(lmerTest)
library(optimx)
library(devtools)
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
md.pattern(dat)
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
val.prob(logit=full$linear.predictor, y=full$y)
##Internal validation
validate(full, B=200) 
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
###Restricted cubic splines
lrm(death ~ rcs(age,4) + gender + rcs(hr,4)+rcs(rr,4), data = dat_c)
###Fractional polynomials
mfp(death ~ fp(age) + gender + fp(hr)+fp(rr), alpha = 0.2, select = 0.2, data = dat_c, family = binomial)
###Generalized additive model
gam(death ~ s(age) + gender+s(hr)+s(rr), data = dat_c, family = binomial)
###LASSO
glmmod <- glmnet(full8$x, y=full8$y, alpha=1, family="binomial")
###Internal validation by bootstrap
B <- 5000
for (i in 1:B) {
  dat_c_b <- dat_c[sample(nrow(dat_c), replace=T),]
  full.gam.gcv.B <- gam(death ~ s(age) + gender+s(hr)+s(rr), data = dat_c, family = binomial)
  val.prob(y=dat_c_b$death, logit=predict(full.gam.gcv.B))
}
###Multilevel model
###Random intercept model
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
mira <- with(dat_mi0, glm(death ~ age + gender + rr + hr,  family = binomial))
result <- summary(pool(mira))
exp(result[2,2])
exp(result[2,3])
###Another way
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
###Multilevel model
###Random intercept model
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