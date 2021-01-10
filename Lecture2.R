###############################################
#USING SIMULATION FOR PSS WHEN DATA ARE GENERATED FROM NORMAL DIST
#PURPOSE: TO ESTIMATE THE SAMPLE SIZE THAT IS NEEDED TO ACHIEVE POWER OF 80% OR 90%
###############################################

B.sim<-1000	#No. of simulations
alpha.sim<-0.05
mu0.sim<-5.5
mu1.sim<-6.0
sd.sim<-1.4	
targ.power.sim<-0.8

N.sim<-30:60	#Sample size to try.

pwr.func<-function(N, B=B.sim, alpha=alpha.sim, mu0=mu0.sim, mu1=mu1.sim, sd=sd.sim)
{
  p.res<-c() #This is to save the p-values from each simulation
  
  for (j in 1:B){	
    set.seed(12300+123*N*j)		#Setting the seed number for random number generation, this is to make the result will be reproducible.
    sim.data<-rnorm(n=N, mean=mu1, sd=sd)
    p.res[j]<-t.test(x=sim.data,  alternative = c( "greater"),  mu = mu0, conf.level = 1-alpha)$"p.value"
  }
  tmp<-ifelse(p.res<alpha,1,0) 
  power.all<-sum(tmp)/B
  return(power.all)
}



res.pwr<-c()
for (i in 1:length(N.sim)){
  res.pwr[i]<-pwr.func(N=N.sim[i])
}

res.pwr1<-cbind(N.sim, res.pwr)
res.pwr1

tmp<-ifelse((res.pwr1[,"res.pwr"]-targ.power.sim)>0, (res.pwr1[,"res.pwr"]-targ.power.sim),100) #This is to identify the N with power closest to target power.
N.need.pwr1<-res.pwr1[which(tmp==min(tmp)),]
N.need.pwr1


plot(x=res.pwr1[,"N.sim"], y=res.pwr1[,"res.pwr"]) #Power curve



###############################################
#USING SIMULATION FOR PSS WHEN DATA ARE GENERATED FROM EXPONENTIAL DIST
###############################################

p.res<-c() #This is to save the p-values from each simulation
B<-1000	#No. of simulations
alpha<-0.05

mu0<-5.5
mu1<-6.0
sd<-1.4	
N<-30	#Sample size

for (i in 1:B){	
  set.seed(12300+123*i)		#Setting the seed number for random number generation, this is to make sure the result will be reproducible.
  sim.data<-rexp(n=N, rate=1/mu1)
  p.res[i]<-t.test(x=sim.data,  alternative = c( "greater"),  mu = mu0, conf.level = 1-alpha)$"p.value"
}


tmp<-ifelse(p.res<alpha,1,0) 
power.res<-sum(tmp)/B

power.res

###############################################
#USING SIMULATION FOR PSS FOR LINEAR RANDOM INTERCEPT MODEL
###############################################
library(MASS)
library(nlme)

#####Initialize parameters 
alpha<-0.05
Bsim<-500
effect<-0.1
corr<-0.4
beta0<-0.6
sd<-0.40

Nsim<-300	#Sample size (i.e total no. of families).
fam.dist<-c(1/3, 1/2, 1/6) #Probability of a family having 1, 2, or 3 kids carrying the genetic trait, respectively.

rmat<-matrix(c(1, corr, corr, corr, 1, corr, corr, corr,1), nrow=3, ncol=3, byrow=TRUE)
p.value<-c()	#Initialize p.value
for (i in 1:Bsim) {
  set.seed<-12300+i*11+Nsim*10
  n.fam<-rmultinom(n=1, size=Nsim, prob=fam.dist) #Numbers of families with 1, 2, or 3 kids carrying the trait, respectively.
  
  #####Generate X the exposure variable
  x.fam1<-rep(c(1,0,0), n.fam[1])	#the x vector for families with 1 kid carrying the trait.
  id.fam1<-rep(1:n.fam[1], each=3) #generate the ids for families with 1 kid carrying the trait.
  
  x.fam2<-rep(c(1,1,0), n.fam[2]) #the x vector for families with 2 kids carrying the trait.
  id.fam2<-rep(1:n.fam[2], each=3)+n.fam[1] #generate the ids for families with 1 kid carrying the trait.
  
  x.fam3<-rep(c(1,1,1), n.fam[3]) #the x vector for families with 3 kids carrying the trait.
  id.fam3<-rep(1:n.fam[3], each=3)+n.fam[1]+n.fam[2] #generate the ids for families with 1 kid carrying the trait.
  
  x.vec<-c(x.fam1, x.fam2, x.fam3)
  id.vec<-c(id.fam1, id.fam2, id.fam3) 
  
  v.mat<-sd*rmat
  
  ###For families with 1 kid carrying the trait
  set.seed<-12300+i*100+Nsim*1
  mean.fam1<-beta0+effect*c(1,0,0)
  y.fam1<-mvrnorm(n.fam[1], mean.fam1, v.mat)
  
  ###For families with 2 kids carrying the trait
  set.seed<-12300+2*i*100+Nsim*2
  mean.fam2<-beta0+effect*c(1,1,0)
  y.fam2<-mvrnorm(n.fam[2], mean.fam2, v.mat)
  
  ###For families with 3 kids carrying the trait
  set.seed<-12300+3*i*100+Nsim*3
  mean.fam3<-beta0+effect*c(1,1,1)
  y.fam3<-mvrnorm(n.fam[3], mean.fam3, v.mat)
  
  ###Run the LME
  y.vec<-c(c(t(y.fam1)), c(t(y.fam2)), c(t(y.fam3)))
  
  dataset<-data.frame(y.vec, x.vec, id.vec)	#Create the dataset
  data.formu<-groupedData(y.vec~x.vec|id.vec, data=dataset)	#Attacheh the formula to the dataset attribute
  res<-lme(data.formu, random=~1)		
  p.value[i]<-summary(res)$tTable["x.vec","p-value"]
}

pwr.func.res<-sum(ifelse(p.value<alpha,1,0))/Bsim
pwr.func.res