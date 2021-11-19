source("Programs/sim_source.R")


#debug(simdat)
set.seed(12)
check1 <- simdat(n=5000, 
                 accru.rate=20, 
                 loss.fu.rate=.02, 
                 marker.pos=.30, 
                 rand.arm=.5,
                 rand.pos=.5,
                 rand.neg=.5,
                 rand=.5,
                 phys.choice.pos=.2,
                 phys.choice.neg=.8,
                 
                 medposB=9,
                 distrposB="weib",
                 shapeposB=3,
                 
                 medposA=6,
                 distrposA="exp",
                 
                 mednegB=9,
                 distrnegB="exp",
                 
                 mednegA=12,
                 distrnegA="exp")

temp <- check1[order(check1$arm),]

x <- survdiff(Surv(surv.timeA, event) ~ arm, data=temp)
str(x)
x$exp

fit <- coxph(Surv(surv.timeA, event) ~ trt.ms, data=temp, cluster=id)
summary(fit)
summary(fit)$waldtest[3]

posB <- check1[check1$marker.stat==1,]
posA <- check1[check1$marker.stat==1,]
negB <- check1[check1$marker.stat==0,]
negA <- check1[check1$marker.stat==0,]

median(posB$surv.timeB)
median(posA$surv.timeA)
median(negB$surv.timeB)
median(negA$surv.timeA)

pseudo <- pseudosurv(check1$surv.time, check1$event, tmax=10)
ipseudo <- ipseudo <- 1-pseudo$pseudo
check2 <- data.frame(cbind(check1, ipseudo))
fit <- geese(formula=ipseudo ~ trt, data=check2, id=id, scale.fix=TRUE,
             family=gaussian, jack=TRUE, mean.link="cloglog", corstr="independence")
sum_fit <- round(cbind(mean=fit$beta, SE=sqrt(diag(fit$vbeta.ajs)), Z=fit$beta/sqrt(diag(fit$vbeta.ajs)),
                       PVal = 2-2*pnorm(abs(fit$beta/sqrt(diag(fit$vbeta.ajs))))),4)

fit1 <- survfit(Surv(surv.time, event) ~ trt, data=check1)
ggsurvplot(fit1, data=check1)

debug(sim.trial)
x <- sim.trial(data=check1,
          n=500,
          num_interim=1,
          int_timing=1,
          alpha=.025,
          low_err=0.1,
          bound.type=c(1,1),
          time.points=20,
          RM.times=20,
          design.type="modstrat",
          test.grp="subgrp"
          )

source("Programs/sim_source.R")
Sys.time()
#debug(eval.scen)
x <- eval.scen(estimand="subgrp", # declaring the estimand of interest; options are:
               # subgrp - estimating the effect within a subgroup
               # clin - estimating clinical utility
               # inter - estimating a differential treatment effect between subgroups
               reps=2, # of trials to simulate
               max_enroll=5000,  # maximum number enrolled over course of trial
               n=500, # total number of events per trial
               num_interim=1, # total number of analyses, including final
               int_timing=1, # vector of proportion of events for each interim analysis timing
               alpha=.025, # one-sided type 1 error rate (upper boundary)
               low_err=0.1,  # lower (futility) boundary error
               bound.type=c(1,1), # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively)
               time.points=15,    # value or vector of time points for which to compare survival probability at time points
               RM.times=15, # value or vector of time points to compare the restricted mean survival probability
               
               accru.rate=20,           # acrrual rate which follows a Poisson distribution with input rate 
               loss.fu.rate=.02,          # lost to follow-up rate; assumed equal in each strata
               marker.pos=0.4,           # proportion of subjects who are marker positive
               rand.arm=.5,             # probability of being randomized to strategy arm
               rand.pos=.5,             # probability of being randomized to trt B in positives
               rand.neg=.5,             # probability of being randomized to trt B in negatives
               rand=0.5,
               phys.choice.pos=.2,        # proportion of strategy-based treatment assignments that agrees with physician's choice in positives
               # strategy is to treat the positives with trt B
               phys.choice.neg=.8,        # proportion of strategy-based treatment assignments that agrees with physician's choice in negatives
               medposB=12,             # desired median for the survival distribution in the positive group, trtmnt B
               
               distrposB="exp",           # distr defines the distribution of survival times; positives, trtmt B
               # exp = exponential survival times
               # weib = weibull survival times
               shapeposB,           # shape parmeter for the weibull distribution; positives, trtmt B
               
               medposA=9,             # desired median for the survival distribution in the positive group, trtmnt A
               distrposA="exp",           # distr defines the distribution of survival times; positives, trtmt A
               shapeposA,           # shape parmeter for the weibull distribution; positives, trtmt A
               
               mednegB=10,             # desired median for the survival distribution in the negative group, trtmnt B
               distrnegB="exp",           # distr defines the distribution of survival times; negatives, trtmt B
               shapenegB,           # shape parmeter for the weibull distribution; negatives, trtmt B
               
               mednegA=15,             # desired median for the survival distribution in the negative group, trtmnt A
               distrnegA="exp",           # distr defines the distribution of survival times; negatives, trtmt A
               shapenegA           # shape parmeter for the weibull distribution; negatives, trtmt A
)
Sys.time()




test1 <- x$inter_strategy


test2 <- x$inter_stratify

filename_enrich <- paste("subgrp", "enrich", smple.size[2], "exp", 
                         "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep=".")
saveRDS(test2, file=paste("Results", "subgrp",filename_enrich, sep="/"))

test4 <- readRDS("Results/inter/test2.rds")

test3 <- x$inter_modstrat


debug(simdat)
set.seed(12)
check1 <- simdat(n=2000, 
                 accru.rate=20, 
                 loss.fu.rate=.01, 
                 marker.pos=.40, 
                 hr.trtB=.7, 
                 hr.pos=1.5, 
                 hr.neg=1.5, 
                 design.type="strategy",
                 rand.arm=.75,
                 rand.pos=.5,
                 rand.neg=.25,
                 rand=.9,
                 phys.choice.pos=.2,
                 phys.choice.neg=.8,
                 med=12,
                 distr="weib",
                 shape=2)
check1 <- check1[order(check1$trt, check1$arm),]
x <- survdiff(Surv(surv.time, event) ~ arm, data=check1)
x
check2 <- check1[check1$marker.stat==1,]
y <- survdiff(Surv(surv.time, event) ~ trt, data=check2)
y

# testing for coxhp
check2 <- check1
for (i in 1:length(check2[,1])){
  if (check2$trt.ms[i]=="B"){check2$surv.time[i]=check2$surv.timeB[i]}
   else {check2$surv.time[i]=check2$surv.timeA[i]}
}
fit <- coxph(Surv(surv.time, event) ~ trt, data=check2, cluster=id)
summary(fit)

check2 <- check1[check1$marker.stat==1,]
check2$surv.time <- check2$surv.timeB
check2$trt <- check2$trt.ms
y2 <- survfit(Surv(surv.time, event) ~ trt, data=check2)
summary(y2)$table["median"]

y1 <- survfit(Surv(surv.time, event) ~ trt, data=check2)


# testing to get the pseudomean

pseudom <- pseudomean(check1$surv.time, check1$event, tmax=12)
test <- cbind(check1, pseudo=pseudo)

fit <- geese(pseudo ~ trt, data=test, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)

sum.fit <- round(cbind(mean=fit$beta, SE=sqrt(diag(fit$vbeta.ajs)), Z=fit$beta/sqrt(diag(fit$vbeta.ajs)),
                       PVal = 2-2*pnorm(abs(fit$beta/sqrt(diag(fit$vbeta.ajs))))),4)
sum.fit

# testing the pseudo observation code for RMST and point-wise survival estimates
M <- 1
check2 <- check1[order(check1$surv.time),]
check2 <- check2[check2$event==1,]
times <- NULL
for (i in 1:M){
  num <- ceiling((i/M)*length(check2$surv.time))-5
  times[i] <- round(check2$surv.time[num])
}
pseudo <- pseudosurv(check1$surv.time, check1$event, tmax=12)
pseudop <- pseudo$pseudo
check3 <- cbind(id=check1$id, time=check1$surv.time, event=check1$event, test)
c <- NULL
for (j in 4:ncol(check3)){
  c <- rbind(c, cbind(check1, pseudo=check3[,j], tpseudo=times[j-3]))
}
c <- c[order(c$id),]
c$tpseudo <- as.factor(c$tpseudo)
c$ipseudo <- 1-c$pseudo

xfit <- geese(formula=test ~ trt, data=c, id=id, scale.fix=TRUE, 
              family=gaussian, jack=TRUE, mean.link="cloglog", corstr="independence")
sum_xfit.f <- round(cbind(mean=xfit$beta, SE=sqrt(diag(xfit$vbeta.fij)), Z=xfit$beta/sqrt(diag(xfit$vbeta.fij)),
                          PVal = 2-2*pnorm(abs(xfit$beta/sqrt(diag(xfit$vbeta.fij))))),4)
sum_xfit.f

sum_xfit.a <- round(cbind(mean=xfit$beta, SE=sqrt(diag(xfit$vbeta.ajs)), Z=xfit$beta/sqrt(diag(xfit$vbeta.ajs)),
                          PVal = 2-2*pnorm(abs(xfit$beta/sqrt(diag(xfit$vbeta.ajs))))),4)
sum_xfit.a


# example from pseudo code paper
data(bmt)
names(bmt)[c(1,3:6,13,20)] <- c("disease", "tdfs", "relapse", "trm", "dfs", "age", "fab")
bmt$disease <- as.factor(bmt$disease)

cutoffs <- c(50,105,170,280,530)
pseudo <- pseudosurv(bmt$tdfs, bmt$dfs, tmax=cutoffs)

pseudo1 <- data.frame(cbind(bmt$tdfs, bmt$dfs, pseudo$pseudo))
names(pseudo1)[c(1,2)] <- c("time", "event")
pseudo1 <- pseudo1[order(pseudo1$time),]

b <- NULL
for (j in 3:ncol(pseudo1)){
  b <- rbind(b, cbind(bmt, pseudo1=pseudo1[,j], tpseudo=cutoffs[j-2], id=1:nrow(bmt)))
}
b <- b[order(b$id),]



time   = check1$surv.time
status = check1$event
arm    = check1$arm


test <- rmst2(time=check1$surv.time, status=check1$event, arm=check1$trtn)


D = rmst2.sample.data()
nrow(D)

time   = D$time
status = D$status
arm    = D$arm
test <- rmst2(time, status, arm, tau=10)

posB <- check1[check1$marker.stat==1 & check1$trt=="B",]
posA <- check1[check1$marker.stat==1 & check1$trt=="A",]
negB <- check1[check1$marker.stat==0 & check1$trt=="B",]
negA <- check1[check1$marker.stat==0 & check1$trt=="A",]

median(posB$surv.time)
median(posA$surv.time)
median(negB$surv.time)
median(negA$surv.time)

check2 <- check1[check1$marker.stat==1 & check1$trt=="B" & check1$phys.trt=="B",]

fit1 <- survfit(Surv(surv.time, event) ~ trt, data=check1)
ggsurvplot(fit1, data=check1)

test <- readRDS(file="Results/clin/clin.modstrat.250.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.9.rds")
test2 <- readRDS(file="Results/clin/clin.stratify.250.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.9.rds")
test3 <- readRDS(file="Results/clin/clin.strategy.250.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.9.rds")






