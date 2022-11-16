###### this is a debug program; is not use for simulation or analysis

### this is a throw away program

source("Programs/sim_source.R")

filename_enrich <- paste("subgrp", "enrich", smple.size[2], "exp", 
                         "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep=".")
saveRDS(test2, file=paste("Results", "subgrp",filename_enrich, sep="/"))

test1 <- readRDS("Results/clin/clin.strategy.500.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.PC.0.rds")
test2 <- readRDS("Results/clin/Run1/clin.strategy.500.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.PC.0.rds")
test3 <- readRDS("Results/subgrp/subgrp.enrich.type1.550.accru.10.exp.Bpos.12.Apos.12.Bneg.12.Aneg.9.M..75.LTFU.0.02.rds")

data <- test3

test3 <- x$inter_modstrat

undebug(analyze.exp)
x <- readRDS("Results/subgrp/subgrp.enrich.1000.accru.25.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.75.LTFU.0.02.rds")
debug(analyze.exp)
check <- analyze.exp(data <- x,
                     estimand="subgrp",
                     design.type="enrich",
                     medianBpos=21,
                     medianBneg=9,
                     medianApos=9,
                     medianAneg=12,
                     marker.pos=.75,
                     phys.choice=0)



undebug(fintable)
fintable(caption="test",
         estimand="inter",
         design.type="stratify",
         smple.size=smple.size <- seq(from=250, to=1000, by=50),
         accru=25,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5)

debug(power.curve)
check <- power.curve(caption="test",
                     estimand="inter",
                     design.type="stratify",
                     smple.size=smple.size <- seq(from=250, to=1000, by=50),
                     accru=25,
                     LTFU=.02,
                     medianBpos=15,
                     medianBneg=9,
                     medianApos=9,
                     medianAneg=12,
                     marker.pos=.25,
                     phys.choice=.25,
                     rpos=0.5            # randomization to treatment b in the positives
)

par(mfrow=c(4,4))
fig9.25.5 <- power.curve(caption="test",
                         estimand="inter",
                         design.type="stratify",
                         smple.size=smple.size <- seq(from=250, to=1000, by=50),
                         accru=25,
                         LTFU=.02,
                         medianBpos=9,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.25,
                         phys.choice=.25,
                         rpos=0.5            # randomization to treatment b in the positives
)
fig12.25.5 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=12,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.5            # randomization to treatment b in the positives
)
fig15.25.5 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=15,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.5            # randomization to treatment b in the positives
)

fig21.25.5 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=21,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.5            # randomization to treatment b in the positives
)

fig9.25.67 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=9,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.67            # randomization to treatment b in the positives
)
fig12.25.67 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=12,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.67            # randomization to treatment b in the positives
)
fig15.25.67 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=15,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.67            # randomization to treatment b in the positives
)
fig21.25.67 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=21,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.67            # randomization to treatment b in the positives
)
fig9.5.5 <- power.curve(caption="test",
                          estimand="inter",
                          design.type="stratify",
                          smple.size=smple.size <- seq(from=250, to=1000, by=50),
                          accru=25,
                          LTFU=.02,
                          medianBpos=9,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.5,
                          phys.choice=.25,
                          rpos=0.5            # randomization to treatment b in the positives
)
fig12.5.5 <- power.curve(caption="test",
                        estimand="inter",
                        design.type="stratify",
                        smple.size=smple.size <- seq(from=250, to=1000, by=50),
                        accru=25,
                        LTFU=.02,
                        medianBpos=12,
                        medianBneg=9,
                        medianApos=9,
                        medianAneg=12,
                        marker.pos=.5,
                        phys.choice=.25,
                        rpos=0.5            # randomization to treatment b in the positives
)
fig15.5.5 <- power.curve(caption="test",
                        estimand="inter",
                        design.type="stratify",
                        smple.size=smple.size <- seq(from=250, to=1000, by=50),
                        accru=25,
                        LTFU=.02,
                        medianBpos=15,
                        medianBneg=9,
                        medianApos=9,
                        medianAneg=12,
                        marker.pos=.5,
                        phys.choice=.25,
                        rpos=0.5            # randomization to treatment b in the positives
)
fig21.5.5 <- power.curve(caption="test",
                        estimand="inter",
                        design.type="stratify",
                        smple.size=smple.size <- seq(from=250, to=1000, by=50),
                        accru=25,
                        LTFU=.02,
                        medianBpos=21,
                        medianBneg=9,
                        medianApos=9,
                        medianAneg=12,
                        marker.pos=.5,
                        phys.choice=.25,
                        rpos=0.5            # randomization to treatment b in the positives
)
fig9.5.67 <- power.curve(caption="test",
                        estimand="inter",
                        design.type="stratify",
                        smple.size=smple.size <- seq(from=250, to=1000, by=50),
                        accru=25,
                        LTFU=.02,
                        medianBpos=9,
                        medianBneg=9,
                        medianApos=9,
                        medianAneg=12,
                        marker.pos=.5,
                        phys.choice=.25,
                        rpos=0.67            # randomization to treatment b in the positives
)
fig12.5.67 <- power.curve(caption="test",
                         estimand="inter",
                         design.type="stratify",
                         smple.size=smple.size <- seq(from=250, to=1000, by=50),
                         accru=25,
                         LTFU=.02,
                         medianBpos=12,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.5,
                         phys.choice=.25,
                         rpos=0.67            # randomization to treatment b in the positives
)
fig15.5.67 <- power.curve(caption="test",
                         estimand="inter",
                         design.type="stratify",
                         smple.size=smple.size <- seq(from=250, to=1000, by=50),
                         accru=25,
                         LTFU=.02,
                         medianBpos=15,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.5,
                         phys.choice=.25,
                         rpos=0.67            # randomization to treatment b in the positives
)
fig21.5.67 <- power.curve(caption="test",
                         estimand="inter",
                         design.type="stratify",
                         smple.size=smple.size <- seq(from=250, to=1000, by=50),
                         accru=25,
                         LTFU=.02,
                         medianBpos=21,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.5,
                         phys.choice=.25,
                         rpos=0.67            # randomization to treatment b in the positives
)

test <- ggarrange(fig9.25.5, fig12.25.5, fig15.25.5, fig21.25.5,
                  fig9.25.67, fig12.25.67, fig15.25.67, fig21.25.67,
                  fig9.5.5, fig12.5.5, fig15.5.5, fig21.5.5,
                  fig9.5.67, fig12.5.67, fig15.5.67, fig21.5.67,
                  ncol=4, nrow=4)

annotate_figure(test, 
                top=text_grob("Median Survival for Positive Patients on Treatment B
        9 months                 12 months                     15 months                  21 months", size=16),
                left=text_grob("Proportion of Biomarker Positive Patients
0.5                                   0.25", size=16, rot=90),
                right=text_grob("Randomization of Positive Patients to B
   0.5                     0.67                 0.5               0.67", size=16, rot=-90))
                
test

test <- plot_grid(NULL,NULL,NULL,NULL,NULL,
                  NULL,fig9.25.5, fig12.25.5, fig15.25.5, fig21.25.5,
                  NULL,fig9.25.67, fig12.25.67, fig15.25.67, fig21.25.67,
                  NULL,fig9.5.5, fig12.5.5, fig15.5.5, fig21.5.5,
                  NULL,fig9.5.67, fig12.5.67, fig15.5.67, fig21.5.67,
          #labels=c("A","B","C","D"),
          ncol=5, nrow=5)

test2 <- test + draw_figure_label(label="Median Survival in Positive Patients on Treatment B


                         9 months                12 months              15 months              21 months", position="top.right", size=16)

test3 <- test2 +  draw_figure_label(label="M=1 Proportion
                    0.5                                                    0.25       
                    Randomization Probability:              Randomization Probability 
                    0.67                0.5                 0.67               0.5", 
                    position="bottom.left", angle=-90, size=16)
test3

test <- c(1,NA,3)

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
                 phys.choice.pos=1,
                 phys.choice.neg=1,
                 
                 medposB=21,
                 distrposB="exp",
                 shapeposB=3,
                 
                 medposA=9,
                 distrposA="exp",
                 
                 mednegB=9,
                 distrnegB="exp",
                 
                 mednegA=12,
                 distrnegA="exp")
clin

debug(clin.util.comp)
test <- clin.util.comp(n=400,
                       accru=25,
                       ltfu=.02,
                       pos=0.5,
                       pc.pos=1,
                       pc.neg=1,
                       medposB=21,
                       medposA=9,
                       mednegB=9,
                       mednegA12,
                       pseudo.diff=12)

test1 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rds")
test2 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")
test3 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.12.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rds")
test4 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.12.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")
test5 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rds")
test6 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")


length(test1[test1$lr.pval.shih<=.05,1])/1000
length(test1[test1$lr.pval.strat<=.05,1])/1000

mean(test1$mean.sd.shih)
mean(test1$mean.sd.strat)

mean(test1$rmst.shih)
mean(test1$rmst.strat) - 1.96*(sd(test1$rmst.strat)/sqrt(length(test1$rmst.strat)))
mean(test1$rmst.strat) + 1.96*(sd(test1$rmst.strat)/sqrt(length(test1$rmst.strat)))

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

test <- sum(posA$surv.timeA>=29.9)/length(posA$surv.timeA)

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

test1 <- readRDS("Results/clin/clin.strategy.500.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.PC.0.rds")
test2 <- readRDS("Results/clin/Run1/clin.strategy.500.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.PC.0.rds")
test3 <- readRDS("Results/subgrp/subgrp.stratify.550.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")

test3 <- x$inter_modstrat

undebug(analyze.exp)
check <- analyze.exp(data <- readRDS("Results/inter/inter.stratify.1000.accru.25.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rpos.0.67.rds"),
                     estimand="inter",
                     design.type="stratify",
                     medianBpos=21,
                     medianBneg=9,
                     medianApos=9,
                     medianAneg=12,
                     marker.pos=.5,
                     phys.choice.pos=0.25,
                     phys.choice.neg=0.25)


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

test1 <- readRDS("Results/inter/inter.stratify.300.accru.25.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rpos.0.5.rds")
test2 <- readRDS("Results/subgrp/subgrp.stratify.600.accru.25.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")
test3 <- readRDS("Results/subgrp/subgrp.enrich.950.accru.25.exp.Bpos.15.Apos.9.Bneg.9.Aneg.12.M.0.75.LTFU.0.02.rds")
test4 <- x$subgrp_enrich

true.rat <- exp(-(log(2)/15)*51.8)/exp(-(log(2)/9)*51.8)



source("Programs/sim_source.R")

# selecting a vecotr of sample sizes
smple.size <- 950



Bpos <- 15
Bneg <- 9
Apos <- 9
Aneg <- 12
M <- .75
LTFU <- .02


for (j in 1:length(Bpos)){
  for (k in 1:length(M)){
    
    # selecting survival distributions and median survival for each of the 4 groups
    distr_Bpos <- "exp"
    medianBpos <- Bpos[j]
    
    distr_Bneg <- "exp"
    medianBneg <- Bneg
    
    distr_Apos <- "exp"
    medianApos <- Apos
    
    distr_Aneg <- "exp"
    medianAneg <- Aneg
    
    Mpos <- M[k]
    
    
    accru <- 25
    
    set.seed(12)
    #getting results for each smple.size based on 1000 replications
    for (i in 1:length(smple.size)){
      debug(eval.scen)
      x <- eval.scen(estimand="subgrp", 
                     reps=100, 
                     max_enroll=10000,  
                     n=smple.size[i], 
                     num_interim=1, 
                     int_timing=1, 
                     alpha=.025, 
                     low_err=0.1,  
                     bound.type=c(1,1), 
                     # time.points=time.point,    
                     # RM.times=time.point, 
                     
                     accru.rate=accru,           
                     loss.fu.rate=LTFU,          
                     marker.pos=Mpos,           
                     rand.arm=.5,             
                     rand.pos=.5,            
                     rand.neg=.5,             
                     rand=0.5,
                     phys.choice.pos=.5,        
                     phys.choice.neg=.5,
                     
                     medposB=medianBpos,             
                     distrposB=distr_Bpos,           
                     shapeposB,         
                     
                     medposA=medianApos,             
                     distrposA=distr_Apos,           
                     shapeposA,           
                     
                     mednegB=medianBneg,             
                     distrnegB=distr_Bneg,           
                     shapenegB,           
                     
                     mednegA=medianAneg,             
                     distrnegA=distr_Aneg,          
                     shapenegA          
      )
      
      temp_enrich <- x$subgrp_enrich
      temp_strat <- x$subgrp_stratify
      # temp_mod <- x$subgrp_modstrat
      
      filename_enrich <- paste("subgrp", "enrich", smple.size[i], "accru", accru, "exp", 
                               "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                               "M", Mpos, "LTFU", LTFU, "rds", sep=".")
      
      filename_strat <- paste("subgrp", "stratify", smple.size[i], "accru", accru,"exp", 
                              "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
                              "M", Mpos, "LTFU", LTFU, "rds", sep=".")
      # filename_mod <- paste("subgrp", "modstrat", smple.size[i], "accru", accru,"exp", 
      #                         "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
      #                       "M", Mpos, "LTFU", LTFU, "rds",sep=".")
      
      saveRDS(temp_enrich, file=paste("Results", "subgrp",filename_enrich, sep="/"))
      saveRDS(temp_strat, file=paste("Results", "subgrp",filename_strat, sep="/"))
      # saveRDS(temp_mod, file=paste("Results", "subgrp",filename_mod, sep="/"))
    }
  }}



