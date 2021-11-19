#==============================================================================
# FILENAME: surv_data_gen.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: generate simulated survival data to evaluate trial designs 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 3.6.1
#==============================================================================
#Notes: 


        


# =============================================================================
#install.packages("gsDesign")
#install.packages("coxed")
#install.packages("survminer")
#install.packages("survival")
#install.packages("lubridate")
#install.packages("ldbounds")
#install.packages("survRM2")
#install.packages("pseudo")
#install.packages("geepack")
#install.packages("xtable")
#install.packages("tidyverse")
#install.packages("cowplot")
library(xtable)
library(geepack)
library(coxed)
library(survminer)
library(gsDesign)
library(survival)
library(lubridate)
library(ldbounds)
library(survRM2)
library(pseudo)
library(tidyverse)
library(cowplot)


##### function to simulate survival times

survt <- function(n,               # desired sample size
                  
                  medposB,             # desired median for the survival distribution in the positive group, trtmnt B
                  distrposB,           # distr defines the distribution of survival times; positives, trtmt B
                                           # exp = exponential survival times
                                           # weib = weibull survival times
                                           # nonpar = non-parametric survival times created using NOT CODED CURRENTY
                                           # the coxed package and sim.survdata function
                  shapeposB=1,           # shape parmeter for the weibull distribution; positives, trtmt B
                  
                  medposA,             # desired median for the survival distribution in the positive group, trtmnt A
                  distrposA,           # distr defines the distribution of survival times; positives, trtmt A
                  shapeposA=1,           # shape parmeter for the weibull distribution; positives, trtmt A
                  
                  mednegB,             # desired median for the survival distribution in the negative group, trtmnt B
                  distrnegB,           # distr defines the distribution of survival times; negatives, trtmt B
                  shapenegB=1,           # shape parmeter for the weibull distribution; negatives, trtmt B
                  
                  mednegA,             # desired median for the survival distribution in the negative group, trtmnt A
                  distrnegA,           # distr defines the distribution of survival times; negatives, trtmt A
                  shapenegA=1           # shape parmeter for the weibull distribution; negatives, trtmt A
                  ){
  
# creating a data frame with specified sample size and number of datasets
surv_time <- data.frame(matrix(nrow=n, ncol=4))  
  
# simulates exponential survival times for positives, trt B
if (distrposB=="exp"){
  # scale parameter for positives, trt B
  lambdaposB <- (log(2)/medposB)
  # simulating survival times for positives, trt B
  surv_time[,1] <- rexp(n, rate=lambdaposB)
}
# simulates weibull survival times for positives, trt B
else if (distrposB=="weib"){
  lambdaposB <- medposB/((log(2))^(1/shapeposB))
  surv_time[,1] <- rweibull(n, shape=shapeposB, scale=lambdaposB)
}

# simulates exponential survival times for positives, trt A
if (distrposA=="exp"){
  # scale parameter for positives, trt A
  lambdaposA <- (log(2)/medposA)
  # simulating survival times for positives, trt A
  surv_time[,2] <- rexp(n, rate=lambdaposA)
}
# simulates weibull survival times for positives, trt A
else if (distrposA=="weib"){
  lambdaposA <- medposA/((log(2))^(1/shapeposA))
  surv_time[,2] <- rweibull(n, shape=shapeposA, scale=lambdaposA)
}


# simulates exponential survival times for negatives, trt B
if (distrnegB=="exp"){
  # scale parameter for negatives, trt B
  lambdanegB <- (log(2)/mednegB)
  # simulating survival times for negatives, trt B
  surv_time[,3] <- rexp(n, rate=lambdanegB)
}
# simulates weibull survival times for negatives, trt B
else if (distrnegB=="weib"){
  lambdanegB <- mednegB/((log(2))^(1/shapenegB))
  surv_time[,3] <- rweibull(n, shape=shapenegB, scale=lambdanegB)
}

# simulates exponential survival times for negatives, trt A
if (distrnegA=="exp"){
  # scale parameter for negatives, trt A
  lambdanegA <- (log(2)/mednegA)
  # simulating survival times for negatives, trt A
  surv_time[,4] <- rexp(n, rate=lambdanegA)
}
# simulates weibull survival times for negatives, trt A
else if (distrnegA=="weib"){
  lambdanegA <- mednegA/((log(2))^(1/shapenegA))
  surv_time[,4] <- rweibull(n, shape=shapenegA, scale=lambdanegA)
}

  # seeting the column names
  colnames(surv_time) <- c("SURV_POS_B", "SURV_POS_A","SURV_NEG_B", "SURV_NEG_A")
  
  #returning survival times
  return(surv_time)
}


simdat <- function(n,                    # total sample size
                   accru.rate,           # acrrual rate which follows a Poisson distribution 
                   loss.fu.rate,          # lost to follow-up rate; assumed equal in each strata
                   marker.pos,           # proportion of subjects who are marker positive
                   rand.arm,              # probability of being randomized to strategy arm                         
                   rand.pos,             # probability of being randomized to trt B in positives
                   rand.neg,             # probability of being randomized to trt B in negatives
                   rand,                  # probability of being randomized to trr B regardless of marker status
                   phys.choice.pos,        # proportion of strategy-based treatment assignments that agrees with physician's choice in positives
                                           # strategy is to treat the positives with trt B
                   phys.choice.neg,        # proportion of strategy-based treatment assignments that agrees with physician's choice in negatives
                                           # strategy is to treat the negtives with trt A
                   medposB,             # desired median for the survival distribution in the positive group, trtmnt B
                   distrposB,           # distr defines the distribution of survival times; positives, trtmt B
                                           # exp = exponential survival times
                                           # weib = weibull survival times
                                           # nonpar = non-parametric survival times created using NOT CODED CURRENTY
                                          # the coxed package and sim.survdata function
                   shapeposB,           # shape parmeter for the weibull distribution; positives, trtmt B
                   
                   medposA,             # desired median for the survival distribution in the positive group, trtmnt A
                   distrposA,           # distr defines the distribution of survival times; positives, trtmt A
                   shapeposA,           # shape parmeter for the weibull distribution; positives, trtmt A
                   
                   mednegB,             # desired median for the survival distribution in the negative group, trtmnt B
                   distrnegB,           # distr defines the distribution of survival times; negatives, trtmt B
                   shapenegB,           # shape parmeter for the weibull distribution; negatives, trtmt B
                   
                   mednegA,             # desired median for the survival distribution in the negative group, trtmnt A
                   distrnegA,           # distr defines the distribution of survival times; negatives, trtmt A
                   shapenegA           # shape parmeter for the weibull distribution; negatives, trtmt A
){
  # generating time between accruing each patient based on the accrual rate
  times <- data.frame(matrix(nrow=n, ncol=1))
  times[,1] <- rexp(n=n, rate=accru.rate)
  # generating study entry times based on time between recruiting each patient
  entry.time <- as.vector(replicate(n=n, 0))
  for (i in 2:n){
    entry.time[i] <- entry.time[(i-1)] + times[i,1]
  }
  
  # generating marker status for each patient, 1=marker positive, 0=marker negative
  marker.stat <- rbinom(n, size=1, prob=marker.pos)
  
  # geeting the character marker variable
  Marker <- NULL
  for (i in 1:n){
    if (marker.stat[i]==1){Marker[i]="+"}
    else {Marker[i]="-"}
  }
  
  # generating subject ids in order of accrual/screening
  id <- seq(from=1, to=n, by=1)
  
  # combining the above into a data frame
  tmp1 <- data.frame(cbind(id, entry.time, marker.stat, Marker))
  tmp1$id <- as.numeric(as.character((tmp1$id)))
  tmp1$entry.time <- as.numeric(as.character((tmp1$entry.time)))
  tmp1$marker.stat <- as.numeric(as.character((tmp1$marker.stat)))
  
  # getting randomized treatment assignment regardless of marker status
  trt <- rbinom(n, size=1, prob=rand)
  for (i in 1:n){
    if (trt[i]==1){tmp1$trt.rnd[i]="B" 
                  tmp1$trtn.rnd[i]=1}
    else {tmp1$trt.rnd[i]="A"
          tmp1$trtn.rnd[i]=0}
    
  }
  
  #separating data frame into positive/negative marker groups
  tmp_pos <- tmp1[tmp1$marker.stat==1,]
  npos <- length(tmp_pos[,1])
  tmp_neg <- tmp1[tmp1$marker.stat==0,]
  nneg <- length(tmp_neg[,1])
  
  # getting the randomized assignment treatment for each marker group
  trt_pos <- rbinom(npos, size=1, prob=rand.pos)
  trt_neg <- rbinom(nneg, size=1, prob=rand.neg)
  for (i in 1:npos){
    if (trt_pos[i]==1){tmp_pos$trt.ms[i]="B" 
                      tmp_pos$trtn.ms[i]=1}
    else {tmp_pos$trt.ms[i]="A" 
          tmp_pos$trtn.ms[i]=0}
  }
  for (i in 1:nneg){
    if (trt_neg[i]==1){tmp_neg$trt.ms[i]="B" 
                      tmp_neg$trtn.ms[i]=1 
                      }
    else {tmp_neg$trt.ms[i]="A" 
          tmp_neg$trtn.ms[i]=0 
          }
  }
  # generating the pysicians choice of treatment based on the phys.choice.pos and phys.choice.neg
  # also generating an indicator phys.choice that the physician choice agrees with the strategy
  phys_pos <- rbinom(npos, size=1, prob=phys.choice.pos)
  phys_neg <- rbinom(nneg, size=1, prob=phys.choice.neg)
  for (i in 1:npos){
    if (phys_pos[i]==1){tmp_pos$phys.trt[i]="B"
                        tmp_pos$phys.trtn[i]=1
                        }
    else {tmp_pos$phys.trt[i]="A"
          tmp_pos$phys.trtn[i]=0
          }
  }
  for (i in 1:nneg){
    if (phys_neg[i]==1){tmp_neg$phys.trt[i]="A" 
                        tmp_neg$phys.trtn[i]=0
                        }
    else {tmp_neg$phys.trt[i]="B"
          tmp_neg$phys.trtn[i]=1
         }
  }
  
  # rejoining the marker groups into single dataset
  simdata <- data.frame(rbind(tmp_pos, tmp_neg))
  
  #generating randomzied arm assignment, either the strategy arm or physician choice/randomized arm 
  trt.arm <- rbinom(n, size=1, prob=rand.arm)
  
  for (i in 1:length(simdata[,1])){
    if (trt.arm[i]==1){simdata$arm[i]="strat"}
    else {simdata$arm[i]="phys"}
  }
  
  # getting survival times
  surv.time <- survt(n=n, medposB=medposB, distrposB = distrposB, shapeposB=shapeposB,
                     medposA=medposA, distrposA = distrposA, shapeposA=shapeposA,
                     mednegB=mednegB, distrnegB = distrnegB, shapenegB=shapenegB,
                     mednegA=mednegA, distrnegA = distrnegA, shapenegA=shapenegA)
  # attaching the surv.times to the dataset
  for (i in 1:n){
    if (simdata$Marker[i]=="+"){simdata$surv.timeB[i] <- surv.time[i,1]
                                simdata$surv.timeA[i] <- surv.time[i,2]}
    else if (simdata$Marker[i]=="-"){simdata$surv.timeB[i] <- surv.time[i,3]
                                    simdata$surv.timeA[i] <- surv.time[i,4]}
  }
  
  
  # simulating loss to follow-up based on loss.f.rate
  simdata$event <- 1
  loss <- rbinom(n, size=1, prob=loss.fu.rate)
  simdata$loss2fu <- NA
  simdata$loss2futA <- NA
  simdata$loss2futB <- NA
  for (i in 1:n){
    if (loss[i]==1){simdata$loss2fu[i]=1
                    simdata$event[i]=0
                    # time when lost to follow-up is distributed uniformly from almst time 0 to their survival time
                    simdata$loss2futA[i]=runif(n=1, min=0.00000000001, max=simdata$surv.timeA[i])
                    simdata$surv.timeA[i]=simdata$loss2futA[i]
                    simdata$loss2futB[i]=runif(n=1, min=0.00000000001, max=simdata$surv.timeB[i])
                    simdata$surv.timeB[i]=simdata$loss2futB[i]}
    else {simdata$loss2fu[i]=0}
  }
  # getting the study time for failure/loss to follow up
    simdata$event.timeA <- simdata$entry.time + simdata$surv.timeA
    simdata$event.timeB <- simdata$entry.time + simdata$surv.timeB
  
  simdata <- simdata[order(simdata$id),]
  
 return(simdata)   
}

  











