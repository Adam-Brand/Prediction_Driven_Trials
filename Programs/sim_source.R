#==============================================================================
# FILENAME: sim_source.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: source functions to simulate trial results 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 3.6.1
#==============================================================================
#Notes: 
### uses functions within gsDesign package in R; some explanations of used functions below




# =============================================================================
source("Programs/surv_data_gen.R")


sim.trial <- function(data, # input survival data from our data generating function
                      n, # total number of events
                      num_interim, # total number of analyses, including final
                      int_timing, # vector of proportion of events for each interim analysis timing
                      alpha=.025, # one-sided type 1 error rate (upper boundary)
                      low_err=0.1,  # lower (futility) boundary error
                      bound.type=c(1,1),  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
                     # time.points,    # value or vector of time points for which to compare survival probability at time points
                     # RM.times, # value or vector of time points to compare the restricted mean survival probability
                      design.type,    # type of design being analyzed choices are:
                                          # enrich
                                          # stratify
                                          # strategy
                                          # modstrat
                      test.grp,       # the group comparison that the logrank test is based on, choices are
                                         # pos - trt B vs A in the positives
                                         # neg - trt B vs A in the negatives
                                         # clin - clinical utility; strategy arm vs phys choice arm
                                         # trtB - pos vs neg for trtB patients
                                         # trtA - pos vs neg for trtA patients
                                         # inter - 4 way comparison for interaction; k=4 logrank test
                     pseudot.ratio    # the set time point to compare survival ratios, set to 90th percentile for lowest survival group
                      ){
  
  # getting the pvalue cutoffs for the analysis times
  x <- bounds(t=int_timing, # the percentage of information times, last one always =1
              iuse=bound.type,             # indicates O'brien fleming boundaries for both upper and lower boundaries
              alpha=c(low_err, alpha))     # the alpas for the lower and upper boundaries, respectively
  # saving the upper pvalue cutoffs for each interim
  upper.bounds <- x$upper.bounds
  
  # saving the upper pvalue cutoffs for each interim
  lower.bounds <- x$lower.bounds
  
  # number of events needed at each timing
  n_events <- ceiling(int_timing*n)
  
  # setting treatment received and event time based on study design, test group and estmiand of interest
  if (design.type=="enrich" | design.type=="stratify"){
    data$trt <- data$trt.ms
    data$trtn <- data$trtn.ms
  }
  
  else if (design.type=="strategy"){
    for (i in 1:length(data[,1])){
      if (data$arm[i]=="phys"){data$trt[i]=data$phys.trt[i]
                              data$trtn[i]=data$phys.trtn[i]}
      else if (data$arm[i]=="strat" & data$marker.stat[i]==1){data$trt[i]="B"
                                                          data$trtn[i]=1}
      else if (data$arm[i]=="strat" & data$marker.stat[i]==0){data$trt[i]="A"
                                                          data$trtn[i]=0}
    }
  }
  
  else if (design.type=="modstrat"){
    for (i in 1:length(data[,1])){
      if (data$arm[i]=="phys"){data$trt[i]=data$trt.rnd[i]
                              data$trtn[i]=data$trtn.rnd[i]}
      else if (data$arm[i]=="strat" & data$marker.stat[i]==1){data$trt[i]="B"
                                                          data$trtn[i]=1}
      else if (data$arm[i]=="strat" & data$marker.stat[i]==0){data$trt[i]="A"
                                                          data$trtn[i]=0}
    }
  }
  
  for (i in 1:length(data[,1])){
    if (data$trt[i]=="B"){data$surv.time[i]=data$surv.timeB[i]
                          data$event.time[i]=data$event.timeB[i]
                          if (data$marker.stat[i]==1) {data$strategy[i]=1}
                          else {data$strategy[i]=0}
                          if (data$phys.trt[i]=="B"){data$phys.choice[i]=1}
                          else {data$phys.choice[i]=0}
    }
    else if (data$trt[i]=="A"){data$surv.time[i]=data$surv.timeA[i]
                              data$event.time[i]=data$event.timeA[i]
                              if (data$marker.stat[i]==1) {data$strategy[i]=0}
                              else {data$strategy[i]=1}
                              if (data$phys.trt[i]=="A"){data$phys.choice[i]=1}
                              else {data$phys.choice[i]=0}
    }
  }
  
  # ordering the data by event time
  trial_data <- data[order(data$event.time),]
  
  # getting the study times for each interim analysis
  stdy_times <- NULL
  
  ### computing stdy_times for different trial types and testing comparisons
  
  if (design.type=="enrich"){
    # keeping only the positive group
    trial_data <- trial_data[trial_data$marker.stat==1,]
    # removing censored observations to get the study time after x events
    event_dat <- trial_data[trial_data$event==1 & trial_data$marker.stat==1,]
    for (i in 1:num_interim){
      stdy_times[i] <- event_dat[n_events[i],"event.time"]
    }
  }
  else if (design.type=="stratify" | design.type=="modstrat"){
    if (test.grp=="pos"){
      event_dat <- trial_data[trial_data$event==1 & trial_data$marker.stat==1,]
        for (i in 1:num_interim){
          stdy_times[i] <- event_dat[n_events[i],"event.time"]
        } 
    }
    else if (test.grp=="neg"){
      event_dat <- trial_data[trial_data$event==1 & trial_data$marker.stat==0,]
      for (i in 1:num_interim){
        stdy_times[i] <- event_dat[n_events[i],"event.time"]
      } 
    }
    else if (test.grp=="clin"){
      event_dat <- trial_data[trial_data$event==1,]
      for (i in 1:num_interim){
        stdy_times[i] <- event_dat[n_events[i],"event.time"]
      } 
    }
    else if (test.grp=="trtB"){
      event_dat <- trial_data[trial_data$event==1 & trial_data$trt=="B",]
      for (i in 1:num_interim){
        stdy_times[i] <- event_dat[n_events[i],"event.time"]
      } 
    }
    else if (test.grp=="trtA"){
      event_dat <- trial_data[trial_data$event==1 & trial_data$trt=="A",]
      for (i in 1:num_interim){
        stdy_times[i] <- event_dat[n_events[i],"event.time"]
      } 
    }
    else if (test.grp=="inter"){
      event_dat <- trial_data[trial_data$event==1,]
      for (i in 1:num_interim){
        stdy_times[i] <- event_dat[n_events[i],"event.time"]
      } 
    }
  }
  
  else if (design.type=="strategy"){
      event_dat <- trial_data[trial_data$event==1,]
      for (i in 1:num_interim){
        stdy_times[i] <- event_dat[n_events[i],"event.time"]
      } 
  }
  
  i <- 1
  cnt <- 1
  while (cnt!=0){
    # subsetting on subjects who entered the trial by the analysis time
    temp <- trial_data[trial_data$entry.time <= stdy_times[i],]
      # changing event indicators to censored for subjects who have not had their event by the interim analysis
      # also changing to observed survival time at the interim time to reflect current information
      for (j in 1:length(temp[,1])){
        if (temp$event.time[j] > stdy_times[i])
        {temp$event[j]=0
          temp$surv.time[j] <- stdy_times[i] - temp$entry.time[j]}
      }
    nevents <- n_events[i]
    stdy_time <- stdy_times[i]
    # getting number enrolled at the study time
    nenroll <- length(temp[unique(temp$id),1]) 
    
    ### testing for the enrichment design
    if (design.type=="enrich"){
      # conducting the logrank test
      temp <- temp[order(temp$trt),]
      x <- survdiff(Surv(surv.time, event) ~ trt, data=temp)
    }
    
    ### testing for the stratify and modstrat design
    else if (design.type=="stratify" | design.type=="modstrat"){
      if (test.grp=="pos"){
        temp <- temp[order(temp$trt),]
        # conducting the logrank test
        x <- survdiff(Surv(surv.time, event) ~ trt, data=temp, subset=marker.stat==1)
      }
      else if (test.grp=="neg"){
        temp <- temp[order(temp$trt),]
        # conducting the logrank test
        x <- survdiff(Surv(surv.time, event) ~ trt, data=temp, subset=marker.stat==0)
      }
      else if (test.grp=="trtA"){
        temp <- temp[order(temp$Marker),]
        # conducting the logrank test
        x <- survdiff(Surv(surv.time, event) ~ Marker, data=temp, subset=trt=="A")
      }
      else if (test.grp=="trtB"){
        temp <- temp[order(temp$Marker),]
        # conducting the logrank test
        x <- survdiff(Surv(surv.time, event) ~ Marker, data=temp, subset=trt=="B")
      }
      else if (test.grp=="inter"){
        temp <- temp[order(temp$trt, temp$Marker),]
        # conducting the logrank test
        x <- survdiff(Surv(surv.time, event) ~ trt + Marker, data=temp)
      }
      else if (test.grp=="clin"){
        # reworking the data to include arm for either phys choice or strategy arm
        temp.phys <- temp[temp$phys.choice==1,]
        temp.phys$arm <- "phys"
        temp.strat <- temp[temp$strategy==1,]
        temp.strat$arm <- "strat"
        temp <- data.frame(rbind(temp.phys, temp.strat))
        temp <- temp[order(temp$arm),]
        # conducting the logrank test
        x <- survdiff(Surv(surv.time, event) ~ arm, data=temp)
      }
    }
    else if (design.type=="strategy"){
        temp <- temp[order(temp$arm),]
        x <- survdiff(Surv(surv.time, event) ~ arm, data=temp)
    }
      # extracting the chi square test stat from logrank test
      chisq <- x$chisq
      ### getting the expected number of deaths for the proper contrast
      ### temp is ordered above to always provide the proper contrast in the 2nd element
      ### for the inter test, directionality has no meaning
      exp_num <- x$exp[2] - x$obs[2]
      # extracting the p-value form the logrank test
      pval <- 1-pchisq(x$chisq, length(x$n) - 1)
      # getting the efficacy and futility boundaries
      upper.bound <- upper.bounds[i]
      lower.bound <- lower.bounds[i]
      # indicator for stopping for efficacy
      if ((design.type=="stratify" | design.type=="modstrat") & test.grp=="inter" & sqrt(chisq) > upper.bound)
          {stop_eff=1
          cnt <- 0}
      else if (sqrt(chisq) > upper.bound & exp_num > 0)
          {stop_eff=1
          cnt <- 0} 
      else {stop_eff=0}
      # indicator for stopping for futility
      if ((design.type=="stratify" | design.type=="modstrat") & test.grp=="inter")
          {stop_fut=0}
      else if (sqrt(chisq) > abs(lower.bound) & exp_num <= 0)
          {stop_fut=1
          cnt <- 0} 
      else {stop_fut=0}
      if (i==num_interim | stop_eff==1 | stop_fut==1)
          {cnt <- 0}
      else {i <- i+1}
  }
  
  # computing the pseudo observations and pseudomean
  pseudo <- pseudosurv(temp$surv.time, temp$event, tmax = pseudot.ratio)
  ipseudo <- pseudo$pseudo[,ncol(pseudo$pseudo)]
  pseudom <- pseudomean(temp$surv.time, temp$event)
  pseudo_time <- max(temp$surv.time)
  temp <- data.frame(cbind(temp, ipseudo=ipseudo, pseudom=pseudom))
  temp_pos <- temp[temp$marker.stat==1,]
  temp_neg <- temp[temp$marker.stat==0,]
  temp_phys <- temp[temp$arm=="phys",]
  temp_strat <- temp[temp$arm=="strat",]
  # declaring vectors for use below
  ratio <- NULL
  ratio_pos <- NULL
  ratio_neg <- NULL
  ratio_int <- NULL
  ratio_arm <- NULL
  RM <- NULL
  RM_pos <- NULL
  RM_neg <- NULL
  RM_int <- NULL
  RM_arm <- NULL
  HR <- NULL
  HR_pos <- NULL
  HR_neg <- NULL
  HR_int <- NULL
  HR_arm <- NULL
  
  ######################### update below here to account for different design types
  if (design.type=="enrich"){
    # extracting median survival for each treatment group
    y <- survfit(Surv(surv.time, event) ~ trt, data=temp)
    medA <- summary(y)$table["trt=A","median"]
    medB <- summary(y)$table["trt=B","median"]
    # getting overall median for positives
    y2 <- survfit(Surv(surv.time, event) ~ 1, data=temp)
    medpos <- summary(y2)$table["median"]
    # getting the ratio of survival probabilities at specific time points
   # M <- length(time.points)
   # for (i in 1:M){
      fit <- geese(formula=ipseudo ~ trt, data=temp, id=id, scale.fix=TRUE, 
                   family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.a <- round(cbind(mean=fit$beta, SE=sqrt(diag(fit$vbeta.ajs)), Z=fit$beta/sqrt(diag(fit$vbeta.ajs)),
                                PVal = 2-2*pnorm(abs(fit$beta/sqrt(diag(fit$vbeta.ajs))))),4)
      ratio_mean <- exp(sum_fit.a["trtB","mean"])
      ratio_lower <- exp(sum_fit.a["trtB","mean"] - (sum_fit.a["trtB","SE"]*qnorm(1-alpha)))
      ratio_upper <- exp(sum_fit.a["trtB","mean"] + (sum_fit.a["trtB","SE"]*qnorm(1-alpha)))
      ratio_pval <- sum_fit.a["trtB","PVal"]
      ratio <- c(ratio, ratio_mean, ratio_lower,ratio_upper, ratio_pval)
   # }
    # getting restricted mean absolute difference, B-A
   # N <- length(RM.times)
   # for (i in 1:N){
      fit.rm <- geese(pseudom ~ trt, data=temp, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm <- round(cbind(mean=fit.rm$beta, SE=sqrt(diag(fit.rm$vbeta.ajs)), Z=fit.rm$beta/sqrt(diag(fit.rm$vbeta.ajs)),
                             PVal = 2-2*pnorm(abs(fit.rm$beta/sqrt(diag(fit.rm$vbeta.ajs))))),4)
      # difference in restricted mean B - A
      RMST <- sum.fit.rm["trtB","mean"]
      RMST_lower <- sum.fit.rm["trtB","mean"] - (sum.fit.rm["trtB","SE"]*qnorm(1-alpha))
      RMST_upper <- sum.fit.rm["trtB","mean"] + (sum.fit.rm["trtB","SE"]*qnorm(1-alpha))
      RMST_pval <- sum.fit.rm["trtB","PVal"]
      RM <- c(RM, RMST, RMST_lower,RMST_upper, RMST_pval)
   # }
    
    # getting the coxph estimated hazard ratio
    fit.hr <- coxph(Surv(surv.time, event) ~ trt, data=temp, cluster=id)
    hr.est <- exp(fit.hr$coefficients["trtB"])
    hr.lower <- exp(fit.hr$coefficients["trtB"] - (summary(fit.hr)$coefficients[4]*qnorm(1-alpha)))
    hr.upper <- exp(fit.hr$coefficients["trtB"] + (summary(fit.hr)$coefficients[4]*qnorm(1-alpha)))
    hr.pval <- summary(fit.hr)$waldtest[3]
    
    HR <- c(hr.est, hr.lower, hr.upper, hr.pval)
      
    # combining results
    result <- c(nenroll, nevents, stdy_time,
                       pval, chisq, upper.bound,lower.bound,
                      exp_num, stop_eff, stop_fut, medA, medB, medpos,
                       HR, ratio, RM, pseudo_time, pseudot.ratio)
  }
  else if (design.type=="stratify" | design.type=="modstrat"){
    # extracting median survival for each group
    y1 <- survfit(Surv(surv.time, event) ~ trt, data=temp)
    medA <- summary(y1)$table["trt=A","median"]
    medB <- summary(y1)$table["trt=B","median"]
    y2 <- survfit(Surv(surv.time, event) ~ trt + Marker, data=temp)
    medAneg <- summary(y2)$table["trt=A, Marker=-","median"]
    medApos <- summary(y2)$table["trt=A, Marker=+","median"]
    medBneg <- summary(y2)$table["trt=B, Marker=-","median"]
    medBpos <- summary(y2)$table["trt=B, Marker=+","median"]
    y3 <- survfit(Surv(surv.time, event) ~ Marker, data=temp)
    medneg <- summary(y3)$table["Marker=-","median"]
    medpos <- summary(y3)$table["Marker=+","median"]
    y4 <- survfit(Surv(surv.time, event) ~ arm, data=temp)
    medphys <- summary(y4)$table["arm=phys","median"]
    medstrat <- summary(y4)$table["arm=strat","median"]
    # getting the ratio of survival probabilities at specific time points
   # M <- length(time.points)
   # for (i in 1:M){
      # getting ratio of B vs A in the positives
      fit_pos <- geese(formula=ipseudo ~ trt, data=temp_pos, id=id, scale.fix=TRUE, 
                   family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.a <- round(cbind(mean=fit_pos$beta, SE=sqrt(diag(fit_pos$vbeta.ajs)), Z=fit_pos$beta/sqrt(diag(fit_pos$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_pos$beta/sqrt(diag(fit_pos$vbeta.ajs))))),4)
      ratio_mean_pos <- exp(sum_fit.a["trtB","mean"])
      ratio_lower_pos <- exp(sum_fit.a["trtB","mean"] - (sum_fit.a["trtB","SE"]*qnorm(1-alpha)))
      ratio_upper_pos <- exp(sum_fit.a["trtB","mean"] + (sum_fit.a["trtB","SE"]*qnorm(1-alpha)))
      ratio_pval_pos <- sum_fit.a["trtB","PVal"]
      ratio_pos <- c(ratio_pos, ratio_mean_pos, ratio_lower_pos,ratio_upper_pos, ratio_pval_pos)
      # getting ratio of B vs A in the negatives
      fit_neg <- geese(formula=ipseudo ~ trt, data=temp_neg, id=id, scale.fix=TRUE, 
                       family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.b <- round(cbind(mean=fit_neg$beta, SE=sqrt(diag(fit_neg$vbeta.ajs)), Z=fit_neg$beta/sqrt(diag(fit_neg$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_neg$beta/sqrt(diag(fit_neg$vbeta.ajs))))),4)
      ratio_mean_neg <- exp(sum_fit.b["trtB","mean"])
      ratio_lower_neg <- exp(sum_fit.b["trtB","mean"] - (sum_fit.b["trtB","SE"]*qnorm(1-alpha)))
      ratio_upper_neg <- exp(sum_fit.b["trtB","mean"] + (sum_fit.b["trtB","SE"]*qnorm(1-alpha)))
      ratio_pval_neg <- sum_fit.b["trtB","PVal"]
      ratio_neg <- c(ratio_neg, ratio_mean_neg, ratio_lower_neg,ratio_upper_neg, ratio_pval_neg)
      # getting ratio of B vs A overall
      fit <- geese(formula=ipseudo ~ trt, data=temp, id=id, scale.fix=TRUE, 
                       family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.c <- round(cbind(mean=fit$beta, SE=sqrt(diag(fit$vbeta.ajs)), Z=fit$beta/sqrt(diag(fit$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit$beta/sqrt(diag(fit$vbeta.ajs))))),4)
      ratio_mean <- exp(sum_fit.c["trtB","mean"])
      ratio_lower <- exp(sum_fit.c["trtB","mean"] - (sum_fit.c["trtB","SE"]*qnorm(1-alpha)))
      ratio_upper <- exp(sum_fit.c["trtB","mean"] + (sum_fit.c["trtB","SE"]*qnorm(1-alpha)))
      ratio_pval <- sum_fit.c["trtB","PVal"]
      ratio <- c(ratio, ratio_mean, ratio_lower,ratio_upper, ratio_pval)
      # getting ratio of interaction term
      fit_int <- geese(formula=ipseudo ~ trt*Marker, data=temp, id=id, scale.fix=TRUE, 
                   family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.d <- round(cbind(mean=fit_int$beta, SE=sqrt(diag(fit_int$vbeta.ajs)), Z=fit_int$beta/sqrt(diag(fit_int$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_int$beta/sqrt(diag(fit_int$vbeta.ajs))))),4)
      ratio_mean_int <- exp(sum_fit.d["trtB:Marker+","mean"])
      ratio_lower_int <- exp(sum_fit.d["trtB:Marker+","mean"] - (sum_fit.d["trtB:Marker+","SE"]*qnorm(1-alpha)))
      ratio_upper_int <- exp(sum_fit.d["trtB:Marker+","mean"] + (sum_fit.d["trtB:Marker+","SE"]*qnorm(1-alpha)))
      ratio_pval_int <- sum_fit.d["trtB:Marker+","PVal"]
      ratio_int <- c(ratio_int, ratio_mean_int, ratio_lower_int,ratio_upper_int, ratio_pval_int)
      # getting ratio between arms
      fit_arm <- geese(formula=ipseudo ~ arm, data=temp, id=id, scale.fix=TRUE, 
                       family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.e <- round(cbind(mean=fit_arm$beta, SE=sqrt(diag(fit_arm$vbeta.ajs)), Z=fit_arm$beta/sqrt(diag(fit_arm$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_arm$beta/sqrt(diag(fit_arm$vbeta.ajs))))),4)
      ratio_mean_arm <- exp(sum_fit.e["armstrat","mean"])
      ratio_lower_arm <- exp(sum_fit.e["armstrat","mean"] - (sum_fit.e["armstrat","SE"]*qnorm(1-alpha)))
      ratio_upper_arm <- exp(sum_fit.e["armstrat","mean"] + (sum_fit.e["armstrat","SE"]*qnorm(1-alpha)))
      ratio_pval_arm <- sum_fit.e["armstrat","PVal"]
      ratio_arm <- c(ratio_arm, ratio_mean_arm, ratio_lower_arm,ratio_upper_arm, ratio_pval_arm)
   # }
    
    # getting restricted means absolute difference
   # N <- length(RM.times)
   # for (i in 1:N){
      # getting absolute difference in B vs A overall
      fit.rm <- geese(pseudom ~ trt, data=temp, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm.a <- round(cbind(mean=fit.rm$beta, SE=sqrt(diag(fit.rm$vbeta.ajs)), Z=fit.rm$beta/sqrt(diag(fit.rm$vbeta.ajs)),
                                PVal = 2-2*pnorm(abs(fit.rm$beta/sqrt(diag(fit.rm$vbeta.ajs))))),4)
      
      # difference in restricted mean B - A
      RMST <- sum.fit.rm.a["trtB","mean"]
      RMST_lower <- sum.fit.rm.a["trtB","mean"] - (sum.fit.rm.a["trtB","SE"]*qnorm(1-alpha))
      RMST_upper <- sum.fit.rm.a["trtB","mean"] + (sum.fit.rm.a["trtB","SE"]*qnorm(1-alpha))
      RMST_pval <- sum.fit.rm.a["trtB","PVal"]
      RM <- c(RM, RMST, RMST_lower,RMST_upper, RMST_pval)
      
      # getting absolute difference in B vs A in positives
      fit_pos.rm <- geese(pseudom ~ trt, data=temp_pos, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm.b <- round(cbind(mean=fit_pos.rm$beta, SE=sqrt(diag(fit_pos.rm$vbeta.ajs)), Z=fit_pos.rm$beta/sqrt(diag(fit_pos.rm$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_pos.rm$beta/sqrt(diag(fit_pos.rm$vbeta.ajs))))),4)
    
      RMST_pos <- sum.fit.rm.b["trtB","mean"]
      RMST_lower_pos <- sum.fit.rm.b["trtB","mean"] - (sum.fit.rm.b["trtB","SE"]*qnorm(1-alpha))
      RMST_upper_pos <- sum.fit.rm.b["trtB","mean"] + (sum.fit.rm.b["trtB","SE"]*qnorm(1-alpha))
      RMST_pval_pos <- sum.fit.rm.b["trtB","PVal"]
      RM_pos <- c(RM_pos, RMST_pos, RMST_lower_pos,RMST_upper_pos, RMST_pval_pos)
      
      # getting absolute difference in B vs A in negatives
      fit_neg.rm <- geese(pseudom ~ trt, data=temp_neg, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm.c <- round(cbind(mean=fit_neg.rm$beta, SE=sqrt(diag(fit_neg.rm$vbeta.ajs)), Z=fit_neg.rm$beta/sqrt(diag(fit_neg.rm$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_neg.rm$beta/sqrt(diag(fit_neg.rm$vbeta.ajs))))),4)
      
      RMST_neg <- sum.fit.rm.c["trtB","mean"]
      RMST_lower_neg <- sum.fit.rm.c["trtB","mean"] - (sum.fit.rm.c["trtB","SE"]*qnorm(1-alpha))
      RMST_upper_neg <- sum.fit.rm.c["trtB","mean"] + (sum.fit.rm.c["trtB","SE"]*qnorm(1-alpha))
      RMST_pval_neg <- sum.fit.rm.c["trtB","PVal"]
      RM_neg <- c(RM_neg, RMST_neg, RMST_lower_neg,RMST_upper_neg, RMST_pval_neg)
      
      # getting absolute difference in B vs A for interaction term
      fit_int.rm <- geese(pseudom ~ trt*Marker, data=temp, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm.d <- round(cbind(mean=fit_int.rm$beta, SE=sqrt(diag(fit_int.rm$vbeta.ajs)), Z=fit_int.rm$beta/sqrt(diag(fit_int.rm$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_int.rm$beta/sqrt(diag(fit_int.rm$vbeta.ajs))))),4)
      
      RMST_int <- sum.fit.rm.d["trtB:Marker+","mean"]
      RMST_lower_int <- sum.fit.rm.d["trtB:Marker+","mean"] - (sum.fit.rm.d["trtB:Marker+","SE"]*qnorm(1-alpha))
      RMST_upper_int <- sum.fit.rm.d["trtB:Marker+","mean"] + (sum.fit.rm.d["trtB:Marker+","SE"]*qnorm(1-alpha))
      RMST_pval_int <- sum.fit.rm.d["trtB:Marker+","PVal"]
      RM_int <- c(RM_int, RMST_int, RMST_lower_int,RMST_upper_int, RMST_pval_int)
      
      # getting absolute difference in B vs A by arm
      fit_arm.rm <- geese(pseudom ~ arm, data=temp, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm.e <- round(cbind(mean=fit_arm.rm$beta, SE=sqrt(diag(fit_arm.rm$vbeta.ajs)), Z=fit_arm.rm$beta/sqrt(diag(fit_arm.rm$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit_arm.rm$beta/sqrt(diag(fit_arm.rm$vbeta.ajs))))),4)
      
      RMST_arm <- sum.fit.rm.e["armstrat","mean"]
      RMST_lower_arm <- sum.fit.rm.e["armstrat","mean"] - (sum.fit.rm.e["armstrat","SE"]*qnorm(1-alpha))
      RMST_upper_arm <- sum.fit.rm.e["armstrat","mean"] + (sum.fit.rm.e["armstrat","SE"]*qnorm(1-alpha))
      RMST_pval_arm <- sum.fit.rm.e["armstrat","PVal"]
      RM_arm <- c(RM_arm, RMST_arm, RMST_lower_arm, RMST_upper_arm, RMST_pval_arm)
   # }
    
    # getting the HRs between groups
    
    # getting HR between B vs A overall
    fit.hr <- coxph(Surv(surv.time, event) ~ trt, data=temp, cluster=id)
    hr.est <- exp(fit.hr$coefficients["trtB"])
    hr.lower <- exp(fit.hr$coefficients["trtB"] - (summary(fit.hr)$coefficients[4]*qnorm(1-alpha)))
    hr.upper <- exp(fit.hr$coefficients["trtB"] + (summary(fit.hr)$coefficients[4]*qnorm(1-alpha)))
    hr.pval <- summary(fit.hr)$waldtest[3]
    HR <- c(hr.est, hr.lower, hr.upper, hr.pval)
    
    # getting HR between B vs A in the positives
    fit.hr.pos <- coxph(Surv(surv.time, event) ~ trt, data=temp_pos, cluster=id)
    hr.est.pos <- exp(fit.hr.pos$coefficients["trtB"])
    hr.lower.pos <- exp(fit.hr.pos$coefficients["trtB"] - (summary(fit.hr.pos)$coefficients[4]*qnorm(1-alpha)))
    hr.upper.pos <- exp(fit.hr.pos$coefficients["trtB"] + (summary(fit.hr.pos)$coefficients[4]*qnorm(1-alpha)))
    hr.pval.pos <- summary(fit.hr.pos)$waldtest[3]
    HR_pos <- c(hr.est.pos, hr.lower.pos, hr.upper.pos, hr.pval.pos)
    
    # getting HR between B vs A in the negatives
    fit.hr.neg <- coxph(Surv(surv.time, event) ~ trt, data=temp_neg, cluster=id)
    hr.est.neg <- exp(fit.hr.neg$coefficients["trtB"])
    hr.lower.neg <- exp(fit.hr.neg$coefficients["trtB"] - (summary(fit.hr.neg)$coefficients[4]*qnorm(1-alpha)))
    hr.upper.neg <- exp(fit.hr.neg$coefficients["trtB"] + (summary(fit.hr.neg)$coefficients[4]*qnorm(1-alpha)))
    hr.pval.neg <- summary(fit.hr.neg)$waldtest[3]
    HR_neg <- c(hr.est.neg, hr.lower.neg, hr.upper.neg, hr.pval.neg)
    
    # getting HR for the interaction term
    fit.hr.int <- coxph(Surv(surv.time, event) ~ trt*Marker, data=temp, cluster=id)
    hr.est.int <- exp(fit.hr.int$coefficients["trtB:Marker+"])
    hr.lower.int <- exp(fit.hr.int$coefficients["trtB:Marker+"] - (summary(fit.hr.int)$coefficients[3,4]*qnorm(1-alpha)))
    hr.upper.int <- exp(fit.hr.int$coefficients["trtB:Marker+"] + (summary(fit.hr.int)$coefficients[3,4]*qnorm(1-alpha)))
    hr.pval.int <- summary(fit.hr.int)$coefficients[3,6]
    HR_int <- c(hr.est.int, hr.lower.int, hr.upper.int, hr.pval.int)
    
    # getting HR between arms
    fit.hr.arm <- coxph(Surv(surv.time, event) ~ arm, data=temp, cluster=id)
    hr.est.arm <- exp(fit.hr.arm$coefficients["armstrat"])
    hr.lower.arm <- exp(fit.hr.arm$coefficients["armstrat"] - (summary(fit.hr.arm)$coefficients[4]*qnorm(1-alpha)))
    hr.upper.arm <- exp(fit.hr.arm$coefficients["armstrat"] + (summary(fit.hr.arm)$coefficients[4]*qnorm(1-alpha)))
    hr.pval.arm <- summary(fit.hr.arm)$waldtest[3]
    HR_arm <- c(hr.est.arm, hr.lower.arm, hr.upper.arm, hr.pval.arm)
    
    # combining results
    result <- c(nenroll, nevents, stdy_time, pval, chisq, upper.bound,lower.bound,
                exp_num, stop_eff, stop_fut,medA, medB,medAneg, medApos, medBneg, medBpos, medneg, medpos,
                medphys, medstrat, HR, HR_pos, HR_neg, HR_int, HR_arm, 
                ratio, ratio_pos, ratio_neg, ratio_int, ratio_arm,
                RM, RM_pos, RM_neg, RM_int, RM_arm, pseudo_time, pseudot.ratio)
    
  }
  else if (design.type=="strategy"){
    # extracting median survival for each treatment group
    y <- survfit(Surv(surv.time, event) ~ arm, data=temp)
    medphys <- summary(y)$table["arm=phys","median"]
    medstrat <- summary(y)$table["arm=strat","median"]
    y2 <- survfit(Surv(surv.time, event) ~ arm + Marker, data=temp)
    medstratpos <- summary(y2)$table["arm=strat, Marker=+","median"]
    medstratneg <- summary(y2)$table["arm=strat, Marker=-","median"]
    medphyspos <- summary(y2)$table["arm=phys, Marker=+","median"]
    medphysneg <- summary(y2)$table["arm=phys, Marker=-","median"]
    y3 <- survfit(Surv(surv.time, event) ~ arm + trt, data=temp)
    medstratB <- summary(y3)$table["arm=strat, trt=B","median"]
    medstratA <- summary(y3)$table["arm=strat, trt=A","median"]
    medphysB <- summary(y3)$table["arm=phys, trt=B","median"]
    medphysA <- summary(y3)$table["arm=phys, trt=A","median"]
    
    # getting the ratio of survival probabilities at specific time points
   # M <- length(time.points)
   # for (i in 1:M){
      fit.rat <- geese(formula=ipseudo ~ arm, data=temp, id=id, scale.fix=TRUE, 
                   family=gaussian, jack=TRUE, mean.link="log", corstr="independence")
      sum_fit.rat <- round(cbind(mean=fit.rat$beta, SE=sqrt(diag(fit.rat$vbeta.ajs)), Z=fit.rat$beta/sqrt(diag(fit.rat$vbeta.ajs)),
                               PVal = 2-2*pnorm(abs(fit.rat$beta/sqrt(diag(fit.rat$vbeta.ajs))))),4)
      ratio_mean_arm <- exp(sum_fit.rat["armstrat","mean"])
      ratio_lower_arm <- exp(sum_fit.rat["armstrat","mean"] - (sum_fit.rat["armstrat","SE"]*qnorm(1-alpha)))
      ratio_upper_arm <- exp(sum_fit.rat["armstrat","mean"] + (sum_fit.rat["armstrat","SE"]*qnorm(1-alpha)))
      ratio_pval_arm <- sum_fit.rat["armstrat","PVal"]
      ratio_arm <- c(ratio_arm, ratio_mean_arm, ratio_lower_arm,ratio_upper_arm, ratio_pval_arm)
   # }
    # getting restricted mean absolute difference, B-A
   # N <- length(RM.times)
   # for (i in 1:N){
      fit.rm <- geese(pseudom ~ arm, data=temp, id=id, jack=TRUE, family=gaussian, corstr="independence", scale.fix=F)
      sum.fit.rm <- round(cbind(mean=fit.rm$beta, SE=sqrt(diag(fit.rm$vbeta.ajs)), Z=fit.rm$beta/sqrt(diag(fit.rm$vbeta.ajs)),
                                PVal = 2-2*pnorm(abs(fit.rm$beta/sqrt(diag(fit.rm$vbeta.ajs))))),4)
      # difference in restricted mean B - A
      RMST_arm <- sum.fit.rm["armstrat","mean"]
      RMST_lower_arm <- sum.fit.rm["armstrat","mean"] - (sum.fit.rm["armstrat","SE"]*qnorm(1-alpha))
      RMST_upper_arm <- sum.fit.rm["armstrat","mean"] + (sum.fit.rm["armstrat","SE"]*qnorm(1-alpha))
      RMST_pval_arm <- sum.fit.rm["armstrat","PVal"]
      RM_arm <- c(RM_arm, RMST_arm, RMST_lower_arm,RMST_upper_arm, RMST_pval_arm)
   # }
    # getting the coxph estimated hazard ratio
    fit.arm <- coxph(Surv(surv.time, event) ~ arm, data=temp, cluster=id)
    hr.est.arm <- exp(fit.arm$coefficients["armstrat"])
    hr.lower.arm <- exp(fit.arm$coefficients["armstrat"] - (summary(fit.arm)$coefficients[4]*qnorm(1-alpha)))
    hr.upper.arm <- exp(fit.arm$coefficients["armstrat"] + (summary(fit.arm)$coefficients[4]*qnorm(1-alpha)))
    hr.pval.arm <- summary(fit.arm)$waldtest[3]
    HR_arm <- c(hr.est.arm, hr.lower.arm, hr.upper.arm, hr.pval.arm)
    
    # combining results
    result <- c(nenroll, nevents, stdy_time, pval, chisq, upper.bound,lower.bound,
                exp_num, stop_eff, stop_fut,  medphys, medstrat, medstratpos, medstratneg, medphyspos, medphysneg,
                medstratB, medstratA, medphysB, medphysA, HR_arm,
                ratio_arm, RM_arm, pseudo_time, pseudot.ratio)
  }
  
    return(result)
}

# function that evaluates the rnichment design over desired number of repetitions
eval.scen <- function(estimand, # declaring the estimand of interest; options are:
                                      # subgrp - estimating the effect within a subgroup
                                      # clin - estimating clinical utility
                                      # inter - estimating a differential treatment effect between subgroups
                       reps, # of trials to simulate
                       max_enroll,  # maximum number enrolled over course of trial
                       n, # total number of events per trial
                       num_interim, # total number of analyses, including final
                       int_timing, # vector of proportion of events for each interim analysis timing
                       alpha=.025, # one-sided type 1 error rate (upper boundary)
                       low_err=0.1,  # lower (futility) boundary error
                       bound.type=c(1,1), # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively)
                      # time.points,    # value or vector of time points for which to compare survival probability at time points
                      # RM.times, # value or vector of time points to compare the restricted mean survival probability
                        
                       accru.rate,           # acrrual rate which follows a Poisson distribution with input rate 
                       loss.fu.rate,          # lost to follow-up rate; assumed equal in each strata
                       marker.pos,           # proportion of subjects who are marker positive
                       rand.arm,             # probability of being randomized to strategy arm
                       rand.pos,             # probability of being randomized to trt B in positives
                       rand.neg,             # probability of being randomized to trt B in negatives
                       rand,                  # probability of being randomized to trr B regardless of marker status
                       phys.choice.pos,        # proportion of strategy-based treatment assignments that agrees with physician's choice in positives
                                              # strategy is to treat the positives with trt B
                       phys.choice.neg,        # proportion of strategy-based treatment assignments that agrees with physician's choice in negatives
                       medposB,             # desired median for the survival distribution in the positive group, trtmnt B
                      
                       distrposB,           # distr defines the distribution of survival times; positives, trtmt B
                                                # exp = exponential survival times
                                                # weib = weibull survival times
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
  
  # creating the result data frame for the enirchment design results
  cols_enrich <- 17+(4*1)+(4*1)+2
  result_enrich <- data.frame(matrix(nrow=reps, ncol = cols_enrich))
  names1_enrich <- NULL
  names2_enrich <- NULL
  for (j in 1:1){
    for (i in 1:4){
      if (i %% 4==1){
        names1_enrich[(i+(4*(j-1)))] <- paste("ratio.mean")
    }
      else if (i %% 4==2){
        names1_enrich[(i+(4*(j-1)))] <- paste("ratio.lower")
    }
      else if (i %% 4==3){
        names1_enrich[(i+(4*(j-1)))] <- paste("ratio.upper")
      }
      else if (i %% 4==0){
        names1_enrich[(i+(4*(j-1)))] <- paste("ratio.pval")
      }
    }
  }
  for (j in 1:1){
    for (i in 1:4){
      if (i %% 4==1){
        names2_enrich[(i+(4*(j-1)))] <- paste("RMST.mean")
    }
      else if (i %% 4==2){
        names2_enrich[(i+(4*(j-1)))] <- paste("RMST.lower")
    }
      else if (i %% 4==3){
        names2_enrich[(i+(4*(j-1)))] <- paste("RMST.upper")
      }
      else if (i %% 4==0){
        names2_enrich[(i+(4*(j-1)))] <- paste("RMST.pval")
      }
    }
  }
  colnames(result_enrich) <- c("enroll", "events", "analysis.time",
                        "pvalue", "chisq", "upper.bound","lower.bound",
                        "exp.diff", "stop.eff", "stop.fut", "medianA", "medianB","medianpos",
                        "hr.mean", "hr.lower", "hr.upper","hr.pval" ,names1_enrich, names2_enrich, "pseudom.time", "pseudor.time")
  
  # creating the result data frame for the stratified design results
  cols_strat <- 40+(20*1)+(20*1)+2
  result_strat <- data.frame(matrix(nrow=reps, ncol = cols_strat))
  names1_strat <- NULL
  names2_strat <- NULL
  for (j in 1:1){
    for (i in 1:4){
      if (i %% 4==1){
        names1_strat[(i+(20*(j-1)))] <- paste("ratio.BvA.mean")
        names1_strat[(4+i+(20*(j-1)))] <- paste("ratio.pos.mean")
        names1_strat[(8+i+(20*(j-1)))] <- paste("ratio.neg.mean")
        names1_strat[(12+i+(20*(j-1)))] <- paste("ratio.int.mean")
        names1_strat[(16+i+(20*(j-1)))] <- paste("ratio.arm.mean")
      }
      else if (i %% 4==2){
        names1_strat[(i+(20*(j-1)))] <- paste("ratio.BvA.lower")
        names1_strat[(4+i+(20*(j-1)))] <- paste("ratio.pos.lower")
        names1_strat[(8+i+(20*(j-1)))] <- paste("ratio.neg.lower")#,time.points[j])
        names1_strat[(12+i+(20*(j-1)))] <- paste("ratio.int.lower")#,time.points[j])
        names1_strat[(16+i+(20*(j-1)))] <- paste("ratio.arm.lower")#,time.points[j])
      }
      else if (i %% 4==3){
        names1_strat[(i+(20*(j-1)))] <- paste("ratio.BvA.upper")#,time.points[j])
        names1_strat[(4+i+(20*(j-1)))] <- paste("ratio.pos.upper")#,time.points[j])
        names1_strat[(8+i+(20*(j-1)))] <- paste("ratio.neg.upper")#,time.points[j])
        names1_strat[(12+i+(20*(j-1)))] <- paste("ratio.int.upper")#,time.points[j])
        names1_strat[(16+i+(20*(j-1)))] <- paste("ratio.arm.upper")#,time.points[j])
      }
      else if (i %% 4==0){
        names1_strat[(i+(20*(j-1)))] <- paste("ratio.BvA.pval")#,time.points[j])
        names1_strat[(4+i+(20*(j-1)))] <- paste("ratio.pos.pval")#,time.points[j])
        names1_strat[(8+i+(20*(j-1)))] <- paste("ratio.neg.pval")#,time.points[j])
        names1_strat[(12+i+(20*(j-1)))] <- paste("ratio.int.pval")#,time.points[j])
        names1_strat[(16+i+(20*(j-1)))] <- paste("ratio.arm.pval")#,time.points[j])
      }
    }
  }
  for (j in 1:1){#length(RM.times)){
    for (i in 1:4){
      if (i %% 4==1){
        names2_strat[(i+(20*(j-1)))] <- paste("RMST.BvA.mean")#,RM.times[j])
        names2_strat[(4+i+(20*(j-1)))] <- paste("RMST.pos.mean")#,RM.times[j])
        names2_strat[(8+i+(20*(j-1)))] <- paste("RMST.neg.mean")#,RM.times[j])
        names2_strat[(12+i+(20*(j-1)))] <- paste("RMST.int.mean")#,RM.times[j])
        names2_strat[(16+i+(20*(j-1)))] <- paste("RMST.arm.mean")#,RM.times[j])
      }
      else if (i %% 4==2){
        names2_strat[(i+(20*(j-1)))] <- paste("RMST.BvA.lower")#,RM.times[j])
        names2_strat[(4+i+(20*(j-1)))] <- paste("RMST.pos.lower")#,RM.times[j])
        names2_strat[(8+i+(20*(j-1)))] <- paste("RMST.neg.lower")#,RM.times[j])
        names2_strat[(12+i+(20*(j-1)))] <- paste("RMST.int.lower")#,RM.times[j])
        names2_strat[(16+i+(20*(j-1)))] <- paste("RMST.arm.lower")#,RM.times[j])
      }
      else if (i %% 4==3){
        names2_strat[(i+(20*(j-1)))] <- paste("RMST.BvA.upper")#,RM.times[j])
        names2_strat[(4+i+(20*(j-1)))] <- paste("RMST.pos.upper")#,RM.times[j])
        names2_strat[(8+i+(20*(j-1)))] <- paste("RMST.neg.upper")#,RM.times[j])
        names2_strat[(12+i+(20*(j-1)))] <- paste("RMST.int.upper")#,RM.times[j])
        names2_strat[(16+i+(20*(j-1)))] <- paste("RMST.arm.upper")#,RM.times[j])
      }
      else if (i %% 4==0){
        names2_strat[(i+(20*(j-1)))] <- paste("RMST.BvA.pval")#,RM.times[j])
        names2_strat[(4+i+(20*(j-1)))] <- paste("RMST.pos.pval")#,RM.times[j])
        names2_strat[(8+i+(20*(j-1)))] <- paste("RMST.neg.pval")#,RM.times[j])
        names2_strat[(12+i+(20*(j-1)))] <- paste("RMST.int.pval")#,RM.times[j])
        names2_strat[(16+i+(20*(j-1)))] <- paste("RMST.arm.pval")#,RM.times[j])
      }
    }
  }
  colnames(result_strat) <- c("enroll", "events", "analysis.time",
                               "pvalue", "chisq", "upper.bound","lower.bound",
                               "exp.diff", "stop.eff", "stop.fut", "medianA", "medianB","medianAneg","medianApos",
                               "medianBneg", "medianBpos", "medianneg", "medianpos", "medianphys", "medianstrat",
                               "hr.BvA.mean","hr.BvA.lower","hr.BvA.upper","hr.BvA.pval",
                               "hr.pos.mean","hr.pos.lower","hr.pos.upper","hr.pos.pval",
                               "hr.neg.mean","hr.neg.lower","hr.neg.upper","hr.neg.pval",
                               "hr.int.mean","hr.int.lower","hr.int.upper","hr.int.pval", 
                               "hr.arm.mean","hr.arm.lower","hr.arm.upper","hr.arm.pval", names1_strat, names2_strat, 
                              "pseudom.time", "pseudor.time")
  
  
  # creating the result data frame for the enirchment design results; same as the stratified results
  result_mod <- result_strat
  
  
  # creating the result data frame for the strategy design results
  cols_clin <- 24+(4*1)+(4*1)+2
  result_clin <- data.frame(matrix(nrow=reps, ncol = cols_clin))
  names1_clin <- NULL
  names2_clin <- NULL
  for (j in 1:1){
    for (i in 1:4){
      if (i %% 4==1){
        names1_clin[(i+(4*(j-1)))] <- paste("ratio.arm.mean")#,time.points[j])
      }
      else if (i %% 4==2){
        names1_clin[(i+(4*(j-1)))] <- paste("ratio.arm.lower")#,time.points[j])
      }
      else if (i %% 4==3){
        names1_clin[(i+(4*(j-1)))] <- paste("ratio.arm.upper")#,time.points[j])
      }
      else if (i %% 4==0){
        names1_clin[(i+(4*(j-1)))] <- paste("ratio.arm.pval")#,time.points[j])
      }
    }
  }
  for (j in 1:1){
    for (i in 1:4){
      if (i %% 4==1){
        names2_clin[(i+(4*(j-1)))] <- paste("RMST.arm.mean")#,RM.times[j])
      }
      else if (i %% 4==2){
        names2_clin[(i+(4*(j-1)))] <- paste("RMST.arm.lower")#,RM.times[j])
      }
      else if (i %% 4==3){
        names2_clin[(i+(4*(j-1)))] <- paste("RMST.arm.upper")#,RM.times[j])
      }
      else if (i %% 4==0){
        names2_clin[(i+(4*(j-1)))] <- paste("RMST.arm.pval")#,RM.times[j])
      }
    }
  }
  colnames(result_clin) <- c("enroll", "events", "analysis.time",
                               "pvalue", "chisq", "upper.bound","lower.bound",
                               "exp.diff", "stop.eff", "stop.fut", "medianphys", "medianstrat", "medianstratpos",
                               "medianstratneg", "medianphyspos", "medianphysneg",
                               "medianstratB", "medianstratA", "medianphysB","medianphysA",
                               "hr.arm.mean", "hr.arm.lower", "hr.arm.upper","hr.arm.pval" ,names1_clin, names2_clin, 
                             "pseudom.time", "pseudor.time")
  
  if (estimand=="subgrp"){
    min.surv <- min(medposB, medposA)
    pseudot <- (-log(0.5))/(log(2)/min.surv)
  
    for (i in 1:reps){
      # simulate data
      x <- simdat(n=max_enroll,                    # total sample size
                   accru.rate=accru.rate,           # acrrual rate which follows a Poisson distribution 
                   loss.fu.rate=loss.fu.rate,          # lost to follow-up rate; assumed equal in each strata
                   marker.pos=marker.pos,           # proportion of subjects who are marker positive
                   rand.arm=rand.arm,              # probability of being randomized to strategy arm                         
                   rand.pos=rand.pos,             # probability of being randomized to trt B in positives
                   rand.neg=rand.neg,             # probability of being randomized to trt B in negatives
                   rand=rand,                  # probability of being randomized to trr B regardless of marker status
                   phys.choice.pos=phys.choice.pos,        # proportion of strategy-based treatment assignments that agrees with physician's choice in positives
                   # strategy is to treat the positives with trt B
                   phys.choice.neg=phys.choice.neg,        # proportion of strategy-based treatment assignments that agrees with physician's choice in negatives
                   # strategy is to treat the negtives with trt A
                   medposB=medposB,             # desired median for the survival distribution in the positive group, trtmnt B
                   distrposB=distrposB,           # distr defines the distribution of survival times; positives, trtmt B
                   # exp = exponential survival times
                   # weib = weibull survival times
                   # nonpar = non-parametric survival times created using NOT CODED CURRENTY
                   # the coxed package and sim.survdata function
                   shapeposB=shapeposB,           # shape parmeter for the weibull distribution; positives, trtmt B
                   
                   medposA=medposA,             # desired median for the survival distribution in the positive group, trtmnt A
                   distrposA=distrposA,           # distr defines the distribution of survival times; positives, trtmt A
                   shapeposA=shapeposA,           # shape parmeter for the weibull distribution; positives, trtmt A
                   
                   mednegB=mednegB,             # desired median for the survival distribution in the negative group, trtmnt B
                   distrnegB=distrnegB,           # distr defines the distribution of survival times; negatives, trtmt B
                   shapenegB=shapenegB,           # shape parmeter for the weibull distribution; negatives, trtmt B
                   
                   mednegA=mednegA,             # desired median for the survival distribution in the negative group, trtmnt A
                   distrnegA=distrnegA,           # distr defines the distribution of survival times; negatives, trtmt A
                   shapenegA=shapenegA           # shape parmeter for the weibull distribution; negatives, trtmt A)
      )
      
      
      # simulate the trial results for the enrichment design
      result_enrich[i,] <- sim.trial(data=x, # input survival data from our data generating function
                                     n=n, # total number of events
                                     num_interim=num_interim, # total number of analyses, including final
                                     int_timing=int_timing, # vector of proportion of events for each interim analysis timing
                                     alpha=alpha, # one-sided type 1 error rate (upper boundary)
                                     low_err=low_err,  # lower (futility) boundary error
                                     bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
                                    # time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
                                    # RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
                                     design.type="enrich",    # type of design being analyzed choices are:
                                     # enrich
                                     # stratify
                                     # strategy
                                     # modstrat
                                     test.grp="pos",       # the group comparison that the logrank test is based on, choices are
                                     # pos - trt B vs A in the positives
                                     # neg - trt B vs A in the negatives
                                     # clin - clinical utility; strategy arm vs phys choice arm
                                     # trtB - pos vs neg for trtB patients
                                     # trtA - pos vs neg for trtA patients
                                     # inter - 4 way comparison for interaction; k=4 logrank test
                                    pseudot.ratio=pseudot
      )
      
      # simulate the trial results for the stratified design
      result_strat[i,] <- sim.trial(data=x, # input survival data from our data generating function
                                     n=n, # total number of events
                                     num_interim=num_interim, # total number of analyses, including final
                                     int_timing=int_timing, # vector of proportion of events for each interim analysis timing
                                     alpha=alpha, # one-sided type 1 error rate (upper boundary)
                                     low_err=low_err,  # lower (futility) boundary error
                                     bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
                                    # time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
                                    # RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
                                     design.type="stratify",    # type of design being analyzed choices are:
                                     # enrich
                                     # stratify
                                     # strategy
                                     # modstrat
                                     test.grp="pos",       # the group comparison that the logrank test is based on, choices are
                                     # pos - trt B vs A in the positives
                                     # neg - trt B vs A in the negatives
                                     # clin - clinical utility; strategy arm vs phys choice arm
                                     # trtB - pos vs neg for trtB patients
                                     # trtA - pos vs neg for trtA patients
                                     # inter - 4 way comparison for interaction; k=4 logrank test
                                    pseudot.ratio=pseudot
      )
      
      # # simulate the trial results for the modified strategy design
      # result_mod[i,] <- sim.trial(data=x, # input survival data from our data generating function
      #                               n=n, # total number of events
      #                               num_interim=num_interim, # total number of analyses, including final
      #                               int_timing=int_timing, # vector of proportion of events for each interim analysis timing
      #                               alpha=alpha, # one-sided type 1 error rate (upper boundary)
      #                               low_err=low_err,  # lower (futility) boundary error
      #                               bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
      #                              # time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
      #                             #  RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
      #                               design.type="modstrat",    # type of design being analyzed choices are:
      #                               # enrich
      #                               # stratify
      #                               # strategy
      #                               # modstrat
      #                               test.grp="pos",       # the group comparison that the logrank test is based on, choices are
      #                               # pos - trt B vs A in the positives
      #                               # neg - trt B vs A in the negatives
      #                               # clin - clinical utility; strategy arm vs phys choice arm
      #                               # trtB - pos vs neg for trtB patients
      #                               # trtA - pos vs neg for trtA patients
      #                               # inter - 4 way comparison for interaction; k=4 logrank test
      #                               pseudot.ratio=pseudot
      # )
    }
    resultList <- list("subgrp_enrich"=result_enrich, "subgrp_stratify"=result_strat
                       #, "subgrp_modstrat"=result_mod
                       )
  }
  
  else if (estimand=="clin"){
    min.surv <- min(medposB, medposA, mednegB, mednegA)
    pseudot <- (-log(0.5))/(log(2)/min.surv)
    for (i in 1:reps){
      # simulate data
      x <- simdat(n=max_enroll,                    # total sample size
                accru.rate=accru.rate,           # acrrual rate which follows a Poisson distribution 
                loss.fu.rate=loss.fu.rate,          # lost to follow-up rate; assumed equal in each strata
                marker.pos=marker.pos,           # proportion of subjects who are marker positive
                rand.arm=rand.arm,              # probability of being randomized to strategy arm                         
                rand.pos=rand.pos,             # probability of being randomized to trt B in positives
                rand.neg=rand.neg,             # probability of being randomized to trt B in negatives
                rand=rand,                  # probability of being randomized to trr B regardless of marker status
                phys.choice.pos=phys.choice.pos,        # proportion of strategy-based treatment assignments that agrees with physician's choice in positives
                # strategy is to treat the positives with trt B
                phys.choice.neg=phys.choice.neg,        # proportion of strategy-based treatment assignments that agrees with physician's choice in negatives
                # strategy is to treat the negtives with trt A
                medposB=medposB,             # desired median for the survival distribution in the positive group, trtmnt B
                distrposB=distrposB,           # distr defines the distribution of survival times; positives, trtmt B
                # exp = exponential survival times
                # weib = weibull survival times
                # nonpar = non-parametric survival times created using NOT CODED CURRENTY
                # the coxed package and sim.survdata function
                shapeposB=shapeposB,           # shape parmeter for the weibull distribution; positives, trtmt B
                
                medposA=medposA,             # desired median for the survival distribution in the positive group, trtmnt A
                distrposA=distrposA,           # distr defines the distribution of survival times; positives, trtmt A
                shapeposA=shapeposA,           # shape parmeter for the weibull distribution; positives, trtmt A
                
                mednegB=mednegB,             # desired median for the survival distribution in the negative group, trtmnt B
                distrnegB=distrnegB,           # distr defines the distribution of survival times; negatives, trtmt B
                shapenegB=shapenegB,           # shape parmeter for the weibull distribution; negatives, trtmt B
                
                mednegA=mednegA,             # desired median for the survival distribution in the negative group, trtmnt A
                distrnegA=distrnegA,           # distr defines the distribution of survival times; negatives, trtmt A
                shapenegA=shapenegA           # shape parmeter for the weibull distribution; negatives, trtmt A)
    )
    
    # simulate the trial results for the enrichment design
    result_clin[i,] <- sim.trial(data=x, # input survival data from our data generating function
                                   n=n, # total number of events
                                   num_interim=num_interim, # total number of analyses, including final
                                   int_timing=int_timing, # vector of proportion of events for each interim analysis timing
                                   alpha=alpha, # one-sided type 1 error rate (upper boundary)
                                   low_err=low_err,  # lower (futility) boundary error
                                   bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
                                  # time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
                                  # RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
                                   design.type="strategy",    # type of design being analyzed choices are:
                                   # enrich
                                   # stratify
                                   # strategy
                                   # modstrat
                                   test.grp="clin",       # the group comparison that the logrank test is based on, choices are
                                   # pos - trt B vs A in the positives
                                   # neg - trt B vs A in the negatives
                                   # clin - clinical utility; strategy arm vs phys choice arm
                                   # trtB - pos vs neg for trtB patients
                                   # trtA - pos vs neg for trtA patients
                                   # inter - 4 way comparison for interaction; k=4 logrank test
                                 pseudot.ratio=pseudot
    )
    
    # simulate the trial results for the stratified design
    # result_strat[i,] <- sim.trial(data=x, # input survival data from our data generating function
    #                               n=n, # total number of events
    #                               num_interim=num_interim, # total number of analyses, including final
    #                               int_timing=int_timing, # vector of proportion of events for each interim analysis timing
    #                               alpha=alpha, # one-sided type 1 error rate (upper boundary)
    #                               low_err=low_err,  # lower (futility) boundary error
    #                               bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
    #                             #  time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
    #                              # RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
    #                               design.type="stratify",    # type of design being analyzed choices are:
    #                               # enrich
    #                               # stratify
    #                               # strategy
    #                               # modstrat
    #                               test.grp="clin",       # the group comparison that the logrank test is based on, choices are
    #                               # pos - trt B vs A in the positives
    #                               # neg - trt B vs A in the negatives
    #                               # clin - clinical utility; strategy arm vs phys choice arm
    #                               # trtB - pos vs neg for trtB patients
    #                               # trtA - pos vs neg for trtA patients
    #                               # inter - 4 way comparison for interaction; k=4 logrank test
    #                                pseudot.ratio=pseudot
    # )
    
    # simulate the trial results for the modified strategy design
    # result_mod[i,] <- sim.trial(data=x, # input survival data from our data generating function
    #                             n=n, # total number of events
    #                             num_interim=num_interim, # total number of analyses, including final
    #                             int_timing=int_timing, # vector of proportion of events for each interim analysis timing
    #                             alpha=alpha, # one-sided type 1 error rate (upper boundary)
    #                             low_err=low_err,  # lower (futility) boundary error
    #                             bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
    #                           #  time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
    #                           #  RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
    #                             design.type="modstrat",    # type of design being analyzed choices are:
    #                             # enrich
    #                             # stratify
    #                             # strategy
    #                             # modstrat
    #                             test.grp="clin",       # the group comparison that the logrank test is based on, choices are
    #                             # pos - trt B vs A in the positives
    #                             # neg - trt B vs A in the negatives
    #                             # clin - clinical utility; strategy arm vs phys choice arm
    #                             # trtB - pos vs neg for trtB patients
    #                             # trtA - pos vs neg for trtA patients
    #                             # inter - 4 way comparison for interaction; k=4 logrank test
    #                               pseudot.ratio=pseudot
    # )
    }
    resultList <- list("clin_strategy"=result_clin
                       #,"clin_stratify"=result_strat, "clin_modstrat"=result_mod
                       )
  }
  
  else if (estimand=="inter"){
    min.surv <- min(medposB, medposA, mednegB, mednegA)
    pseudot <- (-log(0.5))/(log(2)/min.surv)
    for (i in 1:reps){
      # simulate data
      x <- simdat(n=max_enroll,                    # total sample size
                  accru.rate=accru.rate,           # acrrual rate which follows a Poisson distribution 
                  loss.fu.rate=loss.fu.rate,          # lost to follow-up rate; assumed equal in each strata
                  marker.pos=marker.pos,           # proportion of subjects who are marker positive
                  rand.arm=rand.arm,              # probability of being randomized to strategy arm                         
                  rand.pos=rand.pos,             # probability of being randomized to trt B in positives
                  rand.neg=rand.neg,             # probability of being randomized to trt B in negatives
                  rand=rand,                  # probability of being randomized to trr B regardless of marker status
                  phys.choice.pos=phys.choice.pos,        # proportion of strategy-based treatment assignments that agrees with physician's choice in positives
                  # strategy is to treat the positives with trt B
                  phys.choice.neg=phys.choice.neg,        # proportion of strategy-based treatment assignments that agrees with physician's choice in negatives
                  # strategy is to treat the negtives with trt A
                  medposB=medposB,             # desired median for the survival distribution in the positive group, trtmnt B
                  distrposB=distrposB,           # distr defines the distribution of survival times; positives, trtmt B
                  # exp = exponential survival times
                  # weib = weibull survival times
                  # nonpar = non-parametric survival times created using NOT CODED CURRENTY
                  # the coxed package and sim.survdata function
                  shapeposB=shapeposB,           # shape parmeter for the weibull distribution; positives, trtmt B
                  
                  medposA=medposA,             # desired median for the survival distribution in the positive group, trtmnt A
                  distrposA=distrposA,           # distr defines the distribution of survival times; positives, trtmt A
                  shapeposA=shapeposA,           # shape parmeter for the weibull distribution; positives, trtmt A
                  
                  mednegB=mednegB,             # desired median for the survival distribution in the negative group, trtmnt B
                  distrnegB=distrnegB,           # distr defines the distribution of survival times; negatives, trtmt B
                  shapenegB=shapenegB,           # shape parmeter for the weibull distribution; negatives, trtmt B
                  
                  mednegA=mednegA,             # desired median for the survival distribution in the negative group, trtmnt A
                  distrnegA=distrnegA,           # distr defines the distribution of survival times; negatives, trtmt A
                  shapenegA=shapenegA           # shape parmeter for the weibull distribution; negatives, trtmt A)
      )
      
      
      
      # simulate the trial results for the stratified design
      result_strat[i,] <- sim.trial(data=x, # input survival data from our data generating function
                                    n=n, # total number of events
                                    num_interim=num_interim, # total number of analyses, including final
                                    int_timing=int_timing, # vector of proportion of events for each interim analysis timing
                                    alpha=alpha, # one-sided type 1 error rate (upper boundary)
                                    low_err=low_err,  # lower (futility) boundary error
                                    bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
                                  #  time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
                                   # RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
                                    design.type="stratify",    # type of design being analyzed choices are:
                                    # enrich
                                    # stratify
                                    # strategy
                                    # modstrat
                                    test.grp="inter",       # the group comparison that the logrank test is based on, choices are
                                    # pos - trt B vs A in the positives
                                    # neg - trt B vs A in the negatives
                                    # clin - clinical utility; strategy arm vs phys choice arm
                                    # trtB - pos vs neg for trtB patients
                                    # trtA - pos vs neg for trtA patients
                                    # inter - 4 way comparison for interaction; k=4 logrank test
                                  pseudot.ratio=pseudot
      )
      
      # # simulate the trial results for the modified strategy design
      # result_mod[i,] <- sim.trial(data=x, # input survival data from our data generating function
      #                             n=n, # total number of events
      #                             num_interim=num_interim, # total number of analyses, including final
      #                             int_timing=int_timing, # vector of proportion of events for each interim analysis timing
      #                             alpha=alpha, # one-sided type 1 error rate (upper boundary)
      #                             low_err=low_err,  # lower (futility) boundary error
      #                             bound.type=bound.type,  # boundary type 1= OBF, 2=Pockock for lower and upper boundary, respectively
      #                            # time.points=time.points,    # value or vector of time points for which to compare survival probability at time points
      #                           #  RM.times=RM.times, # value or vector of time points to compare the restricted mean survival probability
      #                             design.type="modstrat",    # type of design being analyzed choices are:
      #                             # enrich
      #                             # stratify
      #                             # strategy
      #                             # modstrat
      #                             test.grp="inter",       # the group comparison that the logrank test is based on, choices are
      #                             # pos - trt B vs A in the positives
      #                             # neg - trt B vs A in the negatives
      #                             # clin - clinical utility; strategy arm vs phys choice arm
      #                             # trtB - pos vs neg for trtB patients
      #                             # trtA - pos vs neg for trtA patients
      #                             # inter - 4 way comparison for interaction; k=4 logrank test
      #                             pseudot.ratio=pseudot
      # )
    }
    resultList <- list("inter_stratify"=result_strat
                       #, "inter_modstrat"=result_mod
                       )
  }
  
  return(resultList)
}


