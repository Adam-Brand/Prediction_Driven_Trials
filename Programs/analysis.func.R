#==============================================================================
# FILENAME: analysis.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: functions toread in the result data frames, summarizes results and puts into tables and figures
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 3.6.1
#==============================================================================
#Notes: 
### 




# =============================================================================
## sourcing the source program
source("Programs/sim_source.R")


## function to summarize one results data frame
analyze.exp <- function(data,              # result data frame 
                    estimand,          # subgrp, clin, inter
                    design.type,        # trial design type, i.e, enrich, stratify, strategy, modstrat
                    medianBpos,
                    medianBneg,
                    medianApos,
                    medianAneg,
                    marker.pos,
                    phys.choice,
                    type1 = FALSE
){
  
if (estimand=="subgrp"){  
  if (design.type=="enrich"){
    
    if (type1 == FALSE){
  # getting the power of each testing method
      power.lr <- round(sum(data$stop.eff==1)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.pval<0.05 & data$hr.mean < 1)/length(data$hr.pval), digit=3)
      power.rat <- round(sum(data[,21]<0.05 & data[,18]>1)/length(data[,21]), digits=3)
      power.rmst <- round(sum(data[,25]<0.05 & data[,22]>0)/length(data[,25]), digits=3)
    }
    
    else if (type1==TRUE){
      power.lr <- round(sum(data$pvalue<.05)/length(data$pvalue), digits=3)
      power.hr <- round(sum(data$hr.pval<0.05)/length(data$hr.pval), digit=3)
      power.rat <- round(sum(data[,21]<0.05)/length(data[,21]), digits=3)
      power.rmst <- round(sum(data[,25]<0.05)/length(data[,25]), digits=3)
    }
  # getting the estimated treatment difference using each method
  hr.mean <- mean(data$hr.mean)
  hr.mean.SE <- sd(data$hr.mean)*qnorm(.975)
  
  rat.mean <- mean(data[,18])
  rat.mean.SE <- sd(data[,18])*qnorm(.975)
  
  rmst.mean <- mean(data[,22])
  rmst.mean.SE <- sd(data[,22])*qnorm(.975)
  
  # calculating the true contrasts using exponential distributions and given median survivals
  true.hr <- medianApos/medianBpos
  true.rat <- exp(-(log(2)/medianBpos)*data$pseudor.time[1])/exp(-(log(2)/medianApos)*data$pseudor.time[1])
  integrand1 <- function(x){exp(-(x*(log(2)/medianBpos)))}
  integrand2 <- function(x){exp(-(x*(log(2)/medianApos)))}
  int1 <- NULL
  int2 <- NULL
  true.rmst <- NULL
  cover.rmst <- NULL
  cover.hr <- NULL
  cover.rat <- NULL
  for (i in 1:length(data$pseudom.time)){
    int1[i] <- integrate(integrand1, 0,data$pseudom.time[i])$value
    int2[i] <- integrate(integrand2, 0,data$pseudom.time[i])$value
    true.rmst[i] <- int1[i]-int2[i]
    if ((data$RMST.lower[i]  <= true.rmst[i]) & (data$RMST.upper[i]  >= true.rmst[i])){cover.rmst[i]=1}
    else {cover.rmst[i]=0}
    if ((data$hr.lower[i]  <= true.hr) & (data$hr.upper[i]  >= true.hr)){cover.hr[i]=1}
    else {cover.hr[i]=0}
    if (is.na(data$ratio.lower[i]) | is.na(data$ratio.upper[i])){cover.rat[i]=0}
    else if ((data$ratio.lower[i]  <= true.rat) & (data$ratio.upper[i]  >= true.rat)){cover.rat[i]=1}
    else {cover.rat[i]=0}
  }
  hr.cover <- mean(cover.hr)
  RMST.cover <- mean(cover.rmst)
  rat.cover <- mean(cover.rat)
  
  }
  
  
  else if (design.type=="stratify" | design.type=="modstrat"){
    # getting the power of each testing method for detecting B superior to A in positives
    
    if (type1 == FALSE){
      # getting the power of each testing method
      power.lr <- round(sum(data$stop.eff==1 & data$exp.diff>0)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.pos.pval<0.05 & data$hr.pos.mean < 1)/length(data$hr.pos.pval), digit=3)
      power.rat <- round(sum(data[,48]<0.05 & data[,45]>1)/length(data[,45]), digits=3)
      power.rmst <- round(sum(data[,68]<0.05 & data[,65]>0)/length(data[,65]), digits=3)
    }
    
    else if (type1==TRUE){
      power.lr <- round(sum(data$pvalue<.05)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.pos.pval<0.05)/length(data$hr.pos.pval), digit=3)
      power.rat <- round(sum(data[,48]<0.05)/length(data[,45]), digits=3)
      power.rmst <- round(sum(data[,68]<0.05)/length(data[,65]), digits=3)
    }
    
    
    # getting the estimated treatment difference using each method
    hr.mean <- mean(data$hr.pos.mean)
    hr.mean.SE <- sd(data$hr.pos.mean)*qnorm(.975)
    
    rat.mean <- mean(data[,45])
    rat.mean.SE <- sd(data[,45])*qnorm(.975)
    
    rmst.mean <- mean(data[,65])
    rmst.mean.SE <- sd(data[,65])*qnorm(.975)
    
    # calculating the true contrasts using exponential distributions and given median survivals
    true.hr <- medianApos/medianBpos
    true.rat <- exp(-(log(2)/medianBpos)*data$pseudor.time[1])/exp(-(log(2)/medianApos)*data$pseudor.time[1])
    integrand1 <- function(x){exp(-(x*(log(2)/medianBpos)))}
    integrand2 <- function(x){exp(-(x*(log(2)/medianApos)))}
    int1 <- NULL
    int2 <- NULL
    true.rmst <- NULL
    cover.rmst <- NULL
    cover.hr <- NULL
    cover.rat <- NULL
    for (i in 1:length(data$pseudom.time)){
      int1[i] <- integrate(integrand1, 0,data$pseudom.time[i])$value
      int2[i] <- integrate(integrand2, 0,data$pseudom.time[i])$value
      true.rmst[i] <- int1[i]-int2[i]
      if ((data$RMST.pos.lower[i]  <= true.rmst[i]) & (data$RMST.pos.upper[i]  >= true.rmst[i])){cover.rmst[i]=1}
      else {cover.rmst[i]=0}
      if ((data$hr.pos.lower[i]  <= true.hr) & (data$hr.pos.upper[i]  >= true.hr)){cover.hr[i]=1}
      else {cover.hr[i]=0}
      if (is.na(data$ratio.pos.lower[i]) | is.na(data$ratio.pos.upper[i])){cover.rat[i]=0}
      else if ((data$ratio.pos.lower[i]  <= true.rat) & (data$ratio.pos.upper[i]  >= true.rat)){cover.rat[i]=1}
      else {cover.rat[i]=0}
    }
    hr.cover <- mean(cover.hr)
    RMST.cover <- mean(cover.rmst)
    rat.cover <- mean(cover.rat)
  }
  }
  
else if (estimand=="clin"){
  if (design.type=="stratify" | design.type=="modstrat"){
    # getting the power of each testing method
    
    if (type1 == FALSE){
      # getting the power of each testing method
      power.lr <- round(sum(data$stop.eff==1 & data$exp.diff>0)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.arm.pval<0.05 & data$hr.arm.mean < 1)/length(data$hr.arm.pval), digit=3)
      power.rat <- round(sum(data[,60]<0.05 & data[,57]>1)/length(data[,57]), digits=3)
      power.rmst <- round(sum(data[,80]<0.05 & data[,77]>0)/length(data[,77]), digits=3)
    }
    
    else if (type1==TRUE){
      power.lr <- round(sum(data$pvalue<.05)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.arm.pval<0.05)/length(data$hr.arm.pval), digit=3)
      power.rat <- round(sum(data[,60]<0.05)/length(data[,57]), digits=3)
      power.rmst <- round(sum(data[,80]<0.05)/length(data[,77]), digits=3)
    }
    
    
    # getting the estimated treatment difference using each method
    hr.mean <- mean(data$hr.arm.mean)
    hr.mean.SE <- sd(data$hr.arm.mean)*qnorm(.975)
    
    rat.mean <- mean(data[,57])
    rat.mean.SE <- sd(data[,57])*qnorm(.975)
    
    rmst.mean <- mean(data[,77])
    rmst.mean.SE <- sd(data[,77])*qnorm(.975)
    
  }
  else if (design.type=="strategy"){
    # getting the power of each testing method
    if (type1 == FALSE){
      # getting the power of each testing method
      power.lr <- round(sum(data$stop.eff==1 & data$exp.diff>0)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.arm.pval<0.05 & data$hr.arm.mean < 1)/length(data$hr.arm.pval), digit=3)
      power.rat <- round(sum(data[,28]<0.05 & data[,25]>1)/length(data[,25]), digits=3)
      power.rmst <- round(sum(data[,32]<0.05 & data[,29]>0)/length(data[,29]), digits=3)
    }
    
    else if (type1==TRUE){
      power.lr <- round(sum(data$pvalue<.05)/length(data$stop.eff), digits=3)
      power.hr <- round(sum(data$hr.arm.pval<0.05)/length(data$hr.arm.pval), digit=3)
      power.rat <- round(sum(data[,28]<0.05)/length(data[,25]), digits=3)
      power.rmst <- round(sum(data[,32]<0.05)/length(data[,29]), digits=3)
    }
    
    
    # getting the estimated treatment difference using each method
    hr.mean <- mean(data$hr.arm.mean)
    hr.mean.SE <- sd(data$hr.arm.mean)*qnorm(.975)
    
    rat.mean <- mean(data[,25])
    rat.mean.SE <- sd(data[,25])*qnorm(.975)
    
    rmst.mean <- mean(data[,29])
    rmst.mean.SE <- sd(data[,29])*qnorm(.975)
    
    # calculating the true contrasts using exponential distributions and given median survivals
    
    #probabilities of marker positive and negative patients
    p_pos <- marker.pos
    p_neg <- 1-marker.pos
    # rate parameters for each of the 4 groups
    lam.pos.b <- log(2)/medianBpos
    lam.neg.b <- log(2)/medianBneg
    lam.pos.a <- log(2)/medianApos
    lam.neg.a <- log(2)/medianAneg
    
    # computing the hazard and integrating over x; the proportional hazards assumption does not hold due to different survivals within arms
    haz.integrand.clin <- function(x){((((p_pos)*lam.pos.b*exp(-lam.pos.b*x)) + ((p_neg)*lam.neg.a*exp(-lam.neg.a*x)))/(((p_pos)*exp(-lam.pos.b*x)) + ((p_neg)*exp(-lam.neg.a*x))))}
    haz.integrand.phys <- function(x){((((p_pos*phys.choice)*lam.pos.b*exp(-lam.pos.b*x)) + ((p_pos*(1-phys.choice))*lam.pos.a*exp(-lam.pos.a*x)) +
                    ((p_neg*phys.choice)*lam.neg.a*exp(-lam.neg.a*x)) + ((p_neg*(1-phys.choice))*lam.neg.b*exp(-lam.neg.b*x)))/
                (((p_pos*phys.choice)*exp(-lam.pos.b*x)) + ((p_pos*(1-phys.choice))*exp(-lam.pos.a*x)) +
                    ((p_neg*phys.choice)*exp(-lam.neg.a*x)) + ((p_neg*(1-phys.choice))*exp(-lam.neg.b*x))))}
    
    
    # computing the true ratio probabilities
    true.rat <- (((p_pos)*exp(-lam.pos.b*data$pseudor.time[1])) + ((p_neg)*exp(-lam.neg.a*data$pseudor.time[1])))/
                (((p_pos*phys.choice)*exp(-lam.pos.b*data$pseudor.time[1])) + ((p_pos*(1-phys.choice))*exp(-lam.pos.a*data$pseudor.time[1])) +
                ((p_neg*phys.choice)*exp(-lam.neg.a*data$pseudor.time[1])) + ((p_neg*(1-phys.choice))*exp(-lam.neg.b*data$pseudor.time[1])))
    
    #computing the true RMST
    int.pos.B <- function(x){exp(-(x*(log(2)/medianBpos)))}
    int.pos.A <- function(x){exp(-(x*(log(2)/medianApos)))}
    int.neg.B <- function(x){exp(-(x*(log(2)/medianBneg)))}
    int.neg.A <- function(x){exp(-(x*(log(2)/medianAneg)))}
    
    intposB <- NULL
    intposA <- NULL
    intnegB <- NULL
    intnegA <- NULL
    true.rmst <- NULL
    cover.rmst <- NULL
    cover.hr <- NULL
    cover.rat <- NULL
    haz.clin <- NULL
    haz.phys <- NULL
    true.hr <- NULL
    
    for (i in 1:length(data$pseudom.time)){
      intposB[i] <- integrate(int.pos.B, 0,data$pseudom.time[i])$value
      intposA[i] <- integrate(int.pos.A, 0,data$pseudom.time[i])$value
      intnegB[i] <- integrate(int.neg.B, 0,data$pseudom.time[i])$value
      intnegA[i] <- integrate(int.neg.A, 0,data$pseudom.time[i])$value
      true.rmst[i] <- ((p_pos*intposB[i]) + (p_neg*intnegA[i])) - ((p_pos*phys.choice*intposB[i]) + (p_pos*(1-phys.choice)*intposA[i]) +
                      (p_neg*phys.choice*intnegA[i]) + (p_neg*(1-phys.choice)*intnegB[i]))
      haz.clin[i] <- integrate(haz.integrand.clin, lower=0, upper=data$pseudom.time[i])$value
      haz.phys[i] <- integrate(haz.integrand.phys, lower=0, upper=data$pseudom.time[i])$value
      # average HR up to the pseudo time
      true.hr[i] <- haz.clin[i]/haz.phys[i]
      if ((data$RMST.arm.lower[i]  <= true.rmst[i]) & (data$RMST.arm.upper[i]  >= true.rmst[i])){cover.rmst[i]=1}
      else {cover.rmst[i]=0}
      if ((data$hr.arm.lower[i]  <= true.hr[i]) & (data$hr.arm.upper[i]  >= true.hr[i])){cover.hr[i]=1}
      else {cover.hr[i]=0}
      if (is.na(data$ratio.arm.lower[i]) | is.na(data$ratio.arm.upper[i])){cover.rat[i]=0}
      else if ((data$ratio.arm.lower[i]  <= true.rat) & (data$ratio.arm.upper[i]  >= true.rat)){cover.rat[i]=1}
      else {cover.rat[i]=0}
    }
    hr.cover <- mean(cover.hr)
    RMST.cover <- mean(cover.rmst)
    rat.cover <- mean(cover.rat)
    
  }
}
  
else if (estimand=="inter"){
  # getting the power of each testing method
  
    # getting the power of each testing method
    power.lr <- round(sum(data$pvalue<.05)/length(data$pvalue), digits=3)
    power.hr <- round(sum(data$hr.int.pval<0.05)/length(data$hr.int.pval), digit=3)
    power.rat <- round(sum(data[,56]<0.05)/length(data[,53]), digits=3)
    power.rmst <- round(sum(data[,76]<0.05)/length(data[,73]), digits=3)
  
  # getting the estimated treatment difference using each method
  hr.mean <- mean(data$hr.int.mean)
  hr.mean.SE <- sd(data$hr.int.mean)*qnorm(.975)
  
  rat.mean <- mean(data[,53])
  rat.mean.SE <- sd(data[,53])*qnorm(.975)
  
  rmst.mean <- mean(data[,73])
  rmst.mean.SE <- sd(data[,73])*qnorm(.975)
  
  # calculating the true contrasts using exponential distributions and given median survivals
  true.hr <- (medianApos/medianBpos)/(medianAneg/medianBneg)
  true.rat <- (exp(-(log(2)/medianBpos)*data$pseudor.time[1])/exp(-(log(2)/medianApos)*data$pseudor.time[1]))/(exp(-(log(2)/medianBneg)*data$pseudor.time[1])/exp(-(log(2)/medianAneg)*data$pseudor.time[1]))
  integrand1 <- function(x){exp(-(x*(log(2)/medianBpos)))}
  integrand2 <- function(x){exp(-(x*(log(2)/medianApos)))}
  integrand3 <- function(x){exp(-(x*(log(2)/medianBneg)))}
  integrand4 <- function(x){exp(-(x*(log(2)/medianAneg)))}
  int1 <- NULL
  int2 <- NULL
  int3 <- NULL
  int4 <- NULL
  true.rmst <- NULL
  cover.rmst <- NULL
  cover.hr <- NULL
  cover.rat <- NULL
  for (i in 1:length(data$pseudom.time)){
    int1[i] <- integrate(integrand1, 0,data$pseudom.time[i])$value
    int2[i] <- integrate(integrand2, 0,data$pseudom.time[i])$value
    int3[i] <- integrate(integrand3, 0,data$pseudom.time[i])$value
    int4[i] <- integrate(integrand4, 0,data$pseudom.time[i])$value
    true.rmst[i] <- ((int1[i]-int2[i]) - (int3[i]-int4[i]))
    if ((data$RMST.int.lower[i]  <= true.rmst[i]) & (data$RMST.int.upper[i]  >= true.rmst[i])){cover.rmst[i]=1}
    else {cover.rmst[i]=0}
    if ((data$hr.int.lower[i]  <= true.hr) & (data$hr.int.upper[i]  >= true.hr)){cover.hr[i]=1}
    else {cover.hr[i]=0}
    if (is.na(data$ratio.int.lower[i]) | is.na(data$ratio.int.upper[i])){cover.rat[i]=0}
    else if ((data$ratio.int.lower[i]  <= true.rat) & (data$ratio.int.upper[i]  >= true.rat)){cover.rat[i]=1}
    else {cover.rat[i]=0}
  }
  hr.cover <- mean(cover.hr)
  RMST.cover <- mean(cover.rmst)
  rat.cover <- mean(cover.rat)
}
  
  #getting number of events - same for each obs
  events <- data$events[1]
  # getting the mean number enrolled
  enroll <- ceiling(mean(data$enroll))
  enroll.SE <- round(sd(data$enroll)*qnorm(.975), digits=2)
  # getting mean analysis time
  time <- round(mean(data$analysis.time), digits=1)
  
  # rounding for output
  hr.mean <- round(hr.mean, digits=2)
  hr.mean.SE <- round(hr.mean.SE, digits=2)
  
  rat.mean <- round(rat.mean, digits=2)
  rat.mean.SE <- round(rat.mean.SE, digits=2)
  
  rmst.mean <- round(rmst.mean, digits=2)
  rmst.mean.SE <- round(rmst.mean.SE, digits=2)
  
  # hr <- paste(hr.mean, hr.mean.SE, sep="+")
  # rat <- paste(rat.mean, rat.mean.SE, sep="+")
  # rmst <- paste(rmst.mean, rmst.mean.SE, sep="+")
  hr <- hr.mean
  rat <- rat.mean
  rmst <- rmst.mean
  
  true.rmst <- round(mean(true.rmst), digits=2)
  true.rat <- round(mean(true.rat), digits=2)
  true.hr <- round(mean(true.hr), digits=2)
  
  result <- data.frame(matrix(nrow=1, ncol=16))
  colnames(result) <- c("Events","Size","Duration","Logrank Power",
                        "HR Power","Ratio Power","RMST Power", "HR coverage", "RMST coverage", "Ratio coverage",
                        "Mean HR","Mean Ratio","Mean RMST", "True HR", "True Ratio", "True RMST")
  
  result[1,] <- c(events, enroll, time, power.lr, 
              power.hr,power.rat,power.rmst, hr.cover, RMST.cover, rat.cover,
              hr, rat, rmst, true.hr, true.rat, true.rmst)
  result[,1:16] <- as.numeric(result[,1:16])
  
  return(result)
}


large <- function(x){
  paste0('{\\Large ', x, '}')
}
bold <- function(x){
  paste0('{\\bfseries ', x, '}')
}
italic <- function(x){
  paste0('{\\emph{ ', x, '}}')
}

# function to create a table from multiple result data frames for a given scenario
fintable <- function(caption,
                     estimand,
                     design.type,
                     smple.size,
                     accru,
                     LTFU,
                     medianBpos,
                     medianBneg,
                     medianApos,
                     medianAneg,
                     marker.pos,
                     phys.choice,
                     rpos,              # randomization to treatment b in the positives
                     type1=FALSE
                     ){
                     
    Mpos <- marker.pos
    pc <- phys.choice
    rp <- rpos

    result <- data.frame(matrix(nrow=length(smple.size)*2, ncol=16))
    # colnames(result) <- c("Events","Size","Duration","Logrank Power",
    #                       "HR Power","Ratio Power","RMST Power", "HR coverage", "RMST coverage", "Ratio coverage",
    #                       "Mean HR","Mean Ratio","Mean RMST", "True HR", "True Ratio", "True RMST")
  if (type1=="FALSE"){
    for (i in 1:length(smple.size)){
      if (estimand=="clin"){
        filename <- paste("Results",estimand,
                          paste(estimand, design.type, smple.size[i], "accru", accru, "exp", 
                                "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                                "M", Mpos, "LTFU", LTFU, "PC", pc, "rds", sep="."), sep="/")
      }
      else if (estimand=="subgrp"){
        filename <- paste("Results",estimand,
                          paste(estimand, design.type, smple.size[i], "accru", accru, "exp", 
                                "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                                "M", Mpos, "LTFU", LTFU, "rds", sep="."), sep="/")
      }
      else if (estimand=="inter"){
        filename <- paste("Results",estimand,
                          paste(estimand, design.type, smple.size[i], "accru", accru,"exp", 
                                "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
                                "M", Mpos, "LTFU", LTFU, "rpos",rp, "rds", sep="."), sep="/")
      }
      x <- readRDS(file=filename)
      y <- analyze.exp(data=x, 
                       estimand=estimand, 
                       design.type=design.type,
                       medianBpos=medianBpos,
                       medianBneg=medianBneg,
                       medianApos=medianApos,
                       medianAneg=medianAneg,
                       marker.pos=marker.pos,
                       phys.choice=phys.choice,
                       type1=type1)
      result[1,] <- c("Events","Size","Duration","Logrank Power",
                      "HR Power","Ratio Power","RMST Power", "HR coverage", "RMST coverage", "Ratio coverage",
                      "Mean HR","Mean Ratio","Mean RMST", "True HR", "True Ratio", "True RMST")
      result[(2*i),] <- y
      result[((2*(i+1))-1),] <- NA
    }}
    
    else if (type1==TRUE){
      if (estimand=="clin"){
        filename <- paste("Results",estimand,paste("clin", "strategy","type1",  smple.size, "accru", accru, "exp", 
                          "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                          "M", Mpos, "LTFU", LTFU, "PC", "1", "rds", sep="."), sep="/")
      }
      else if (estimand=="subgrp"){
        filename <- paste("Results",estimand,paste("subgrp", "enrich","type1",  smple.size, "accru", accru, "exp", 
                          "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                          "M", ".75", "LTFU", LTFU, "rds", sep="."), sep="/")
      }
      else if (estimand=="inter"){
        filename <- paste("Results",estimand,paste("inter", "stratify","tpe1",  smple.size, "accru", accru,"exp", 
                          "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
                          "M", Mpos, "LTFU", LTFU, "rpos",rp, "rds", sep="."), sep="/")
      }
      x <- readRDS(file=filename)
      y <- analyze.exp(data=x, 
                       estimand=estimand, 
                       design.type=design.type,
                       medianBpos=medianBpos,
                       medianBneg=medianBneg,
                       medianApos=medianApos,
                       medianAneg=medianAneg,
                       marker.pos=marker.pos,
                       phys.choice=phys.choice,
                       type1=type1)
      result[1,] <- c("Events","Size","Duration","Logrank Power",
                      "HR Power","Ratio Power","RMST Power", "HR coverage", "RMST coverage", "Ratio coverage",
                      "Mean HR","Mean Ratio","Mean RMST", "True HR", "True Ratio", "True RMST")
      result[3,] <- y
      result[2,] <- NA
      
    }
    
    result <- data.frame(sapply(result, as.character))
    
    final <- xtable(result, include.rownames=FALSE, caption=caption)
    align(final) <- c("rr|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|", rep("p{0.4in}",16))
    
    return(print(final, include.rownames=FALSE, caption.placement="top"))
}




power.curve <- function(estimand,
                        design.type,
                        smple.size,
                        accru,
                        LTFU,
                        medianBpos,
                        medianBneg,
                        medianApos,
                        medianAneg,
                        marker.pos,
                        phys.choice,
                        rpos              # randomization to treatment b in the positives
){
  
  power.lr <- matrix(nrow=length(smple.size), ncol=4)
  power.hr <- matrix(nrow=length(smple.size), ncol=4)
  power.rat <- matrix(nrow=length(smple.size), ncol=4)
  power.rmst <- matrix(nrow=length(smple.size), ncol=4)
  
  Mpos <- marker.pos
  pc <- phys.choice
  rp <- rpos
  for (i in 1:length(smple.size)){
    if (estimand=="clin"){
      filename <- paste("Results",estimand,
                        paste(estimand, design.type, smple.size[i], "accru", accru, "exp", 
                              "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                              "M", Mpos, "LTFU", LTFU, "PC", pc, "rds", sep="."), sep="/")
    }
    else if (estimand=="subgrp"){
      filename <- paste("Results",estimand,
                        paste(estimand, design.type, smple.size[i], "accru", accru, "exp", 
                              "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, 
                              "M", Mpos, "LTFU", LTFU, "rds", sep="."), sep="/")
    }
    else if (estimand=="inter"){
      filename <- paste("Results",estimand,
                        paste(estimand, design.type, smple.size[i], "accru", accru,"exp", 
                              "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
                              "M", Mpos, "LTFU", LTFU, "rpos",rp, "rds", sep="."), sep="/")
    }
    
    x <- readRDS(file=filename)
    result <- analyze.exp(data=x,              # result data frame 
                         estimand=estimand,          # subgrp, clin, inter
                         design.type=design.type,        # trial design type, i.e, enrich, stratify, strategy, modstrat
                         medianBpos=medianBpos,
                         medianBneg=medianBneg,
                         medianApos=medianApos,
                         medianAneg=medianAneg,
                         marker.pos=marker.pos,
                         phys.choice=phys.choice)
    
    power.lr[i,] <- c(result[1,1],design.type,"power.type"="lr", result[1,4])
    power.hr[i,] <- c(result[1,1],design.type,"power.type"="hr", result[1,5])
    power.rat[i,] <- c(result[1,1],design.type,"power.type"="rat", result[1,6])
    power.rmst[i,] <- c(result[1,1],design.type,"power.type"="rmst", result[1,7])
    if (estimand=="inter"){
      temp <- rbind(power.hr, power.rat, power.rmst)
    }
    else {temp <- rbind(power.lr, power.hr, power.rat, power.rmst)}
    
}
  power <- data.frame(temp)
  colnames(power) <- c("events","design","power.type" ,
                       "power")
  
  power[,c(1,4)] <- apply(power[ , c(1,4)], 2,            # Specify own function within apply
                          function(x) as.numeric(as.character(x)))
  
    xscale <- scale_x_continuous(name="Events", breaks=seq(from=smple.size[1], 
                                  to=smple.size[length(smple.size)], by=(smple.size[length(smple.size)]-smple.size[1])/5), 
                                 limits <- c(smple.size[1],smple.size[length(smple.size)]))
    expand <- expand_limits(x=smple.size[1], y=c(0.2, 1))
  
  yscale <- scale_y_continuous(name="Power", breaks=seq(from=0.2, to=1, by=.2), limits <- c(0.2,1))
  
  fig1 <- ggplot(data=power, mapping=aes(x=events, y=power, color=power.type, group=power.type)) + 
    geom_point() + geom_line() + 
    theme(legend.position = "none") + 
    expand + xscale + yscale + theme(axis.text.x = element_text(angle=90))
  
  return(fig1)
}








