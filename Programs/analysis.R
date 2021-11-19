#==============================================================================
# FILENAME: analysis.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: reads in the result data frames, summarizes results and puts into tables and figures
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
analyze <- function(data,              # result data frame 
                    estimand,          # subgrp, clin, inter
                    design.type        # trial design type, i.e, enrich, stratify, strategy, modstrat
){
  
if (estimand=="subgrp"){  
  if (design.type=="enrich"){
  # getting the power of each testing method
  power.lr <- round(sum(data$stop.eff==1)/length(data$stop.eff), digits=3)
  power.hr <- round(sum(data$hr.pval<0.05 & data$hr.mean < 1)/length(data$hr.pval), digit=3)
  power.rat <- round(sum(data[,21]<0.05 & data[,18]<1)/length(data[,21]), digits=3)
  power.rmst <- round(sum(data[,25]<0.05 & data[,22]>0)/length(data[,25]), digits=3)
  
  # getting the estimated treatment difference using each method
  hr.mean <- round(mean(data$hr.mean), digits=2)
  hr.mean.SE <- round(sd(data$hr.mean)*qnorm(.975), digits=2)
  
  rat.mean <- round(mean(data[,18]), digits=2)
  rat.mean.SE <- round(sd(data[,18])*qnorm(.975), digits=2)
  
  rmst.mean <- round(mean(data[,22]), digits=2)
  rmst.mean.SE <- round(sd(data[,22])*qnorm(.975), digits=2)
  }
  
  
  else if (design.type=="stratify" | design.type=="modstrat"){
    # getting the power of each testing method
    power.lr <- round(sum(data$stop.eff==1 & data$exp.diff>0)/length(data$stop.eff), digits=3)
    power.hr <- round(sum(data$hr.pos.pval<0.05 & data$hr.pos.mean < 1)/length(data$hr.pos.pval), digit=3)
    power.rat <- round(sum(data[,48]<0.05 & data[,45]<1)/length(data[,45]), digits=3)
    power.rmst <- round(sum(data[,68]<0.05 & data[,65]>0)/length(data[,65]), digits=3)
    
    # getting the estimated treatment difference using each method
    hr.mean <- round(mean(data$hr.pos.mean), digits=2)
    hr.mean.SE <- round(sd(data$hr.pos.mean)*qnorm(.975), digits=2)
    
    rat.mean <- round(mean(data[,45]), digits=2)
    rat.mean.SE <- round(sd(data[,45])*qnorm(.975), digits=2)
    
    rmst.mean <- round(mean(data[,65]), digits=2)
    rmst.mean.SE <- round(sd(data[,65])*qnorm(.975), digits=2)
  }
    
}
  
else if (estimand=="clin"){
  if (design.type=="stratify" | design.type=="modstrat"){
    # getting the power of each testing method
    power.lr <- round(sum(data$stop.eff==1 & data$exp.diff>0)/length(data$stop.eff), digits=3)
    power.hr <- round(sum(data$hr.arm.pval<0.05 & data$hr.arm.mean < 1)/length(data$hr.arm.pval), digit=3)
    power.rat <- round(sum(data[,60]<0.05 & data[,57]<1)/length(data[,57]), digits=3)
    power.rmst <- round(sum(data[,80]<0.05 & data[,77]>0)/length(data[,77]), digits=3)
    
    # getting the estimated treatment difference using each method
    hr.mean <- round(mean(data$hr.arm.mean), digits=2)
    hr.mean.SE <- round(sd(data$hr.arm.mean)*qnorm(.975), digits=2)
    
    rat.mean <- round(mean(data[,57]), digits=2)
    rat.mean.SE <- round(sd(data[,57])*qnorm(.975), digits=2)
    
    rmst.mean <- round(mean(data[,77]), digits=2)
    rmst.mean.SE <- round(sd(data[,77])*qnorm(.975), digits=2)
  }
  else if (design.type=="strategy"){
    # getting the power of each testing method
    power.lr <- round(sum(data$stop.eff==1 & data$exp.diff>0)/length(data$stop.eff), digits=3)
    power.hr <- round(sum(data$hr.arm.pval<0.05 & data$hr.arm.mean < 1)/length(data$hr.arm.pval), digit=3)
    power.rat <- round(sum(data[,28]<0.05 & data[,25]<1)/length(data[,25]), digits=3)
    power.rmst <- round(sum(data[,32]<0.05 & data[,29]>0)/length(data[,29]), digits=3)
    
    # getting the estimated treatment difference using each method
    hr.mean <- round(mean(data$hr.arm.mean), digits=2)
    hr.mean.SE <- round(sd(data$hr.arm.mean)*qnorm(.975), digits=2)
    
    rat.mean <- round(mean(data[,25]), digits=2)
    rat.mean.SE <- round(sd(data[,25])*qnorm(.975), digits=2)
    
    rmst.mean <- round(mean(data[,29]), digits=2)
    rmst.mean.SE <- round(sd(data[,29])*qnorm(.975), digits=2)
  }
}
  
else if (estimand=="inter"){
  # getting the power of each testing method
  power.lr <- round(sum(data$pvalue<.05)/length(data$pvalue), digits=3)
  power.hr <- round(sum(data$hr.int.pval<0.05)/length(data$hr.int.pval), digit=3)
  power.rat <- round(sum(data[,56]<0.05)/length(data[,53]), digits=3)
  power.rmst <- round(sum(data[,76]<0.05)/length(data[,73]), digits=3)
  
  # getting the estimated treatment difference using each method
  hr.mean <- round(mean(data$hr.int.mean), digits=2)
  hr.mean.SE <- round(sd(data$hr.int.mean)*qnorm(.975), digits=2)
  
  rat.mean <- round(mean(data[,53]), digits=2)
  rat.mean.SE <- round(sd(data[,53])*qnorm(.975), digits=2)
  
  rmst.mean <- round(mean(data[,73]), digits=2)
  rmst.mean.SE <- round(sd(data[,73])*qnorm(.975), digits=2) 
}
  
  #getting number of events - same for each obs
  events <- data$events[1]
  # getting the mean number enrolled
  enroll <- ceiling(mean(data$enroll))
  # getting mean analysis time
  time <- round(mean(data$analysis.time), digits=1)
  
  hr <- paste(hr.mean, hr.mean.SE, sep="+")
  rat <- paste(rat.mean, rat.mean.SE, sep="+")
  rmst <- paste(rmst.mean, rmst.mean.SE, sep="+")
  
  result <- data.frame(matrix(nrow=1, ncol=10))
  colnames(result) <- c("Events","Size","Duration","Logrank Power","HR Power",
                        "Ratio Power","RMST Power","Mean HR","Mean Ratio","Mean RMST")
  
  result[1,] <- c(events, enroll, time, power.lr, 
              power.hr,power.rat,power.rmst,
              hr, rat, rmst)
  result[,1:7] <- as.numeric(result[,1:7])
  
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
                     time.point,
                     accru,
                     medianBpos,
                     medianBneg,
                     medianApos,
                     medianAneg){

    result <- data.frame(matrix(nrow=length(smple.size)*2, ncol=10))
    colnames(result) <- c("Events","Size","Duration","Logrank Power","HR Power","Ratio Power","RMST Power",
                          "Mean HR", "Mean Ratio", "Mean RMST")
    for (i in 1:length(smple.size)){
      filename <- paste("Results",estimand,
                        paste(estimand, design.type, smple.size[i], "accru", accru,"time",time.point, "exp", 
                              "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
      x <- readRDS(file=filename)
      y <- analyze(data=x, estimand=estimand, design.type=design.type)
      result[((2*i)-1),] <- y
      result[(2*i),] <- NA
    }
    
    result <- data.frame(sapply(result, as.character))
    
    final <- xtable(result, include.rownames=FALSE, caption=caption)
    align(final) <- "rr|c|c|c|c|c|c|c|c|c"
    
    return(print(final, include.rownames=FALSE, sanitize.colnames.function = bold,
                 caption.placement="top"))
}




power.curve <- function(caption,
                        estimand,
                        smple.size,
                        time.point,
                        accru,
                        medianBpos,
                        medianBneg,
                        medianApos,
                        medianAneg){
  
  ## reading the data and extracting the power for the different designs at different sample sizes
if (estimand=="subgrp"){
  power <- data.frame(matrix(nrow=(length(smple.size)*3), ncol=6))
  colnames(power) <- c("size", "design","lr","hr","rat","rmst")
  for (i in 1:length(smple.size)){
    filename_en <- paste("Results",estimand,
                      paste(estimand, "enrich", smple.size[i], "accru", accru,"time",time.point, "exp", 
                            "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
    filename_strat <- paste("Results",estimand,
                         paste(estimand, "stratify", smple.size[i], "accru", accru,"time",time.point, "exp", 
                               "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
    filename_mod <- paste("Results",estimand,
                            paste(estimand, "modstrat", smple.size[i], "accru", accru,"time",time.point, "exp", 
                                  "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
    en <- readRDS(file=filename_en)
    strat<- readRDS(file=filename_strat)
    mod <- readRDS(file=filename_mod)
    
    pow_en <- analyze(data=en, estimand=estimand, design.type="enrich")
    pow_strat <- analyze(data=strat, estimand=estimand, design.type="stratify")
    pow_mod <- analyze(data=mod, estimand=estimand, design.type="modstrat")
    
    power[((3*i)-2),] <- c("size"=pow_en[2],"design"="enrich", "lr"=pow_en[4], 
                         "hr"=pow_en[5],"rat"=pow_en[6],"rmst"=pow_en[7])
    power[((3*i)-1),] <- c("size"=pow_strat[2],"design"="stratify", "lr"=pow_strat[4], 
                         "hr"=pow_strat[5],"rat"=pow_strat[6],"rmst"=pow_strat[7])
    power[(3*i),] <- c("size"=pow_mod[2],"design"="modstrat", "lr"=pow_mod[4], 
                         "hr"=pow_mod[5],"rat"=pow_mod[6],"rmst"=pow_mod[7])
    power$size <- as.numeric(power$size)
    
  }
}
  
 else if (estimand=="clin"){
    power <- data.frame(matrix(nrow=(length(smple.size)*3), ncol=6))
    colnames(power) <- c("size", "design","lr","hr","rat","rmst")
    for (i in 1:length(smple.size)){
      filename_clin <- paste("Results",estimand,
                           paste(estimand, "strategy", smple.size[i], "accru", accru,"time",time.point, "exp", 
                                 "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
      filename_strat <- paste("Results",estimand,
                              paste(estimand, "stratify", smple.size[i], "accru", accru,"time",time.point, "exp", 
                                    "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
      filename_mod <- paste("Results",estimand,
                            paste(estimand, "modstrat", smple.size[i], "accru", accru,"time",time.point, "exp", 
                                  "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
      clin <- readRDS(file=filename_clin)
      strat<- readRDS(file=filename_strat)
      mod <- readRDS(file=filename_mod)
      
      pow_clin <- analyze(data=clin, estimand=estimand, design.type="strategy")
      pow_strat <- analyze(data=strat, estimand=estimand, design.type="stratify")
      pow_mod <- analyze(data=mod, estimand=estimand, design.type="modstrat")
      
      power[((3*i)-2),] <- c("size"=pow_clin[2],"design"="strategy", "lr"=pow_clin[4], 
                             "hr"=pow_clin[5],"rat"=pow_clin[6],"rmst"=pow_clin[7])
      power[((3*i)-1),] <- c("size"=pow_strat[2],"design"="stratify", "lr"=pow_strat[4], 
                             "hr"=pow_strat[5],"rat"=pow_strat[6],"rmst"=pow_strat[7])
      power[(3*i),] <- c("size"=pow_mod[2],"design"="modstrat", "lr"=pow_mod[4], 
                         "hr"=pow_mod[5],"rat"=pow_mod[6],"rmst"=pow_mod[7])
      power$size <- as.numeric(power$size)
    }
 }
  else if (estimand=="inter"){
    power <- data.frame(matrix(nrow=(length(smple.size)*2), ncol=6))
    colnames(power) <- c("size", "design","lr","hr","rat","rmst")
    for (i in 1:length(smple.size)){
      
      filename_strat <- paste("Results",estimand,
                              paste(estimand, "stratify", smple.size[i], "accru", accru,"time",time.point, "exp", 
                                    "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
      filename_mod <- paste("Results",estimand,
                            paste(estimand, "modstrat", smple.size[i], "accru", accru,"time",time.point, "exp", 
                                  "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep="."), sep="/")
      
      strat<- readRDS(file=filename_strat)
      mod <- readRDS(file=filename_mod)
      
      
      pow_strat <- analyze(data=strat, estimand=estimand, design.type="stratify")
      pow_mod <- analyze(data=mod, estimand=estimand, design.type="modstrat")
      
      
      power[((2*i)-1),] <- c("size"=pow_strat[2],"design"="stratify", "lr"=pow_strat[4], 
                             "hr"=pow_strat[5],"rat"=pow_strat[6],"rmst"=pow_strat[7])
      power[(2*i),] <- c("size"=pow_mod[2],"design"="modstrat", "lr"=pow_mod[4], 
                         "hr"=pow_mod[5],"rat"=pow_mod[6],"rmst"=pow_mod[7])
      power$size <- as.numeric(power$size)
    }
  }
  
  return(power)
}

# selecting a vecotr of sample sizes
smple.size <- seq(from=250, to=1000, by=50)

# selecting survival distributions and median survival for each of the 4 groups
distr_Bpos <- "exp"
medianBpos <- 12

distr_Bneg <- "exp"
medianBneg <- 9

distr_Apos <- "exp"
medianApos <- 9

distr_Aneg <- "exp"
medianAneg <- 12

time.point <- 12

accru <- 25

fintable(caption="Stratified Design. 1/3 increase in Median survival for Positives on Treatment B (12v9) and for Negatives on Treatment A (9v12). Accrual Rate 25. Spot Comparisons at time 12. Exponential Distribution. 0.4 positive proportion. 1:1 randomization. .2 physician agreement in positives, and .8 in negatives.",
                 estimand="inter",
                 design.type="stratify",
                 smple.size=smple.size,
                 time.point=time.point,
                 accru=accru,
                 medianBpos=medianBpos,
                 medianBneg=medianBneg,
                 medianApos=medianApos,
                 medianAneg=medianAneg)

#debug(power.curve)
test <- power.curve(caption="Power",
                 estimand="inter",
                 smple.size=smple.size,
                 time.point=time.point,
                 accru=accru,
                 medianBpos=medianBpos,
                 medianBneg=medianBneg,
                 medianApos=medianApos,
                 medianAneg=medianAneg)

fig1 <- ggplot(data=test, mapping=aes(x=size, y=lr, color=design, group=design)) + geom_point() + geom_line()
fig1

fig2 <- ggplot(data=test, mapping=aes(x=size, y=hr, color=design, group=design)) + geom_point() + geom_line()
fig2

fig3 <- ggplot(data=test, mapping=aes(x=size, y=rat, color=design, group=design)) + geom_point() + geom_line()
fig3

fig4 <- ggplot(data=test, mapping=aes(x=size, y=rmst, color=design, group=design)) + geom_point() + geom_line()
fig4

plot_grid(fig1, fig2, fig3, fig4,
          labels=c("A","B","C","D"),
          ncol=2, nrow=2)

x <- readRDS(file="Results/subgrp/subgrp.enrich.250.accru.25.time.12.exp.Bpos.12.Apos.9.Bneg.9.Aneg.12.rds")

undebug(analyze)
y <- analyze(data=x,
             estimand="subgrp",
             design.type="enrich")

# selecting a vecotr of sample sizes
smple.size <- seq(from=250, to=300, by=50)

# selecting survival distributions and median survival for each of the 4 groups
distr_Bpos <- "exp"
medianBpos <- 12

distr_Bneg <- "exp"
medianBneg <- 9

distr_Apos <- "exp"
medianApos <- 9

distr_Aneg <- "exp"
medianAneg <- 12

time.point <- 12

accru <- 25

debug(fintable)
test <- fintable(caption="Optimal Treatment for Positives, Enrichment Design. 1/3 increase in Median survival for Positives on Treatment B (12v9). 1/4 Decrease in survival for Negatives on Treatment B (9v12). Accrual Rate 25. Spot Comparisons at 12.",
                 estimand="subgrp",
                 design.type="enrich",
                 smple.size=smple.size,
                 time.point=time.point,
                 accru=accru,
                 medianBpos=medianBpos,
                 medianBneg=medianBneg,
                 medianApos=medianApos,
                 medianAneg=medianAneg)







