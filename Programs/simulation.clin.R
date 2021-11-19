#==============================================================================
# FILENAME: simulation.clin.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: calls the simulation functions to simulate the different scenarios and ave the result datasets in Results folder 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 3.6.1
#==============================================================================
#Notes: 
### this program simulates different scenarios for the estimand for clinical utility




# =============================================================================

## sourcing the source program
source("Programs/sim_source.R")

# selecting a vecotr of sample sizes
smple.size <- seq(from=250, to=1000, by=50)

Bpos <- c(9,12)
Bneg <- c(9,12)
Apos <- c(9,12)
Aneg <- c(9,12)


 for (j in 1:length(Bpos)){
   for (k in 1:length(Bneg)){
     for (s in 1:length(Apos)){
       for (t in 1:length(Aneg)){
        
        # selecting survival distributions and median survival for each of the 4 groups
        distr_Bpos <- "exp"
        medianBpos <- Bpos[j]
        
        distr_Bneg <- "exp"
        medianBneg <- Bneg[k]
        
        distr_Apos <- "exp"
        medianApos <- Apos[s]
        
        distr_Aneg <- "exp"
        medianAneg <- Aneg[t]
        

accru <- 25



set.seed(12)
#getting results for each smple.size based on 1000 replications
for (i in 1:length(smple.size)){
  #undebug(eval.scen)
  x <- eval.scen(estimand="clin", 
                 reps=100, 
                 max_enroll=4000,  
                 n=smple.size[i], 
                 num_interim=1, 
                 int_timing=1, 
                 alpha=.025, 
                 low_err=0.1,  
                 bound.type=c(1,1), 
               #  time.points=time.point,    
                # RM.times=time.point, 
                 
                 accru.rate=accru,           
                 loss.fu.rate=.02,          
                 marker.pos=0.4,           
                 rand.arm=.5,             
                 rand.pos=.5,            
                 rand.neg=.5,             
                 rand=0.5,
                 phys.choice.pos=.2,        
                 phys.choice.neg=.8,
                 
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
  
  temp_strategy <- x$clin_strategy
  temp_stratify <- x$clin_stratify
  temp_modstrat <- x$clin_modstrat
  
  filename_strategy <- paste("clin", "strategy", smple.size[i], "accru", accru, "exp", 
                           "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds", sep=".")
  
  filename_stratify <- paste("clin", "stratify", smple.size[i], "accru", accru,"exp", 
                          "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,"rds", sep=".")
  filename_modsrtat <- paste("clin", "modstrat", smple.size[i], "accru", accru,"exp", 
                        "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg, "rds",sep=".")
  
  saveRDS(temp_strategy, file=paste("Results", "clin",filename_strategy, sep="/"))
  saveRDS(temp_stratify, file=paste("Results", "clin",filename_stratify, sep="/"))
  saveRDS(temp_modstrat, file=paste("Results", "clin",filename_modsrtat, sep="/"))
  
  
}
 }}}}