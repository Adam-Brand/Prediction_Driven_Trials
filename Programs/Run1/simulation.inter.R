#==============================================================================
# FILENAME: simulation.inter.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: calls the simulation functions to simulate the different scenarios and ave the result datasets in Results folder 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 3.6.1
#==============================================================================
#Notes: 
### this program simulates different scenarios for the estimand for interaction term




# =============================================================================

## sourcing the source program
source("Programs/sim_source.R")

# selecting a vecotr of sample sizes
smple.size <- seq(from=250, to=1000, by=50)

Bpos <- c(9,12,15,21)
Bneg <- 9
Apos <- 9
Aneg <- 12
M <- c(0.25,.5)
LTFU <- .02
rpos <- c(0.5, 0.67)



for (j in 1:length(Bpos)){
  for (k in 1:length(M)){
    for (s in 1:length(rpos)){
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
    
    rp <- rpos[s]
    rn <- 1-rp


accru <- 25


set.seed(12)
#getting results for each smple.size based on 1000 replications
for (i in 1:length(smple.size)){
  #undebug(eval.scen)
  x <- eval.scen(estimand="inter", 
                 reps=100, 
                 max_enroll=4000,  
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
                 rand.pos=rp,            
                 rand.neg=rn,             
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
  
  
  temp_stratify <- x$inter_stratify
 # temp_modstrat <- x$inter_modstrat
  
  
  filename_stratify <- paste("inter", "stratify", smple.size[i], "accru", accru,"exp", 
                             "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
                             "M", Mpos, "LTFU", LTFU, "rpos",rp, "rds", sep=".")
  # filename_modsrtat <- paste("inter", "modstrat", smple.size[i], "accru", accru,"exp", 
  #                            "Bpos",medianBpos, "Apos", medianApos, "Bneg",medianBneg, "Aneg",  medianAneg,
  #                            "M", Mpos, "LTFU", LTFU, "rds",sep=".")

  
  saveRDS(temp_stratify, file=paste("Results", "inter",filename_stratify, sep="/"))
  #saveRDS(temp_modstrat, file=paste("Results", "inter",filename_modsrtat, sep="/"))
  
}
}}}