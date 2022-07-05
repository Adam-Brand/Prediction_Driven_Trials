#==============================================================================
# FILENAME: simulation.inter.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: calls the simulation functions to simulate the different scenarios and ave the result datasets in Results folder 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 4.0.5
#==============================================================================
#Notes: 
### this program simulates different scenarios for the estimand for interaction term




# =============================================================================

## sourcing the source program
source("Programs/sim_source.R")

# selecting a vector of sample sizes
smple.size <- seq(from=250, to=1000, by=50)

Bpos <- c(9,12,15)
Bneg <- 9
Apos <- 9
Aneg <- 12
M <- c(0.25,.5)
LTFU <- .02
rpos <- 0.5



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


accru <- 10


set.seed(22456785)
#getting results for each smple.size based on 1000 replications
for (i in 1:length(smple.size)){
  #undebug(eval.scen)
  x <- eval.scen(estimand="inter", 
                 reps=1000, 
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

  
  saveRDS(temp_stratify, file=paste("Results","inter",filename_stratify, sep="/"))
  #saveRDS(temp_modstrat, file=paste("Results", "inter",filename_modsrtat, sep="/"))
  
}
}}}






# increasing sample sizes to get 90% power for one scenario
smple.size <- seq(from=1050, to=2000, by=50)

Bpos <- 9
Bneg <- 9
Apos <- 9
Aneg <- 12
M <- c(0.25,.5)
LTFU <- .02
rpos <- 0.5



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
      
      
      accru <- 10
      
      
      set.seed(22356785)
      #getting results for each smple.size based on 1000 replications
      for (i in 1:length(smple.size)){
        #undebug(eval.scen)
        x <- eval.scen(estimand="inter", 
                       reps=1000, 
                       max_enroll=6000,  
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
        
        
        saveRDS(temp_stratify, file=paste("Results","inter",filename_stratify, sep="/"))
        #saveRDS(temp_modstrat, file=paste("Results", "inter",filename_modsrtat, sep="/"))
        
      }
    }}}


### function calls to get type 1 error

# selecting a vecotr of sample sizes
smple.size <- c(2000,700,400,2000,500,300)

Bpos <- c(9,12,15,9,12,15)
Bneg <- Bpos
Apos <- c(12,9,12,12,9,12)
Aneg <- Apos
M <- c(.25,.25,.25,.5,.5,.5)
LTFU <- .02
rpos <- 0.5


rp <- rpos
rn <- 1-rp

set.seed(22456785)
for (i in 1:length(smple.size)){ 
  Mpos <- M[i]
  type1_int <- eval.scen(estimand="inter", 
                         reps=1000, 
                         max_enroll=6000,  
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
                         
                         medposB=Bpos[i],             
                         distrposB="exp",           
                         shapeposB=1,         
                         
                         medposA=Apos[i],             
                         distrposA="exp",           
                         shapeposA=1,           
                         
                         mednegB=Bneg[i],             
                         distrnegB="exp",           
                         shapenegB=1,           
                         
                         mednegA=Aneg[i],             
                         distrnegA="exp",          
                         shapenegA=1          
  )
  
  temp_stratify_tp1 <- type1_int$inter_stratify
  
  filename_type1 <- paste("inter", "stratify","tpe1",  smple.size[i], "accru", accru,"exp", 
                          "Bpos",Bpos[i], "Apos", Apos[i], "Bneg",Bneg[i], "Aneg",  Aneg[i],
                          "M", Mpos, "LTFU", LTFU, "rpos",rp, "rds", sep=".")
  
  
  saveRDS(temp_stratify_tp1, file=paste("Results","inter",filename_type1, sep="/"))
}



