#==============================================================================
# FILENAME: simulation.clin.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: Program to simulate a comparison of the Shih estimand and the Paper's estimand 
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 4.0.5
#==============================================================================
#Notes: 
### this program simulates different scenarios for the estimand for clinical utility




# =============================================================================

## sourcing the source program
source("Programs/sim_source.R")

# selecting a vector of sample sizes
reps <- 1000
smple.size <- 400
Bpos <- c(9,12,21)
LTFU <- .02
PC <- 1
accru <- 10
pos <- c(0.25,0.5)


for (j in 1:length(Bpos)){
  for (i in 1:length(pos)){
    
    result <- data.frame(matrix(nrow=reps, ncol=16))
    colnames(result) <- c("events","stdy.time","lr.pval.shih","lr.pval.strat",
                          "mean.sd.shih","pval.sd.shih","rmst.shih",
                          "pval.rmst.shih","hr.shih","ht.pval.shih",
                          "mean.sd.strat","pval.sd.strat","rmst.strat",
                          "pval.rmst.strat","hr.strat","ht.pval.strat")

    for (k in 1:reps){
      result[k,] <- clin.util.comp(n=smple.size,
                         accru=accru,
                         ltfu=LTFU,
                         pos=pos[i],
                         pc.pos=PC,
                         pc.neg=PC,
                         medposB=Bpos[j],
                         medposA=9,
                         mednegB=9,
                         mednegA=12,
                         pseudo.diff=12)
    }
    
  


  filename <- paste("comp", "pc1.", "n", smple.size, "accru", accru, "exp", 
                           "Bpos",Bpos[j], "Apos", "9", "Bneg","9", "Aneg",  "12", 
                           "M", pos[i], "LTFU", LTFU, "rds", sep=".")


  saveRDS(result, file=paste("Results", "comp",filename, sep="/"))
  }
}