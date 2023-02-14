#==============================================================================
# FILENAME: analysis.R
# PROJECT: 	Prediction Driven Trial Designs, Phase III
# PURPOSE: calls the analysis.func.R functions to produce the resulting figures and table
#          
# AUTHOR: Adam Brand


# OUTPUT: 
#        

# R VERSION: 4.0.5
#==============================================================================
#Notes: 
### 




# =============================================================================
## sourcing the source program
source("Programs/analysis.func.R")


#### this generates the simulation results table in the paper

comp9.25 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")
comp9.5 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.9.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rds")
comp12.25 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.12.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")
comp12.5 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.12.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rds")
comp21.25 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.25.LTFU.0.02.rds")
comp21.5 <- readRDS("Results/comp/comp.pc1..n.400.accru.10.exp.Bpos.21.Apos.9.Bneg.9.Aneg.12.M.0.5.LTFU.0.02.rds")

### getting the output for table 2 in the paper
tab2.r1 <- round(abs(analyze.comp(comp9.25)), digits=4)
tab2.r2 <- round(abs(analyze.comp(comp9.5)), digits=4)
tab2.r3 <- round(abs(analyze.comp(comp12.25)), digits=4)
tab2.r4 <- round(abs(analyze.comp(comp12.5)), digits=4)
tab2.r5 <- round(abs(analyze.comp(comp21.25)), digits=4)
tab2.r6 <- round(abs(analyze.comp(comp21.5)), digits=4)

### creating table 2
table2 <- matrix(nrow=6, ncol=16)
table2[,1] <- c(9,9,12,12,21,21)
table2[,2] <- c(0.25,0.50)

for (i in 3:16){
  table2[,i] <- c(tab2.r1[(i-2)],tab2.r2[(i-2)],tab2.r3[(i-2)],tab2.r4[(i-2)],tab2.r5[(i-2)],tab2.r6[(i-2)])
}

table2 <- data.frame(table2)

colnames(table2) <- c("Median","Mark.pos","LR.exp","LR.comp","HR.rej.exp","HR.rej.comp","HR.mean.exp",
                      "HR.mean.comp","rmst.rej.exp","rmst.rej.comp","rmst.mean.exp","rmst.mean.comp",
                      "sd.rej.exp","sd.rej.comp","sd.mean.exp","sd.mean.comp")

print(table2)

##### NOTE, THE TABLE OUTPUTS MORE DECIMALS THAN IN THE TABLE IN THE PAPER. THIS IS BECAUSE THE ROUNDING
##### FUNCTION IN R DOES NOT ROUND A 5 DECIMAL CONSISTENTLY. THE TABLE IN THE PAPER IS ROUNDED MANUALLY.


####### BELOW HERE IS ANALYSIS IN THE SUPPLEMENT


## creating the figure for the power curves for the differential treatment effect estimand


fig9.25.5_2 <- power.curve(estimand="inter",
                         design.type="stratify",
                         smple.size= seq(from=250, to=2000, by=150),
                         accru=10,
                         LTFU=.02,
                         medianBpos=9,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.25,
                         phys.choice=.25,
                         rpos=0.5            # randomization to treatment b in the positives
)


fig12.25.5 <- power.curve(estimand="inter",
                          design.type="stratify",
                          smple.size= seq(from=250, to=1000, by=50),
                          accru=10,
                          LTFU=.02,
                          medianBpos=12,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.5            # randomization to treatment b in the positives
)
fig15.25.5 <- power.curve(estimand="inter",
                          design.type="stratify",
                          smple.size= seq(from=250, to=1000, by=50),
                          accru=10,
                          LTFU=.02,
                          medianBpos=15,
                          medianBneg=9,
                          medianApos=9,
                          medianAneg=12,
                          marker.pos=.25,
                          phys.choice=.25,
                          rpos=0.5            # randomization to treatment b in the positives
)




fig9.5.5_2 <- power.curve(estimand="inter",
                        design.type="stratify",
                        smple.size= seq(from=250, to=2000, by=150),
                        accru=10,
                        LTFU=.02,
                        medianBpos=9,
                        medianBneg=9,
                        medianApos=9,
                        medianAneg=12,
                        marker.pos=.5,
                        phys.choice=.25,
                        rpos=0.5            # randomization to treatment b in the positives
)

fig12.5.5 <- power.curve(estimand="inter",
                         design.type="stratify",
                         smple.size= seq(from=250, to=1000, by=50),
                         accru=10,
                         LTFU=.02,
                         medianBpos=12,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.5,
                         phys.choice=.25,
                         rpos=0.5            # randomization to treatment b in the positives
)
fig15.5.5 <- power.curve(estimand="inter",
                         design.type="stratify",
                         smple.size= seq(from=250, to=1000, by=50),
                         accru=10,
                         LTFU=.02,
                         medianBpos=15,
                         medianBneg=9,
                         medianApos=9,
                         medianAneg=12,
                         marker.pos=.5,
                         phys.choice=.25,
                         rpos=0.5            # randomization to treatment b in the positives
)


diff.fig <- ggarrange(fig9.25.5_2, fig12.25.5, fig15.25.5,
                  fig9.5.5_2, fig12.5.5, fig15.5.5, 
                  ncol=3, nrow=2)

annotate_figure(diff.fig, 
                top=text_grob("Median Survival for Positive Patients on Treatment B
                    9 months                                               12 months                                             15 months          ", size=12),
                left=text_grob("Proportion of Biomarker Positive Patients
                        0.5                                          0.25", size=12, rot=90)
                 ,bottom=text_grob("Figure S2: Power of the Biomarker Stratified Design 
 to Detect Differential Treatment Effect by Number of Events: 
 Red=HR, Blue=RMST, Green=SD", size=15)
                )


#### creating the figure for single subgroup effect

fig12.5 <- power.curve(estimand="subgrp",
                      design.type="enrich",
                      smple.size= seq(from=50, to=800, by=50),
                      accru=10,
                      LTFU=.02,
                      medianBpos=12,
                      medianBneg=9,
                      medianApos=9,
                      medianAneg=12,
                      marker.pos=.9,
                      phys.choice=.25,
                      rpos=0.67            # randomization to treatment b in the positives
)
fig15.5 <- power.curve(estimand="subgrp",
                      design.type="enrich",
                      smple.size= seq(from=50, to=800, by=50),
                      accru=10,
                      LTFU=.02,
                      medianBpos=15,
                      medianBneg=9,
                      medianApos=9,
                      medianAneg=12,
                      marker.pos=.9,
                      phys.choice=.25,
                      rpos=0.67            # randomization to treatment b in the positives
)

subgrp.fig <- ggarrange(fig12.5, fig15.5,
                      ncol=2, nrow=1)

annotate_figure(subgrp.fig, 
                top=text_grob("Median Survival for Positive Patients on Treatment B
                          12 months                                                     15 months                  ", size=12)
                 , bottom=text_grob("Figure S1: Power of the Enrichment Design 
 to Detect Treatment Effect for Biomarker 
 Positive Patients by Number of Events: 
 Red=HR, Blue=SD, Green=LR, Purple=RMST", size=15)
                )

### creating the power figure for clinical utility


fig9.5.25.clin <- power.curve(estimand="clin",
                             design.type="strategy",
                             smple.size= seq(from=500, to=3000, by=250),
                             accru=10,
                             LTFU=.02,
                             medianBpos=9,
                             medianBneg=9,
                             medianApos=9,
                             medianAneg=12,
                             marker.pos=.5,
                             phys.choice=0.25,
                             rpos=0.67            # randomization to treatment b in the positives
)
fig12.5.25.clin <- power.curve(estimand="clin",
                              design.type="strategy",
                              smple.size= seq(from=500, to=3000, by=250),
                              accru=10,
                              LTFU=.02,
                              medianBpos=12,
                              medianBneg=9,
                              medianApos=9,
                              medianAneg=12,
                              marker.pos=.5,
                              phys.choice=0.25,
                              rpos=0.67            # randomization to treatment b in the positives
)
fig21.5.25.clin <- power.curve(estimand="clin",
                              design.type="strategy",
                              smple.size= seq(from=500, to=3000, by=250),
                              accru=10,
                              LTFU=.02,
                              medianBpos=21,
                              medianBneg=9,
                              medianApos=9,
                              medianAneg=12,
                              marker.pos=.5,
                              phys.choice=0.25,
                              rpos=0.67            # randomization to treatment b in the positives
)
fig9.5.5.clin <- power.curve(estimand="clin",
                              design.type="strategy",
                              smple.size= seq(from=500, to=3000, by=250),
                              accru=10,
                              LTFU=.02,
                              medianBpos=9,
                              medianBneg=9,
                              medianApos=9,
                              medianAneg=12,
                              marker.pos=.5,
                              phys.choice=0.5,
                              rpos=0.67            # randomization to treatment b in the positives
)
fig12.5.5.clin <- power.curve(estimand="clin",
                             design.type="strategy",
                             smple.size= seq(from=500, to=3000, by=250),
                             accru=10,
                             LTFU=.02,
                             medianBpos=12,
                             medianBneg=9,
                             medianApos=9,
                             medianAneg=12,
                             marker.pos=.5,
                             phys.choice=0.5,
                             rpos=0.67            # randomization to treatment b in the positives
)
fig21.5.5.clin <- power.curve(estimand="clin",
                             design.type="strategy",
                             smple.size= seq(from=500, to=3000, by=250),
                             accru=10,
                             LTFU=.02,
                             medianBpos=21,
                             medianBneg=9,
                             medianApos=9,
                             medianAneg=12,
                             marker.pos=.5,
                             phys.choice=0.5,
                             rpos=0.67            # randomization to treatment b in the positives
)
fig9.5.75.clin <- power.curve(estimand="clin",
                             design.type="strategy",
                             smple.size= seq(from=500, to=3000, by=250),
                             accru=10,
                             LTFU=.02,
                             medianBpos=9,
                             medianBneg=9,
                             medianApos=9,
                             medianAneg=12,
                             marker.pos=.5,
                             phys.choice=0.75,
                             rpos=0.67            # randomization to treatment b in the positives
)
fig12.5.75.clin <- power.curve(estimand="clin",
                              design.type="strategy",
                              smple.size= seq(from=500, to=3000, by=250),
                              accru=10,
                              LTFU=.02,
                              medianBpos=12,
                              medianBneg=9,
                              medianApos=9,
                              medianAneg=12,
                              marker.pos=.5,
                              phys.choice=0.75,
                              rpos=0.67            # randomization to treatment b in the positives
)
fig21.5.75.clin <- power.curve(estimand="clin",
                              design.type="strategy",
                              smple.size= seq(from=500, to=3000, by=250),
                              accru=10,
                              LTFU=.02,
                              medianBpos=21,
                              medianBneg=9,
                              medianApos=9,
                              medianAneg=12,
                              marker.pos=.5,
                              phys.choice=0.75,
                              rpos=0.67            # randomization to treatment b in the positives
)

clin.fig <- ggarrange( fig9.5.25.clin, fig9.5.5.clin, fig9.5.75.clin,
                       fig12.5.25.clin, fig12.5.5.clin, fig12.5.75.clin,
                       fig21.5.25.clin, fig21.5.5.clin, fig21.5.75.clin,
                      ncol=3, nrow=3)

annotate_figure(clin.fig, 
                top=text_grob("Proportion of Physician's Choice Equalling Biomarker-directed Strategy
         0.25                                                      0.5                                                     0.75", size=12),
                left=text_grob("Median Survival in Positive Patients on Treatment B
                      21 months                       12 months                         9 months", size=12, rot=90)
                 , bottom=text_grob("Figure S3: Power of the Strategy Design 
 to Detect Clinical Utility for by Number of Events: 
 Red=HR, Blue=SD, Green=LR, Purple=RMST", size=15)
                )

#### creating individual summary tables for each scenario


##### differential treatment effect scenarios

# type 1 error rate tables
fintable(caption="Differential Treatment Effect Type 1 Error Rate: 
         Median Survival in M=1 Patients and X=B is 9 months, 
          Proportion of Positive Patients in 0.25, and Sample size is 2000",
         estimand="inter",
         design.type="stratify",
         smple.size= 2000,
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=12,
         medianAneg=12,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Differential Treatment Effect Type 1 Error Rate: 
         Median Survival in M=1 Patients and X=B is 9 months, 
          Proportion of Positive Patients in 0.5, and Sample size is 2000",
         estimand="inter",
         design.type="stratify",
         smple.size= 2000,
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=12,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Differential Treatment Effect Type 1 Error Rate: 
         Median Survival in M=1 Patients and X=B is 15 months, 
          Proportion of Positive Patients in 0.5, and Sample size is 300",
         estimand="inter",
         design.type="stratify",
         smple.size= 300,
         accru=10,
         LTFU=.02,
         medianBpos=15,
         medianBneg=15,
         medianApos=12,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Differential Treatment Effect Type 1 Error Rate: 
         Median Survival in M=1 Patients and X=B is 15 months, 
          Proportion of Positive Patients in 0.25, and Sample size is 400",
         estimand="inter",
         design.type="stratify",
         smple.size= 400,
         accru=10,
         LTFU=.02,
         medianBpos=15,
         medianBneg=15,
         medianApos=12,
         medianAneg=12,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Differential Treatment Effect Type 1 Error Rate: 
         Median Survival in M=1 Patients and X=B is 12 months, 
          Proportion of Positive Patients in 0.25, and Sample size is 700",
         estimand="inter",
         design.type="stratify",
         smple.size= 700,
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=12,
         medianApos=9,
         medianAneg=9,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Differential Treatment Effect Type 1 Error Rate: 
         Median Survival in M=1 Patients and X=B is 12 months, 
          Proportion of Positive Patients in 0.5, and Sample size is 500",
         estimand="inter",
         design.type="stratify",
         smple.size= 500,
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=12,
         medianApos=9,
         medianAneg=9,
         marker.pos=.5,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

## operating characteristic tables
fintable(caption="Differential Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 9, Proportion of 
         Positive Patients in 0.25, and Probability of randomizing Positive Patients to 
         Treatment B is 0.5",
         estimand="inter",
         design.type="stratify",
         smple.size= seq(from=250, to=2000, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5)


fintable(caption="Differential Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 12, Proportion of 
         Positive Patients in 0.25, and Probability of randomizing Positive Patients to 
         Treatment B is 0.5",
         estimand="inter",
         design.type="stratify",
         smple.size= seq(from=250, to=1000, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5)
fintable(caption="Differential Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 15, Proportion of 
         Positive Patients in 0.25, and Probability of randomizing Positive Patients to 
         Treatment B is 0.5",
         estimand="inter",
         design.type="stratify",
         smple.size= seq(from=250, to=1000, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=15,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.25,
         phys.choice=.25,
         rpos=0.5)

fintable(caption="Differential Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 9, Proportion of 
         Positive Patients in 0.5, and Probability of randomizing Positive Patients to 
         Treatment B is 0.5",
         estimand="inter",
         design.type="stratify",
         smple.size= seq(from=250, to=2000, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=.25,
         rpos=0.5)


fintable(caption="Differential Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 12, Proportion of 
         Positive Patients in 0.5, and Probability of randomizing Positive Patients to 
         Treatment B is 0.5",
         estimand="inter",
         design.type="stratify",
         smple.size= seq(from=250, to=1000, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=.25,
         rpos=0.5)
fintable(caption="Differential Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 15, Proportion of 
         Positive Patients in 0.5, and Probability of randomizing Positive Patients to 
         Treatment B is 0.5",
         estimand="inter",
         design.type="stratify",
         smple.size= seq(from=250, to=1000, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=15,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=.25,
         rpos=0.5)


##### clinical utility scenarios

### type 1 error rate tables

fintable(caption="Clinical Utility Type 1 Error Rate: 
         Median Survival in M=1, X=B Patients is 9 months, and
         Sample Size is 3000",
         estimand="clin",
         design.type="strategy",
         smple.size= 3000,
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=1,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Clinical Utility Type 1 Error Rate: 
         Median Survival in M=1, X=B Patients is 12 months and
         Sample Size is 1000",
         estimand="clin",
         design.type="strategy",
         smple.size= 1000,
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=1,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Clinical Utility Type 1 Error Rate: 
         Median Survival in M=1, X=B Patients is 12 months and
         Sample Size is 2000",
         estimand="clin",
         design.type="strategy",
         smple.size= 2000,
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=1,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Clinical Utility Type 1 Error Rate: 
         Median Survival in M=1, X=B Patients is 21 months and
         Sample Size is 500",
         estimand="clin",
         design.type="strategy",
         smple.size= 500,
         accru=10,
         LTFU=.02,
         medianBpos=21,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=1,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Clinical Utility Type 1 Error Rate: 
         Median Survival in M=1, X=B Patients is 21 months and
         Sample Size is 750",
         estimand="clin",
         design.type="strategy",
         smple.size= 750,
         accru=10,
         LTFU=.02,
         medianBpos=21,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=1,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Clinical Utility Type 1 Error Rate: 
         Median Survival in M=1, X=B Patients is 21 months and
         Sample Size is 3000",
         estimand="clin",
         design.type="strategy",
         smple.size= 3000,
         accru=10,
         LTFU=.02,
         medianBpos=21,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=1,
         rpos=0.5,
         type1=TRUE)

## operating characteristic tables


fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 9, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.25",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.25,
         rpos=0.67)
fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 9, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.5",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.5,
         rpos=0.67)
fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 9, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.75",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=9,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.75,
         rpos=0.67)


fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 12, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.25",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.25,
         rpos=0.67)
fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 12, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.5",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.5,
         rpos=0.67)
fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 12, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.75",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.75,
         rpos=0.67)

fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 21, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.25",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=21,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.25,
         rpos=0.67)
fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 21, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.5",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=21,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.5,
         rpos=0.67)
fintable(caption="Clinical Utility Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 21, Proportion of 
         Positive Patients in 0.5, and Proportion Biomarker-directed Treatment Strategy Equalling
         Physician's Choice of Treatment is 0.75",
         estimand="clin",
         design.type="strategy",
         smple.size= seq(from=500, to=3000, by=250),
         accru=10,
         LTFU=.02,
         medianBpos=21,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.5,
         phys.choice=0.75,
         rpos=0.67)

### creating summary tables for the single subgroup results

fintable(caption="Single Subgroup Treatment Effect Type 1 Error Rate:
         Median Survival in M=1, X=B is 12 months and Sample size is 550",
         estimand="subgrp",
         design.type="enrich",
         smple.size= 550,
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=12,
         medianApos=12,
         medianAneg=9,
         marker.pos=.75,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

fintable(caption="Single Subgroup Treatment Effect Type 1 Error Rate:
         Median Survival in M=1, X=B is 15 months and Sample size is 150",
         estimand="subgrp",
         design.type="enrich",
         smple.size= 150,
         accru=10,
         LTFU=.02,
         medianBpos=15,
         medianBneg=12,
         medianApos=15,
         medianAneg=9,
         marker.pos=.75,
         phys.choice=.25,
         rpos=0.5,
         type1=TRUE)

### getting the operating characteristics summaries

fintable(caption="Single Subgroup Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 12 and Proportion of 
         Positive Patients is 0.25",
         estimand="subgrp",
         design.type="enrich",
         smple.size= seq(from=50, to=800, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=12,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.9,
         phys.choice=0.5,
         rpos=0.5)

fintable(caption="Single Subgroup Treatment Effect Operating Characteristics Summary: 
         Median Survival in Positive Patients on Treatment B is 15 and Proportion of 
         Positive Patients is 0.25",
         estimand="subgrp",
         design.type="enrich",
         smple.size= seq(from=50, to=800, by=50),
         accru=10,
         LTFU=.02,
         medianBpos=15,
         medianBneg=9,
         medianApos=9,
         medianAneg=12,
         marker.pos=.9,
         phys.choice=0.5,
         rpos=0.5)



