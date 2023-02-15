Prediction_Driven_Trials



This project contains all of the code to produce the results in the paper titled "Confirmatory Prediction-driven RCTs in Comparative Effectiveness Settings for Cancer Treatment" published in British Journal of Cancer.

Below are descriptions of programs and the files in the project.

In order to replicate the simulation results:


1. run analysis.R to produce Table 2 in the paper as well as the supplemental tables and figures in the supplemental materials attached to the paper. First run all of analysis.R. Then:

  - line 57 will print Table 2 in the paper in the R readout. This can then be manually input into Overleaf as Table 2, as was done for the paper.

  - lines 156-168 will print figure S2

  - lines 200-210 will print figure S1

  - lines 333-346 will print figure S3

  - to reproduce table S1

    - run lines 819-883. This will output Latex code for 4 tables. Copy and paste that code into an Overleaf page to create the 4 tables. You may need to make the Overleaf document landscape and the font of the tables small in order to view th full tables. Table S1 is a summary taking outputs from each of those tables. Rows 1 and 2, except for the T1E columns, are taken from Tables 3 and 4, respectively, where SD=Ratio. The T1E columns are taken from the 'Power' columns from Tables 1 and 2.

  - to produce table S2

    - run lines 354-578. Again copy the output into an Overleaf document, and adjust layout and font size. Table S2 gets the output from 7-12 for rows 1-6, respectively, except for the T1E columns, where SD=Ratio. The T1E columns are taken from the 'Power' columns in tables 1-6.

  - to produce table S3

   - run lines 585-797. Again copy the output into an Overleaf document, and adjust layout and font size. Table S3 gets the output from 7-12 for rows 1-6, respectively, except for the T1E columns, where SD=Ratio. The T1E columns are taken from the 'Power' columns in tables 1-6.

2. In order to produce the results datasets in the Results folder, run simulation.clin.R, simulation.inter.R, simulation.shih.R and simulation.subgrp.R. If running on a PC, these simulations will take multiple days each.
