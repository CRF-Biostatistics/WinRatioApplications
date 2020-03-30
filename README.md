# WinRatioApplications
R code to reproduce sample size calculations in the "The Win Ratio Approach for Composite Endpoints: Practical Guidance Based on Previous Experience" manuscript by Redfors et. al.

The R scripts 'win ratio functions.R', 'simulateJF.R', 'power_WR.R', and 'power_CoxFrailty.R' contain R functions that perform the Win Ratio test, simulate terminal and recurrent event data under the joint frailty model, and calculate sample sizes for joint frailty data using the Win Ratio and Cox proportional hazards regression model.

Additionally, the R scripts '01-runSim1.R', '02-runSim2.R', and '03-runSim3.R' are examples which reproduce the results summarized in Table 4 of the manuscript referenced above. Please note that these simulations were run in batch mode and may take up to 44 hours to finish. Try a smaller number of simulations (n.sims) to get started.
