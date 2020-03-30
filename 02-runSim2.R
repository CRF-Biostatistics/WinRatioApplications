# include utility functions
source("./win ratio functions.R")
source("./simulateJF.R")
source("./power_WR.R")
source("./power_CoxFrailty.R")

set.seed(153207293)

scenario2 <- WRPower(n.start = 1000,
                     fu.start = 1,
                     power = 0.80,
                     alpha = 0.05,
                     n.sim = 10000,
                     n.iter = 5,
                     p.active = 0.50,
                     hr.death = 0.80,
                     hr.hfh = 0.80,
                     rate.control.d = 0.25,
                     rate.control.hfh = 0.90,
                     frailty.scale = 1,
                     frailty.exponent = 1,
                     criteria = c("death", "nHFH"))

save(scenario2, file = "scenario2.Rdata")
