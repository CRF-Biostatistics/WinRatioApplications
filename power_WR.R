# Requires functions in 'win ratio functions.R', 'simulateJF.R', and the
# "binom" package
# install.packages("binom", dependencies = TRUE)
require("binom")

###############################################################################

# The following functions estimate the statistical power for a hierachical
# composite endpoint under the joint frailty model to be analyzed using the Win
# Ratio. There is no closed solution so the function iteratively calculates the
# power using Monte Carlo simulation. A mathematical formulation is described
# in the Statistical Appendix in the 'The Win Ratio Approach for Composite
# Endpoints: Practical Guidance Based on Previous Experience' manuscript by
# Redfors et. al.

###############################################################################

## Analyzes a 3-level hierarchical composite endpoint by the Win Ratio
evaluate_jf <-
  # jf_data = a 'frailtySim' object (see 'simulateJF' above)
  # criteria = a numeric vector of the number of levels in the endpoint
  # KCCQ_control = vector of KCCQ scores for control patients
  # KCCQ_treatment = vector of KCCQ scores for treatment patients
  function(jf_data,
           criteria,
           KCCQ_control = NULL,
           KCCQ_treatment = NULL) {
    fu_vector1 <-
      slot(jf_data, "followup")[slot(jf_data, "treatment") == 0]
    fu_vector2 <-
      slot(jf_data, "followup")[slot(jf_data, "treatment") == 1]
    
    t_death1 <-
      slot(jf_data, "t_death")[slot(jf_data, "treatment") == 0]
    t_death2 <-
      slot(jf_data, "t_death")[slot(jf_data, "treatment") == 1]
    
    winner_death <-
      win_matrix_firstEvent(t_death1, t_death2, fu_vector1, fu_vector2)
    
    if (length(criteria) == 1) {
      if (criteria == "death") {
        winner_list <- list()
        winner_list[[1]] <- winner_death
      }
    }
    
    winner_list <- list()
    
    if (length(criteria) == 2 | length(criteria) == 3) {
      if (criteria[1] == "death" & criteria[2] == "firstHFH") {
        t_hfh1 <- slot(jf_data, "t_hfh")[slot(jf_data, "treatment") == 0, 1]
        t_hfh2 <-
          slot(jf_data, "t_hfh")[slot(jf_data, "treatment") == 1, 1]
        winner_hfh <-
          win_matrix_firstEvent(t_hfh1, t_hfh2, fu_vector1, fu_vector2)
      }
      
      if (criteria[1] == "death" & criteria[2] == "nHFH") {
        t_hfh1 <- slot(jf_data, "t_hfh")[slot(jf_data, "treatment") == 0, ]
        t_hfh2 <-
          slot(jf_data, "t_hfh")[slot(jf_data, "treatment") == 1, ]
        winner_hfh <-
          winner_nEvents(t_hfh1, t_hfh2, fu_vector1, fu_vector2)
      }
      winner_list[[1]] <- winner_death
      winner_list[[2]] <- winner_hfh
    }
    
    if (length(criteria) == 3) {
      winner_KCCQ <-
        win_matrix_continuous(vector1 = KCCQ_control, vector2 = KCCQ_treatment)
      winner_list[[3]] <- winner_KCCQ
    }
    
    return(winner_list)
    
  }

## Create a class to store the win ratio power calculation results
setClass(
  "WRPower",
  representation(
    required_n = "numeric",
    meanWR = "numeric",
    perc_decisions = "vector",
    logWR = "vector",
    logWR_se = "vector",
    fu.start = "vector",
    power = "numeric",
    n.sim = "numeric",
    n.iter = "numeric",
    hr.death = "numeric",
    hr.hfh = "numeric",
    p.active = "numeric",
    rate.control.d = "numeric",
    rate.control.hfh = "numeric",
    frailty.scale = "numeric",
    frailty.exponent = "numeric",
    criteria = "vector",
    cont_parameters = "vector"
  )
)

## Main function to perform power calculations for the Win Ratio by performing
# Monte Carlo simulations from joint frailty data
WRPower <-
  # n.start = number of patients in the initial simulation
  # fu.start = numeric vector of follow-up times for each patient in the initial simulation
  # power = target power desired
  # alpha = type I error rate as a numeric decimal (i.e. 0.05 for 5%)
  # n.sim = number of simulations to perform for each iteration
  # n.iter = number of iterations
  # p.active = proportion of treatment vs. control patients (entered as a decimal; see stats::rbinom)
  # hr.death = hazard ratio between treatment vs. control for death
  # hr.hfh = hazard ratio between treatment vs. control for recurrent heart failure hospitalizations
  # rate.control.d = control hazard rate for death
  # rate.control.hfh = control hazard rate for heart failure hospitalizations
# frailty.scale = scale parameter for the joint frailty link
# frailty.exponent = exponent of the joint frailty link
# criteria = a numeric vector of the number of levels in the endpoint
# cont_parms =  optional vector of length 4, containing the mean and
#               standard deviation for the control and treatment patients,
#               respectively
function(n.start = NULL,
         fu.start,
         power = 0.8,
         alpha = 0.05,
         n.sim = 1000,
         n.iter = 3,
         p.active = 0.5,
         hr.death,
         hr.hfh,
         rate.control.d,
         rate.control.hfh,
         frailty.scale = 0,
         frailty.exponent = 0,
         criteria,
         cont_params = "vector") {
  #Choosing a sensible starting place number of patients that is reasonably quick to run (if none is specified)
  if (is.null(n.start)) {
    n.start <- 500
  }
  
  #Generating z-statistics from supplied alpha and power
  z_power <- qnorm(power)
  z_alpha_desired <- qnorm(1 - alpha / 2)
  
  #Setting starting number of observations, and censoring distribution for first iteration
  n.obs <- n.start
  length.fu <- fu.start
  
  #Looping through iterative process to determine required sample size
  for (iter in 1:n.iter) {
    print(paste("Iteration", iter))
    
    logWR <- vector(length = n.sim)
    logWR_se <- vector(length = n.sim)
    is.sig <- vector(length = n.sim)
    
    #Looping through n.sim simulated datasets within each iteration
    for (sim in 1:n.sim) {
      winner_list <- NULL
      
      #Simulating data in each iteration and posting results to vectors
      jf_data <-
        simulateJF(
          n.obs = n.obs,
          p.active = p.active,
          hr.death = hr.death,
          hr.hfh = hr.hfh,
          rate.control.d = rate.control.d,
          rate.control.hfh = rate.control.hfh,
          length.fu = length.fu,
          frailty.scale = frailty.scale,
          frailty.exponent = frailty.exponent
        )
      
      #If there are three tiers then assuming KCCQ data is supplied (will update to generalize later)
      if (length(criteria) == 3) {
        if (length(cont_params) != 4) {
          print("Parameters for KCCQ not correctly specified")
        }
        n_obs1 <-
          length(slot(jf_data, "treatment")[slot(jf_data, "treatment") == 0])
        n_obs2 <-
          length(slot(jf_data, "treatment")[slot(jf_data, "treatment") == 1])
        KCCQ1 <-
          rnorm(n_obs1, mean = cont_params[1], sd = cont_params[2])
        KCCQ2 <-
          rnorm(n_obs2, mean = cont_params[3], sd = cont_params[4])
        winner_list <-
          evaluate_jf(
            jf_data,
            criteria = criteria,
            KCCQ_control = KCCQ1,
            KCCQ_treatment = KCCQ2
          )
      } else {
        winner_list <- evaluate_jf(jf_data, criteria = criteria)
      }
      
      #Calculating win ratio statistics and storing them
      wrstats <- getWRStats(winner_list, alpha)
      logWR[sim] <- slot(wrstats, "logWR")
      logWR_se[sim] <- slot(wrstats, "logWR_se")
      is.sig[sim] <- slot(wrstats, "is.sig")
      
      #Displaying progress of simulations
      if (sim %% 10 == 0) {
        print(paste("Simulation", sim))
      }
    }
    
    #size is proportional to the square of [z(alpha/2) + z(beta)]-adapting sample size on this basis
    z <- logWR / logWR_se
    beta <- 1 - power
    #Idenfitying which alpha, based upon the current sample size, would lead to the required level of statistical power
    z_alpha_obs <- quantile(abs(z), beta)
    
    #Estimating required N based on required power, observed alpha above, and formula above
    required_n <-
      n.obs * (z_alpha_desired + z_power) ^ 2 / (z_power + z_alpha_obs) ^ 2
    required_n <- round(required_n)
    
    #Displaying adaptation of sample size
    current_power <- mean(is.sig)
    current_z <- qnorm(current_power)
    print(
      paste(
        "Iteration",
        iter,
        "Power:" ,
        current_power,
        "N",
        n.obs,
        "Estimated required:",
        required_n
      )
    )
    
    #Redefining inputs for joint frailty simulations based on updated sapmle size
    n.obs_old <- n.obs
    n.obs <- required_n
    length.fu <- sample(fu.start, size = n.obs, replace = TRUE)
    
    #Simulations become very slow once N gets too large
    if (required_n > 10000) {
      stop("Required sample size too large to run simulations")
    }
    
    #Displaying final estimated sample size and results from final iteration
    if (iter == n.iter) {
      power_ci <- binom.confint(sum(is.sig), n.sim, method = "exact")
      power_ci <-
        binom.confint(sum(is.sig), n.sim, method = "exact")
      print(paste(
        "95% CI for power with n=",
        n.obs_old ,
        ":",
        power_ci$lower,
        power_ci$upper
      ))
      print(paste("Refined estimate (not simulated)", n.obs))
    }
    
  }
  
  #Preparing output
  mean_WR <- exp(mean(logWR))
  #Percentage decisions, just being lazy and using last simulation
  perc_decisions <- slot(wrstats, "perc_decisions")
  
  #Creating object that contains the results and the input
  output <-
    new(
      "WRPower",
      required_n = required_n,
      meanWR = mean_WR,
      perc_decisions = perc_decisions,
      logWR = logWR,
      logWR_se = logWR_se,
      p.active = p.active,
      hr.death = hr.death,
      hr.hfh = hr.hfh,
      rate.control.d = rate.control.d,
      rate.control.hfh = rate.control.hfh,
      frailty.scale = frailty.scale,
      frailty.exponent = frailty.exponent,
      criteria = criteria,
      cont_parameters = cont_params,
      fu.start = fu.start,
      power = power,
      n.sim = n.sim,
      n.iter = n.iter
    )
  
}