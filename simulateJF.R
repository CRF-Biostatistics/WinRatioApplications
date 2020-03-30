# Requires functions in 'win ratio functions.R'

###############################################################################

# The following functions simulate patient level data according to the joint
# frailty model for death and heart failure hospitalizations. Parameterization
# follows the "Joint frailty models for recurring events and death using
# maximum penalized likelihood estimation: application on cancer events"
# paper by Rondeau et. al.

###############################################################################

## Returns the number of finite observations in a vector 
length_notInf <-
  # vector = a numeric vector
  function(vector) {
    output <- length(vector[vector != Inf])
    return(output)
  }

## Create a class to store the simulated terminal and recurrent event data
setClass(
  "frailtySim",
  representation(
    indiv_hfh = "vector",
    death = "vector",
    t_death = "vector",
    t_hfh = "matrix",
    treatment = "vector",
    followup = "vector",
    n = "numeric",
    n_hfh = "numeric",
    n_death = "numeric",
    timeAtRisk = "numeric"
  )
)

## Main function to simulate data with a joint frailty link between non-fatal and fatal events
simulateJF <-
  # n.obs = number of patients
  # p.active = proportion of treatment vs. control patients (entered as a decimal; see stats::rbinom)
  # hr.death = hazard ratio between treatment vs. control for death
  # hr.hfh = hazard ratio between treatment vs. control for recurrent heart failure hospitalizations
  # rate.control.d = control hazard rate for death
  # rate.control.hfh = control hazard rate for heart failure hospitalizations
  # length.fu = numeric vector of follow-up times for each patient
  # frailty.scale = scale parameter for the joint frailty link
  # frailty.exponent = exponent of the joint frailty link
  function(n.obs,
           p.active,
           hr.death,
           hr.hfh,
           rate.control.d,
           rate.control.hfh,
           length.fu,
           frailty.scale,
           frailty.exponent) {
    #Checking sanity for censoring distribution and frailty distribution
    if (length(length.fu) != 1 & length(length.fu) != n.obs) {
      print("Breaking function, #items in censoring distribution neither 1 nor N")
      break
    }
    
    if (frailty.scale < 0) {
      print("Frailty scale must be non-negative")
      break
    }
    
    #If frailty scale set to 0, then all patients have the same risk
    if (frailty.scale == 0) {
      frailty.hfh <- 1
      frailty.d <- 1
    }
    
    #Creating frailty term. Link is the same as in a standard joint frailty model
    #Can set frailty.hfh to 1 for no variation in patient frailty
    #Shape=1/scale in order to maintain mean frailty=1
    if (frailty.scale != 0) {
      frailty.shape <- 1 / frailty.scale
      frailty.hfh <-
        rgamma(n = n.obs,
               scale = frailty.scale,
               shape = frailty.shape)
      frailty.d <- frailty.hfh ^ frailty.exponent
    }
    
    #Creating active treatment rates based on simulation parameters
    rate.active.d <- rate.control.d * hr.death
    rate.active.hfh <- rate.control.hfh * hr.hfh
    
    #Randomly generating treatment
    treatment <- rbinom(n = n.obs, size = 1, prob = p.active)
    
    #Generating rates dependent upon treatment
    #First splitting into active/control
    #Then multiplying by frailty term
    rate.d <- vector(length = n.obs)
    rate.d[treatment == 1] <- rate.active.d
    rate.d[treatment == 0] <- rate.control.d
    indiv.rate.d <- rate.d * frailty.d
    
    rate.hfh <- vector(length = n.obs)
    rate.hfh[treatment == 1] <- rate.active.hfh
    rate.hfh[treatment == 0] <- rate.control.hfh
    indiv.rate.hfh <- rate.hfh * frailty.hfh
    
    #Creating time to death
    t.death.matrix <- rexp(n.obs, rate = indiv.rate.d)
    
    #Creating first one column of HFH data
    t.hfh.matrix <- matrix(nrow = n.obs, ncol = 1)
    t.hfh.matrix[, 1] <- rexp(n.obs, rate = indiv.rate.hfh)
    #Sum HFH is used to define the timing of the last HFH
    sum_hfh <- apply(t.hfh.matrix, 1, sum)
    
    #Keep adding more columns the earliest HFH is after the end of follow up
    while (sum(sum_hfh > length.fu) < n.obs) {
      next.hfh.col <- rexp(n.obs, rate = indiv.rate.hfh)
      t.hfh.matrix <- cbind(t.hfh.matrix, next.hfh.col)
      sum_hfh <- apply(t.hfh.matrix, 1, sum)
    }
    t.hfh.matrix <- t(apply(t.hfh.matrix, 1, cumsum))
    
    #Censoring death end of follow up
    t.death.matrix[t.death.matrix > length.fu] <- Inf
    death <- rep(0, times = n.obs)
    death[t.death.matrix < length.fu] <- 1
    
    #Censoring HFH at death or end of follow up
    t.hfh.matrix[t.hfh.matrix > length.fu] <- Inf
    t.hfh.matrix[t.hfh.matrix > t.death.matrix] <- Inf
    
    #Creating a follow up time based on death and administrative censoring at end of study
    follow.up <- cbind(length.fu, t.death.matrix)
    follow.up <- apply(follow.up, 1, min)
    
    #Calculating some summary statistics
    n_hfh <- length(t.hfh.matrix[t.hfh.matrix != Inf])
    n_death <- length(t.death.matrix[t.death.matrix != Inf])
    n_patients_hfh <-
      length(t.hfh.matrix[, 1][t.hfh.matrix[, 1] != Inf])
    indiv_hfh <- apply(t.hfh.matrix, 1, length_notInf)
    
    #Outputting information
    output <-
      new(
        "frailtySim",
        t_death = t.death.matrix,
        death = death,
        t_hfh = t.hfh.matrix,
        treatment = treatment,
        n = length(treatment),
        followup = follow.up,
        n_hfh = n_hfh,
        n_death = n_death,
        timeAtRisk = sum(pmin(length.fu, t.death.matrix)),
        indiv_hfh = indiv_hfh
      )
    return(output)
    
  }

## Transforms result from the 'simulateJF' function (a 'frailtySim' object) into
## an analyzable format for functions in the 'frailtypack' package
## see frailtypack::frailtyPenal
transform_jf <-
  function(jf_data) {
    # jf_data = a 'frailtySim' object (see 'simulateJF' above)
    n.obs <- slot(jf_data, "n")
    
    #Formatting data to be compliant with joint frailty programmes
    t_hfh <-
      data.frame(cbind(matrix(1:n.obs, nrow = n.obs), slot(jf_data, "t_hfh")))
    colnames(t_hfh)[2:3] <- c("next.hfh.col.0", "next_hfh.col.100")
    
    #Reshaping HFH data into long format
    hfh_long <-
      reshape(
        t_hfh,
        idvar = "subject",
        ids = t_hfh$V1,
        times = colnames(t_hfh)[2:dim(t_hfh)[2]],
        timevar = "rep",
        direction = "long",
        varying = list(colnames(t_hfh)[2:dim(t_hfh)[2]])
      )
    
    hfh_long <-
      data.frame(cbind(hfh_long$subject, hfh_long$next.hfh.col))
    names(hfh_long) <- c("subject", "time")
    
    #Taking only the HFHs (i.e. discaring empty rows where patient does not experience a HFH)
    hfh_long <- hfh_long[hfh_long$time != Inf, ]
    hfh_long$hfh <- 1
    hfh_long$death <- 0
    
    #Combining death data
    death <-
      data.frame(cbind(1:n.obs, slot(jf_data, "t_death"), rep(0, times = n.obs)))
    names(death) <- c("subject", "time", "death")
    death$death[death$time != Inf] <- 1
    death$time[death$time == Inf] <-
      slot(jf_data, "followup")[death$time == Inf]
    death$hfh <- 0
    
    #Putting HFH and death data together, then merging in treatment codes
    treatment <-
      data.frame(cbind(1:n.obs, slot(jf_data, "treatment")))
    names(treatment) <- c("subject", "treatment")
    jfSim <- rbind(hfh_long, death)
    jfSim <- merge(jfSim, treatment)
    
    #Ordering data by Subject ID and then by time
    jfSim <- jfSim[order(jfSim$subject, jfSim$time), ]
    
    #Where a patient has more than one row, adding in the start time in the second row, which is the time at which
    #the first event occurs. This process repeats for third, fourth row etc.
    jfSim$start <- lag(jfSim$time)
    jfSim$prevSubj <- lag(jfSim$subject)
    jfSim$start[jfSim$prevSubj != jfSim$subject] <- 0
    jfSim$start[1] <- 0
    
    #Returning formatted database
    return(jfSim)
  }
