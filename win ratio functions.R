# Requires the 'WWR' package
# install.packages("WWR", dependencies = TRUE)
require("WWR")

###############################################################################

# The following functions compare all pairs of treatment and control patients
# as described in "The win ratio: a new approach to the analysis of composite
# endpoints in clinical trials based on clinical priorities" by Pocock et. al.

# The result of 'win_matrix_firstEvent', 'win_matrix_continuous',
# 'winner_nEvents', and 'multi_criteria_winner' are a numeric matrix of win
# loss indicators. This is compatable with the WWR::genwr function.

# The 'getWRStats' consolidates a list of win loss matrices into one, where
# the 'winners' are picked hierarchically.

###############################################################################

## Perform pairwise comparisons for time-to-first event data
# Patient who is known to be event-free longer "wins"
win_matrix_firstEvent <-
  # t_vector1 = vector of binary event flags for control patients
  # t_vector2 = vector of binary event flags for treatment patients
  # fu_vector1 = vector of follow-up times for control patients
  # fu_vector2 = vector of follow-up times for treatment patients
  function(t_vector1,
           t_vector2,
           fu_vector1,
           fu_vector2) {
    #Collecting number of treatment and control patients and forming matrix
    n_obs1 <- length(t_vector1)
    n_obs2 <- length(t_vector2)
    
    win_matrix_firstEvent <- matrix(nrow = n_obs1, ncol = n_obs2)
    
    #Forming matrices of repeated values of time to first event for each group
    #They are set into matrices of the same dimension for easy direct comparison
    bigMatrix1 <-
      matrix(rep(t_vector1, times = n_obs2), ncol = n_obs2)
    bigMatrix2 <-
      matrix(rep(t_vector2, times = n_obs1),
             nrow = n_obs1,
             byrow = T)
    diffMatrix <- bigMatrix1 - bigMatrix2
    
    #Forming matrices of repeated values of follow up for each group
    #They are set into matrices of the same dimension for easy direct comparison
    sharedMatrix1 <-
      matrix(rep(fu_vector1, times = n_obs2), ncol = n_obs2)
    sharedMatrix2 <-
      matrix(rep(fu_vector2, times = n_obs1),
             nrow = n_obs1,
             byrow = T)
    sharedMatrix <-
      matrix(pmin(sharedMatrix1, sharedMatrix2),
             nrow = n_obs1,
             ncol = n_obs2)
    
    #Defining the winner
    winner <- matrix(nrow = n_obs1, ncol = n_obs2)
    winner[diffMatrix < 0 & bigMatrix1 <= sharedMatrix] <- 1
    winner[diffMatrix > 0 & bigMatrix2 <= sharedMatrix] <- -1
    winner[diffMatrix == 0] <- 0
    winner[sharedMatrix < bigMatrix1 &
             sharedMatrix < bigMatrix2] <- 0
    
    #Returning output
    return(winner)
  }

## Perform pairwise comparisons for continuous or categorical data
# Patient with a larger value "wins"
win_matrix_continuous <- function(vector1, vector2) {
  # vector1 = vector of values for control patients
  # vector2 = vector of values for treatment patients
  n_obs1 <- length(vector1)
  n_obs2 <- length(vector2)
  
  #Forming matrices of repeated values of continous variable for each group
  bigMatrix1 <- matrix(rep(vector1, times = n_obs2), ncol = n_obs2)
  bigMatrix2 <-
    matrix(rep(vector2, times = n_obs1),
           nrow = n_obs1,
           byrow = T)
  
  #Comparing them and then returning win matrix
  diffMatrix <- bigMatrix2 - bigMatrix1
  winner <- sign(diffMatrix)
  return(winner)
}

## Perform pairwise comparisons for recurrent event data
# Patient who is known to be event-free longer "wins"
winner_nEvents <-
  # timings1 = matrix of event time(s) for control patients
  # timings2 = vector of event time(s) for treatment patients
  ## organized with one row per patient and J columns for each event
  ## where J = max number of events across patients
  # fu_vector1 = vector of follow-up times for control patients
  # fu_vector2 = vector of follow-up times for treatment patients
  function(timings1,
           timings2,
           fu_vector1,
           fu_vector2) {
    #Calculate number observations
    n_obs1 <- length(fu_vector1)
    n_obs2 <- length(fu_vector2)
    
    #Forming matrices of repeated values of follow up for each group
    #They are set into matrices of the same dimension for easy direct comparison
    sharedMatrix1 <-
      matrix(rep(fu_vector1, times = n_obs2), ncol = n_obs2)
    sharedMatrix2 <-
      matrix(rep(fu_vector2, times = n_obs1),
             nrow = n_obs1,
             byrow = T)
    sharedMatrix <-
      matrix(pmin(sharedMatrix1, sharedMatrix2),
             nrow = n_obs1,
             ncol = n_obs2)
    
    #For each row of HFH matrix need to identify the number of elements that are
    #less than the shared follow-up from shared follow-up matrix
    #This is going to require an array of matrices
    max_hfh1 <- ncol(timings1)
    max_hfh2 <- ncol(timings2)
    
    hfh_array1 <- array(dim = c(n_obs1, n_obs2, max_hfh1))
    shared_array1 <- array(dim = c(n_obs1, n_obs2, max_hfh1))
    hfh_array2 <- array(dim = c(n_obs1, n_obs2, max_hfh2))
    shared_array2 <- array(dim = c(n_obs1, n_obs2, max_hfh2))
    
    for (i in 1:max_hfh1) {
      hfh_array1[, , i] <-
        matrix(rep(timings1[, i], times = n_obs2),
               ncol = n_obs2,
               byrow = FALSE)
      shared_array1[, , i] <- sharedMatrix
    }
    for (i in 1:max_hfh2) {
      hfh_array2[, , i] <-
        matrix(rep(timings2[, i], times = n_obs1),
               nrow = n_obs1,
               byrow = TRUE)
      shared_array2[, , i] <- sharedMatrix
    }
    
    nhfh_array1 <- shared_array1 - hfh_array1
    nhfh_array2 <- shared_array2 - hfh_array2
    
    #For each matrix creating an indicator variable for whether the event was before
    #the end of shared follow up
    nhfh_array1 <- sign(nhfh_array1)
    nhfh_array1[nhfh_array1 == 0] <- 1
    nhfh_array1[nhfh_array1 < 0] <- 0
    
    #Now summing these indicators across the arrays (i.e. summing across columns of input
    #time to HFH matrices)
    nhfh_matrix1 <- apply(nhfh_array1, 1:2, sum)
    
    #Repeating for the treatment group
    nhfh_array2 <- sign(nhfh_array2)
    nhfh_array2[nhfh_array2 == 0] <- 1
    nhfh_array2[nhfh_array2 < 0] <- 0
    nhfh_matrix2 <- apply(nhfh_array2, 1:2, sum)
    
    #Now comparing the number of HFH during shared follow up in treatment and control grou
    #Fewer hospitalisations is better with treatment
    winner <- sign(nhfh_matrix1 - nhfh_matrix2)
    return(winner)
  }

## Picks out the "winner" based on a hierarchical set of win ratio matrices
setClass(
  "multiWinner",
  representation(
    win_matrix = "matrix",
    decisions = "vector",
    perc_decisions = "vector"
  )
)

multi_criteria_winner <-
  # winner_list = a list of numeric matrices of win loss indicators
  ## first level is highest in the hierarchy and last level is lowest
  function(winner_list) {
    n_criteria <- length(winner_list)
    decisions <- vector(length = n_criteria)
    
    n_obs1 <- nrow(winner_list[[1]])
    n_obs2 <- ncol(winner_list[[1]])
    
    #Simple idea is to loop through each criteria (most important first), replacing any ties
    multi_winner <- matrix(0, nrow = n_obs1, ncol = n_obs2)
    
    for (i in 1:n_criteria) {
      multi_winner[multi_winner == 0] <-
        winner_list[[i]][multi_winner == 0]
      decisions[i] <- sum(abs(multi_winner))
    }
    
    #These decision statistics help to show how many of the winners were decided by which tier
    n_decisions <- max(decisions)
    decisions <- c(decisions[1], diff(decisions))
    perc_decisions <- 100 * decisions / n_decisions
    
    #Output data
    output <-
      new(
        "multiWinner",
        win_matrix = multi_winner,
        decisions = decisions,
        perc_decisions = perc_decisions
      )
    return(output)
  }

## Calculates the Win Ratio and associated statistics based on one win loss matrix
setClass(
  "WRStats",
  representation(
    winratio = "numeric",
    is.sig = "logical",
    logWR = "numeric",
    logWR_se = "numeric",
    winmatrix = "matrix",
    decisions = "vector",
    perc_decisions = "vector"
  )
)

getWRStats <-
  # winner_list = a numeric matrix of win loss indicators (see WWR:: genwr)
  # alpha = type I error rate as a numeric decimal (i.e. 0.05 for 5%)
  function(winner_list, alpha) {
    z_alpha <- qnorm(1 - alpha / 2)
    
    MCW <- multi_criteria_winner(winner_list)
    
    wrstats <- genwr(slot(MCW, "win_matrix"))
    winratio <- wrstats$wr
    logWR <- log(wrstats$wr)
    logWR_se <- logWR / wrstats$tr
    is.sig <- (logWR - z_alpha * logWR_se) > 0 &
      (logWR + z_alpha * logWR_se) > 0
    
    output <- new(
      "WRStats",
      winratio = winratio,
      is.sig = is.sig,
      logWR = logWR,
      logWR_se = logWR_se,
      winmatrix = slot(MCW, "win_matrix"),
      decisions = slot(MCW, "decisions"),
      perc_decisions = slot(MCW, "perc_decisions")
    )
  }