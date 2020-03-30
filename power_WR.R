#Required packages 
require("binom") 
library("binom")



#Defining object to hold win ratio power calculation results
setClass("WRPower", representation(required_n= "numeric", meanWR="numeric", perc_decisions="vector", logWR="vector",
                                   logWR_se="vector", 
                                   fu.start="vector", power="numeric", n.sim="numeric", n.iter="numeric", hr.death="numeric", hr.hfh="numeric", 
                                   p.active="numeric",  rate.control.d="numeric", rate.control.hfh="numeric",
                                   frailty.scale="numeric", frailty.exponent="numeric", criteria="vector", cont_parameters="vector"))


#This function performs power calculations for win ratio on data simulated from joint frailty data
WRPower<-function(n.start=NULL,fu.start, power=0.8,alpha=0.05, n.sim=1000, n.iter=3, p.active=0.5, hr.death, hr.hfh, rate.control.d, rate.control.hfh, 
                      frailty.scale=0, frailty.exponent=0, criteria, cont_params=NULL) {


#Choosing a sensible starting place number of patients that is reasonably quick to run (if none is specified)
if(is.null(n.start)) {
  n.start<-500
}
    
  
  
#Generating z-statistics from supplied alpha and power
z_power<-qnorm(power)
z_alpha_desired<-qnorm(1-alpha/2)

#Setting starting number of observations, and censoring distribution for first iteration
n.obs<-n.start
length.fu<-fu.start 

#Looping through iterative process to determine required sample size
  for(iter in 1:n.iter) {
  print(paste("Iteration", iter))
  
  logWR<-vector(length=n.sim)
  logWR_se<-vector(length=n.sim)
  is.sig<-vector(length=n.sim) 
  

#Looping through n.sim simulated datasets within each iteration
  	for(sim in 1:n.sim) {
  	winner_list<-NULL
  	
  	#Simulating data in each iteration and posting results to vectors 
  	jf_data<-simulateJF(n.obs=n.obs, p.active=p.active, hr.death=hr.death, hr.hfh=hr.hfh, 
  	rate.control.d=rate.control.d, rate.control.hfh=rate.control.hfh, 
  	length.fu=length.fu, frailty.scale=frailty.scale, frailty.exponent=frailty.exponent)
 	
	
	#If there are three tiers then assuming KCCQ data is supplied (will update to generalize later)	
        if(length(criteria)==3) {
      	 	if(length(cont_params)!=4) {
        	print("Parameters for KCCQ not correctly specified")
        	}
      	n_obs1<-length(slot(jf_data, "treatment")[slot(jf_data, "treatment")==0]) 
      	n_obs2<-length(slot(jf_data, "treatment")[slot(jf_data, "treatment")==1])
      	KCCQ1<-rnorm(n_obs1, mean=cont_params[1], sd=cont_params[2])
      	KCCQ2<-rnorm(n_obs2, mean=cont_params[3], sd=cont_params[4])
      #	print("KCCQ loop")
      #	print(KCCQ1)
      	winner_list<-evaluate_jf(jf_data, criteria=criteria, 	KCCQ_control=KCCQ1, KCCQ_treatment=KCCQ2)
      	}

        else {
        winner_list<-evaluate_jf(jf_data, criteria=criteria)  
        }
  	
  	#Calculating win ratio statistics and storing them
  	wrstats<-getWRStats(winner_list, alpha)
  	logWR[sim]<-slot(wrstats, "logWR")
  	logWR_se[sim]<-slot(wrstats, "logWR_se")
  	is.sig[sim]<-slot(wrstats, "is.sig")
	
	#Displaying progress of simulations
  	  if(sim%%10==0){
  	  print(paste("Simulation", sim))
  	  }
  	}
  
  #size is proportional to the square of [z(alpha/2) + z(beta)]-adapting sample size on this basis
  z<-logWR/logWR_se
  beta<-1-power
  #Idenfitying which alpha, based upon the current sample size, would lead to the required level of statistical power
  z_alpha_obs<-quantile(abs(z), beta)
  
  #Estimating required N based on required power, observed alpha above, and formula above 
  required_n<-n.obs*(z_alpha_desired+z_power)^2 / (z_power+z_alpha_obs)^2
  required_n<-round(required_n)
  
  #Displaying adaptation of sample size 
  current_power<-mean(is.sig)
  current_z<-qnorm(current_power)
  print(paste("Iteration", iter, "Power:" , current_power, "N", n.obs, "Estimated required:", required_n))
  
  #Redefining inputs for joint frailty simulations based on updated sapmle size
  n.obs_old<-n.obs
  n.obs<-required_n
  length.fu<-sample(fu.start, size=n.obs, replace=TRUE)
  
  	#Simulations become very slow once N gets too large 
  	if(required_n>10000) {
  	  stop("Required sample size too large to run simulations")
  	}
  
  	#Displaying final estimated sample size and results from final iteration
  	if(iter==n.iter) {
  	power_ci<-binom.confint(sum(is.sig),n.sim, method="exact")
  	power_ci<-binom.confint(sum(is.sig),n.sim, method="exact")
  	print(paste("95% CI for power with n=", n.obs_old ,":", power_ci$lower, power_ci$upper))
  	print(paste("Refined estimate (not simulated)", n.obs))
  	}
  
  }

#Preparing output
mean_WR<-exp(mean(logWR))
#Percentage decisions, just being lazy and using last simulation
perc_decisions<-slot(wrstats, "perc_decisions")

#Creating object that contains the results and the input
output<-new("WRPower", required_n=required_n, meanWR=mean_WR, perc_decisions=perc_decisions,
            logWR=logWR, logWR_se=logWR_se, 
            p.active=p.active, hr.death=hr.death, hr.hfh=hr.hfh, 
            rate.control.d=rate.control.d, rate.control.hfh=rate.control.hfh, frailty.scale=frailty.scale, 
            frailty.exponent=frailty.exponent, criteria=criteria, cont_parameters=cont_params, fu.start=fu.start,
            power=power, n.sim=n.sim, n.iter=n.iter) 
}



