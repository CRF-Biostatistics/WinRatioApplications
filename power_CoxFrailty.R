require("survival")
require("binom")
library(survival)
library("binom")

#Creating a function to simulate power from using a Cox proportional hazards model on joint frailty data
power_CoxFrailty<-function(n.start=NULL,fu.start, power=0.8,alpha=0.05, n.sim=1000, n.iter=3, 
                           p.active=0.5, hr.death, hr.hfh, rate.control.d, rate.control.hfh, 
                  frailty.scale=0, frailty.exponent=0, composite=TRUE) {
  
  #Generating sensible starting values where none are given
   if(is.null(n.start)) {
    n.start<-500
  }
  
  #Generating z-statistics
  z_power<-qnorm(power)
  z_alpha_desired<-qnorm(1-alpha/2)
  
  #Setting starting number of observations, and censoring distribution for first iteration
  n.obs<-n.start
  length.fu<-fu.start 
  
   
  for(iter in 1:n.iter) {
    print(paste("Iteration", iter))
    logHR<-vector(length=n.sim)
    logHR_se<-vector(length=n.sim)
    is.sig<-vector(length=n.sim) 
    z<-vector(length=n.sim) 
    p<-vector(length=n.sim) 
  	
	  for(i in 1:n.sim) {
  	  jf_data<-simulateJF(n.obs=n.obs, p.active=p.active, hr.death=hr.death, hr.hfh=hr.hfh, 
                          rate.control.d=rate.control.d, rate.control.hfh=rate.control.hfh, 
                          length.fu=length.fu, frailty.scale=frailty.scale, frailty.exponent=frailty.exponent)
      
  	  #Grabbing HFH data and time to censoring
  	  t_hfh<-slot(jf_data, "t_hfh")[,1]
  	  t_hfh[t_hfh==Inf]<-slot(jf_data, "t_death")[t_hfh==Inf]
  	  t_hfh[t_hfh==Inf]<-slot(jf_data, "followup")[t_hfh==Inf]
  	  hfh<-vector(length=n.obs)
  	  hfh<-slot(jf_data, "indiv_hfh")
  	  hfh[hfh>1]<-1

	  #Make it a composite outcome if specified
	  	if(composite==TRUE) {
		hfh[slot(jf_data, "death")==1]<-1
		}	 
  	  
	  #Running Cox model and grabbing the relevant statistics 
  	  cox<-summary(coxph(Surv(time=t_hfh, event=hfh)~slot(jf_data, "treatment")))
  	  p[i]<-cox$coefficients[5]
  	  logHR[i]<-cox$coefficients[1]
  	  logHR_se[i]<-cox$coefficients[3]
  	  is.sig[i]<-(p[i]<0.05)
  	  z[i]<-cox$coefficients[4]
  	  }
    
	
    beta<-1-power

    #Idenfitying observed alpha at required level of statistical power
    z_alpha_obs<-quantile(abs(z), beta)

    #Estimating required N based on required power, observed alpha above, and formula above 
    print(paste("Z-statistic at 1-power (", beta, ")=",round(z_alpha_obs, digits=2)))
    print(paste("Z- desired=", round(z_alpha_desired, digits=2)))
    required_n<-n.obs*(z_alpha_desired+z_power)^2 / (z_power+z_alpha_obs)^2
    required_n<-round(required_n)
    
    #Relaying adaptation of sample size back to user
    current_power<-mean(is.sig)
    current_z<-qnorm(current_power)
    print(paste("Iteration", iter, "Power:" , current_power, "N", n.obs, "Estimated required:", required_n))
    
    #Redefining inputs for joint frailty simulations based on updated sapmle size
    n.obs_old<-n.obs
    n.obs<-required_n
    length.fu<-sample(fu.start, size=n.obs, replace=TRUE)
    
    if(required_n>30000) {
      stop("Required sample size too large to run simulations")
    }
    
    #Displaing final information for user 
    if(iter==n.iter) {
      power_ci<-binom.confint(sum(is.sig),n.sim, method="exact")
      power_ci<-binom.confint(sum(is.sig),n.sim, method="exact")
      print(paste("95% CI for power with n=", n.obs_old ,":", power_ci$lower, power_ci$upper))
      print(paste("Refined estimate (not simulated)", n.obs))
    }
    
  }
  
  #Preparing output, simply required N here
  mean_HR<-exp(mean(logHR))
  print(paste("Mean marginal HR", round(mean_HR, digits=2) ))
  output<-required_n 
  return(output)
}







