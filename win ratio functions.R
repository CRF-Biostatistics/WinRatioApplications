#Requirements WWR package
#Loading required libraries
#install.packages("WWR")
library("WWR")

############################################################################
#Win ratio for time to first event
#Inputs are vectors of time to first event for control (1) and treament (2)  groups , and 
#follow-up times for control (1) and treatment (2) groups
win_matrix_firstEvent<-function(t_vector1, t_vector2, fu_vector1, fu_vector2) {

#Collecting number of treatment and control patients and forming matrix
n_obs1<-length(t_vector1) 
n_obs2<-length(t_vector2) 
win_matrix_firstEvent<-matrix(nrow=n_obs1, ncol=n_obs2)

#Forming matrices of repeated values of time to first event for each group
#They are set into matrices of the same dimension for easy direct comparison
bigMatrix1<-matrix(rep(t_vector1, times=n_obs2), ncol=n_obs2)
bigMatrix2<-matrix(rep(t_vector2, times=n_obs1), nrow=n_obs1, byrow=T)
diffMatrix<- bigMatrix1-bigMatrix2	

#Forming matrices of repeated values of follow up for each group
#They are set into matrices of the same dimension for easy direct comparison
sharedMatrix1<-matrix(rep(fu_vector1, times=n_obs2), ncol=n_obs2)
sharedMatrix2<-matrix(rep(fu_vector2, times=n_obs1), nrow=n_obs1, byrow=T)
sharedMatrix<-matrix(pmin(sharedMatrix1, sharedMatrix2), nrow=n_obs1, ncol=n_obs2)

#Defining the winner
winner<-matrix(nrow=n_obs1, ncol=n_obs2)
winner[diffMatrix<0 & bigMatrix1<=sharedMatrix]<- 1
winner[diffMatrix>0 & bigMatrix2<=sharedMatrix]<- -1
winner[diffMatrix==0]<-0
winner[sharedMatrix<bigMatrix1 & sharedMatrix<bigMatrix2]<-0

#Returning output
return(winner)
}
#############################################################################



#############################################################################
#Comparison for a continuous variables (easy)
#Inputs are vectors of values for treatment and control group (larger values are better,invert as necessary)
win_matrix_continuous<-function(vector1, vector2) {
n_obs1<-length(vector1)
n_obs2<-length(vector2)

#Forming matrices of repeated values of continous variable for each group
bigMatrix1<-matrix(rep(vector1, times=n_obs2), ncol=n_obs2)
bigMatrix2<-matrix(rep(vector2, times=n_obs1), nrow=n_obs1, byrow=T)

#Comparing them and then returning win matrix 
diffMatrix<- bigMatrix2-bigMatrix1	
winner<-sign(diffMatrix)
return(winner)
}
##############################################################################


########################################################################
#Win ratio based on number of hospitalisations in shared follow up 
#Input is two matrices of control (1) and treatment (2) times to HFH, organised with one row per patient 
#and J columns where J is the maximum number of events per patient

#These were data used for testing 
timings1<-matrix(1:9, nrow=3, ncol=3)
timings2<-matrix(20:1, nrow=5, ncol=4)
fu_vector1<-c(6,5,15)
fu_vector2<-c(5,4,6, 7, 20)

winner_nEvents<-function(timings1,timings2, fu_vector1, fu_vector2){

#Calculate number observations
n_obs1<-length(fu_vector1)
n_obs2<-length(fu_vector2)


#Forming matrices of repeated values of follow up for each group
#They are set into matrices of the same dimension for easy direct comparison
sharedMatrix1<-matrix(rep(fu_vector1, times=n_obs2), ncol=n_obs2)
sharedMatrix2<-matrix(rep(fu_vector2, times=n_obs1), nrow=n_obs1, byrow=T)
sharedMatrix<-matrix(pmin(sharedMatrix1, sharedMatrix2), nrow=n_obs1, ncol=n_obs2)


#For each row of HFH matrix need to identify the number of elements that are 
#less than the shared follow-up from shared follow-up matrix 
#This is going to require an array of matrices

max_hfh1<-ncol(timings1)
max_hfh2<-ncol(timings2)

hfh_array1<-array(dim=c(n_obs1, n_obs2,max_hfh1))
shared_array1<-array(dim=c(n_obs1, n_obs2,max_hfh1))
hfh_array2<-array(dim=c(n_obs1, n_obs2,max_hfh2))
shared_array2<-array(dim=c(n_obs1, n_obs2,max_hfh2))


	for(i in 1:max_hfh1) {
	hfh_array1[,,i]<-matrix(rep(timings1[,i], times=n_obs2), ncol=n_obs2, byrow=FALSE)
	shared_array1[,,i]<-sharedMatrix
	}
	for(i in 1:max_hfh2) {
	hfh_array2[,,i]<-matrix(rep(timings2[,i], times=n_obs1), nrow=n_obs1, byrow=TRUE)
	shared_array2[,,i]<-sharedMatrix
	}

nhfh_array1<-shared_array1-hfh_array1
nhfh_array2<-shared_array2-hfh_array2

#For each matrix creating an indicator variable for whether the event was before
#the end of shared follow up
nhfh_array1<-sign(nhfh_array1)
nhfh_array1[nhfh_array1==0] <-1
nhfh_array1[nhfh_array1<0] <-0
#Now summing these indicators across the arrays (i.e. summing across columns of input 
#time to HFH matrices)
nhfh_matrix1<-apply(nhfh_array1, 1:2, sum)

#Repeating for the treatment group
nhfh_array2<-sign(nhfh_array2)
nhfh_array2[nhfh_array2==0] <-1
nhfh_array2[nhfh_array2<0] <-0
nhfh_matrix2<-apply(nhfh_array2, 1:2, sum)

#Now comparing the number of HFH during shared follow up in treatment and control grou
#Fewer hospitalisations is better with treatment
winner<-sign(nhfh_matrix1-nhfh_matrix2) 
return(winner)

}
############################################################################



################################################################################
##This function decides the "winner" based on a set of win ratio matrices that are given in a list
##and chucks out some summary statistics 
setClass("multiWinner", representation(win_matrix= "matrix", decisions="vector", perc_decisions="vector"))

multi_criteria_winner<-function(winner_list) {
n_criteria<-length(winner_list)
decisions<-vector(length=n_criteria)

n_obs1<-nrow(winner_list[[1]])
n_obs2<-ncol(winner_list[[1]])

#Simple idea is to loop through each criteria (most important first), replacing any ties
multi_winner<-matrix(0, nrow=n_obs1, ncol=n_obs2)

	for(i in 1:n_criteria) {
	multi_winner[multi_winner==0] <-winner_list[[i]][multi_winner==0]	
	decisions[i]<-sum(abs(multi_winner))
	}

#These decision statistics help to show how many of the winners were decided by which tier 
n_decisions<-max(decisions)
decisions<-c(decisions[1], diff(decisions))
perc_decisions<-100*decisions/n_decisions

#Output data
output<-new("multiWinner",win_matrix=multi_winner,decisions=decisions, perc_decisions=perc_decisions) 
return(output)
}
#########################################################################################


##This function piggybacks existing functions (genwr) to give win ratio statistics once
##a single matrix of winners and losers exists 
setClass("WRStats", representation(winratio= "numeric", is.sig="logical", 
logWR="numeric", logWR_se="numeric", winmatrix="matrix",decisions="vector", perc_decisions="vector"))

getWRStats<-function(winner_list, alpha) {
z_alpha<-qnorm(1-alpha/2)

MCW<-multi_criteria_winner(winner_list)

wrstats<-genwr(slot(MCW, "win_matrix"))
winratio<-wrstats$wr
logWR<-log(wrstats$wr)
logWR_se<-logWR/wrstats$tr
is.sig<-(logWR-z_alpha*logWR_se)>0 & (logWR+z_alpha*logWR_se)>0


output<-new("WRStats", winratio=winratio, is.sig=is.sig,
logWR=logWR, logWR_se=logWR_se, winmatrix=slot(MCW, "win_matrix"), 
decisions=slot(MCW, "decisions"), perc_decisions=slot(MCW, "perc_decisions") )

}


#############################################################
#Short wrapper function to get win list

evaluate_jf<-function(jf_data, criteria, KCCQ_control=NULL, KCCQ_treatment=NULL){

fu_vector1<-slot(jf_data, "followup")[slot(jf_data, "treatment")==0]
fu_vector2<-slot(jf_data, "followup")[slot(jf_data, "treatment")==1]

t_death1<-slot(jf_data, "t_death")[slot(jf_data, "treatment")==0]
t_death2<-slot(jf_data, "t_death")[slot(jf_data, "treatment")==1]

winner_death<-win_matrix_firstEvent(t_death1,t_death2, fu_vector1, fu_vector2)

  if(length(criteria)==1) {
    if(criteria=="death"){
    winner_list<-list()
    winner_list[[1]]<-winner_death
    }
  }

winner_list<-list()

  if(length(criteria)==2 | length(criteria)==3) {
    
  	if(criteria[1]=="death" & criteria[2]=="firstHFH") {
  	t_hfh1<-slot(jf_data, "t_hfh")[slot(jf_data, "treatment")==0,1]
  	t_hfh2<-slot(jf_data, "t_hfh")[slot(jf_data, "treatment")==1,1]
  	winner_hfh<-win_matrix_firstEvent(t_hfh1,t_hfh2, fu_vector1, fu_vector2)
  	}
  	
  	if(criteria[1]=="death" & criteria[2]=="nHFH") {
  	t_hfh1<-slot(jf_data, "t_hfh")[slot(jf_data, "treatment")==0,]
  	t_hfh2<-slot(jf_data, "t_hfh")[slot(jf_data, "treatment")==1,]
  	winner_hfh<-winner_nEvents(t_hfh1,t_hfh2, fu_vector1, fu_vector2)
  	}
  winner_list[[1]]<-winner_death
  winner_list[[2]]<-winner_hfh
  }

	if(length(criteria)==3) {
	winner_KCCQ<-win_matrix_continuous(vector1=KCCQ_control,vector2=KCCQ_treatment)
	winner_list[[3]]<-winner_KCCQ
	}

return(winner_list)
}
#####################################################################




