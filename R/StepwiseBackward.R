#---FUNCTION-------------------------
#Custom function for backward stepwise.
StepwiseBackward <- function(X,y){
	
	Xwhole = X				#X matrix with all values.
	
	#Loop through possible numbers of predictors to add.
	for (i in 1:ncol(X)){ 
		
		cp_xnew = myCp(Xwhole,y)	#Benchmark cp for xnew model.
		cp = rep(0,ncol(Xwhole))	#Vector to hold CP values for remaining predictors.
		vars = names(Xwhole)		#Vector of variable names still remaining to add.
		
		#Loop through each column of remaining X variables.
		for (j in 1:ncol(Xwhole)){
			#Calculate cp for each potential new model.
			cp[j] = myCp(cbind(Xwhole[,-j]),y)
		}
		
		#If none of the new models have lower cp than the old model,
		#break and return the old model.
		if (min(cp) >= cp_xnew){
			break;
		}
			
		#Pick the Xwhole column with the loweset cp to remove.
		#Remove it from Xwhole and add it to Xnew.
		removeVar = vars[which(cp==min(cp))]	#Pick the variable to move.

		Xwhole = Xwhole[,-which(cp==min(cp)),drop=F]	
	}
	
	final_cp = myCp(Xwhole,y)
	return(list(X=Xwhole,vars=colnames(Xwhole),cp = final_cp))
		
} #End backward stepwise function.
#---END FUNCTION---------------------