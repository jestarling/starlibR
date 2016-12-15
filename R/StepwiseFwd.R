StepwiseFwd <- function(X,y){
	
	Xwhole = X				#X matrix with all values.
	Xnew = NULL				#Start with null model.  (myLinearReg adds intercept)
	keeps = character()		#Empty char vector to hold names of kept vars.
	
	#Loop through possible numbers of predictors to add.
	for (i in 1:ncol(X)){ 
		
		cp_xnew = myCp(Xnew,y)		#Benchmark cp for xnew model.
		cp = rep(0,ncol(Xwhole))	#Vector to hold CP values for remaining predictors.
		vars = names(Xwhole)		#Vector of variable names still remaining to add.
		
		
		#Loop through each column of remaining X variables.
		for (j in 1:ncol(Xwhole)){
			#Calculate cp for each potential new model.
			cp[j] = myCp(cbind(Xnew,Xwhole[,j]),y)
		}
		
		#If none of the new models have lower cp than the old model,
		#break and return the old model.
		if (min(cp) >= cp_xnew){
			break;
		}
			
		#Pick the Xwhole column with the loweset cp.
		#Remove it from Xwhole and add it to Xnew.
		moveVar = vars[which(cp==min(cp))]	#Pick the variable to move.
		Xnew = cbind(Xnew,Xwhole[,moveVar])
		Xwhole = Xwhole[,-which(cp==min(cp)),drop=F]	
		keeps = c(keeps,moveVar)	#Keep track of vars added to model.
	}
	
	colnames(Xnew) = keeps
	final_cp = myCp(Xnew,y)
	return(list(X=Xnew,vars=keeps,cp = final_cp))	
} #End fwd stepwise function.
#---END FUNCTION---------------------