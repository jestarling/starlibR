#Runs the EM algorithm for the Bayesian Document Model.  (Not LDA.)

rm(list=ls())

#------------------------------------------------------------
#DATA LOADING:

#Read in data.
setwd("/Users/jennstarling/UTAustin/starlib/data")

#List of documents.
docs = read.table(file="./bbcsport.docs")

words = read.table(file='./bbcsport.terms')
colnames(words)='words'

classes = read.table(file='./bbcsport.classes',skip=4)
colnames(classes) = c('doc','class')
classes$doc = classes$doc + 1	#Add 1 so docs range from 1 to 737, to match 'docs' numbering.

mtx = read.table(file='./bbcsport.mtx',skip=2)
colnames(mtx) = c('word','doc','count')
#mtx: First column are terms.  Second column is documents.  Third is count of times each word appears in the doc in col 2.

#Set up key for labels.  Don't really need this.
class.names = c('athletics','cricket','football','rugby','tennis')

#Set up some useful indices.
D = nrow(docs)						#Total number of docs.
V = nrow(words)						#Total number of words in vocab.
K = length(unique(classes$class))	#Total number of topics (classes).

#Calculate true parameter values.
data=mtx
data$class = classes$class[match(data$doc,classes$doc)]	#Add true classes to data.

pi_k.true = table(classes$class) / nrow(classes)

beta_kw.true = matrix(0,nrow=K,ncol=V)

for (k in 1:K){
	data_k = data[data$class==(k-1),]
	w_idx = unique(data_k$word)
	words_aggregated = aggregate(count ~ word, data=data_k, FUN = sum)
	beta_kw.true[k,w_idx] = words_aggregated$count / sum(data_k$count)
}

#------------------------------------------------------------
#PARAMETERS SETUP:

#Initialize pi_k values.  (First parameter to be estimated.)
set.seed(123)
pi_0 = runif(K)
pi_0 = pi_0 / sum(pi_0)

#Initialize beta_kw values. (Each row is a vector of beta_w's for class k.)
set.seed(123)
beta_0 = matrix(runif(K*V),ncol=V)
beta_0 = beta_0 / rowSums(beta_0)

#------------------------------------------------------------
#Run EM algorithm.

results = em_bayesian.doc.model(data=data,D,K,V,pi_0,beta_0,maxiter=1000,tol=1E-14)

#Check convergence.
results$converged
results$iter

#Check pi_k's.
results$tracking_pi_k[1:results$iter,]
round(results$pi_0,4)
round(pi_k.true,4)
round(results$pi_k,4)

#Check betas.
round(beta_kw.true[,1:10],3)
round(results$beta_kw[,1:10],3)

#Plot loglikelihood.
jpeg(file='/Users/jennstarling/UTAustin/2016_Fall_SDS 383C_Statistical Modeling 1/Homework/HW 05/Images/logl.jpg')
plot(results$ll[10:results$iter],type='l',col='blue',main='EM Loglikelihood',xlab='iteration',ylab='log-ll')
dev.off()

#-----------------------
#2-a-ii & iiii: Find words with top 10 beta_kw true and estimated proportions for each topic.

#Set up true beta_kw's.
beta_kw.true = matrix(0,nrow=K,ncol=V)

for (k in 1:K){
	data_k = data[data$class==(k-1),]
	w_idx = unique(data_k$word)
	words_aggregated = aggregate(count ~ word, data=data_k, FUN = sum)
	beta_kw.true[k,w_idx] = words_aggregated$count / sum(data_k$count)
}

#Show top true and estimated words.
word.rank.true = data.frame(matrix("a",nrow=10,ncol=K))
word.rank.est = data.frame(matrix("a",nrow=10,ncol=K))
colnames(word.rank.true) = colnames(word.rank.est) = paste('k',(1:K)-1,sep='')

for (k in 1:K){
	idx.true = order(beta_kw.true[k,], decreasing = T)[1:10]
	idx.est = order(results$beta_kw[k,], decreasing = T)[1:10]
	
	word.rank.true[,k] = as.character(words[idx.true,])
	word.rank.est[,k] = as.character(words[idx.est,])
}

#Display results.
word.rank.true
word.rank.est

#Calculate proportion of top 10 words in estimated that were in true.
props.match = rep(0,K)
names(props.match) = paste('k',(1:K)-1,sep='')
for(k in 1:K){
	props.match[k] = sum(word.rank.est[,k] %in% word.rank.true[,k]) / 10
}
props.match 	#Display result.

#------------------------------------------------------------
#LOGLIKELIHOOD FUNCTION
#Set up expected value of the log-likelihood function for plotting.
logl = function(data,pi_k,beta_kw,pkd){
	
	ll = 0	#Initialize loglikelihood value.
	
	for (d in 1:D){ #Loop over docs.
		data_doc = data[data$doc==d, ]			#Subset of data for just that document.

		for (k in 1:K){ #Loop over topics
			
			betas_doc = beta_kw[k,data_doc$word]		#Ordered by word-order in doc.
			nwd_doc = data_doc$count				#Counts for each word in doc.
			w_idx = data_doc$word					#Words indices for words in doc.
			
			temp = nwd_doc * log(beta_kw[k,w_idx])	#Handling infinite values caused by 0 probabilities.
			p1 = sum(temp[is.finite(temp)])
			p2 = log(pi_k[k]) * pkd[k,d]
			ll = ll + p1 + p2
		} #end topics loop
	} #end docs loop
	
	return(ll)
}	

#------------------------------------------------------------
#EM ALGORITHM:

em_bayesian.doc.model = function(data,D,K,V,pi_0,beta_0,maxiter=100,tol=1E-6){
	
	#INITIALIZATIONS:
	
	#Store log-likelihood for each iteration.
	ll = rep(0,maxiter)
	
	tracking_pi_k = matrix(0,nrow=maxiter,ncol=K)
	tracking_pi_k[1,] = pi_0
		
	beta_kw = beta_0	#Initialize beta_kw to value provided.
	pi_k = pi_0			#Initialize pi_k to value provided.
	
	#Set up p(z_d = k | w, pi, beta)
	pkd = matrix(0,nrow=K,ncol=D)	#One value per doc per class.
	
	for (i in 1:maxiter){
		
		#CALCULATE LOGLIKELIHOOD:
		ll[i] = logl(data,pi_k,beta_kw,pkd)
		
		#-------------------------
		#EXPECTATION STEP:
		for (d in 1:D){	#Loop over docs.
		
			data_doc = data[data$doc==d, ]			#Subset of data for just that document.
			word_idx = data_doc$word				#Indices of words in doc d.
			
			for (k in 1:K){	#Loop over topics.
				
				betas_doc = beta_kw[k,data_doc$word]		#Ordered by word-order in doc.
				nwd_doc = data_doc$count				#Counts for each word in doc.
			
				#Update in log-scale.
				pkd[k,d] = sum(nwd_doc * log(betas_doc)) + log(pi_k[k])
		
		}	#end topics loop.
		
			#Normalize each doc column via log-exponentiation trick.
			c <- max(pkd[,d])
			denom <- sum(exp(pkd[,d] - c))
			pkd[,d] = exp(pkd[,d] - c) / denom
		
		} #end docs loop.
			
		#-------------------------
		#MAXIMIZATION STEP:
	
		#1. Update pi's.
		pi_k = rowSums(pkd) / sum(pkd)
		tracking_pi_k[i,] = pi_k
	
		#2. Update beta's.
		
		#Reset beta_kw to zero for new iteration of algorithm.
		beta_kw = matrix(0,nrow=K,ncol=V)
		
		for (d in 1:D){ #Loop over docs.
			
			data_doc = data[data$doc==d,]
			nwd_doc = data_doc$count
			w_idx = data_doc$word
			
			for (k in 1:K){ #Loop over topics.
				beta_kw[k,w_idx] = beta_kw[k,w_idx] + pkd[k,d] * nwd_doc
			} #End topics loop.
		} #End docs loop.
		
		#Normalize beta_kw
		beta_kw = beta_kw / rowSums(beta_kw)
	
		#-------------------------
		#CONVERGENCE CHECK:
		converged = 0
	
			if(i>2){
				if(abs(ll[i] - ll[i-1])/abs(ll[i]) < tol){
					converged = 1
					break
				}
			}
	
		#-------------------------
	
	} #end main for loop.
	
	iter = i
	
	return(list(pi_k=pi_k, 
		beta_kw = beta_kw, 
		ll=ll, 
		converged=converged, 
		iter=iter, 
		pi_0=pi_0, 
		beta_0=beta_0, 
		tracking_pi_k = tracking_pi_k))
	
} #end function.

#------------------------------------------------------------




