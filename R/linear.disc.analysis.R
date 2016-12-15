#------------------------------------------------------------------------
#MY CUSTOM LDA FUNCTION FUNCTION:

linear.disc.analysis = function(train,test){
	k = length(unique(train$y))	#Number of classes.
	p = ncol(train)-1			#Number of predictors.
	n = nrow(train)				#Number of obs.

	#1. Prior probabilities of groups, using training sample proportions.
	prior.probs = table(train$y)/n

	#2. Calculate centroids (class-dependent means) and pooled sample cov matrix W.
	M = matrix(0, nrow=k, ncol=p)
	W = matrix(0,nrow=p,ncol=p)

	for (i in 1:k){
		class.idx = train$y == i	#Identify class indices.
		X.temp = X.train[class.idx,]
	
		M[i,] = apply(X.temp,2,mean) #Centroids for class i.
		W = W + cov(X.temp)			 #Running sum of W cov matrix.
	}
	
	#3. Calculate eigen-decomposition of W.
	W.eigen <- eigen(W,symmetric=T)
	W.vecs <- W.eigen$vectors
	W.neghalf <- W.vecs %*% diag(1 / sqrt(W.eigen$values)) %*% t(W.vecs)

	#4. Calculate projection matrix V.  (See pg 114 of ESL for steps.)
	Mstar = M %*% W.neghalf

	Bstar = cov(Mstar)
	Bstar.eigen = eigen(Bstar)

	Vstar = Bstar.eigen$vectors 	#Using negative to match ESL plot.
	Db = Bstar.eigen$values

	V = W.neghalf %*% Vstar		#Projection matrix.

	#Project data.
	projected.data = X.train %*% V
	
	return(V,projected.data)
}
