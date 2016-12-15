#BAYESIAN DOCUMENT MODEL - TOY EXAMPLE

#(Full doc model, not LDA)
#Docs can only have 1 topic.

###############################
###        TOY EXAMPLE:     ###
###############################

#3 docs: doc1 (true class sports), doc2 (true class entertainment), doc3 (true class entertainment)
#11 words in vocabulary: sports, hockey, soccer, football, tennis, tv, radio, actress, movie, fashion, netflix
#2 topics: sports (pi_k1) and entertainment (pi_k2).

doc1 = c('sports', 'sports', 'sports', 'hockey', 'soccer', 'football', 'tennis', 'tv', 'radio')
doc2 = c('tv','hockey','tv','actress','movie','movie','movie','fashion','netflix')
doc3 = c('tv','actress','movie','football')

words = cbind.data.frame(1:length(unique(c(doc1,doc2))), unique(c(doc1,doc2)))
colnames(words) = c('word','full.word')

topics = c('sports','entertainment')
topics = cbind.data.frame(classes=c(1,2),topics=topics)

dat.doc = c(rep(1,length(unique(doc1))), rep(2,length(unique(doc2))), rep(3,length(unique(doc3))))
dat.words.full = c(unique(doc1),unique(doc2),unique(doc3))
dat.word = words[match(dat.words.full,words$full.word),1]

dat.count = c(3,1,1,1,1,1,1,2,1,1,3,1,1,1,1,1,1)
dat.class = c(rep(1,length(unique(doc1))), rep(2,length(unique(doc2))), rep(2,length(unique(doc3))))

dat = cbind.data.frame(words.full = dat.words.full, word=dat.word, doc=dat.doc, count=dat.count, class=dat.class)

#Set up some useful indices.
D = 3			#Total number of docs.
V = 11			#Total number of words in vocab.
K = 2			#Total number of topics (classes).

#Calculate true parameter values.
pi_k.true = c(1/3,2/3) 		#Proportion of docs in each class.

beta_kw.true = matrix(0,nrow=K,ncol=V) #Probability of each word in each class.

for (k in 1:K){
	data_k = dat[dat$class==k,]
	w_idx = unique(data_k$word)
	words_aggregated = aggregate(count ~ word, data=data_k, FUN = sum)
	beta_kw.true[k,w_idx] = words_aggregated$count / sum(data_k$count)
}

#------------------------------------------------------------
#Initialize values for EM algorithm.

#Initialize pi_0.
pi_0 = c(.5,.5)

#Initialize beta_0. (Each row is a vector of beta_w's for class k.)
set.seed(100)
beta_0 = matrix(runif(K*V),ncol=V)
beta_0 = beta_0 / rowSums(beta_0)

#------------------------------------------------------------
#Run EM algorithm.

results.toy = em_bayesian.doc.model(data=dat,D,K,V,pi_0,beta_0,maxiter=20,tol=1E-14)

#Check convergence.
results.toy$converged
results.toy$iter

#Check pi_k's.
round(pi_k.true,4)
round(results.toy$pi_k,4)

results.toy$pi_0

results.toy$tracking_pi_k[1:results.toy$iter,]

#Check betas.
round(beta_kw.true,4)
round(results.toy$beta_kw,4)

par(mfrow=c(1,2))
plot(results.toy$ll[1:results.toy$iter],type='l',col='blue',main='Loglikelihood',xlab='iteration',ylab='log-ll')
plot(exp(results.toy$ll[1:results.toy$iter]),type='l',col='blue',main='likelihood',xlab='iteration',ylab='log-ll')

#Show top true and estimated words.
word.rank.true = data.frame(matrix("a",nrow=10,ncol=K))
word.rank.est = data.frame(matrix("a",nrow=10,ncol=K))
colnames(word.rank.true) = colnames(word.rank.est) = paste('k',1:K,sep='')

for (k in 1:K){
	idx.true = order(beta_kw.true[k,], decreasing = T)[1:10]
	idx.est = order(results.toy$beta_kw[k,], decreasing = T)[1:10]
	
	word.rank.true[,k] = as.character(words[match(idx.true,words$word),2])
	word.rank.est[,k] = as.character(words[match(idx.est,words$word),2])
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