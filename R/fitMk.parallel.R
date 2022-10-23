## parallelized version of fitMk using optimParallel

fitMk.parallel<-function(tree,x,model="SYM",ncores=1,...){
	## compute states
	ss<-sort(unique(x))
	m<-length(ss)
	## set pi
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-"equal"
	if(is.numeric(pi)) root.prior<-"given"
	if(pi[1]=="equal"){ 
		pi<-setNames(rep(1/m,m),ss)
		root.prior<-"flat"
	} else if(pi[1]=="fitzjohn") root.prior<-"nuisance"
	if(is.numeric(pi)){ 
		pi<-pi/sum(pi)
		if(is.null(names(pi))) pi<-setNames(pi,ss)
		pi<-pi[ss]
	} 
	## create object of class "fitMk"
	unfitted<-fitMk(tree,x,model=model,pi=pi,opt.method="none")
	## get initial values for optimization
	if(hasArg(q.init)) q.init<-list(...)$q.init
	else q.init<-rep(m/sum(tree$edge.length),
		max(unfitted$index.matrix,na.rm=TRUE))
	## create likelihood function
	loglik<-function(par,lik,index.matrix){
		Q<-makeQ(nrow(index.matrix),exp(par),index.matrix)
		-lik(Q)
	}
	## create cluster
	cl<-makeCluster(ncores)
	## optimize model
	fit.parallel<-optimParallel(
		log(q.init),
		loglik,lik=unfitted$lik,
		index.matrix=unfitted$index.matrix,
		lower=log(1e-12),
		upper=log(max(nodeHeights(tree))*100),
		parallel=list(cl=cl,forward=FALSE,
		loginfo=TRUE)
	)
	## stop cluster
	## setDefaultCluster(cl=NULL)
	stopCluster(cl)
	## create object
	estQ<-makeQ(nrow(unfitted$index.matrix),
		exp(fit.parallel$par),
		unfitted$index.matrix)
	colnames(estQ)<-rownames(estQ)<-ss
	temp<-fitMk(tree,x,fixedQ=estQ,pi=pi)
	object<-list()
	object$logLik<-(-fit.parallel$value[1])
	object$rates<-exp(fit.parallel$par)
	object$index.matrix<-unfitted$index.matrix
	object$states<-unfitted$states
	object$pi<-temp$pi
	object$method<-"optimParallel"
	object$root.prior<-temp$root.prior
	object$lik<-temp$lik
	class(object)<-"fitMk"
	object
}
