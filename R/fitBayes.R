# this function calls various MCMC methods
# written by Liam J. Revell 2011, 2015

fitBayes<-function(tree,x,ngen=10000,model="BM",method="reduced",control=list()){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(model=="BM"){
		if(method=="reduced")
			X<-mcmcBM(tree,x,ngen,control)
		else if(method=="full")
			X<-mcmcBM.full(tree,x,ngen,control)
		else
			stop("This is not a recognized method for model='BM'.")
	} else if(model=="lambda"){
		if(method=="reduced")
			X<-mcmcLambda(tree,x,ngen,control)
		else if(method=="full")
			stop("This method is not yet implmented for model='lambda'.")
		else
			stop("This is not a recognized method for model='lambda'.")
	} else
		stop("This is not a recognized model.")

	return(X)

}
