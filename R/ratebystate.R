## simulation based test for a correlation between the state of x & the rate of y
## written by Liam J. Revell 2013, 2017, 2019, 2021

ratebystate<-function(tree,x,y,nsim=100,corr=c("pearson","spearman"),...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	corr<-corr[1]
	if(hasArg(sim.method)) sim.method<-list(...)$sim.method
	else sim.method<-"sim.corrs"
	if(hasArg(method)) method<-list(...)$method
	else method<-"by.node"
	if(hasArg(message)) message<-list(...)$message
	else message<-TRUE
	if(hasArg(logarithm)) logarithm<-list(...)$logarithm
	else logarithm<-FALSE
	if(!is.binary(tree)) tree<-multi2di(tree,random=FALSE)
	V<-phyl.vcv(cbind(x[tree$tip.label],y[tree$tip.label]),vcv(tree),lambda=1)$R
	if(method=="by.branch"){
		aa<-c(x[tree$tip.label],fastAnc(tree,x))
		names(aa)[1:length(tree$tip)]<-1:length(tree$tip)
		aa<-rowMeans(matrix(aa[tree$edge],nrow(tree$edge),2))
		a<-vector()
		for(i in 1:tree$Nnode+length(tree$tip)){
			j<-which(tree$edge[,1]==i)
			a[i-length(tree$tip)]<-sum(aa[j]*tree$edge.length[j])/sum(tree$edge.length[j])
		}
		names(a)<-1:tree$Nnode+length(tree$tip)
	}	
	else a<-fastAnc(tree,x)
	if(logarithm) a<-exp(a)
	b<-pic(y,tree)[names(a)]^2
	r<-cor(a,b,method=corr)
	beta<-setNames(lm(b~a)$coefficients[2],NULL)
	foo<-function(tree,V){
		if(sim.method=="fastBM") XY<-fastBM(tree,nsim=2)%*%sqrt(diag(diag(V)))
		else if(sim.method=="sim.corrs") XY<-sim.corrs(tree,V)
		a<-fastAnc(tree,XY[,1])
		b<-pic(XY[,2],tree)[names(a)]^2
		r<-cor(a,b,method=corr)
		return(r)
	}
	r.null<-c(r,replicate(nsim-1,foo(tree,V)))
	P<-mean(abs(r.null)>=abs(r))
	obj<-list(beta=beta,r=r,P=P,corr=corr,method=method)
	class(obj)<-"ratebystate"
	obj
}

# function simulates rate by state evolution for x & y
# written by Liam J. Revell 2013

sim.ratebystate<-function(tree,sig2x=1,sig2y=1,beta=c(0,1),...){
	if(hasArg(method)) method<-list(...)$method
	else method<-"by.node"
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(logarithm)) logarithm<-list(...)$logarithm
	else logarithm<-FALSE
	x<-fastBM(tree,a=if(logarithm) beta[1] else 0,sig2=sig2x,internal=TRUE)
	if(method=="by.node") ss<-x[1:tree$Nnode+length(tree$tip.label)]
	else if(method=="by.branch") ss<-rowMeans(matrix(x[tree$edge],nrow(tree$edge),2))
	zz<-tree
	if(!logarithm) zz$edge.length<-beta[2]*zz$edge.length*(beta[1]+ss-min(ss))
	else zz$edge.length<-beta[2]*zz$edge.length*exp(ss)
	y<-fastBM(zz,sig2=sig2y)
	if(plot) phenogram(zz,x,type="b",colors="blue",ftype="off",
		xlab="expected variance",ylab="independent variable (x)")
	x<-x[tree$tip.label]
	return(cbind(x,y))
}
	 
## S3 print method
print.ratebystate<-function(x,digits=6,...){
	cat("\nObject of class \"ratebystate\".\n")
	cat("\nSummary of object:\n")
	cat(paste("  beta[1] = ",round(x$beta,digits),"\n",sep=""))
	cat(paste("  ",if(x$corr=="pearson") "Pearson " else
		"Spearman ","correlation (r) = ",round(x$r,digits),
		"\n",sep=""))
	cat(paste("  P-value (from simulation) = ",round(x$P,digits),
		"\n\n",sep=""))
	cat(paste("Analysis was conducted using \"",x$method,
		"\" method.\n\n",sep=""))
}
