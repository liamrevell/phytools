# this function calls various MCMC methods
# written by Liam J. Revell 2011, 2015, 2017

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
	obj<-list(mcmc=as.data.frame(X),model=model,method=method,tree=tree)
	class(obj)<-"fitBayes"
	obj
}

## S3 methods for the object class

print.fitBayes<-function(x,digits=6,...){
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-0.2*max(x$mcmc$gen)
	ii<-which(((x$mcmc$gen-burnin)^2)==min((x$mcmc$gen-burnin)^2))+1
	cat("\nObject of class \"fitBayes\" which consists of a posterior sample")
	cat("\n(in component \'mcmc\') from a Bayesian MCMC of the model presented")
	cat("\nin Revell & Reynolds (2012; Evolution).\n\n")
	cat("Object summary:\n")
	cat(paste("\tgenerations of MCMC: ",max(x$mcmc$gen),".\n",sep=""))
	cat(paste("\tsample interval: ",x$mcmc$gen[3]-x$mcmc$gen[2],".\n",sep=""))
	cat(paste("\tmean sigma^2 from posterior sample: ",
		round(mean(x$mcmc$sig2[ii:nrow(x$mcmc)]),6),".\n",sep=""))
	if(x$model=="lambda") cat(paste("\tmean lambda from posterior sample: ",
		round(mean(x$mcmc$lambda[ii:nrow(x$mcmc)]),6),".\n",sep=""))
	cat(paste("\nCalculations based on burn-in of",burnin,"generations.\n"))
	cat("\n")
}

plot.fitBayes<-function(x,...){
	args<-list(...)
	if(is.null(args$what)) what<-"logLik"
	else {
		what<-args$what
		args$what<-NULL
	}
	if(is.null(args$burnin)) burnin<-0.2*max(x$mcmc$gen)
	else {
		burnin<-args$burnin
		args$burnin<-NULL
	}
	if(what=="logLik"){
		args$x<-x$mcmc$gen
		args$y<-x$mcmc$logLik
		if(is.null(args$xlab)) args$xlab<-"generation"
		if(is.null(args$ylab)) args$ylab<-"log(L)"
		if(is.null(args$type)) args$type<-"l"
		if(is.null(args$col)) args$col<-make.transparent("blue",0.5)
		do.call(plot,args)
	} else if(what=="sig2"){
		ii<-which(((x$mcmc$gen-burnin)^2)==min((x$mcmc$gen-burnin)^2))+1
		if(is.null(args$bw)) bw<-0.05*diff(range(x$mcmc$sig2[ii:nrow(x$mcmc)]))
		else {
			bw<-args$bw
			args$bw<-NULL
		}
		d<-density(x$mcmc$sig2[ii:nrow(x$mcmc)],bw=bw)
		args$x<-d$x
		args$y<-d$y
		if(is.null(args$xlab)) args$xlab<-expression(paste("posterior distribution of ",sigma^2))
		if(is.null(args$ylab)) args$ylab<-"density"
		if(is.null(args$main)) args$main<-""
		if(is.null(args$type)) args$type<-"l"
		if(is.null(args$col)) args$col<-"blue"
		do.call(plot,args)
		polygon(x=c(min(d$x),d$x,max(d$x)),y=c(0,d$y,0),
			col=make.transparent(args$col,0.2),border=NA)
	} else if(what=="lambda"){
		if(x$model!="lambda") stop("Model of \"fitBayes\" object is not \"lambda\".")
		else {
			ii<-which(((x$mcmc$gen-burnin)^2)==min((x$mcmc$gen-burnin)^2))+1
			if(is.null(args$bw)) bw<-0.05*diff(range(x$mcmc$lambda[ii:nrow(x$mcmc)]))
			else {
				bw<-args$bw
				args$bw<-NULL
			}
			d<-density(x$mcmc$lambda[ii:nrow(x$mcmc)],bw=bw)
			args$x<-d$x
			args$y<-d$y
			if(is.null(args$xlab)) args$xlab<-expression(paste("posterior distribution of ",lambda))
			if(is.null(args$ylab)) args$ylab<-"density"
			if(is.null(args$main)) args$main<-""
			if(is.null(args$type)) args$type<-"l"
			if(is.null(args$col)) args$col<-"blue"
			do.call(plot,args)
			polygon(x=c(min(d$x),d$x,max(d$x)),y=c(0,d$y,0),
				col=make.transparent(args$col,0.2),border=NA)
		}
	}
}

