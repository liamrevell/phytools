## written by Liam J. Revell

mcmcMk<-function(tree,x,model="ER",ngen=10000,...){
	log.prior<-function(x,prior) sum(dexp(x,prior,log=TRUE))
	proposal<-function(q,pv) abs(q+rnorm(n=length(q),sd=sqrt(pv)))
	makeQ<-function(m,q,index.matrix){
		Q<-matrix(0,m,m)
		Q[]<-c(0,q)[index.matrix+1]
		diag(Q)<-0
		diag(Q)<--rowSums(Q)
		Q
	}		
	if(hasArg(q)) q<-list(...)$q
	else q<-if(is.matrix(x)) nrow(x)/sum(tree$edge.length) else
		length(unique(x))/sum(tree$edge.length)
	if(hasArg(prop.var)) prop.var<-list(...)$prop.var
	else prop.var<-if(is.matrix(x)) sqrt(0.01*ncow(x)/sum(tree$edge.length)) else
		sqrt(0.01*length(unique(x))/sum(tree$edge.length))
	if(hasArg(prior.rate)) prior.rate<-list(...)$prior.rate
	else prior.rate<-if(is.matrix(x)) ncol(x)/sum(tree$edge.length) else
		length(unique(x))/sum(tree$edge.length)
	if(hasArg(print)) print<-list(...)$print
	else print<-100
	if(is.matrix(x)){
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	} else {
		x<-to.matrix(x,sort(unique(x)))
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	}
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-"equal"
	if(pi[1]=="equal") pi<-setNames(rep(1/m,m),states)
	else pi<-pi/sum(pi)
	if(is.character(model)){
		rate<-matrix(NA,m,m)
		if(model=="ER"){ 
			k<-rate[]<-1
			diag(rate)<-NA
		} else if(model=="ARD"){
			k<-m*(m-1)
			rate[col(rate)!=row(rate)]<-1:k
		} else if(model=="SYM"){
			k<-m*(m-1)/2
			ii<-col(rate)<row(rate)
			rate[ii]<-1:k
			rate<-t(rate)
			rate[ii]<-1:k
		}
	} else {
		if(ncol(model)!=nrow(model)) 
			stop("model is not a square matrix")
		if(ncol(model)!=ncol(x)) 
			stop("model does not have the right number of columns")
		rate<-model
		k<-max(rate)
	}
	index.matrix<-rate
	if(length(q)!=k) q<-rep(q[1],k)
	if(length(prop.var)!=k) prop.var<-rep(prop.var[1],k)
	if(length(prior.rate)!=k) prior.rate<-rep(prior.rate[1],k)
	likQ<-logLik(fitMk(tree,x,model=model,fixedQ=makeQ(m,q,index.matrix),pi=pi))
	nn<-vector(length=k,mode="character")
	for(i in 1:k) nn[i]<-paste("[",paste(states[which(rate==i,arr.ind=TRUE)[1,]],
		collapse=","),"]",sep="")
	PS<-matrix(NA,ngen,k+2,dimnames=list(1:ngen,c("gen",nn,"logLik")))
	PS[1,]<-c(1,q,likQ)
	cat("Running MCMC....\n")
	if(print){
		cat(paste(colnames(PS),collapse=" \t"))
		cat("\n")
		i<-1
		cat(paste(round(PS[i,1]),paste(round(PS[i,1:k+1],4),collapse="\t"),
			round(PS[i,ncol(PS)],4),sep="\t"))
		cat("\n")
	}
	flush.console()
	qp<-q
	for(i in 2:ngen){
		qp[i%%k+1]<-proposal(q[i%%k+1],prop.var[i%%k+1])
		likQp<-logLik(fitMk(tree,x,model=model,fixedQ=makeQ(m,qp,index.matrix),pi=pi))
		por<-exp(likQp-likQ+log.prior(qp,prior.rate)-log.prior(q,prior.rate))
		if(por>runif(n=1)){
			q<-qp
			likQ<-likQp
		}
		PS[i,]<-c(i,q,likQ)
		if(print) if(i%%print==1){
			cat(paste(round(PS[i,1]),paste(round(PS[i,1:k+1],4),collapse="\t"),
				round(PS[i,ncol(PS)],4),sep="\t"))
			cat("\n")
			flush.console()
		}
	}
	cat("Done.\n")
	class(PS)<-"mcmcMk"
	attr(PS,"model")<-model
	attr(PS,"index.matrix")<-index.matrix
	attr(PS,"states")<-states
	PS
}	

print.mcmcMk<-function(x,...){
	cat("\nPosterior sample from mcmcMk consisting of a posterior sample obtained using\n")
	cat("Bayesian MCMC in the form of a matrix.\n\n")
	cat("1. plot(\'object_name\') will create a likelihood profile plot.\n")
	cat("2. summary(\'object_name\') will compute a summary of the MCMC.\n")
	cat("3. density(\'object_name\') will calculate a posterior density from the sample.\n")
	cat("4. Finally, plot(density(\'object_name\')) will plot the posterior density and\n")
	cat("   and high probability density intervals.\n\n")
	cat("To work best, we recommend users install the package \'code\'.\n\n")
}

plot.mcmcMk<-function(x,...){
	plot(x[,"gen"],x[,"logLik"],type="s",xlab="generation",ylab="log(L)",
		col="grey",bty="l",main="Likelihood profile from MCMC",
		font.main=1)
}

summary.mcmcMk<-function(object,...){
	makeQ<-function(m,q,index.matrix){
		Q<-matrix(0,m,m)
		Q[]<-c(0,q)[index.matrix+1]
		diag(Q)<-0
		diag(Q)<--rowSums(Q)
		Q
	}	
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else { 
		burnin<-floor(0.2*nrow(object))
		cat("Assuming 20% burn-in as no burn-in was specified....\n\n")
	}
	ii<-which((object[,"gen"]-burnin)^2==min((object[,"gen"]-burnin)^2))
	PD<-object[(ii+1):nrow(object),]
	## compute the average value of Q
	Q<-makeQ(length(attr(object,"states")),
		colMeans(PD[,2:(ncol(object)-1),drop=FALSE]),
		attr(object,"index.matrix"))
	colnames(Q)<-rownames(Q)<-attr(object,"states")
	## print Q
	cat("Mean value of Q from the post burn-in posterior sample:\n")
	print(Q)
	if(.check.pkg("coda")){
			hpd95<-function(x) HPDinterval(as.mcmc(x))
	} else {
		cat("  HPDinterval requires package coda.\n")
		cat("  Computing 95% interval from samples only.\n\n")
		hpd95<-function(x){
			obj<-setNames(c(sort(x)[round(0.025*length(x))],
				sort(x)[round(0.975*length(x))]),
				c("lower","upper"))
			attr(obj,"Probability")<-0.95
			obj
		}
	}
	HPD<-t(apply(PD[,2:(ncol(object)-1),drop=FALSE],2,hpd95))
	colnames(HPD)<-c("lower","upper")
	cat("\n95% HPD interval computed either from the post burn-in\nsamples or using \'coda\':\n")
	print(HPD)
	cat("\n")
	object<-list(mean.Q=Q,HPD95=HPD)
	class(object)<-"summary.mcmcMk"
	invisible(object)
}

density.mcmcMk<-function(x,...){
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else { 
		burnin<-floor(0.2*nrow(x))
		cat("Assuming 20% burn-in as no burn-in was specified....\n")
	}
	ii<-which((x[,"gen"]-burnin)^2==min((x[,"gen"]-burnin)^2))
	PD<-x[(ii+1):nrow(x),]
	if(hasArg(bw)) bw<-list(...)$bw
	else bw<-rep("nrd0",ncol(x)-2)
	if(length(bw)==1) bw<-rep(bw,ncol(x)-2)
	d<-list()
	for(i in 2:(ncol(PD)-1))
		d[[i-1]]<-density(PD[,i],bw=bw[i-1])
	class(d)<-"density.mcmcMk"
	names(d)<-paste("Density",colnames(x)[2:(ncol(x)-1)])
	attr(d,"model")<-attr(x,"model")
	attr(d,"index.matrix")<-attr(x,"index.matrix")
	attr(d,"states")<-attr(x,"states")
	nulo<-capture.output(attr(d,"summary")<-summary(x,...))
	d
}

print.density.mcmcMk<-function(x, digits=NULL, ...){
	for(i in 1:length(x)){
		cat(paste("\n",names(x)[1],"\n\n"))
		print(summary(as.data.frame(x[[1]][c("x", "y")])), 
			digits = digits, ...)
	}
	cat("\n")
	cat("To plot enter plot(\'object_name\') at the command line interface.\n\n")
	invisible(x)
}

plot.density.mcmcMk<-function(x,...){
	if(hasArg(show.matrix)) show.matrix<-list(...)$show.matrix
	else show.matrix<-FALSE
	if(length(x)>1){
		xlim<-range(sapply(x,function(x) x$x))
		ylim<-c(0,1.1*max(sapply(x,function(x) x$y)))
	}
	if(length(x)==1){
		plot(x[[1]],main="estimated posterior density for q",
			bty="l",font.main=1,ylim=c(0,1.05*max(x[[1]]$y)),
			ylab="q")
		polygon(x[[1]],col=make.transparent("blue",0.5))
		lines(x=attr(x,"summary")$HPD95[1,],y=rep(1.01*max(x[[1]]$y),2))
		text(x=mean(attr(x,"summary")$HPD95[1,]),
			y=1.01*max(x[[1]]$y),"95% HPD",pos=3)
	} else if(length(x)==2&&show.matrix==FALSE){
		plot(x[[1]],xlim=xlim,ylim=ylim,
			main=expression(paste("estimated posterior density for ",
			Q[ij])),xlab="q",bty="l")
		polygon(x[[1]],col=make.transparent("blue",0.25))
		lines(x=attr(x,"summary")$HPD95[1,],
			y=rep(max(x[[1]]$y),2)+0.01*diff(ylim))
		text(x=mean(attr(x,"summary")$HPD95[1,]),
			y=max(x[[1]]$y)+0.01*diff(ylim),paste("95% HPD",
			rownames(attr(x,"summary")$HPD95)[1]),pos=3)
		lines(x[[2]])
		polygon(x[[2]],col=make.transparent("red",0.25))
		lines(x=attr(x,"summary")$HPD95[2,],
			y=rep(max(x[[2]]$y),2)+0.01*diff(ylim))
		text(x=mean(attr(x,"summary")$HPD95[2,]),
			y=max(x[[2]]$y)+0.01*diff(ylim),paste("95% HPD",
			rownames(attr(x,"summary")$HPD95)[2]),pos=3)
		nam.legend<-sapply(strsplit(names(x),"Density "),function(x) x[2])
		legend(x="topright",legend=nam.legend,bty="n",pt.cex=2,pch=22,
			pt.bg=c(make.transparent("blue",0.25),make.transparent("red",0.25)))
	} else {
		k<-length(attr(d,"states"))
		par(mfrow=c(k,k))
		NAMES<-sapply(strsplit(names(x),"Density "),function(x) x[2])
		NAMES<-sapply(strsplit(NAMES,""),function(x) paste(x[2:(length(x)-1)],
			collapse=""))
		s.i<-sapply(strsplit(NAMES,","),function(x) x[1])
		s.j<-sapply(strsplit(NAMES,","),function(x) x[2])
		for(i in 1:k){
			for(j in 1:k){
				ii<-intersect(which(s.i==attr(x,"states")[i]),
					which(s.j==attr(x,"states")[j]))
				if(length(ii)==0) { 
					plot(c(0,1),c(1,0),type="l",lty="dotted",xaxt="n",yaxt="n",
						xlab="",ylab="")
				} else {
					plot(x[[ii]],xlim=xlim,ylim=ylim,
						main=paste("estimated posterior density for Q[",
						s.i[ii],",",s.j[ii],"]",sep=""),xlab="q",bty="l",
						font.main=1,cex.main=1,mar=c(4.1,4.1,3.1,1.1))
					polygon(x[[ii]],col=make.transparent("blue",0.25))
					lines(x=attr(x,"summary")$HPD95[ii,],
						y=rep(max(x[[ii]]$y),2)+0.01*diff(ylim))
					text(x=mean(attr(x,"summary")$HPD95[ii,]),
						y=max(x[[ii]]$y)+0.01*diff(ylim),"95% HPD",pos=3)
				}
			}
		}
	}
}