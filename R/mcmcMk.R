## written by Liam J. Revell

mtom<-function(model){
	model<-t(model)
	im<-matrix(NA,nrow(model),ncol(model))
	indices<-setNames(1:length(unique(as.vector(model)))-1,
		unique(as.vector(model)))
	for(i in 1:length(as.vector(model)))
		im[i]<-indices[as.character(model[i])]
	t(im)
}

minChanges<-function(tree,x){
	if(is.matrix(x)) 
		x<-as.factor(apply(x,1,function(x) names(which(x==max(x)))[1]))
	parsimony(tree,as.phyDat(x))
}

makeq<-function(Q,index.matrix){
	q<-vector(length=max(index.matrix,na.rm=TRUE),mode="numeric")
	for(i in 1:length(q)) q[i]<-Q[match(i,index.matrix)]
	q
}
	
mcmcMk<-function(tree,x,model="ER",ngen=10000,...){
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	log.prior<-function(x,prior) sum(dgamma(x,shape=prior$alpha,rate=prior$beta,log=TRUE))
	proposal<-function(q,pv) abs(q+rnorm(n=length(q),sd=sqrt(pv)))	
	if(hasArg(q)) q<-list(...)$q
	else q<-minChanges(tree,x)/sum(tree$edge.length)
	if(hasArg(prop.var)) prop.var<-list(...)$prop.var
	else prop.var<-1/max(nodeHeights(tree))
	if(hasArg(auto.tune)) auto.tune<-list(...)$auto.tune
	else auto.tune<-TRUE
	if(is.numeric(auto.tune)){ 
		target.accept<-auto.tune
		auto.tune<-TRUE
	} else target.accept<-0.5
	if(hasArg(prior)) prior<-list(...)$prior
	else prior<-NULL
	if(hasArg(print)) print<-list(...)$print
	else print<-100
	if(hasArg(likelihood)) likelihood<-list(...)$likelihood
	else {
		if(.check.pkg("geiger")) likelihood<-"fitDiscrete"
		else likelihood<-"fitMk"
	}
	if(likelihood=="fitDiscrete"){
		if(.check.pkg("geiger")==FALSE){
			cat("geiger is not installed. Setting likelihood method to \"fitMk\".\n\n")
			likelihood<-"fitMk"
			fitDiscrete<-function(...) NULL
		} else if(is.matrix(x)){
			cat("likelihood=\"fitDiscrete\" doesn't work for data input as matrix.\n")
			cat("Setting likelihood method to \"fitMk\".\n\n")
			fitDiscrete<-function(...) NULL
		}
	}
	if(is.matrix(x)&&likelihood=="fitDiscrete") likelihood<-"fitMk"
	if(is.matrix(x)){
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	} else {
		y<-x
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
		rate<-t(rate)
	} else {
		if(ncol(model)!=nrow(model)) 
			stop("model is not a square matrix")
		if(ncol(model)!=ncol(x)) 
			stop("model does not have the right number of columns")
		model<-mtom(model)
		rate<-model
		k<-max(rate)
	}
	index.matrix<-rate
	if(length(q)!=k) q<-rep(q[1],k)
	if(length(prop.var)!=k) prop.var<-rep(prop.var[1],k)
	default.alpha<-0.1
	if(is.null(prior)){
		prior<-list(
			alpha=rep(default.alpha,k),
			beta=rep(default.alpha*sum(tree$edge.length)/minChanges(tree,x),k)
		)
	} else if(is.list(prior)){
		if(is.null(prior$alpha)) prior$alpha<-rep(default.alpha,k)
		else if(is.null(prior$beta)) 
			prior$beta<-rep(prior$alpha*sum(tree$edge.length)/minChanges(tree,x),k)
	} else if(is.numeric(prior)){
		if(length(prior)==1)
			prior<-list(alpha=rep(prior,k),
				beta=rep(prior*sum(tree$edge.length)/minChanges(tree,x),k))
		else if(length(prior)==2)
			prior<-list(alpha=rep(prior[1],k),beta=rep(prior[2],k))
	}
	if(likelihood=="fitDiscrete"){
		LIK<-fitDiscrete(tree,y,model=model,niter=1)$lik
		lik.func<-function(Q) LIK(makeq(Q,index.matrix),root="given",root.p=pi)
	} else lik.func<-fitMk(tree,x,fixedQ=makeQ(m,q,index.matrix),pi=pi)$lik
	likQ<-lik.func(makeQ(m,q,index.matrix))
	nn<-vector(length=k,mode="character")
	for(i in 1:k) nn[i]<-paste("[",paste(states[which(rate==i,arr.ind=TRUE)[1,]],
		collapse=","),"]",sep="")
	PS<-matrix(NA,ngen,k+2,dimnames=list(1:ngen,c("gen",nn,"logLik")))
	PS[1,]<-c(1,q,likQ)
	cat("Running MCMC....\n")
	accept<-0
	if(print){
		cat(paste(paste(colnames(PS),collapse=" \t"),"\taccept\n",sep=""))
		i<-1
		cat(paste(round(PS[i,1]),paste(round(PS[i,1:k+1],4),collapse="\t"),
			round(PS[i,ncol(PS)],4),sep="\t"))
		cat("\n")
	}
	flush.console()
	qp<-q
	for(i in 2:ngen){
		qp[i%%k+1]<-proposal(q[i%%k+1],prop.var[i%%k+1])
		likQp<-lik.func(makeQ(m,qp,index.matrix))
		por<-exp(likQp-likQ+log.prior(qp,prior)-log.prior(q,prior))
		if(por>runif(n=1)){
			q<-qp
			likQ<-likQp
			accept<-accept+1/100
		}
		PS[i,]<-c(i,q,likQ)
		if(print) if(i%%print==1){
			cat(paste(round(PS[i,1]),paste(round(PS[i,1:k+1],4),collapse="\t"),
				round(PS[i,ncol(PS)],4),round(accept,3),sep="\t"))
			cat("\n")
			flush.console()
		} 
		if(i%%100==1){
			if(auto.tune)
				prop.var<-if(accept>target.accept) 1.1*prop.var else prop.var/1.1
			accept<-0
		}
		if(plot){
			dev.hold()
			par(mfrow=c(2,1),mar=c(5.1,4.1,2.1,1.1))
			plot(1:i,PS[1:i,"logLik"],col="darkgrey",xlab="",
				ylab="log(L)",xlim=c(0,ngen),type="l",bty="l")
			mtext("a) log-likelihood profile plot",side=3,line=1,cex=1,
				at=0,outer=FALSE,adj=0)
			plot(1:i,PS[1:i,2],col=Palette(1),
				xlab="generation",xlim=c(0,ngen),
				ylim=c(0,max(PS[1:i,2:(ncol(PS)-1)])),
				ylab="q",type="l",bty="l")
			abline(h=mean(PS[1:i,2]),lty="dotted",col=Palette(1))
			text(x=ngen,y=mean(PS[1:i,2]),colnames(PS)[2],cex=0.5,
				col=Palette(1),pos=3)
			if(length(q)>1){
				for(j in 2:length(q)){ 
					lines(1:i,PS[1:i,j+1],col=Palette(j))
					abline(h=mean(PS[1:i,j+1]),lty="dotted",col=Palette(j))
					text(x=ngen,y=mean(PS[1:i,j+1]),colnames(PS)[j+1],
						cex=0.5,col=Palette(j),pos=3)
				}
			}
			mtext("b) rates",side=3,line=1,cex=1,at=0,outer=FALSE,adj=0)
			dev.flush()
		}
	}
	cat("Done.\n")
	class(PS)<-"mcmcMk"
	attr(PS,"model")<-model
	attr(PS,"index.matrix")<-index.matrix
	attr(PS,"states")<-states
	PS
}

Palette<-function(i){
	if(!.check.pkg("RColorBrewer")){
		brewer.pal<-function(...) NULL
		COLOR<-rep(palette(),ceiling(i/8))[i]
		COLOR<-if(COLOR=="black") "darkgrey"
	} else COLOR<-rep(brewer.pal(8,"Accent"),ceiling(i/8))[i]
	COLOR
}

print.mcmcMk<-function(x,...){
	cat("\nPosterior sample from mcmcMk consisting of a posterior sample obtained using\n")
	cat("Bayesian MCMC in the form of a matrix.\n\n")
	cat("1. plot(\'object_name\') will create a likelihood profile plot.\n")
	cat("2. summary(\'object_name\') will compute a summary of the MCMC.\n")
	cat("3. density(\'object_name\') will calculate a posterior density from the sample.\n")
	cat("4. Finally, plot(density(\'object_name\')) will plot the posterior density and\n")
	cat("   and high probability density intervals.\n\n")
	cat("To work best, we recommend users install the package \'coda\'.\n\n")
}

plot.mcmcMk<-function(x,...){
	par(mfrow=c(2,1),mar=c(5.1,4.1,2.1,1.1))
	plot(x[,"gen"],x[,"logLik"],col="darkgrey",xlab="",
		ylab="log(L)",type="l",bty="l")
	mtext("a) log-likelihood profile plot",side=3,line=1,cex=1,
		at=0,outer=FALSE,adj=0)
	plot(x[,"gen"],x[,2],col=Palette(1),xlab="generation",
		ylim=c(0,max(x[,2:(ncol(x)-1)])),ylab="q",type="l",
		bty="l")
	abline(h=mean(x[,2]),lty="dotted",col=Palette(1))
	text(x=max(x[,"gen"]),y=mean(x[,2]),colnames(x)[2],cex=0.5,
		col=Palette(1),pos=3)
	if(ncol(x)>3){
		for(j in 2:(ncol(x)-2)){ 
			lines(x[,"gen"],x[,j+1],col=Palette(j))
			abline(h=mean(x[,j+1]),lty="dotted",col=Palette(j))
			text(x=max(x[,"gen"]),y=mean(x[,j+1]),colnames(x)[j+1],
				cex=0.5,col=Palette(j),pos=3)
		}
	}
	mtext("b) rates",side=3,line=1,cex=1,at=0,outer=FALSE,adj=0)
}

summary.mcmcMk<-function(object,...){
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
	## compute the median value of Q
	median.Q<-makeQ(length(attr(object,"states")),
		apply(PD[,2:(ncol(object)-1),drop=FALSE],2,median),
		attr(object,"index.matrix"))
	colnames(median.Q)<-rownames(median.Q)<-attr(object,"states")
	## print Q
	cat("\nMedian value of Q from the post burn-in posterior sample:\n")
	print(median.Q)		
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
	object<-list(mean.Q=Q,median.Q=median.Q,HPD95=HPD)
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
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-range(sapply(x,function(x) x$x))
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(0,1.1*max(sapply(x,function(x) x$y)))
	if(length(x)==1){
		if(hasArg(main)) main<-list(...)$main
		else main<-"estimated posterior density for q"
		plot(x[[1]],main=main,
			bty="l",font.main=1,xlim=xlim,ylim=ylim,
			xlab="q",bty="l")
		polygon(x[[1]],col=make.transparent("blue",0.5))
		lines(x=attr(x,"summary")$HPD95[1,],y=rep(1.01*max(x[[1]]$y),2))
		text(x=mean(attr(x,"summary")$HPD95[1,]),
			y=1.01*max(x[[1]]$y),"95% HPD",pos=3)
	} else if(length(x)==2&&show.matrix==FALSE){
		if(hasArg(main)) main<-list(...)$main
		else main<-expression(paste("estimated posterior density for ",
			Q[ij]))
		plot(x[[1]],xlim=xlim,ylim=ylim,
			main=main,xlab="q",bty="l")
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
			pt.bg=c(make.transparent("blue",0.25),make.transparent("red",
			0.25)))
	} else {
		k<-length(attr(x,"states"))
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
