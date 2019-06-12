# 95% CI on ltts
# written by Liam J. Revell 2013, 2014, 2015, 2019

ltt95<-function(trees,alpha=0.05,log=FALSE,method=c("lineages","times"),mode=c("median","mean"),...){
	if(!inherits(trees,"multiPhylo")) stop("trees should be an object of class \"multiPhylo\".")
	method<-method[1]
	mode<-mode[1]
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	X<-ltt(trees,plot=FALSE,gamma=FALSE)
	if(method=="times"){
		N<-length(X)
		tt<-sapply(X,function(x) max(x$times))
		zz<-max(tt)-tt
		for(i in 1:N) X[[i]]$times<-X[[i]]$times+zz[i]
		n<-sapply(X,function(x) max(x$ltt))
		if(all(n==max(n))) n<-max(n) 
		else stop("for method=\"times\" all trees must contain the same number of lineages")
		LL<-sapply(X,function(x) x$times[1:length(x$times)])
		ii<-max(floor(alpha/2*N)+1,1)
		jj<-min(ceiling((1-alpha/2)*N),N)
		low<-apply(LL,1,function(x) sort(x)[ii])
		high<-apply(LL,1,function(x) sort(x)[jj])
		ll<-if(mode=="median") apply(LL,1,function(x) median(x)[1]) else rowMeans(LL)
		obj<-cbind(c(1:n,n),low,ll,high)
		colnames(obj)<-c("lineages","low(time)","time","high(time)")
		rownames(obj)<-NULL
	} else if(method=="lineages"){
		N<-length(X)
		tt<-sapply(X,function(x) max(x$times))
		zz<-max(tt)-tt
		for(i in 1:N) X[[i]]$times<-X[[i]]$times+zz[i]
		tt<-0:res*max(tt)/res
		ll<-low<-high<-vector()
		for(i in 1:length(tt)){
			ss<-vector()
			for(j in 1:N){
				ii<-2
				while(tt[i]>X[[j]]$times[ii]&&ii<length(X[[j]]$times)) ii<-ii+1
				ss[j]<-X[[j]]$ltt[ii-1]
			}
			ll[i]<-if(mode=="median") median(ss) else mean(ss)
			low[i]<-sort(ss)[max(floor(alpha/2*N)+1,1)]
			high[i]<-sort(ss)[min(ceiling((1-alpha/2)*N),N)]
		}
		obj<-cbind(tt,low,ll,high)
		colnames(obj)<-c("time","low(lineages)","lineages","high(lineages)")
		rownames(obj)<-NULL
	}
	attr(obj,"class")<-"ltt95"
	attr(obj,"alpha")<-alpha
	attr(obj,"method")<-method
	attr(obj,"mode")<-mode
	attr(obj,"log")<-log
	if(plot) plot(obj,...)
	invisible(obj)
}

## S3 plotting method for objects of class 'ltt95'
## written by Liam J. Revell 2014, 2019
plot.ltt95<-function(x,...){
	if(hasArg(lend)) lend<-list(...)$lend
	else lend<-3
	old.lend<-par()$lend
	par(lend=lend)
	if(hasArg(log)) log<-list(...)$log
	else log<-attr(x,"log")
	if(hasArg(xaxis)) xaxis<-list(...)$xaxis
	else xaxis<-"standard"
	if(hasArg(shaded)) shaded<-list(...)$shaded
	else shaded<-FALSE
	if(shaded) bg<-if(hasArg(bg)) list(...)$bg else rgb(0,0,1,0.25)
	if(attr(x,"method")=="times"){
		n<-max(x[,1])
		if(xaxis=="negative"){ 
			x[,2:4]<-x[,2:4]-max(x[,2:4])
		}
		if(xaxis=="flipped"){
			x[,2:4]<-max(x[,2:4])-x[,2:4]
			x.lim<-c(max(x[,2:4]),min(x[,2:4]))
		} else x.lim<-range(x[,2:4])
		plot(x[,3],x[,1],lwd=2,xlim=x.lim,
			type=if(attr(x,"mode")=="median") "s" else "l",main=NULL,
			xlab=if(xaxis=="standard") "time from the oldest root" else 
			if(xaxis=="negative") "time from the present" else 
			if(xaxis=="flipped") "time before the present day",
			ylab="lineages",log=if(log) "y" else "")
		if(!shaded){
			lines(x[,2],x[,1],lty="dashed",type=if(attr(x,"mode")=="median") "s" else "l")
			lines(x[,4],x[,1],lty="dashed",type=if(attr(x,"mode")=="median") "s" else "l")
		} else { 
			xx<-c(x[1,2],rbind(x[2:nrow(x),2],x[2:nrow(x),2]),
				rbind(x[nrow(x):2,4],x[nrow(x):2,4]),x[1,4])
			yy<-c(rbind(x[1:(nrow(x)-1),1],x[1:(nrow(x)-1),1]),x[nrow(x),1],
				x[nrow(x),1],rbind(x[(nrow(x)-1):1,1],x[(nrow(x)-1):1,1]))
			polygon(xx,yy,border=NA,col=bg)
			lines(x[,3],x[,1],lwd=2,
				type=if(attr(x,"mode")=="median") "s" else "l")
		}
	} else if(attr(x,"method")=="lineages"){
		if(xaxis=="negative") x[,1]<-x[,1]-max(x[,1])
		if(xaxis=="flipped"){
			x[,1]<-max(x[,1])-x[,1]
			x.lim<-c(max(x[,1]),min(x[,1]))
		} else x.lim<-range(x[,1])
		plot(x[,1],x[,3],xlim=x.lim,ylim=c(min(x[,2]),max(x[,4])),lwd=2,
			type=if(attr(x,"mode")=="median") "s" else "l",main=NULL,
			xlab=if(xaxis=="standard") "time from the oldest root" else 
			if(xaxis=="negative") "time from the present" else 
			if(xaxis=="flipped") "time before the present day",
			ylab="lineages",log=if(log) "y" else "")
		if(!shaded){
			lines(x[,1],x[,2],lty="dashed",type="s")
			lines(x[,1],x[,4],lty="dashed",type="s")
		} else { 
			xx<-c(x[1,1],rbind(x[2:nrow(x),1],x[2:nrow(x),1]),
				rbind(x[nrow(x):2,1],x[nrow(x):2,1]),x[1,1])
			yy<-c(rbind(x[1:(nrow(x)-1),2],x[1:(nrow(x)-1),2]),x[nrow(x),2],
				x[nrow(x),4],rbind(x[(nrow(x)-1):1,4],x[(nrow(x)-1):1,4]))
			polygon(xx,yy,border=NA,col=bg)
			lines(x[,1],x[,3],lwd=2,
				type=if(attr(x,"mode")=="median") "s" else "l")
		}
	}
	par(lend=old.lend)
}

## S3 print method for object of class 'ltt95'
## written by Liam J. Revell 2014
print.ltt95<-function(x,...){
	cat("Object of class \"ltt95\".\n")
	cat(paste("  alpha:\t",round(attr(x,"alpha"),3),"\n",sep=""))
	cat(paste("  method:\t",attr(x,"method"),"\n",sep=""))
	cat(paste("  mode:  \t",attr(x,"mode"),"\n",sep=""))
	n<-if(attr(x,"method")=="times") max(x[,1]) else range(apply(x[,2:4],2,max))
	if(length(n)==2&&n[1]==n[2]) n<-n[1]
	if(length(n)==1) cat(paste("  N(lineages):\t",n,"\n\n",sep=""))
	else cat(paste("  N(lineages):\t[",n[1],",",n[2],"]\n\n",sep=""))
}
