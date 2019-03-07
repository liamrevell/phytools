## fitpolyMk 
## written by Liam J. Revell 2019

fitpolyMk<-function(tree,x,model="SYM",ordered=FALSE,...){
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(is.factor(x)) x<-setNames(as.character(x),names(x))
	X<-strsplit(x,"+",fixed=TRUE)
	ns<-sapply(X,length)
	if(all(ns==1)){
		cat("No polymorphic species found. Use fitMk.\n\n")
		object<-NULL
	} else {
		## get the states
		states<-sort(unique(unlist(X)))
		if(ordered){
			ss<-vector()
			for(i in 1:(length(states)-1)) 
				ss<-c(ss,states[i],paste(states[i],states[i+1],sep="+"))
			ss<-c(ss,states[i+1])
			tmodel<-matrix(0,length(ss),length(ss),dimnames=list(ss,ss))
		} else {
			ss<-vector()
			for(i in 1:length(states))
				ss<-c(ss,apply(combinations(length(states),i,states),
					1,paste,collapse="+"))
			tmodel<-matrix(0,length(ss),length(ss),dimnames=list(ss,ss))
		}
		poly<-strsplit(ss,"+",fixed=TRUE)
		index<-0
		for(i in 1:nrow(tmodel)){
			for(j in i:ncol(tmodel)){
				INT<-intersect(poly[[i]],poly[[j]])
				SDij<-setdiff(poly[[i]],poly[[j]])
				SDji<-setdiff(poly[[j]],poly[[i]])
				if(length(INT)>0
					&&(0%in%c(length(SDij),length(SDji)))
					&&(1%in%c(length(SDij),length(SDji)))){					
					if(model=="ER"){ 
						tmodel[i,j]<-tmodel[j,i]<-1
					} else if(model%in%c("ARD","SYM")){
						index<-index+1
						tmodel[i,j]<-index
						if(model=="SYM") tmodel[j,i]<-index
						else {
							tmodel[j,i]<-index+1
							index<-index+1
						}
					} else if(model=="transient") {
						if(length(poly[[i]])>length(poly[[j]])){ 
							tmodel[i,j]<-1
						 	tmodel[j,i]<-2
						} else {
							tmodel[i,j]<-2
							tmodel[j,i]<-1
						}
					}
				}
			}
		}
		if(!quiet){
			cat("\nThis is the design matrix of the fitted model. Does it make sense?\n\n")
			print(tmodel)
			cat("\n")
		}
		X<-to.matrix(x,ss)
		object<-fitMk(tree,X,model=tmodel,...)
	}
	object$model<-model
	object$ordered<-ordered
	class(object)<-"fitpolyMk"
	object
}

## print method for objects of class "fitpolyMk"
print.fitpolyMk<-function(x,digits=6,...){
	cat("Object of class \"fitpolyMk\".\n\n")
	if(x$ordered){ 
		cat("Evolution was modeled as \'ordered\' (i.e., transitions are assumed\n")
		cat(paste("to occur ",x$states[1]," <-> ",x$states[2]," <-> ",x$states[3],
			" and so on) ",sep=""))
		cat(paste("using the \"",x$model,"\" model.\n\n",sep=""))
	} else {
		cat("Evolution was modeled as \'unordered\' ")
		cat(paste("using the \"",x$model,"\" model.\n\n",sep=""))
	}
	cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-x$states
	print(round(Q,digits))
	cat("\nFitted (or set) value of pi:\n")
	print(x$pi)
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",sep=""))
}

## logLik method for objects of class "fitpolyMk"
logLik.fitpolyMk<-function(object,...){ 
	lik<-object$logLik
	attr(lik,"df")<-length(object$rates)
	lik
}

## AIC method
AIC.fitpolyMk<-function(object,...,k=2){
	np<-length(object$rates)
	-2*logLik(object)+np*k
}

## S3 plot method for objects of class "fitpolyMk"
plot.fitpolyMk<-function(x,...){
	if(hasArg(signif)) signif<-list(...)$signif
	else signif<-3
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(cex.main)) cex.main<-list(...)$cex.main
	else cex.main<-1.2
	if(hasArg(cex.traits)) cex.traits<-list(...)$cex.traits
	else cex.traits<-1
	if(hasArg(cex.rates)) cex.rates<-list(...)$cex.rates
	else cex.rates<-0.6
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(1.1,1.1,3.1,1.1)
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	Q<-matrix(NA,length(x$states),length(x$states))
    	Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	spacer<-0.1
	plot.new()
	par(mar=mar)
	xylim<-c(-1.2,1.2)
	plot.window(xlim=xylim,ylim=xylim,asp=1)
	if(!is.null(main)) title(main=main,cex.main=cex.main)
	nstates<-length(x$states)
	if(x$ordered){
		step<-360/nstates
		angles<-seq(-floor(nstates/2)*step,360-ceiling(nstates/2)*step,by=step)/180*pi
		if(nstates==2) angles<-angles+pi/2
		v.x<-sin(angles)
		v.y<-cos(angles)
		for(i in 1:nstates) for(j in 1:nstates)
			if(if(!isSymmetric(Q)) i!=j else i>j){
				dx<-v.x[j]-v.x[i]
				dy<-v.y[j]-v.y[i]
				slope<-abs(dy/dx)
				shift.x<-0.02*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
				shift.y<-0.02*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
				s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
					if(isSymmetric(Q)) 0 else shift.x,
					v.y[i]+spacer*sin(atan(slope))*sign(dy)+
					if(isSymmetric(Q)) 0 else shift.y)
				e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
					if(isSymmetric(Q)) 0 else shift.x,
					v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
					if(isSymmetric(Q)) 0 else shift.y)
				if(x$index.matrix[i,j]!=0){
					if(abs(diff(c(i,j)))==1||abs(diff(c(i,j)))==(nstates-1))
						text(mean(c(s[1],e[1]))+1.5*shift.x,
							mean(c(s[2],e[2]))+1.5*shift.y,
							round(Q[i,j],signif),cex=cex.rates,
							srt=atan(dy/dx)*180/pi)
					else
						text(mean(c(s[1],e[1]))+0.3*diff(c(s[1],e[1]))+
							1.5*shift.x,
							mean(c(s[2],e[2]))+0.3*diff(c(s[2],e[2]))+
							1.5*shift.y,
							round(Q[i,j],signif),cex=cex.rates,
							srt=atan(dy/dx)*180/pi)
					arrows(s[1],s[2],e[1],e[2],length=0.05,
						code=if(isSymmetric(Q)) 3 else 2,lwd=lwd)
				}
			}
		text(v.x,v.y,x$states,cex=cex.traits,
			col=make.transparent("black",0.7))
	} else {
		Ns<-inv.ncombn(nstates)
		step.y<-2/(Ns-1)
		v.x<-v.y<-vector()
		for(i in 1:Ns){
			nc<-ncombn(Ns,i)
			v.x<-c(v.x,if(nc>1) seq(-1,1,by=2/(nc-1)) else 0)
			v.y<-c(v.y,rep(1-rep((i-1)*step.y,nc)))
		}
		for(i in 1:nstates) for(j in 1:nstates)
			if(if(!isSymmetric(Q)) i!=j else i>j){
				dx<-v.x[j]-v.x[i]
				dy<-v.y[j]-v.y[i]
				slope<-abs(dy/dx)
				shift.x<-0.02*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
				shift.y<-0.02*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
				s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
					if(isSymmetric(Q)) 0 else shift.x,
					v.y[i]+spacer*sin(atan(slope))*sign(dy)+
					if(isSymmetric(Q)) 0 else shift.y)
				e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
					if(isSymmetric(Q)) 0 else shift.x,
					v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
					if(isSymmetric(Q)) 0 else shift.y)
				if(x$index.matrix[i,j]!=0){
					if(abs(diff(c(i,j)))==1||abs(diff(c(i,j)))==(nstates-1))
						text(mean(c(s[1],e[1]))+1.5*shift.x,
						mean(c(s[2],e[2]))+1.5*shift.y,
						round(Q[i,j],signif),cex=cex.rates,
						srt=atan(dy/dx)*180/pi)
				else
					text(mean(c(s[1],e[1]))+0.3*diff(c(s[1],e[1]))+
						1.5*shift.x,
						mean(c(s[2],e[2]))+0.3*diff(c(s[2],e[2]))+
						1.5*shift.y,
						round(Q[i,j],signif),cex=cex.rates,
						srt=atan(dy/dx)*180/pi)
				arrows(s[1],s[2],e[1],e[2],length=0.05,
					code=if(isSymmetric(Q)) 3 else 2,lwd=lwd)
			}
		}
		text(v.x,v.y,x$states,cex=cex.traits,
			col=make.transparent("black",0.7))

	}
}

## internally used functions to either calculate the number of combinations
## of r elements of n or calculate the number of elements n from the sum of
## all combinations from 1:r of n elements
ncombn<-function(n,r) factorial(n)/(factorial(n-r)*factorial(r))
inv.ncombn<-function(N){
	n<-Nc<-1
	while(Nc!=N){
		Nc<-0
		for(r in 1:n) Nc<-Nc+factorial(n)/(factorial(n-r)*factorial(r))
		n<-n+1
	}
	return(n-1)
}
