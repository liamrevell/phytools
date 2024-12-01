## fitpolyMk 
## fits several polymorphic discrete character evolution models
## written by Liam J. Revell 2019, 2020, 2022, 2023, 2024

anova.fitpolyMk<-function(object,...) anova.fitMk(object,...)

as.Qmatrix.fitpolyMk<-function(x,...){
	class(x)<-"fitMk"
	as.Qmatrix(x,...)
}

Combinations<-function(n,r,v=1:n){
	if(n!=length(v)) stop("n and v should have the same length")
	else return(t(combn(v,r)))
}

fitpolyMk<-function(tree,x,model="SYM",ordered=FALSE,...){
	if(hasArg(return_matrix)) return_matrix<-list(...)$return_matrix
	else return_matrix<-FALSE
	if(return_matrix) quiet<-TRUE
	else { 
		if(hasArg(quiet)) quiet<-list(...)$quiet
		else quiet<-FALSE
	}
	if(is.factor(x)) x<-setNames(as.character(x),names(x))
	if(is.matrix(x)) X<-strsplit(colnames(x),"+",fixed=TRUE)
	else X<-strsplit(x,"+",fixed=TRUE)
	ns<-sapply(X,length)
	## get the states
	states<-sort(unique(unlist(X)))
	## check if ordered
	if(ordered){
		if(hasArg(max.poly)) max.poly<-list(...)$max.poly
		else if(hasArg(max.states)){ 
			max.states<-list(...)$max.states
			max.poly<-max.states
		} else max.poly<-max(ns)
		## get any user-supplied ordering
		if(hasArg(order)) order<-list(...)$order
		else order<-NULL
		if(!is.null(order)){
			if(setequal(order,states)) states<-order
			else cat("order & states do not match. using alphabetical order.\n")
		}
	} else {
		if(hasArg(max.poly)) max.poly<-list(...)$max.poly
		else max.poly<-length(states)
	}
	if(all(ns==1)&&return_matrix==FALSE){
		cat("No polymorphic species found. Use fitMk.\n\n")
		object<-NULL
	} else {
		## fix the order of the input data
		if(is.matrix(x)){
			Levs<-sapply(X,function(x) paste(sort(x),collapse="+"))
			colnames(x)<-Levs
		} else x<-sapply(X,function(x) paste(sort(x),collapse="+"))
		if(ordered){
			ss<-vector()
			for(i in 1:(length(states)-1)){
				ss<-c(ss,states[i])
				for(j in (i+1):length(states)) 
					if((j-i)<max.poly) ss<-c(ss,paste(sort(states[i:j]),collapse="+"))
			}
			ss<-c(ss,states[i+1])	
			tmodel<-matrix(0,length(ss),length(ss),dimnames=list(ss,ss))
		} else {
			ss<-vector()
			for(i in 1:max.poly)
				ss<-c(ss,apply(Combinations(length(states),i,states),
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
			cat("\nThis is the design matrix of the fitted model.\nDoes it make sense?\n\n")
			print(tmodel)
			cat("\n")
			flush.console()
		}
		if(is.matrix(x)){
			X<-matrix(0,nrow(x),length(ss),dimnames=list(rownames(x),ss))
			X[rownames(x),colnames(x)]<-x
		} else X<-to.matrix(x,ss)
		if(return_matrix) return(X)
		object<-fitMk(tree,X,model=tmodel,...)
	}
	object$model<-model
	object$ordered<-ordered
	attr(object$ordered,"max.poly")<-max.poly
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
	print(round(x$pi,digits))
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",sep=""))
	if(x$opt_results$convergence==0) 
		cat("R thinks it has found the ML solution.\n\n")
	else cat("R thinks optimization may not have converged.\n\n")
}

## logLik method for objects of class "fitpolyMk"
logLik.fitpolyMk<-function(object,...){ 
	lik<-object$logLik
	attr(lik,"df")<-length(object$rates)
	lik
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
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-1
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-1.2,1.2)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(-1.2,1.2)
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-1.5
	if(hasArg(spacer)) spacer<-list(...)$spacer
	else spacer<-0.1
	Q<-matrix(NA,length(x$states),length(x$states))
    Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	if(!add) plot.new()
	par(mar=mar)
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	if(!is.null(main)) title(main=main,cex.main=cex.main)
	nstates<-length(x$states)
	if(x$ordered){
		if(attr(x$ordered,"max.poly")==2){
			step<-360/nstates
			angles<-seq(-floor(nstates/2)*step,360-ceiling(nstates/2)*step,by=step)/180*pi
			if(nstates==2) angles<-angles+pi/2
			v.x<-sin(angles)
			v.y<-cos(angles)
			for(i in 1:nstates) for(j in 1:nstates){
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
							text(mean(c(s[1],e[1]))+offset*shift.x,
								mean(c(s[2],e[2]))+offset*shift.y,
								round(Q[i,j],signif),cex=cex.rates,
								srt=atan(asp*dy/dx)*180/pi)
						else
							text(mean(c(s[1],e[1]))+0.3*diff(c(s[1],e[1]))+
								offset*shift.x,
								mean(c(s[2],e[2]))+0.3*diff(c(s[2],e[2]))+
								offset*shift.y,
								round(Q[i,j],signif),cex=cex.rates,
								srt=atan(asp*dy/dx)*180/pi)
						arrows(s[1],s[2],e[1],e[2],length=0.05,
							code=if(isSymmetric(Q)) 3 else 2,lwd=lwd)
					}
				}
			}
			text(v.x,v.y,x$states,cex=cex.traits,
					col=make.transparent(par()$fg,0.7))
		} else {
			nlevs<-attr(x$ordered,"max.poly")
			Ns<-inv.ncombn2(nstates,nlevs)
			v.x<-v.y<-vector()
			xx<-seq(-1,1,length.out=Ns)
			for(i in 1:Ns){
				for(j in 1:min(nlevs,Ns-i+1)){
					v.x<-c(v.x,mean(xx[i:(i+j-1)]))
					v.y<-c(v.y,1-2*(j-1)/(nlevs-1))
				}
			}
			for(i in 1:nstates) for(j in 1:nstates){
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
						text(mean(c(s[1],e[1]))+offset*shift.x,
							mean(c(s[2],e[2]))+offset*shift.y,
							round(Q[i,j],signif),cex=cex.rates,
							srt=atan(asp*dy/dx)*180/pi)
						arrows(s[1],s[2],e[1],e[2],length=0.05,
							code=if(isSymmetric(Q)) 3 else 2,lwd=lwd)
					}
				}
			}
			text(v.x,v.y,x$states,cex=cex.traits,
				col=make.transparent(par()$fg,0.7))
		}
	} else {
		nlevs<-attr(x$ordered,"max.poly")
		Ns<-inv.ncombn(nstates,nlevs)
		step.y<-2/(min(Ns,nlevs)-1)
		v.x<-v.y<-vector()
		for(i in 1:min(Ns,nlevs)){
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
						text(mean(c(s[1],e[1]))+offset*shift.x,
						mean(c(s[2],e[2]))+offset*shift.y,
						round(Q[i,j],signif),cex=cex.rates,
						srt=atan(asp*dy/dx)*180/pi)
				else
					text(mean(c(s[1],e[1]))+0.3*diff(c(s[1],e[1]))+
						offset*shift.x,
						mean(c(s[2],e[2]))+0.3*diff(c(s[2],e[2]))+
						offset*shift.y,
						round(Q[i,j],signif),cex=cex.rates,
						srt=atan(asp*dy/dx)*180/pi)
				arrows(s[1],s[2],e[1],e[2],length=0.05,
					code=if(isSymmetric(Q)) 3 else 2,lwd=lwd)
			}
		}
		text(v.x,v.y,x$states,cex=cex.traits,
			col=make.transparent(par()$fg,0.7))

	}
}

## internally used functions to either calculate the number of combinations
## of r elements of n or calculate the number of elements n from the sum of
## all combinations from 1:r of n elements
ncombn<-function(n,r) factorial(n)/(factorial(n-r)*factorial(r))
inv.ncombn<-function(N,m){
	Nc<-n<-0
	while(Nc!=N){
 		n<-n+1
		tmp<-sapply(1:n,ncombn,n=n)
		tmp<-if(length(tmp)<m) tmp else tmp[1:m]
		Nc<-sum(tmp)
	}
	n
}
inv.ncombn2<-function(N,m){
	n<-2
	Nc<-0
	while(Nc!=N){
		Nc<-sum((n:1)[1:min(m,n)])
		n<-n+1
	}
	return(n-1)
}

## maybe this should be plot.mkModel or something?
## then you create a model class & plot it?

graph.polyMk<-function(k=2,model="SYM",ordered=FALSE,...){
	if(hasArg(states)) states<-list(...)$states
	else states<-0:(k-1)
	if(hasArg(max.poly)) max.poly<-list(...)$max.poly
	else max.poly<-k
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(cex.main)) cex.main<-list(...)$cex.main
	else cex.main<-1.2
	if(hasArg(cex.traits)) cex.traits<-list(...)$cex.traits
	else cex.traits<-1
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(1.1,1.1,3.1,1.1)
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-1.2,1.2)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(-1.2,1.2)
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-1
	if(ordered){
		ss<-vector()
		for(i in 1:(length(states)-1)){
			ss<-c(ss,states[i])
			for(j in (i+1):length(states)) 
				if((j-i)<max.poly) ss<-c(ss,paste(states[i:j],collapse="+"))
		}
		ss<-c(ss,states[i+1])	
		tmodel<-matrix(0,length(ss),length(ss),dimnames=list(ss,ss))
	} else {
		ss<-vector()
		for(i in 1:length(states))
			ss<-c(ss,apply(Combinations(length(states),i,states),
				1,paste,collapse="+"))
		tmodel<-matrix(0,length(ss),length(ss),dimnames=list(ss,ss))
	}
	nstates<-length(ss)
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
				} else if(model=="transient"){
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
	if(plot){
		spacer<-if(hasArg(spacer)) list(...)$spacer else 0.1
		plot.new()
		par(mar=mar)
		plot.window(xlim=xlim,ylim=ylim,asp=asp)
		if(!is.null(main)) title(main=main,cex.main=cex.main)
		if(ordered){
			if(max.poly==2){
				step<-360/nstates
				angles<-seq(-floor(nstates/2)*step,360-ceiling(nstates/2)*step,by=step)/180*pi
				if(k==2) angles<-angles+pi/2
				v.x<-sin(angles)
				v.y<-cos(angles)
				for(i in 1:nstates) for(j in 1:nstates){
					if(if(!isSymmetric(tmodel)) i!=j else i>j){
						dx<-v.x[j]-v.x[i]
						dy<-v.y[j]-v.y[i]
						slope<-abs(dy/dx)
						shift.x<-0.02*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
						shift.y<-0.02*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
						s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
							if(isSymmetric(tmodel)) 0 else shift.x,
							v.y[i]+spacer*sin(atan(slope))*sign(dy)+
							if(isSymmetric(tmodel)) 0 else shift.y)
						e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
							if(isSymmetric(tmodel)) 0 else shift.x,
							v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
							if(isSymmetric(tmodel)) 0 else shift.y)
						if(tmodel[i,j]!=0){
							arrows(s[1],s[2],e[1],e[2],length=0.05,
								code=if(isSymmetric(tmodel)) 3 else 2,lwd=lwd)
						}
					}
				}
				text(v.x,v.y,colnames(tmodel),cex=cex.traits,
					col=make.transparent(par()$fg,0.7))
			} else {
				nlevs<-max.poly
				Ns<-inv.ncombn2(nstates,nlevs)
				v.x<-v.y<-vector()
				xx<-seq(-1,1,length.out=Ns)
				for(i in 1:k){
					for(j in 1:min(nlevs,Ns-i+1)){
						v.x<-c(v.x,mean(xx[i:(i+j-1)]))
						v.y<-c(v.y,1-2*(j-1)/(nlevs-1))
					}
				}
				for(i in 1:nstates) for(j in 1:nstates){
					if(if(!isSymmetric(tmodel)) i!=j else i>j){
						dx<-v.x[j]-v.x[i]
						dy<-v.y[j]-v.y[i]
						slope<-abs(dy/dx)
						shift.x<-0.02*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
						shift.y<-0.02*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
						s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
							if(isSymmetric(tmodel)) 0 else shift.x,
							v.y[i]+spacer*sin(atan(slope))*sign(dy)+
							if(isSymmetric(tmodel)) 0 else shift.y)
						e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
							if(isSymmetric(tmodel)) 0 else shift.x,
							v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
							if(isSymmetric(tmodel)) 0 else shift.y)
						if(tmodel[i,j]!=0){
							arrows(s[1],s[2],e[1],e[2],length=0.05,
								code=if(isSymmetric(tmodel)) 3 else 2,lwd=lwd)
						}
					}
				}
				text(v.x,v.y,rownames(tmodel),cex=cex.traits,
					col=make.transparent(par()$fg,0.7))
			}
		} else {
			step.y<-2/(k-1)
			v.x<-v.y<-vector()
			for(i in 1:k){
				nc<-ncombn(k,i)
				v.x<-c(v.x,if(nc>1) seq(-1,1,by=2/(nc-1)) else 0)
				v.y<-c(v.y,rep(1-rep((i-1)*step.y,nc)))
			}
			for(i in 1:ncol(tmodel)) for(j in 1:ncol(tmodel)){
				if(if(!isSymmetric(tmodel)) i!=j else i>j){
					dx<-v.x[j]-v.x[i]
					dy<-v.y[j]-v.y[i]
					slope<-abs(dy/dx)
					shift.x<-0.02*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
					shift.y<-0.02*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
					s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
						if(isSymmetric(tmodel)) 0 else shift.x,
						v.y[i]+spacer*sin(atan(slope))*sign(dy)+
						if(isSymmetric(tmodel)) 0 else shift.y)
					e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
						if(isSymmetric(tmodel)) 0 else shift.x,
						v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
						if(isSymmetric(tmodel)) 0 else shift.y)
					if(tmodel[i,j]!=0){
						arrows(s[1],s[2],e[1],e[2],length=0.05,
							code=if(isSymmetric(tmodel)) 3 else 2,lwd=lwd)
					}
				}
			}
			text(v.x,v.y,rownames(tmodel),cex=cex.traits,
				col=make.transparent(par()$fg,0.7))
		}
	}
	invisible(tmodel)
}
