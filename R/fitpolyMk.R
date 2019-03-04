## fitpolyMk written by Liam J. Revell 2019

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
				if(length(intersect(poly[[i]],poly[[j]]))>0&&
					(length(poly[[i]])!=length(poly[[j]])&&
					(length(setdiff(poly[[i]],poly[[j]]))==1||
					length(setdiff(poly[[j]],poly[[i]]))==1))){
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
