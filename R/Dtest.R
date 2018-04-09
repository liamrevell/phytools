calcD<-function(t1,t2){
	Do<-mapply(Map.Overlap,t1,t2,SIMPLIFY=FALSE)
	foo<-function(M){
		m1<-rowSums(M)
		m2<-colSums(M)
		as.matrix(m1)%*%t(m2)
	}
	De<-lapply(Do,foo)
	Dij<-mapply("-",Do,De,SIMPLIFY=FALSE)
	D<-sapply(Dij,function(x) sum(abs(x)))
	E_DX<-mean(D)
	E_Dij<-Reduce("+",Dij)/length(t1)
	list(E_DX=E_DX,E_Dij=E_Dij)
}

Dtest<-function(t1,t2,nsim=100,...){
	levs1<-sort(unique(as.vector(getStates(t1,"tips"))))
	k1<-length(levs1)
	levs2<-sort(unique(as.vector(getStates(t2,"tips"))))
	k2<-length(levs2)
	nrep<-length(t1)
	obj<-calcD(t1,t2)
	E_DX<-obj$E_DX
	E_Dij<-obj$E_Dij
	## posterior prediction
	PD<-0
	Pdij<-matrix(0,k1,k2,dimnames=list(levs1,levs2))
	cat("|")
	for(i in 1:nrep){
		x<-to.matrix(sim.Mk(t1[[i]],t1[[i]]$Q),levs1)
		y<-to.matrix(sim.Mk(t2[[i]],t2[[i]]$Q),levs2)
		T1<-make.simmap(t1[[i]],x,nsim=nsim,message=FALSE)
		T2<-make.simmap(t2[[i]],y,nsim=nsim,message=FALSE)
		tmp<-calcD(T1,T2)
		PD<-PD+(tmp$E_DX>=E_DX)/nrep
		Pdij<-Pdij+(tmp$E_Dij>=E_Dij)/nrep
		cat(".")
		if(i%%10==0) cat("\n")
		dev.flush()
	}
	list("E(D|X)"=E_DX,"P(D)"=PD,"E(Dij)"=E_Dij,
		"P(Dij)"=Pdij)
}
