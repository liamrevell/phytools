## cospeciation method
## written by Liam J. Revell 2016

cospeciation<-function(t1,t2,distance=c("RF","SPR"),
	method=c("simulation","permutation"),assoc=NULL,
	nsim=100,...){
	distance<-distance[1]
	if(!distance%in%c("RF","SPR")) distance<-"RF"
	method<-method[1]
	if(!method%in%c("simulation","permutation")) 
		method<-"simulation"
	if(is.null(assoc)){
		## assume exact match
		tips<-intersect(t1$tip.label,t2$tip.label)
		assoc<-cbind(tips,tips)
	}
	if(any(!t1$tip.label%in%assoc[,1])) 
		t1<-drop.tip(t1,setdiff(t1$tip.label,assoc[,1]))
	if(any(!assoc[,1]%in%t1$tip.label))
		assoc<-assoc[assoc[,1]%in%t1$tip.label,]
	if(any(!t2$tip.label%in%assoc[,2])) 
		t2<-drop.tip(t2,setdiff(t2$tip.label,assoc[,2]))
	if(any(!assoc[,2]%in%t2$tip.label))
		assoc<-assoc[assoc[,2]%in%t2$tip.label,]	
	if(method=="permutation"){
		perm.labels<-function(tree){
			tree$tip.label<-sample(tree$tip.label)
			tree
		}
		tt1<-replicate(nsim,perm.labels(t1),simplify=FALSE)
		swap.t2<-t2
		swap.t2$tip.label<-sapply(t2$tip.label,function(x,y)
			y[which(y[,2]==x),1],y=assoc)
		tt2<-replicate(nsim,perm.labels(swap.t2),simplify=FALSE)
	} else {
		tt1<-pbtree(n=Ntip(t1),tip.label=t1$tip.label,
			nsim=nsim)
		swap.t2<-t2
		swap.t2$tip.label<-sapply(t2$tip.label,function(x,y)
			y[which(y[,2]==x),1],y=assoc)
		tt2<-pbtree(n=Ntip(t2),tip.label=swap.t2$tip.label,
			nsim=nsim)
	
	}
	if(distance=="SPR"){
		d.null<-mapply(SPR.dist,tt1,tt2)
		d<-SPR.dist(t1,swap.t2)
		P.val<-mean(c(d,d.null)<=d)
	} else {
		d.null<-mapply(RF.dist,tt1,tt2)
		d<-RF.dist(t1,swap.t2)
		P.val<-mean(c(d,d.null)<=d)
	}
	obj<-list(d=d,d.null=d.null,P.val=P.val,
		distance=if(distance=="SPR") "SPR" else "RF",
		method=if(method=="simulation") "simulation" else
		"permutation")
	class(obj)<-"cospeciation"
	obj
}

print.cospeciation<-function(x,...){
	cat(paste("\nCo-speciation test based on",x$distance,
		"distance.\n"))
	cat(paste("P-value obtained via",x$method,".\n\n"))
	if(x$distance=="SPR")
		cat(paste("   SPR distance:",x$d,"\n"))
	else
		cat(paste("   RF distance:",x$d,"\n"))
	cat(paste("   Mean(SD) from null: ",
		round(mean(x$d.null),1),"(",
		round(sd(x$d.null),1),")\n",sep=""))
	cat(paste("   P-value:",round(x$P.val,6),"\n\n"))
}

plot.cospeciation<-function(x,...){
	if(x$distance=="RF")
		p<-hist(x$d.null,breaks=seq(min(c(x$d,x$d.null))-3,
			max(c(x$d,x$d.null))+3,2),plot=FALSE)
	else if(x$distance=="SPR")
		p<-hist(x$d.null,breaks=seq(min(c(x$d,x$d.null))-1.5,
			max(c(x$d,x$d.null))+1.5,1),plot=FALSE)
	plot(p$mids,p$density,xlim=c(min(c(x$d,x$d.null))-2,
		max(c(x$d,x$d.null))+1),ylim=c(0,1.2*max(p$density)),
		type="n",
		xlab=paste(x$distance," distance (null by ",
		x$method,")",sep=""),ylab="relative frequency")
	y2<-rep(p$density,each=2)
	y2<-y2[-length(y2)]
	x2<-rep(p$breaks[2:length(p$breaks)-1],each=2)[-1]
	x3<-c(min(x2),x2,max(x2))
	y3<-c(0,y2,0)
	polygon(x3,y3,col=make.transparent("blue",0.3),
		border=FALSE)
	lines(p$breaks[2:length(p$breaks)-1],p$density,type="s")
	arrows(x$d,max(c(0.2*max(p$density),
		1.1*p$density[which(p$mids==x$d)])),
		x$d,0,lend="round",
		length=0.15,lwd=2,col="black")
	text(x$d-diff(par()$usr[1:2])/40,
		1.1*max(c(0.2*max(p$density),
		1.1*p$density[which(p$mids==x$d)])),
		"observed distance",
		srt=60,pos=4)
}