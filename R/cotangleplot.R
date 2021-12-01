cotangleplot<-function(tr1,tr2,type=c("cladogram","phylogram"),
	use.edge.length=TRUE,tangle=c("both","tree1","tree2"),...){
	tr1<-untangle(tr1,"read.tree")
	tr2<-untangle(tr2,"read.tree")
	type<-type[1]
	tangle<-tangle[1]
	if(!use.edge.length){
		tr1<-compute.brlen(tr1)
		tr2<-compute.brlen(tr2)
	}
	if(hasArg(layout)) layout<-list(...)$layout
	else layout<-c(0.45,0.1,0.45)
	if(hasArg(nodes)) nodes<-list(...)$nodes
	else nodes<-"centered"
	if(hasArg(lty)) lty<-list(...)$lty
	else lty<-if(type=="phylogram") c("dotted","solid") else
		if(type=="cladogram") "solid"
	if(type=="phylogram") if(length(lty)==1) lty<-rep(lty,2)
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-1
	if(hasArg(color)) color<-list(...)$color
	else color<-palette()[4]
	if(tangle=="both"){
		capture.output(tmp<-cophylo(tr1,tr2))
		tips.tr1<-setNames(1:Ntip(tr1),tmp$trees[[1]]$tip.label)
		tips.tr2<-setNames(1:Ntip(tr2),tmp$trees[[2]]$tip.label)
		tips<-sort(rowMeans(cbind(tips.tr1,tips.tr2[names(tips.tr1)])))
		tips[]<-1:length(tips)
	} else if(tangle=="tree1"){
		tips<-setNames(1:Ntip(tr2),tr2$tip.label)
	} else if(tangle=="tree2"){
		tips<-setNames(1:Ntip(tr1),tr1$tip.label)
	}
	layout(matrix(c(1,2,3),1,3),widths=layout)
	plotTree(tr1,color="transparent",ftype="off",tips=tips,type=type,
		nodes=nodes)
	h<-par()$usr[2]
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	for(i in 1:Ntip(tr1)) lines(c(pp$xx[i],h),rep(pp$yy[i],2),
		lty="dotted")
	if(type=="phylogram"){
		for(i in 1:nrow(tr1$edge))
			lines(pp$xx[tr1$edge[i,]],rep(pp$yy[tr1$edge[i,2]],2),
				lwd=lwd,col=color,lty=lty[2])
		for(i in Ntip(tr1)+1:tr1$Nnode){
			dd<-Children(tr1,i)
			lines(rep(pp$xx[i],length(dd)),pp$yy[dd],lwd=lwd,
				col=color,lty=lty[1])	
		}
	} else if(type=="cladogram"){
		par(ljoin=2)
		for(i in 1:Ntip(tr1)){
			AA<-c(i,Ancestors(tr1,i))
			ii<-sapply(AA[1:(length(AA)-1)],function(x,y) 
				which(y==x), y=tr1$edge[,2])
			lines(c(pp$xx[tr1$edge[ii,2]],pp$xx[Ntip(tr1)+1]),
				c(pp$yy[tr1$edge[ii,2]],pp$yy[Ntip(tr1)+1]),
				lwd=lwd+2,
				col=if(par()$bg=="transparent") "white" else par()$bg,
				lty="solid")
			lines(c(pp$xx[tr1$edge[ii,2]],pp$xx[Ntip(tr1)+1]),
				c(pp$yy[tr1$edge[ii,2]],pp$yy[Ntip(tr1)+1]),
				lwd=lwd,col=color,lty=lty)
		}
	}	
	plot(NA,xlim=c(-1,1),ylim=pp$y.lim,axes=FALSE,xlab="",ylab="")
	text(rep(0,Ntip(tr1)),tips,gsub("_"," ",names(tips)),cex=cex,font=3)
	plotTree(tr2,color="transparent",ftype="off",direction="leftwards",
		tips=tips,type=type)
	h<-par()$usr[1]
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	for(i in 1:Ntip(tr2)) lines(c(pp$xx[i],h),rep(pp$yy[i],2),
		lty="dotted")
	if(type=="phylogram"){
		for(i in 1:nrow(tr2$edge))
			lines(pp$xx[tr2$edge[i,]],rep(pp$yy[tr2$edge[i,2]],2),
				lwd=lwd,col=color,lty=lty[2])
		for(i in Ntip(tr2)+1:tr2$Nnode){
			dd<-Children(tr2,i)
			lines(rep(pp$xx[i],length(dd)),sort(pp$yy[dd]),
				lwd=lwd,col=color,lty=lty[1])
		}
	} else if(type=="cladogram"){
		for(i in 1:Ntip(tr2)){
			AA<-c(i,Ancestors(tr2,i))
			ii<-sapply(AA[1:(length(AA)-1)],function(x,y) 
				which(y==x), y=tr2$edge[,2])
			lines(c(pp$xx[tr2$edge[ii,2]],pp$xx[Ntip(tr2)+1]),
				c(pp$yy[tr2$edge[ii,2]],pp$yy[Ntip(tr2)+1]),
				lwd=lwd+2,
				col=if(par()$bg=="transparent") "white" else par()$bg,
				lty="solid")
			lines(c(pp$xx[tr2$edge[ii,2]],pp$xx[Ntip(tr2)+1]),
				c(pp$yy[tr2$edge[ii,2]],pp$yy[Ntip(tr2)+1]),
				lwd=lwd,col=color,lty=lty)
		}
	}	
}	