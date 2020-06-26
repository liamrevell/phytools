## function creates a phenogram (i.e., 'traitgram')
## written by Liam J. Revell 2011, 2012, 2013, 2014, 2015, 2016, 2020

phenogram<-function(tree,x,fsize=1.0,ftype="reg",colors=NULL,axes=list(),add=FALSE,...){
	## get optional arguments
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(log)) log<-list(...)$log
	else log<-""
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(sub)) sub<-list(...)$sub
	else sub<-NULL
	if(hasArg(xlab)) xlab<-list(...)$xlab
	else xlab<-"time"
	if(hasArg(ylab)) ylab<-list(...)$ylab
	else ylab<-"phenotype"
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-NA
	if(hasArg(type)) type<-list(...)$type
	else type<-"l"
	if(hasArg(lty)) lty<-list(...)$lty
	else lty<-1
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-0.2
	if(hasArg(offsetFudge)) offsetFudge<-list(...)$offsetFudge
	else offsetFudge<-1.37
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-2
	if(hasArg(nticks)) nticks<-list(...)$nticks
	else nticks<-5
	if(hasArg(spread.labels)) spread.labels<-list(...)$spread.labels
	else spread.labels<-TRUE
	if(ftype=="off") spread.labels<-FALSE
	if(hasArg(spread.cost)) spread.cost<-list(...)$spread.cost
	else spread.cost<-c(1,0.4)
	if(hasArg(spread.range)) spread.range<-list(...)$spread.range
	else spread.range<-range(x)
	if(hasArg(link)) link<-list(...)$link
	else link<-if(spread.labels) 0.1*max(nodeHeights(tree)) else 0
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(label.pos)) label.pos<-list(...)$label.pos
	else label.pos<-NULL
	## end optional arguments
	# check tree
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# check font
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	if(!ftype&&!add) fsize=0 
	H<-nodeHeights(tree)
	if(length(x)<(length(tree$tip)+tree$Nnode))
		x<-c(x,fastAnc(tree,x))
	else
		x<-c(x[tree$tip.label],x[as.character(length(tree$tip)+1:tree$Nnode)])
	x[1:length(tree$tip)]<-x[tree$tip.label]
	names(x)[1:length(tree$tip)]<-1:length(tree$tip)
	X<-matrix(x[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
	# legacy 'axes' argument trumps ylim & xlim from optional (...)
	if(is.null(axes$trait)&&is.null(ylim)) ylim<-c(min(x),max(x))
	else if(!is.null(axes$trait)) ylim<-axes$trait
	if(!is.null(axes$time)) xlim<-axes$time
	if(!add&&is.null(xlim)){
		pp<-par("pin")[1]
		sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
			offsetFudge*offset*fsize*strwidth("W",units="inches")
		alp<-optimize(function(a,H,link,sw,pp) (a*1.04*(max(H)+link)+sw-pp)^2,H=H,
			link=link,sw=sw,pp=pp,interval=c(0,1e6))$minimum
		xlim<-c(min(H),max(H)+link+sw/alp)
	}
	if(!quiet&&Ntip(tree)>=40&&spread.labels){ 
		cat("Optimizing the positions of the tip labels...\n")
		flush.console()
	}
	## matrix for tip coordinates
	tip.coords<-matrix(NA,Ntip(tree),2,dimnames=list(tree$tip.label,c("x","y")))
	if(hold) null<-dev.hold()
	if(is.null(tree$maps)){
		if(is.null(colors)) colors<-"black"
		if(!add){ 
			plot(H[1,],X[1,],type=type,lwd=lwd,lty=lty,col=colors,xlim=xlim,ylim=ylim,
				log=log,asp=asp,xlab="",ylab="",frame=FALSE, axes=FALSE)
			if(spread.labels) tt<-spreadlabels(tree,x,fsize=fsize,cost=spread.cost,
				range=spread.range,label.pos=label.pos) else tt<-x[1:length(tree$tip)]
			if(tree$edge[1,2]<=length(tree$tip)){
				if(fsize&&!add){
					text(gsub("_"," ",tree$tip.label[tree$edge[1,2]]),x=H[1,2]+link,
						y=tt[tree$edge[1,2]],cex=fsize,font=ftype,pos=4,offset=offset)
					tip.coords[tree$tip.label[tree$edge[1,2]],]<-c(H[1,2]+link,
						tt[tree$edge[1,2]])
					if(link>0) lines(x=c(H[1,2],H[1,2]+link),y=c(X[1,2],
						tt[tree$edge[1,2]]),lty=3)
				}
			}
			s<-2
		} else s<-1
		for(i in s:nrow(H)){ 
			lines(H[i,],X[i,],type=type,lwd=lwd,lty=lty,col=colors)
			if(tree$edge[i,2]<=length(tree$tip)){
				if(fsize&&!add){ 
					text(gsub("_"," ",tree$tip.label[tree$edge[i,2]]),x=H[i,2]+link,
						y=tt[tree$edge[i,2]],cex=fsize,font=ftype,pos=4,offset=offset)
					tip.coords[tree$tip.label[tree$edge[i,2]],]<-c(H[i,2]+link,
						tt[tree$edge[i,2]])
					if(link>0) lines(x=c(H[i,2],H[i,2]+link),y=c(X[i,2],tt[tree$edge[i,2]]),
						lty=3)
				}
			}
		}
	} else {
		if(is.null(colors)){
			nn<-sort(unique(c(getStates(tree,"tips"),getStates(tree,"nodes"))))
			colors<-setNames(palette()[1:length(nn)],nn)
		}
		for(i in 1:nrow(H)){
			y<-H[i,1]
			m<-diff(X[i,])/diff(H[i,])
			for(j in 1:length(tree$maps[[i]])){
				a<-c(y,y+tree$maps[[i]][j])
				b<-m*(a-H[i,1])+X[i,1]
				if(i==1&&j==1&&!add) {
					plot(a,b,col=colors[names(tree$maps[[i]])[j]],type=type,lwd=lwd,
						lty=lty,xlim=xlim,ylim=ylim,log=log,asp=asp,axes=FALSE,xlab="",
						ylab="")
					if(spread.labels) tt<-spreadlabels(tree,x[1:length(tree$tip)],
						fsize=fsize,cost=spread.cost,range=spread.range) else 
						tt<-x[1:length(tree$tip)]
				} else lines(a,b,col=colors[names(tree$maps[[i]])[j]],lwd=lwd,lty=lty,
					type=type)
				y<-a[2]
			}
			if(tree$edge[i,2]<=length(tree$tip)){
				if(fsize&&!add){ 
					text(gsub("_"," ",tree$tip.label[tree$edge[i,2]]),x=H[i,2]+link,
						y=tt[tree$edge[i,2]],cex=fsize,font=ftype,pos=4,offset=offset)
					tip.coords[tree$tip.label[tree$edge[i,2]],]<-c(H[i,2]+link,
						tt[tree$edge[i,2]])
					if(link>0) lines(x=c(H[i,2],H[i,2]+link),y=c(X[i,2],
						tt[tree$edge[i,2]]),lty=3)
				}
			}
		}
	}
	if(!add){
		at<-round(0:(nticks-1)*max(H)/(nticks-1),digits)
		axis(1,at=at); axis(2); title(xlab=xlab,ylab=ylab,main=main,sub=sub)
	}
	if(hold) null<-dev.flush()
	xx<-setNames(c(H[1,1],H[,2]),c(tree$edge[1,1],tree$edge[,2]))
	xx<-xx[order(as.numeric(names(xx)))]
	yy<-setNames(c(X[1,1],X[,2]),c(tree$edge[1,1],tree$edge[,2]))
	yy<-yy[order(as.numeric(names(yy)))]
	PP<-list(type="phenogram",use.edge.length=TRUE,node.pos=1,
		show.tip.label=if(ftype!="off") TRUE else FALSE,show.node.label=FALSE,
		font=ftype,cex=fsize,adj=0,srt=NULL,no.margin=FALSE,label.offset=offset,
		x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
		direction=NULL,tip.color="black",Ntip=Ntip(tree),Nnode=tree$Nnode,
		edge=tree$edge,xx=xx,yy=yy)
	assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	invisible(tip.coords)
}

## function to spread labels
## written by Liam J. Revell 2013, 2014, 2016
spreadlabels<-function(tree,x,fsize=1,cost=c(1,1),range=NULL,label.pos=NULL){
	if(!is.null(label.pos)) return(label.pos[tree$tip.label])
	else {
		if(is.null(range)) range<-range(x)
		yy<-x[1:Ntip(tree)]
		zz<-setNames((rank(yy,ties.method="random")-1)/(length(yy)-1)*diff(range(yy))+
			range(yy)[1],names(yy))
		mm<-max(fsize*strheight(tree$tip.label))
		ff<-function(zz,yy,cost,mo=1,ms=1){
			ZZ<-cbind(zz-mm/2,zz+mm/2)
			ZZ<-ZZ[order(zz),]
			oo<-0
			for(i in 2:nrow(ZZ)) 
				oo<-if(ZZ[i-1,2]>ZZ[i,1]) oo<-oo+ZZ[i-1,2]-ZZ[i,1] else oo<-oo
			pp<-sum((zz-yy)^2)
			oo<-if(oo<(1e-6*diff(par()$usr[3:4]))) 0 else oo
			pp<-if(pp<(1e-6*diff(par()$usr[3:4]))) 0 else pp
			oo/mo*cost[1]+pp/ms*cost[2]
		}
		mo<-ff(yy,zz,cost=c(1,0))
		ms<-ff(yy,zz,cost=c(0,1))
		if(mo==0&&ms==0) return(yy)
		else {
			rr<-optim(zz,ff,yy=yy,mo=mo,ms=ms,cost=cost,method="L-BFGS-B",
				lower=rep(range[1],length(yy)),upper=rep(range[2],length(yy)))
			return(rr$par)
		}
	}
}

