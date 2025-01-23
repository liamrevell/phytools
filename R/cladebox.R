## function draws a box over a clade specified by node

cladebox<-function(node,col="#0000FF40",...){
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	phy<-list(edge=pp$edge,Nnode=pp$Nnode,
		tip.label=1:pp$Ntip)
	class(phy)<-"phylo"
	dd<-getDescendants(phy,node)
	dd<-dd[which(dd<=Ntip(phy))]
	if(hasArg(tip_buffer)) tip_buffer<-list(...)$tip_buffer
	else tip_buffer<-NULL
	if(hasArg(edge_buffer)) edge_buffer<-list(...)$edge_buffer
	else edge_buffer<-0.4
	if(pp$type%in%c("fan","arc")){
		if(node==(Ntip(phy)+1)) stem<-Inf
		else {
			parent<-phy$edge[which(phy$edge[,2]==node),1]
			stem<-sqrt(pp$xx[node]^2+pp$yy[node]^2)-
				sqrt(pp$xx[parent]^2+pp$yy[parent]^2)
		}
		if(hasArg(npoints)) npoints<-list(...)$npoints
		else npoints<-30
		if(hasArg(part)) part<-list(...)$part
		else {
			if(all(pp$xx>=0)&&all(pp$yy>=0)) part<-0.25
			else if(all(pp$yy>=0)) part<-0.5
			else part<-1
		}
		theta<-atan(pp$yy[dd]/pp$xx[dd])
		theta[pp$xx[dd]<0]<-pi+theta[pp$xx[dd]<0]
		theta[pp$xx[dd]>0&pp$yy[dd]<0]<-
			2*pi+theta[pp$xx[dd]>0&pp$yy[dd]<0]
		max_theta<-max(theta)+edge_buffer*2*pi/Ntip(phy)*
			max(c(part,0.5))
		min_theta<-min(theta)-edge_buffer*2*pi/Ntip(phy)*
			max(c(part,0.5))
		max_h<-max(sqrt(pp$xx^2+pp$yy^2))
		ff<-if(part<=0.25) 0.6 else 1.2
		if(is.null(tip_buffer)){
			if(pp$show.tip.label){
				tip_buffer<-ff*(pp$x.lim[2]-max_h)
			} else tip_buffer<-0.02*max_h
		}
		h1<-sqrt((pp$xx[node]^2+pp$yy[node]^2))-
			min(c(0.5*stem,0.02*max_h))
		h2<-max(sqrt(pp$xx[dd]^2+pp$yy[dd]^2))+tip_buffer
		tt<-seq(min_theta,max_theta,length.out=npoints)
		out_arc.x<-h2*cos(tt)
		out_arc.y<-h2*sin(tt)
		in_arc.x<-h1*cos(tt)
		in_arc.y<-h1*sin(tt)
		xx<-c(out_arc.x,in_arc.x[npoints:1])
		yy<-c(out_arc.y,in_arc.y[npoints:1])
		polygon(xx,yy,lwd=5,col=col,border=FALSE)
	} else if(pp$type%in%c("phylogram","cladogram")){
		if(pp$direction=="rightwards"){
			if(node==(Ntip(phy)+1)) stem<-0
			else {
				parent<-phy$edge[which(phy$edge[,2]==node),1]
				stem<-pp$xx[node]-pp$xx[parent]
			}
			max_h<-max(pp$xx)
			if(is.null(tip_buffer)){
				tip_buffer<-
					if(pp$show.tip.label) pp$x.lim[2]-max_h else 
						0.02*max_h
			}
			h1<-pp$xx[node]-min(c(0.5*stem,0.02*max_h))
			h2<-max(pp$xx[dd])+tip_buffer
			xx<-c(h1,h2,h2,h1)
			yy<-c(rep(min(pp$yy[dd])-edge_buffer,2),
				rep(max(pp$yy[dd])+edge_buffer,2))
			polygon(xx,yy,border=FALSE,col=col)
		} else if(pp$direction=="leftwards"){
			if(node==(Ntip(phy)+1)) stem<-0
			else {
				parent<-phy$edge[which(phy$edge[,2]==node),1]
				stem<-pp$xx[parent]-pp$xx[node]
			}
			max_h<-max(pp$xx)
			if(is.null(tip_buffer)){
				tip_buffer<-if(pp$show.tip.label) pp$x.lim[1] else 
					-0.02*max_h
			}
			h1<-pp$xx[node]+min(c(0.02*max_h,0.5*stem))
			h2<-max(pp$xx[dd])+tip_buffer
			xx<-c(h1,h2,h2,h1)
			yy<-c(rep(min(pp$yy[dd])-edge_buffer,2),
				rep(max(pp$yy[dd])+edge_buffer,2))
			polygon(xx,yy,border=FALSE,col=col)
		} else if(pp$direction=="upwards"){
			if(node==(Ntip(phy)+1)) stem<-0
			else {
				parent<-phy$edge[which(phy$edge[,2]==node),1]
				stem<-pp$yy[node]-pp$yy[parent]
			}
			max_h<-max(pp$yy)
			if(is.null(tip_buffer)){
				tip_buffer<-
					if(pp$show.tip.label) pp$y.lim[2]-max_h else 
						0.02*max_h
			}
			xx<-c(rep(min(pp$xx[dd])-edge_buffer,2),
				rep(max(pp$xx[dd])+edge_buffer,2))
			h1<-pp$yy[node]-min(c(0.02*max_h,0.5*stem))
			h2<-max(pp$yy[dd])+tip_buffer
			yy<-c(h1,h2,h2,h1)
			polygon(xx,yy,border=FALSE,col=col)
		} else if(pp$direction=="downwards"){
			if(node==(Ntip(phy)+1)) stem<-0
			else {
				parent<-phy$edge[which(phy$edge[,2]==node),1]
				stem<-pp$yy[parent]-pp$yy[node]
			}
			max_h<-max(pp$yy)
			if(is.null(tip_buffer)){
				tip_buffer<-
					if(pp$show.tip.label) pp$y.lim[1] else 
						-0.02*max_h
			}
			xx<-c(rep(min(pp$xx[dd])-edge_buffer,2),
				rep(max(pp$xx[dd])+edge_buffer,2))
			h1<-pp$yy[node]+min(c(0.02*max_h,0.5*stem))
			h2<-max(pp$yy[dd])+tip_buffer
			yy<-c(h1,h2,h2,h1)
			polygon(xx,yy,border=FALSE,col=col)
		}
	}
	invisible(list(x=xx,y=yy))
}
