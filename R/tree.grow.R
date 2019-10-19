## function to 'grow' a birth-death tree from left-to-right or from upwards
## written by Liam J. Revell 2019

tree.grow<-function(...,res=200,direction="rightwards",ladderize=TRUE){
	tree<-if(ladderize) ladderize(pbtree(...),FALSE) else pbtree(...)
	h<-max(nodeHeights(tree))
	for(i in 1:res){
		dev.hold()
		if(direction=="upwards"){
			if(i<res){ 
				tt<-make.era.map(tree,c(0,i*h/res))
				plot(tt,colors=setNames(c('blue','transparent'),1:2),
					ftype="off",lwd=2,direction="downwards",
					ylim=c(h,0),mar=c(1.1,4.1,1.1,1.1))
			} else { 
				tt<-paintSubTree(tree,Ntip(tree)+1,"1","2")
				tt$mapped.edge<-cbind(tt$mapped.edge,
					rep(0,nrow(tt$mapped.edge)))
				plot(tt,colors=setNames('blue',1),
					ftype="off",lwd=2,direction="downwards",
					ylim=c(h,0),mar=c(1.1,4.1,1.1,1.1))
			}
			axis(2)
			title(ylab="time before present")
		} else {
			if(i<res){ 
				tt<-make.era.map(tree,c(0,i*h/res))
				plot(tt,colors=setNames(c('blue','transparent'),1:2),
					ftype="off",lwd=2,direction="leftwards",
					xlim=c(h,0),mar=c(4.1,1.1,1.1,1.1))
			} else { 
				tt<-paintSubTree(tree,Ntip(tree)+1,"1","2")
				tt$mapped.edge<-cbind(tt$mapped.edge,
					rep(0,nrow(tt$mapped.edge)))
				plot(tt,colors=setNames('blue',1),
					ftype="off",lwd=2,direction="leftwards",
					xlim=c(h,0),mar=c(4.1,1.1,1.1,1.1))
			}
			axis(1)
			title(xlab="time before present")
		}
		dev.flush()
		Sys.sleep(1/res)
	}
	invisible(tree)
}
