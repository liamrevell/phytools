## function to draw a tree from an image
## written by Liam J. Revell 2017

tree.drawer<-function(img,file=NULL){
	par(fg=make.transparent("grey",0.8))
	## img<-readJPEG(file) ## file option disabled
	plot.new()
	par(mar=rep(0.1,4))
	plot.window(xlim=c(0,10),ylim=c(0,10))
	rasterImage(img,0,0,10,10)
	cat("  Click the position of the GLOBAL ROOT.\n")
	flush.console()
	root<-unlist(locator(1))
	cat("  Enter the name of a tip RIGHT of the root. > ")
	flush.console()
	right<-readLines(n=1)
	cat(paste("  Click on the position of ",right,".\n",sep=""))
	flush.console()
	right.xy<-unlist(locator(1))
	cat("  Enter the name of a tip LEFT of the root. > ")	
	flush.console()
	left<-readLines(n=1)
	cat(paste("  Click on the position of ",left,".\n",sep=""))
	flush.console()
	left.xy<-unlist(locator(n=1))
	left.xy
	tree<-list(edge=matrix(c(3,3,1,2),2,2),
		edge.length=c(right.xy[1]-root[1],left.xy[1]-root[1]),
		Nnode=1,tip.label=c(right,left))
	class(tree)<-"phylo"
	tips<-setNames(c(right.xy[2],left.xy[2]),tree$tip.label)
	names(tips)<-gsub(" ","_",names(tips))
	plotTree(tree,add=TRUE,tips=tips,xlim=c(0,10)-root[1],ylim=c(0,10),
		color=make.transparent("blue",0.4),lwd=4)
	tip<-0
	cat("  Enter the name of tip to add (or press ENTER). > ")
	flush.console()
	tip<-readLines(n=1)
	while(tip!=""){
		cat(paste("  Click on the position of ",tip,".\n",sep=""))
		flush.console()
		xy<-unlist(locator(1))
		cat("  Click on the position of its MRCA in the built tree.\n")
		flush.console()
		obj<-get.treepos(message=FALSE)
		tree<-bind.tip(tree,tip,edge.length=xy[1]-(nodeheight(tree,obj$where)-
		obj$pos),where=obj$where,position=obj$pos)
		tips<-c(tips,setNames(xy[2],tip))
		names(tips)<-gsub(" ","_",names(tips))
		plot.new()
		par(mar=rep(0.1,4))
		plot.window(xlim=c(0,10),ylim=c(0,10))
		rasterImage(img,0,0,10,10)
		plotTree(tree,add=TRUE,tips=tips,xlim=c(0,10)-root[1],ylim=c(0,10),
			color=make.transparent("blue",0.4),lwd=4)
		old<-tip
		cat("  Enter the name of tip to add (or press ENTER). > ")
		flush.console()
		tip<-readLines(n=1)
		while(tip=="GOBACK"){
			cat(paste("  Dropping ",old,".\n",sep=""))
			tree<-drop.tip(tree,gsub(" ","_",old))
			plot.new()
			par(mar=rep(0.1,4))
			plot.window(xlim=c(0,10),ylim=c(0,10))
			rasterImage(img,0,0,10,10)
			plotTree(tree,add=TRUE,tips=tips,xlim=c(0,10)-root[1],ylim=c(0,10),
			color=make.transparent("blue",0.4),lwd=4)
			cat("  Enter the name of tip to add (or press ENTER). > ")
			flush.console()
			tip<-readLines(n=1)
		}	
	}
	par(fg="black")
	tree
}
