## function animates branching random diffusion
## written by Liam Revell 2011, 2013, 2015, 2020

branching.diffusion<-function(sig2=1,b=0.0023,time.stop=1000,ylim=NULL,smooth=TRUE,
	pause=0.02,record=NULL,path=NULL,...){
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"l"
	N<-1
	Y<-matrix(0,1,N)
	if(is.null(ylim)) ylim<-c(-2*sqrt(sig2*time.stop),2*sqrt(sig2*time.stop))
	par(bg="white")
	plot(0,0,xlim=c(0,time.stop),ylim=ylim,xlab="time",ylab="phenotype",
		main="branching diffusion",font.main=3,bty=bty)
	chk<-.check.pkg("animation")
	if(!chk){
		cat("  record != NULL requires the package \"animation\"\n")
		cat("  Animation will play but not record\n\n")
		record<-NULL
	}
	if(!is.null(record)){
		if(is.null(path)) path="C:/Program Files/ffmpeg/bin/ffmpeg.exe"
		tmp<-strsplit(path,"/")[[1]]
		ll<-list.files(paste(tmp[2:length(tmp)-1],collapse="/"))
		if(length(grep(tmp[length(tmp)],ll))>0){
			ani.options(interval=0.01,ffmpeg=path,outdir=getwd())
			ani.record(reset=TRUE)
		} else {
			cat("  record != NULL requires correct path supplied to video renderer\n")
			cat("  Animation will play but not record\n\n")
			record<-NULL
		}
	}
	for(i in 2:time.stop){
		time<-1:i
		Y<-rbind(Y,Y[i-1,])
		for(j in 1:N){
			if(b>runif(n=1)){
				Y<-cbind(Y,Y[,j])
				N<-N+1
				Y[i,N]<-Y[i,j]+rnorm(n=1,sd=sqrt(sig2))
			}
			Y[i,j]<-Y[i,j]+rnorm(n=1,sd=sqrt(sig2))
		}
		dev.hold()
		plot(0,0,xlim=c(0,time.stop),ylim=ylim,xlab="time",ylab="phenotype",
			main="branching diffusion",font.main=3,bty=bty)
		apply(Y,2,lines)
		dev.flush()
		if(!is.null(record)) ani.record()
		Sys.sleep(pause)
	}
	if(!is.null(record)) saveVideo(ani.replay(),video.name=record,other.opts="-b 300k",
		clean=TRUE)
}
