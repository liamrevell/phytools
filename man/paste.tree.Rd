\name{paste.tree}
\alias{paste.tree}
\title{Paste two trees together}
\usage{
paste.tree(tr1, tr2)
}
\arguments{
	\item{tr1}{receptor tree.}
	\item{tr2}{donor clade.}
}
\description{
	Internal function for \code{\link{posterior.evolrate}}.
}
\details{
	Primarily designed as an internal function for \code{\link{posterior.evolrate}}; however, can be used to graft a clade onto a receptor tree at the "sticky tip" labeled with \code{"NA"}.

	The donor clade needs to have a root edge, even if it is zero length.
}
\value{
	A tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## create two random trees
tr1<-rtree(10)
tr2<-rtree(10)

## make one tip "sticky" by renaming it "NA"
tr1$tip.label[sample(1:10,1)]<-"NA"
tr2$root.edge<-0

## paste the two trees together
tr3<-paste.tree(tr1,tr2)

## plot the results
layout(matrix(c(1,1,2,2,4,3,3,4),2,4,byrow=TRUE))
plotTree(tr1,direction="upwards",
    mar=c(1.1,1.1,3.1,1.1))
add.arrow(tree=tr1,tip=which(tr1$tip.label=="NA"),
    hedl=0.2,arrl=0.2,offset=1)
mtext("a) recipient tree",adj=0)
plotTree(tr2,direction="upwards",plot=FALSE)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
plotTree(tr2,direction="upwards",
    mar=c(1.1,1.1,3.1,1.1),
    ylim=c(-0.2,pp$y.lim[2]),add=TRUE)
arrows(pp$xx[11],pp$yy[11],pp$xx[11],-0.2,
    length=0.1,lwd=2)
mtext("b) donor tree",adj=0)
plotTree(tr3,direction="upwards",
    mar=c(1.1,1.1,3.1,1.1))
mtext("c) grafted tree",adj=0)

## reset par to default
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{utilities}
