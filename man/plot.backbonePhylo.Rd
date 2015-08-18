\name{plot.backbonePhylo}
\alias{plot.backbonePhylo}
\title{Plots backbone tree with triangles as clades}
\usage{
\method{plot}{backbonePhylo}(x, ...)
}
\arguments{
	\item{x}{an object of class \code{"backbonePhylo"}.}
	\item{...}{optional arguments. So far includes only \code{vscale}, which is used to rescale the vertical dimension in plotting.}
}
\description{
	Function plots a backbone tree (stored as an object of class \code{"backbonePhylo"}) with triangles as subtrees.
}
\value{
	Plots a tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.toBackbone}}
}
\examples{
## first create our backbone tree with
## random subtree diversities
tree<-phytools:::lambdaTree(pbtree(n=10),lambda=0.5)
## for old versions of geiger, use lambdaTree
## create a translation table
## leaving a couple of single-taxon clades for fun
tip.label<-sample(tree$tip.label,8)
clade.label<-LETTERS[1:8]
N<-ceiling(runif(n=8,min=1,max=20))
## set crown node depth to 1/2 the maximum depth
depth<-sapply(tip.label,function(x,y) 
	0.5*y$edge.length[which(tree$edge[,2]==
	which(y$tip.label==x))],y=tree)
trans<-data.frame(tip.label,clade.label,N,depth)
rownames(trans)<-NULL
rm(tip.label,clade.label,N,depth)
## here's what trans looks like
trans
## convert
tt<-phylo.toBackbone(tree,trans)
## plot
plot(tt)
}
\keyword{phylogenetics}
\keyword{plotting}

