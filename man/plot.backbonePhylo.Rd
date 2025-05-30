\name{plot.backbonePhylo}
\alias{plot.backbonePhylo}
\title{Plots backbone tree with triangles as clades}
\usage{
\method{plot}{backbonePhylo}(x, ...)
}
\arguments{
	\item{x}{an object of class \code{"backbonePhylo"}.}
	\item{...}{optional arguments. Includes \code{vscale} (to rescale the vertical dimension in plotting), \code{fixed.height} (logical value to fix the plotted height of subtree triangles), \code{print.clade.size} (logical), \code{fixed.n1} (logical value indicating whether or not to use the fixed triangle height for subtrees containing only one taxon, or to plot as a leaf - defaults to \code{FALSE}), and \code{col} (a single value, or a vector with names, giving the clade colors), as well as \code{xlim}, \code{ylim}, and \code{lwd} (as well as perhaps other standard plotting arguments).}
}
\description{
	Plots a backbone tree (stored as an object of class \code{"backbonePhylo"}) with triangles as subtrees.
}
\value{
	Plots a tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.toBackbone}}
}
\examples{
## first create our backbone tree with
## random subtree diversities
tree<-phytools:::lambdaTree(pbtree(n=10),lambda=0.5)

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

## here's what trans looks like
print(trans)

## convert
obj<-phylo.toBackbone(tree,trans)

## plot
plot(obj)

## reset par to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}

