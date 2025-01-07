\name{bind.tip}
\alias{bind.tip}
\title{Attaches a new tip to a tree}
\usage{
bind.tip(tree, tip.label, edge.length=NULL, where=NULL, position=0,
   interactive=FALSE, ...)
}
\arguments{
	\item{tree}{receptor tree.}
	\item{tip.label}{a string containing the species name for the new tip.}
	\item{edge.length}{edge length for the new tip (a scalar).}
	\item{where}{node number to attach new tip. If \code{position>0} then then tip will be attached \emph{rootward} of the specified node. Node numbers can also be tips, in which case the new tip will be added along the terminal edge. To find out the tip number for given species with name \emph{"species"} type: \code{which(tree$tip.label=="species")}.}
	\item{position}{distance \emph{rootward} of the node to add the new tip.}
	\item{interactive}{logical value indicating whether or not the species should be added interactively. (Defaults to \code{FALSE}.)}
	\item{...}{arguments to be passed to \code{plotTree} (for \code{interactive=TRUE}.)}
}
\description{
	Adds a new tip to the tree.
}
\details{
	 If the tree is ultrametric and no branch length is specified, then \code{edge.length} is scaled so that the tree remains ultrametric after the new tip is added.

	Wrapper function for \pkg{ape} \code{\link[ape]{bind.tree}}. Note that \code{interactive=TRUE} works only for right-facing phylograms.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## set seed
set.seed(123)

## subdivide plotting device
par(mfrow=c(3,1))

## generate tree
tree<-pbtree(b=0.1, n=10)

## plot original tree
plotTree(tree,mar=c(2.1,0.1,2.1,0.1))
axisPhylo()
mtext("a)",adj=0)

## add an extant tip ("t_extant") sister to taxon 't5'
## with divergence time of 4.5 Ma
node<-which(tree$tip.label=="t5")
tree<-bind.tip(tree, tip.label="t_extant",
    where=node, position=4.5)

## re-plot to see the result
plotTree(tree,mar=c(2.1,0.1,2.1,0.1))
axisPhylo()
mtext("b)",adj=0)

## add an extinct tip ("t_extinct") sister to 't2' with
## divergence time of 7.8 Ma and duration (edge length) of
## 3.3 Ma
node<-which(tree$tip.label=="t2")
tree<-bind.tip(tree, tip.label="t_extinct", where=node,
    position=7.8, edge.length=3.3)

## re-plot to see the result
plotTree(tree,mar=c(2.1,0.1,2.1,0.1))
axisPhylo()
mtext("c)",adj=0)

## reset par to default
par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
}
\keyword{phylogenetics}
\keyword{utilities}
