\name{make.era.map}
\alias{make.era.map}
\title{Create "era" map on a phylogenetic tree}
\usage{
make.era.map(tree, limits, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{limits}{a vector containing the temporal limits, in time since the root node of the tree, for the mappings.}
	\item{...}{optional arguments.}
}
\description{
	This function creates a temporal map on the tree based on \code{limits} provided by the user.
}
\value{
	A modified phylogenetic tree of class \code{"phylo"} with the following additional elements:
	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\examples{
tree<-pbtree(n=1000,scale=100)
tree<-make.era.map(tree,c(0,25,50,75))
plotSimmap(tree,pts=FALSE,ftype="off")
}
\keyword{phylogenetics}
\keyword{comparative method}


