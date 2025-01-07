\name{make.era.map}
\alias{make.era.map}
\title{Create "era" map on a phylogenetic tree}
\usage{
make.era.map(tree, limits, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{limits}{a vector containing the temporal limits, in time since the root node of the tree, for the mappings. The first number should be \code{0}, and each subsequent number should be the start of each subsequent regime or era to be mapped on the tree.}
	\item{...}{optional arguments.}
}
\description{
	Creates a temporal map on the tree based on \code{limits} provided by the user.
}
\value{
	An object of class \code{"simmap"} with the specified eras mapped as different regimes.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\examples{
## simulate tree with total depth of 100
tree<-pbtree(n=1000,scale=100)

## subdivide into mapped eras
tree<-make.era.map(tree,c(0,25,50,75))

## create graph of mapped tree
plot(tree,ftype="off",lwd=1)

## reset margins to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{comparative method}


