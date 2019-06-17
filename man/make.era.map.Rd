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
	This function creates a temporal map on the tree based on \code{limits} provided by the user.
}
\value{
	An object of class \code{"simmap"} with the specified eras mapped as different regimes.
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
plot(tree,ftype="off",lwd=1)
}
\keyword{phylogenetics}
\keyword{comparative method}


