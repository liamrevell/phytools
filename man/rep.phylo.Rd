\name{rep.phylo}
\alias{rep.phylo}
\alias{rep.multiPhylo}
\alias{repPhylo}
\title{Replicate a tree or set of trees}
\usage{
\method{rep}{phylo}(x, ...)
\method{rep}{multiPhylo}(x, ...)
repPhylo(tree, times)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{times}{number of times to replicate tree.}
	\item{x}{for S3 method an object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{...}{other arguments for \code{rep} (specifically, \code{times}).}
}
\description{
	\code{\link{rep}} method for object of class \code{"phylo"} or \code{"multiPhylo"}.
}
\details{
	\code{repPhylo} is just an alias for \code{rep.phylo} and \code{rep.multiPhylo}.
}
\value{
	An object of class \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
tree<-pbtree(n=100)
trees<-rep(tree,100)
}
\seealso{
	\code{\link[ape]{c.phylo}}, \code{\link{rep}}
}
\keyword{phylogenetics}
\keyword{utilities}
