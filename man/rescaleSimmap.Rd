\name{rescale.simmap}
\alias{rescale.simmap}
\alias{rescale.multiSimmap}
\alias{rescaleSimmap}
\title{Rescale object of class \code{"simmap"}}
\usage{
\method{rescale}{simmap}(x, model="depth", ...)
\method{rescale}{multiSimmap}(x, model="depth", ...)
rescaleSimmap(tree, ...)
}
\arguments{
	\item{x}{object of class \code{"simmap"} or \code{"multiSimmap"} to be rescaled.}
	\item{model}{model to use to rescale the tree. Currently the only option is \code{"depth"}.}
	\item{...}{parameter of the model to use in rescaling. Currently the only parameter is \code{depth} for \code{model="depth"}.}
	\item{tree}{for \code{rescaleSimmap}, object of class \code{"simmap"} to be rescaled.}
}
\description{
	Scales a tree with a mapped discrete character (\code{"simmap"} object), or a set of such trees, to an arbitrary total height, preserving the relative time spent in each state along each edge.
}
\details{
	Replaces \code{rescaleTree} (now \code{rescale.phylo}) in the \pkg{geiger} package for the \code{"simmap"} object class. \code{rescaleSimmap} is now a redundant alias for the method \code{rescale.simmap}.
}
\value{
	An object of class \code{"simmap"} or \code{"multiSimmap"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\examples{
## load anoletree
data(anoletree)
## rescale to have total depth of 50
rescaled_anoletree<-rescale(anoletree,depth=50)
## plot rescaled tree
plot(rescaled_anoletree,ftype="i",fsize=0.6,
	mar=c(5.1,1.1,1.1,1.1))
axis(1,at=seq(0,50,by=10))
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margin to default
}
\keyword{phylogenetics}
\keyword{utilities}
