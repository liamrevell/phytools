\name{drop.tip.multiSimmap}
\alias{drop.tip.multiSimmap}
\alias{keep.tip.multiSimmap}
\title{Drop or keep tip or tips from an object of class \code{"multiSimmap"}}
\usage{
\method{drop.tip}{multiSimmap}(phy, tip, ...)
\method{keep.tip}{multiSimmap}(phy, tip, ...)
}
\arguments{
	\item{phy}{an object of class \code{"multiSimmap"}.}
	\item{tip}{name or names of species to be dropped, in a vector.}
	\item{...}{optional arguments to be passed to \code{\link{drop.tip.simmap}}. Most optional arguments work, with the exception of \code{interactive=TRUE} which will return an error.}
}
\description{
	This function drops or keeps one or multiple tips from all the trees of an object of class \code{"multiSimmap"}.	
}
\details{
	Equivalent to \code{\link{drop.tip}} and \code{\link{keep.tip}} in \pkg{ape}.

	This function merely wraps \code{\link{drop.tip.simmap}}.
}
\value{
	An object of class \code{"multiSimmap"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{drop.tip}}, \code{\link{drop.tip.multiPhylo}}, \code{\link{drop.tip.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
