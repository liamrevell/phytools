\name{drop.tip.multiPhylo}
\alias{drop.tip.multiPhylo}
\alias{drop.tip.multiSimmap}
\title{Drop tip or tips from an object of class "multiPhylo" or "multiSimmap"}
\usage{
drop.tip.multiPhylo(phy, tip, ...)
drop.tip.multiSimmap(phy, tip)
}
\arguments{
	\item{phy}{an object of class \code{"multiPhylo"} or \code{"multiSimmap"}.}
	\item{tip}{name or names of species to be dropped.}
	\item{...}{optional arguments to be passed to \code{\link{drop.tip}}. Most optional arguments work, with the exception of \code{interactive=TRUE} which will return an error.}
}
\description{
	This function drops one or multiple tips from all the trees of an object of class \code{"multiPhylo"} or \code{"multiSimmap"}.	
}
\details{
	This function merely wraps \code{\link{drop.tip}} and \code{\link{drop.tip.simmap}}. Note that \code{drop.tip.multiSimmap} is merely just an alias of \code{drop.tip.multiPhylo}.
}
\value{
	An object of class \code{"multiPhylo"} or \code{"multiSimmap"}, depending on the input object class.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{drop.tip}}, \code{\link{drop.tip.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
