\name{rescale}
\alias{rescale}
\title{Rescale phylogenetic objects of different types}
\usage{
rescale(x, ...)
}
\arguments{
	\item{x}{phylogenetic tree object to be rescaled: e.g., object of class \code{"phylo"} or \code{"simmap"}.}
	\item{...}{other arguments to be used in rescaling, depending on the object class. (E.g., see \code{\link[geiger]{rescale.phylo}} in \pkg{geiger} and \code{\link{rescale.simmap}}.)}
}
\description{
	Generic method for rescaling different types of phylogenetic trees.
}
\details{
	See \code{\link[geiger]{rescale.phylo}} in \pkg{geiger} and \code{\link{rescale.simmap}} for details.
}
\value{
	A rescaled phylogenetic tree object.
}
\references{
	Pennell, M.W., J. M. Eastman, G. J. Slater, J. W. Brown, J. C. Uyeda, R. G. FitzJohn, M. E. Alfaro, and L. J. Harmon (2014) geiger v2.0: an expanded suite of methods for fitting macroevolutionary models to phylogenetic trees. \emph{Bioinformatics}, \bold{30}, 2216-2218.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[geiger]{rescale.phylo}}, \code{\link{rescale.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
