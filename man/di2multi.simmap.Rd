\name{di2multi.simmap}
\alias{di2multi.simmap}
\title{Collapse branches of zero length to polytomy in stochastic map style tree}
\usage{
di2multi.simmap(phy, ...)
\method{di2multi}{simmap}(phy, ...)
}
\arguments{
	\item{tree}{modified object of class \code{"phylo"} containing a stochastically mapped discrete character.}
	\item{...}{optional argument \code{tol}, length below which edges should be treated as having zero length.}
}
\description{
	This function collapses branches of zero length (or, more specifically, branches with length shorter than \code{tol}) to create a polytomy in a stochastic-map style tree.
}
\details{
	This function should theoretically perform similarly to \code{\link{di2multi}} in ape.
}
\value{
	A tree with a stochastically mapped discrete character
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{di2multi}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{utilities}
