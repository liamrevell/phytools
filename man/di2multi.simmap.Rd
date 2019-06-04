\name{di2multi.simmap}
\alias{di2multi.simmap}
\alias{multi2di.simmap}
\alias{di2multi.multiSimmap}
\alias{multi2di.multiSimmap}
\title{Collapse or resolve polytomies in a tree with a character painted on the edges}
\usage{
\method{di2multi}{simmap}(phy, ...)
\method{multi2di}{simmap}(phy, ...)
\method{di2multi}{multiSimmap}(phy, ...)
\method{multi2di}{multiSimmap}(phy, ...)
}
\arguments{
	\item{phy}{object of class \code{"simmap"} or \code{"multiSimmap"} containing a stochastically mapped discrete character.}
	\item{...}{optional arguments: \code{tol}, length below which edges should be treated as having zero length; and \code{random}, specifying whether to resolve polytomies randomly (if \code{TRUE}) or in the order in which they are encountered.}
}
\description{
	The method \code{di2multi} collapses branches of zero length (or, more specifically, branches with length shorter than \code{tol}) to create a polytomy in a stochastic-map style tree or set of trees.
	The method \code{multi2di} resolves polytomies by adding branches of zero length (while preserving the mappings) in a stochastic-map style tree or set of trees.
}
\details{
	This methods should theoretically behave similarly to \code{\link{di2multi}} and \code{\link{multi2di}} from the \emph{ape} package.
}
\value{
	An object of class \code{"simmap"} or \code{"multiSimmap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{di2multi}}, \code{\link{make.simmap}}, \code{\link{multi2di}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{utilities}
