\name{di2multi.simmap}
\alias{di2multi.simmap}
\alias{multi2di.simmap}
\alias{di2multi.multiSimmap}
\alias{multi2di.multiSimmap}
\alias{di2multi.contMap}
\alias{multi2di.contMap}
\alias{di2multi.densityMap}
\alias{multi2di.densityMap}
\title{Collapse or resolve polytomies in a tree with a character painted on the edges}
\usage{
\method{di2multi}{simmap}(phy, ...)
\method{multi2di}{simmap}(phy, ...)
\method{di2multi}{multiSimmap}(phy, ...)
\method{multi2di}{multiSimmap}(phy, ...)
\method{di2multi}{contMap}(phy, ...)
\method{multi2di}{contMap}(phy, ...)
\method{di2multi}{densityMap}(phy, ...)
\method{multi2di}{densityMap}(phy, ...)
}
\arguments{
	\item{phy}{object of class \code{"simmap"}, \code{"multiSimmap"}, \code{"contMap"}, or \code{"densityMap"} containing a character mapped onto the edges of a tree or set of trees.}
	\item{...}{optional arguments: \code{tol}, length below which edges should be treated as having zero length; and \code{random}, specifying whether to resolve polytomies randomly (if \code{TRUE}) or in the order in which they are encountered.}
}
\description{
	Collapses or resolves polytomies in special types of phylogenetic trees. (See \code{\link[ape]{di2multi}} and \code{\link[ape]{multi2di}}.)
}
\details{
	The method \code{di2multi} collapses branches of zero length (or, more specifically, branches with length shorter than \code{tol}) to create a polytomy in a tree or set of trees.
	
	The method \code{multi2di} resolves polytomies by adding branches of zero length (while preserving the mappings) in a tree or set of trees.

	This methods should theoretically behave similarly to \code{\link[ape]{di2multi}} and \code{\link[ape]{multi2di}} from the \pkg{ape} package.
}
\value{
	An object of class \code{"simmap"}, \code{"multiSimmap"}, \code{"contMap"}, or \code{"densityMap"}, depending on the class of \code{phy}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link[ape]{di2multi}}, \code{\link{make.simmap}}, \code{\link[ape]{multi2di}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{utilities}
