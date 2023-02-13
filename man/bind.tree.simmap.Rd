\name{bind.tree.simmap}
\alias{bind.tree.simmap}
\title{Binds two trees of class \code{"simmap"}}
\usage{
bind.tree.simmap(x, y, where="root")
}
\arguments{
	\item{x}{an object of class \code{"simmap"}. (The receptor tree.)}
	\item{y}{an object of class \code{"simmap"}. (The tree being grafted.)}
	\item{where}{node number to attach new tip, or the root node if \code{where="root"}.}
}
\description{
	This function grafts tree \code{y} onto tree \code{x} at node \code{where}.
}
\details{
	This function wraps around \code{\link{bind.tree}} for objects of class \code{"simmap"}; however it presently only allows \code{y} to be grafted at a node of \code{x} and it does not allow \code{y} to possess a root edge.
}
\value{
	An object of class \code{"simmap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
