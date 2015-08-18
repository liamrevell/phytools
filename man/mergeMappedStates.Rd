\name{mergeMappedStates}
\alias{mergeMappedStates}
\title{Merge two or more mapped states into one state}
\usage{
mergeMappedStates(tree, old.states, new.state)
}
\arguments{
	\item{tree}{a modified object of class \code{"phylo"} or \code{"multiPhylo"} with a mapped discrete character.}
	\item{old.states}{states to merge.}
	\item{new.state}{name for new state.}
}
\description{
	This function merges two or mapped states on the tree to get one new state. For instance, one could merge the states \code{"C"}, \code{"G"}, and \code{"T"} to get the state \code{"not-A"}.
}
\value{
	A modified object of class \code{"phylo"} or \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{bayesian}
\keyword{utilities}
