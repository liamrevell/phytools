\name{mergeMappedStates}
\alias{mergeMappedStates}
\title{Merge two or more mapped states into one state}
\usage{
mergeMappedStates(tree, old.states, new.state)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"} or \code{"multiSimmap"} containing one or more phylogenetic trees with a mapped discrete character.}
	\item{old.states}{states to merge.}
	\item{new.state}{name for new state.}
}
\description{
	Merges two or mapped states on the tree to get one new state.
}
\details{
	\code{mergeMappedStates} can be used to merge two or more mapped states into a single, new state. For instance, one could merge the states \code{"C"}, \code{"G"}, and \code{"T"} and define the new the state \code{"not-A"}.
}
\value{
	An object of class \code{"simmap"} or \code{"multiSimmap"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{bayesian}
\keyword{utilities}
\keyword{discrete character}
