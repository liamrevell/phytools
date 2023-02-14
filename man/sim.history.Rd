\name{sim.history}
\alias{sim.history}
\alias{sim.Mk}
\alias{sim.multiMk}
\title{Simulate character history or a discrete character at the tips of the tree under some model}
\usage{
sim.history(tree, Q, anc=NULL, nsim=1, direction=c("column_to_row",
	"row_to_column"), ...)
sim.Mk(tree, Q, anc=NULL, nsim=1, ...)
sim.multiMk(tree, Q, anc=NULL, nsim=1, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}. For the case of \code{sim.multiMk} \code{tree} should be an object of class \code{"simmap"} in which the regimes for simulation have been mapped onto the tree.}
	\item{Q}{a matrix containing the instantaneous transition rates between states. Note that for \code{sim.history} \emph{by default} (i.e., when \code{direction="column_to_row"}, see below) normally this is the \emph{transpose} of the matrix produced by \code{fitDiscrete} in the geiger package or \code{\link{make.simmap}} in phytools; that is to say the transition rate from \code{i -> j} should be given by \code{Q[j,i]}. However, if your matrix is properly conformed (i.e., rows \emph{or} columns sum to 0), then \code{sim.history} will attempt to transpose your matrix correctly & will return an informative message (if \code{message=TRUE}, see below). For \code{sim.Mk} and \code{sim.multiMk} this matrix has the same conformation as in \code{fitContinuous} and \code{make.simmap}. For \code{sim.multiMk} \code{Q} should be a list of transition matrices with names that correspond to the states mapped onto the tree.}
	\item{anc}{an optional value for the state at the root node; if \code{NULL}, a random state will be assigned. \code{anc} can be a vector of states, in which one of the states will be chosen randomly for each simulation. For \code{sim.history} \code{anc} can be a vector of probabilities with names, in which case a state will be chosen in proportion to the given probabilities.}
	\item{nsim}{number of simulations.}
	\item{direction}{row/column direction of the input transition matrix, \code{Q}. \code{direction="column_to_row"} indicates that the transition rate from \code{i -> j} should be given by \code{Q[j,i]}, while \code{direction="row_to_column"} indicates the converse.}
	\item{...}{other optional arguments. Currently only \code{message}, a logical value indicating whether or not to turn on informational messages (defaults to \code{message=TRUE}).}
}
\description{
	Simulates discrete character evolution on a phylogenetic tree.
}
\details{
	The function \code{sim.history} simulates a stochastic character history for a discretely valued character trait on the tree. The resultant tree is stored as a modified \code{"phylo"} object in stochastic character map (e.g., \code{\link{make.simmap}}) format.
	
	The function \code{sim.Mk} simulates the states for a discrete character at the tips of the tree only.
	
	Finally, the function \code{sim.multiMk} is the same as \code{sim.Mk} except that it permits the user to simulate under different values of \code{Q} in different parts of the tree.
}
\value{
	\code{sim.history} returns an object of class \code{"simmap"} (a tree with a mapped discrete character) or \code{"multiSimmap"} for <code>nsim</code> greater than one.
	
	\code{sim.Mk} and \code{sim.multiMk} return a factor with the states of our discrete character at the tips of the tree only.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}, \code{\link{fitmultiMk}}, \code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}, \code{\link{sim.rates}}
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{discrete character}
