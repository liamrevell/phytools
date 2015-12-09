\name{pbtree}
\alias{pbtree}
\title{Simulate pure-birth or birth-death stochastic tree or trees}
\usage{
pbtree(b=1, d=0, n=NULL, t=NULL, scale=NULL, nsim=1,
	type=c("continuous","discrete"), ...)
}
\arguments{
	\item{b}{birth rate or speciation rate for \code{type="continuous"}; the probability of speciating per time-step for \code{type="discrete"}.}
	\item{d}{death rate or extinction rate for \code{type="continuous"}; the probability of going extinct per time-step for \code{type="discrete"}.}
	\item{n}{desired number of species (i.e., taxa-stop criterion).}
	\item{t}{total time for simulation (i.e., time-stop criterion).}
	\item{scale}{if set, rescales tree to have total length \code{scale}.}
	\item{nsim}{number of simulated trees to return.}
	\item{type}{string to indicate whether to simulate trees in continuous or discrete time. If the former, then wait times between speciation events are drawn from an exponential distribution; whereas if the latter then wait times comes from a geometric distribution.}
	\item{...}{optional arguments including \code{ape}, a logical value indicating whether to return nodes in a 'ape' compatible ordering (default is \code{TRUE}); \code{extant.only} a logical value indicating whether or not to return only extant species (defaults to \code{FALSE}); \code{max.count} a numeric value indicating the maximum number of iterations to run is sampling conditioned on both \code{n} and \code{t} (defaults to \code{1e5}); \code{method} gives the method used for simultaneously conditioning on \code{n} and \code{t} - options are \code{"rejection"} and \code{"direct"}; \code{tip.label}, a vector of tip labels (only works for \code{n!=NULL}); and, finally, \code{quiet}, a logical value indicating whether or not to suppress certain message (defaults to \code{FALSE}).}
}
\description{
	This function simulates stochastic birth-death trees. Simulation can be performed conditioning on \code{n}, on \code{t}, or on both simultaneously. If the both, then (for optional argument \code{method="rejection"}) rejection sampling is performed whereby trees are simulated given \code{b} and \code{t} until a tree containing \code{n} taxa is found. The giving-up point can be set using the optional argument \code{max.count}. Simulations can also be performed in continuous time (the default) or discrete time; the difference being that wait times in the continous-time simulation come from the exponential distribution; whereas waiting times in discrete-time simulations come from the geometric distribution. In addition, discrete-time simulations allow for the possibility that multiple speciation events can occur at (exactly) the same time, so long as they are on separate branches. Finally, sometimes for stopping criterion \code{n} in discrete-time there will be a number of tips different from \code{n}. This indicates that the last event contained more than one speciation event, and a warning is printed.

	\code{method="direct"} is presently experimental. It does not really perform direct sampling; however waiting times & birth or death events are sampled first - with only wait-times consistent with \code{n} and \code{t} being retained. This rejection sampling occurs one layer earlier than for \code{method="rejection"}. This results in a significant (several-fold) speed-up of the code and enables sampling conditioned on \code{n} and \code{t} simultaneously for much higher \code{b} and \code{d}. At the present time, \code{extant.only=TRUE} does not work for this mode, nor does \code{type="discrete"}.

	Note that if \code{ape=FALSE}, then the function will run faster, and the tree is theoretically compatible with the ape \code{"phylo"} standard; however some downstream errors with functions such as \code{\link{bind.tree}} have been observed.
}
\value{
	A tree or list of trees as an object of class \code{"phylo"} or \code{"multiPhylo"}, respectively.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
# simulate a pure-birth tree with 400 tips, scaled to a length of 1.0
tree<-pbtree(n=400,scale=1)
# simulate a pure-birth tree conditioning on n & t
tt<-log(50)-log(2)
tree<-pbtree(n=50,t=tt)
}
\keyword{phylogenetics}
\keyword{simulation}
