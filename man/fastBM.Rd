\name{fastBM}
\alias{fastBM}
\title{(Reasonably) fast quantitative trait simulation on phylogenies}
\usage{
fastBM(tree, a=0, mu=0, sig2=1, bounds=c(-Inf,Inf), internal=FALSE, nsim=1, 
   ...)
}
\arguments{
	\item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
	\item{a}{a value for ancestral state at the root node.}
	\item{mu}{an optional value for the mean of random normal changes along branches of the tree - can be used to simulate a trend if \code{mu!=0}.}
	\item{sig2}{instantaneous variance of the BM process, \eqn{\sigma^2}.}
	\item{bounds}{a vector with the lower and upper bounds (respectively) for bounded Brownian simulation - by default simulation is unbounded.}
	\item{internal}{logical value indicating whether or not to return states for internal nodes.}
	\item{nsim}{number of simulations.}
	\item{...}{optional arguments \code{alpha} and \code{theta} used for OU simulation. If \code{alpha} is set then \code{mu} and \code{bounds} are ignored with a warning.}
}
\description{
	Simulates one or multiple continuous traits on the tree under various evolutionary models.
}
\details{
	This function conducts (reasonably) fast quantitative trait simulation on a phylogeny under several different models: Brownian motion (default), BM with a trend (for \code{mu!=0}), bounds (for \code{bounds!=c(-Inf,} \code{Inf)}), and OU. 
}
\value{
	A vector (for \code{nsim=1}) or matrix containing the tip states for the \code{n} species in the tree, and (optionally) the ancestral states for internal nodes.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{sim.corrs}}
}
\examples{
## simulate 10 characters on the Anolis tree 
## under Brownian motion
data(anoletree)
X<-fastBM(anoletree,nsim=10)
head(X)
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{continuous character}
