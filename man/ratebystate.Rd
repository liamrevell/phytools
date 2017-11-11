\name{ratebystate}
\alias{ratebystate}
\title{Method for investigating the rate of one trait as a function of the state of another}
\usage{
ratebystate(tree, x, y, nsim=100, corr=c("pearson","spearman"), ...)
}
\arguments{
	\item{tree}{phylogenetic tree.}
	\item{x}{a continuous character - the dependent variable in the model.}
	\item{y}{a second continuous trait - the response variable.}
	\item{nsim}{number of simulations for hypothesis testing.}
	\item{corr}{correlation method to use. Same as in \code{\link{cor}}.}
	\item{...}{optional arguments which include \code{sim.method} (\code{"fastBM"} or \code{"sim.corrs"}; see \code{\link{fastBM}} and \code{\link{sim.corrs}}); \code{method} (\code{"by.node"} or \code{"by.branch"} indicating whether to assume the rate varies as a function of the node state or the mean branch state); \code{message} - a logical value indicating whether or not to return \code{corr} and \code{method}; finally \code{logarithm} - indicating whether or not to fit a model in which the variance of Brownian evolution in \code{y} changes as a multiplicative function of \code{x}. The default is \code{logarithm=FALSE}.}
}
\description{
	This function attempts to ask if the rate of a continuous character, \code{y}, depends on the state of a separate continuous trait, \code{x}. This is accomplished by regressing the squared contrasts in \code{y} on the branch or node ancestral estimates of \code{x}.
}
\value{
	This function returns an object of class \code{"ratebystate"} with up to the following four elements:
	\item{beta}{value of the regression coefficient for square of the contrasts in \code{y} regressed on the ancestral or branch-wise estimated states for \code{x}.}
	\item{r}{correlation coefficient for \code{corr=corr}.}
	\item{corr}{string giving the value of \code{corr}.}
	\item{method}{string giving the value of \code{method}.}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fastAnc}}, \code{\link{pic}}
}
\keyword{phylogenetics}
\keyword{comparative method}
