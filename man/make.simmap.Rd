\name{make.simmap}
\alias{make.simmap}
\alias{simmap}
\title{Simulate stochastic character maps on a phylogenetic tree or trees}
\usage{
make.simmap(tree, x, model="SYM", nsim=1, ...)
simmap(object, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}, or a list of trees as an object of class \code{"multiPhylo"}.}
	\item{x}{a vector containing the tip states for a discretely valued character, or a matrix containing the prior probabilities of tip states in rows and character states as column names. The names (if \code{x} is a vector) or row names (if \code{x} is a matrix) should match the tip labels of the tree. The vector can be of class \code{"factor"}, \code{"character"}, or \code{"numeric"} (although in the lattermost case its content should obviously be only integer values).}
	\item{model}{a character string containing the model or a transition model specified in the form of a matrix. See \code{\link[ape]{ace}} for more details.}
	\item{nsim}{number of simulations. If \code{tree} is an object of class \code{"multiPhylo"}, then \code{nsim} simulations will be conducted \emph{per} input tree.}
	\item{...}{optional arguments. So far, \code{pi} gives the prior distribution on the root node of the tree. Acceptable values for \code{pi} are \code{"equal"}, \code{"estimated"}, or a vector with the frequencies. If \code{pi="estimated"} then the stationary distribution is estimated by numerically solving \code{pi*Q=0} for \code{pi}, and this is used as a prior on the root. If \code{pi="fitzjohn"}, then the Fitzjohn et al. (2009) root prior is used. Finally, if \code{pi} is a numeric vector then the root state will be sampled from this vector. The function defaults to \code{pi="equal"} which results in the root node being sampled from the conditional scaled likelihood distribution at the root. \code{message} tells whether or not to print a message containing the rate matrix, \emph{Q} and state frequencies. \code{message} defaults to \code{TRUE}. For optional argument \code{Q="mcmc"} (see below) the mean value of \code{Q} from the posterior sample is printed. \code{tol} gives the tolerance for zero elements in \code{Q}. (Elements less then \code{tol} will be reset to \code{tol}). Optional argument \code{Q} can be a string (\code{"empirical"} or \code{"mcmc"}), or a fixed value of the transition matrix, \emph{Q}. If \code{"empirical"} than a single value of \emph{Q}, the most likely value, is used for all simulations. If \code{"mcmc"}, then \code{nsim} values of \emph{Q} are first obtained from the posterior distribution for \emph{Q} using Bayesian MCMC, then a simulated stochastic character map is generated for each sampled value of \emph{Q}. Optional argument \code{vQ} can consist of a single numeric value or a vector containing the variances of the (normal) proposal distributions for the MCMC. The order of \code{vQ} is assumed to be in the order of the \code{index.matrix} in \code{\link[ape]{ace}} for the chosen model. \code{prior} is a list containing \code{alpha} and \code{beta} parameters for the \eqn{\Gamma} prior distribution on the transition rates in \emph{Q}. Note that \code{alpha} and \code{beta} can be single values or vectors, if different priors are desired for each value in the transition matrix \emph{Q}. As for \code{vQ}, the order of \code{prior} is assumed to correspond with the order of \code{index.matrix} as in \code{\link[ape]{ace}}. \code{prior} can also be given the optional logical value \code{use.empirical} which tells the function whether or not to give the prior distribution the empirical mean for \emph{Q}. If \code{TRUE} then only \code{prior$beta} is used and \code{prior$alpha} is set equal to \code{prior$beta} times the empirical mean of \emph{Q}. \code{burnin} and \code{samplefreq} are burn-in and sample frequency for the MCMC, respectively.}
	\item{object}{for generic \code{simmap} method, object of various classes: for instance, an object of class \code{"fitMk"} from \code{\link{fitMk}}.}
}
\description{
	Performs stochastic character mapping (Huelsenbeck et al., 2003) using several different alternative methods.
}
\details{
	For \code{Q="empirical"}, \code{make.simmap} first fits a continuous-time reversible Markov model for the evolution of \code{x} and then simulates stochastic character histories using that model and the tip states on the tree. This is the same procedure that is described in Bollback (2006), except that simulation is performed using a fixed value of the transition matrix, \emph{Q}, instead of by sampling \emph{Q} from its posterior distribution.

	For \code{Q="mcmc"}, \code{make.simmap} first samples \emph{Q} \code{nsim} times from the posterior probability distribution of \emph{Q} using MCMC, then it simulates \code{nsim} stochastic maps conditioned on each sampled value of \emph{Q}.

	For \code{Q} set to a matrix, \code{make.simmap} samples stochastic mappings conditioned on the fixed input matrix.

	\code{make.simmap} uses code that has been adapted from \pkg{ape}'s function \code{\link[ape]{ace}} (by Paradis et al.) to perform Felsenstein's pruning algorithm to compute the likelihood.

	As of \pkg{phytools} >= 0.2-33 \code{x} can be a vector of states or a matrix containing the prior probabilities of tip states in rows. In this case the column names of \code{x} should contain the states, and the row names should contain the tip names.

	Note that there was a small (but potentially significant) bug in how node states were simulated by \code{make.simmap} in versions of \pkg{phytools} <= 0.2-26. Between \pkg{phytools} 0.2-26 and 0.2-36 there was also a bug for asymmetric models of character change (e.g., \code{model="ARD"}). Finally, between \pkg{phytools} 0.2-33 and \pkg{phytools} 0.2-47 there was an error in use of the conditional likelihoods for the root node, which caused the root node of the tree to be sampled incorrectly. Giorgio Bianchini pointed out that in \pkg{phytools} 1.0-1 (and probably prior recent versions) there was an error sampling the state at the root node of the tree based on the input prior (\code{pi}) supplied by a user -- except for \code{pi="equal"} (a flat prior, the default) or for a prior distribution in which one or another state was known to be the global root state (e.g., \code{pi=c(1,0)}, \code{pi=c(0,1)}, etc.). All of these issues should be fixed in the current and all later versions.

	If \code{tree} is an object of class \code{"multiPhylo"} then \code{nsim} stochastic maps are generated for each input tree.
}
\value{
 	A object of class \code{"simmap"} or \code{"multiSimmap"} which consists of an object of class \code{"phylo"} (or a list of such objects with class \code{"multiPhylo"}), with the following additional elements:
 	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
 	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
	\item{Q}{the assumed or sampled value of \code{Q}.}
	\item{logL}{the log-likelihood of the assumed or sampled \code{Q}.}
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.
	
	FitzJohn, R. G., W. P. Maddison, and S. P. Otto (2009) Estimating trait-dependent speciation and extinction rates from incompletely resolved phylogenies. \emph{Systematic Biology}, \bold{58}, 595-611.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

 	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
	
	Revell, L. J. and L. J. Harmon (2022) \emph{Phylogenetic Comparative Methods in R}. Princeton University Press.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancr}}, \code{\link{brownie.lite}}, \code{\link{brownieREML}}, \code{\link{countSimmap}}, \code{\link{describe.simmap}}, \code{\link{evol.vcv}}, \code{\link{plotSimmap}}, \code{\link{read.simmap}}, \code{\link{write.simmap}}
}
\examples{
\dontrun{
## load tree and data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)

## extract discrete character (feeding mode)
fmode<-setNames(sunfish.data$feeding.mode,
    rownames(sunfish.data))
  
## fit model
er_model<-fitMk(sunfish.tree,fmode,model="ER",
    pi="fitzjohn")
	
## do stochastic mapping
sunfish_smap<-simmap(er_model)

## print a summary of the stochastic mapping
summary(sunfish_smap)

## plot a posterior probabilities of ancestral states
cols<-setNames(c("blue","red"),levels(fmode))
plot(summary(sunfish_smap),colors=cols,ftype="i")
legend("topleft",c("non-piscivorous","piscivorous"),
    pch=21,pt.bg=cols,pt.cex=2)
par(mar=c(5.1,4.1,4.1,2.1),las=1)

## plot posterior density on the number of changes
plot(density(sunfish_smap),bty="l")
title(main="Posterior distribution of changes of each type",
    font.main=3)
## reset par to default
par(las=0)}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{bayesian}
\keyword{discrete character}
