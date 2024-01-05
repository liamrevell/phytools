\name{anc.ML}
\alias{anc.ML}
\title{Ancestral character estimation using likelihood}
\usage{
anc.ML(tree, x, maxit=2000, model=c("BM","OU","EB"), ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization.}
	\item{model}{model of continuous character evolution on the tree. It's possible that only \code{model="BM"} & \code{model="EB"} work in the present version as \code{model="OU"} has not be thoroughly tested & some bugs were reported for an earlier version.}
	\item{...}{optional arguments.}
}
\description{
	This function estimates the evolutionary parameters and ancestral states for Brownian evolution using likelihood. It is also possible (for \code{model="BM"}) to allow for missing data for some tip taxa.
}
\details{
	Because this function relies on a high dimensional numerical optimization of the likelihood function, \code{\link{fastAnc}} should probably be preferred for most purposes. If using \code{\link{anc.ML}}, users should be cautious to ensure convergence. 
	
	This issue has been ameliorated in \pkg{phytools} >= 0.2-48 by seeding the ML optimization with the result from \code{\link{fastAnc}}.
	
	For \code{model="EB"} this should also not be a problem as the numerical optimization is performed for only \code{sig2} and \code{r}, while the ML values of the ancestral states are obtained during every iteration of the optimization algorithmically using the re-rooting method.
}
\value{
	An object of class \code{"anc.ML"} with at least the following four elements (if not more, depending on \code{model}):
	\item{sig2}{the variance of the BM process.}
	\item{ace}{a vector with the ancestral states.}
	\item{logLik}{the log-likelihood.}
	\item{convergence}{the value of \code{convergence} returned by \code{\link{optim}} (0 is good).}
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Schluter, D., Price, T., Mooers, A. O., and Ludwig, D. (1997) Likelihood of ancestor states in adaptive radiation. \emph{Evolution} \bold{51}, 1699-1711.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{anc.Bayes}}, \code{\link{fastAnc}}, \code{\link{optim}}
}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## extract character of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
    rownames(mammal.data)))
## estimate ancestral state under BM model
fit.BM<-anc.ML(mammal.tree,ln.bodyMass)
print(fit.BM)
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
