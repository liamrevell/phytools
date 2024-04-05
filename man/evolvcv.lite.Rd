\name{evolvcv.lite}
\alias{evolvcv.lite}
\title{Likelihood test for a shift in the evolutionary correlation between traits}
\usage{
evolvcv.lite(tree, X, maxit=2000, tol=1e-10, ...)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}. If \code{tree} is an object of class \code{"phylo"} then a simple multivariate Brownian motion model will be fit to the data in \code{X}.}
	\item{X}{an \code{n} x \code{m} matrix of tip values for \code{m} continuously valued traits in \code{n} species - row names should be species names. If \code{X} is supplied as a data frame it will be coerced into a matrix without warning.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization - may need to be increased for large trees.}
	\item{tol}{tolerance value for \code{"L-BFGS-B"} optimization.}
	\item{...}{other optional arguments. The most important optional argument is probably \code{models} which species the models to be fit. See \emph{Description} for more information. A second useful argument is \code{error_vcv} which should be supplied as a list of matrices of sampling \emph{variances} and covariances for (and between) the means of each species. The sampling variance for the mean is just the square of the sampling error. Sampling covariances will tend to be zero (or close to zero) if error for different traits is uncorrelated, for instance, because different specimens were used to estimate the means for different traits, and non-zero otherwise.}
}
\description{
	Fits a multi-regime model for differences in the evolutionary variance-covariance structure or correlation between two continuous traits, following Revell & Collar (2009) and Revell et al. (2022).
}
\details{
	This function takes an object of class \code{"simmap"} with a mapped binary or multi-state trait and data for two and only two continuously valued character. It then fits (by default) four different evolutionary models: common rates and correlation; different rates, common correlation; different correlations, common rates; no common structure.
	
	In addition to the four default models specified above, \code{evolvcv.lite} now fits an additional four additional models.
	
	The set of models to be fit can be specified using the optional argument \code{models} in multiple ways.
	
	First, if left unspecified, then the four models listed above will be fit.
	
	Second, if \code{models} is set to \code{"all models"} than eight models will be fit.
	
	Lastly, one or more (up to all eight) models can be fit by encoding the models to be fit into a single vector containing a subset or all of the following elements: \code{"1"}, \code{"2"}, \code{"2b"}, \code{"3"}, \code{"3b"}, \code{"3c"}, and \code{"4"}. These codes correspond to the following eight models: 1. common rates, common correlation; 2. different rates, common correlation; 2b. different rates for trait 1 only, common correlation; 2c. different rates for trait 2 only, common correlation; 3. common rates, different correlations; 3b. different rates for trait 1 only, different correlations; 3c. different rates for trait 2 only, different correlation; and 4. no common structure.
}
\value{
	A list with the results summarized for each model.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.
	
	Revell, L. J., K. S. Toyama, and D. L. Mahler (2022) A simple hierarchical model for heterogeneity in the evolutionary correlation on a phylogenetic tree. \emph{PeerJ}, \bold{10}, e13910.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{evol.vcv}}
}
\examples{
## load data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)
## fit heirarchical common-structure models
sunfish.fit<-evolvcv.lite(sunfish.tree,
    sunfish.data[,2:3],models="all models")
## print fitted models
print(sunfish.fit)
## compare models
anova(sunfish.fit)
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
