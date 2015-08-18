\name{pgls.Ives}
\alias{pgls.Ives}
\title{Phylogenetic regression with intraspecific sampling error}
\usage{
pgls.Ives(tree, X, y, Vx=NULL, Vy=NULL, Cxy=NULL, lower=c(1e-8,1e-8))
}
\arguments{
	\item{tree}{a phylogeny as an object of class \code{"phylo"}.}
	\item{X}{a named vector containing a \emph{single} independent variable (multiple independent variables to be added in future). \code{X} can contain the species means, or a single long vector containing the sample of values for each species. In the latter case the \code{names(X)} will be repeating - all samples from the same species should have the same name.}
	\item{y}{vector the dependent variable. Can be species means or individual values, as for \code{X}.}
	\item{Vx}{sampling variances for \code{X}. If \code{NULL}, then the within-species variance is computed from the data assuming that individual samples, not species means, have been provided in \code{X}.}
	\item{Vy}{sampling variances for \code{y}. If \code{NULL}, then the within-species variance is computed from the data assuming that individual samples, not species means, have been provided in \code{y}.}
	\item{Cxy}{sampling covariances between \code{X} and \code{y}. This will also be computed from the data if \code{Cxy==NULL}. Note than in this case - but not for the calculation of \code{Vx} and \code{Vy}, the same number of observations and the same ordering must be provided for \code{X} and \code{y}. If this is not the case, then it is assumed that different individuals have been sampled for \code{X} and \code{y} and thus \code{Cxy} is assumed to be zero for all species.}
	\item{lower}{vector specifying the lower bounds for estimation for \code{sig2x} and \code{sig2y}, respectively. (Must be >0.)}
}
\description{
	This function fits the phylogenetic regression model with within-species sampling error following Ives et al. (2007).
}
\details{
	Presently only the bivariate regression model is implemented. Uses \code{optim(...,method="L-BFGS-B")} for optimization. Note that some problems have been reported with the optimization algorithm for this model, which is simple and thus may fail to find the ML solution.
}
\value{
	A list with the following components:
	\item{beta}{a vector or matrix of regression coefficients.}
	\item{sig2x}{fitted BM rate for \code{X}.}
	\item{sig2y}{fitted BM rate for \code{y}.}
	\item{a}{fitted ancestral states for \code{X} and \code{y}.}
	\item{logL}{log-likelihood.}
	\item{convergence}{a value for convergence. \code{convergence=0} is good; see \code{\link{optim}} for more details.}
	\item{message}{a message for convergence.}
}
\references{
	Ives, A. R., P. E. Midford, and T. Garland Jr. (2007) Within-species measurement error in phylogenetic comparative methods. \emph{Systematic Biology}, \bold{56}, 252-270.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{phylosig}}, \code{\link{phyl.resid}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}
\keyword{maximum likelihood}
