\name{phylo.impute}
\alias{phylo.impute}
\title{Phylogenetic imputation for multivariate continuous character data}
\usage{
phylo.impute(tree, X, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{data matrix with species names as row labels. Missing data to be imputed should be coded \code{NA}.}
	\item{...}{optional arguments.}
}
\description{
	This function performs phylogenetic multiple imputation using maximum likelihood.
}
\details{
	This function performs phylogenetic imputation in which the evolution of the characters in \code{X} is assumed to have occurred by correlation multivariate Brownian motion.
	
	Missing values are imputed by maximizing their likelihood jointly with the parameters of the Brownian model. The function \code{\link{evol.vcv}} is used internally to compute the likelihood.
	
	Note that the \pkg{Rphylopars} package also does phylogenetic imputation for multivariate trait data and it seems to be much faster.
}
\value{
	An object of class \code{"phylo.impute"} consisting of a complete data frame with missing values imputed.
	
	Since optimization is performed numerically using likelihood, a summary of the optimization can be seen by evaluating \code{attr(object,"optim")}, in which \code{object} is of class \code{"phylo.impute"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evol.vcv}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{continuous character}
\keyword{statistics}
