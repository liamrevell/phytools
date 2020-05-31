\name{phyl.pairedttest}
\alias{phyl.pairedttest}
\title{Phylogenetic paired t-test}
\usage{
phyl.pairedttest(tree, x1, x2=NULL, se1=NULL, se2=NULL, lambda=1.0, h0=0.0,
   fixed=FALSE, ...)
}
\arguments{
	\item{tree}{a phylogeny as an object of class \code{"phylo"}.}
	\item{x1}{data vector for first trait, or matrix with two traits in columns.}
	\item{x2}{data vector for second trait (or null if \code{x1} is a matrix).}
	\item{se1}{standard errors for \code{x1}.}
	\item{se2}{standard errors for \code{x2}.}
	\item{lambda}{starting value for Pagel's \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} (or fixed value, if \code{fixed=TRUE}).}
	\item{h0}{null hypothesis (to be tested) for the mean difference between \code{x1} and \code{x2}.}
	\item{fixed}{logical value specifying whether or not to optimize \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}.}
	\item{...}{optional arguments.}
}
\description{
	This function conducts a phylogenetic paired t-test, roughly following Lindenfors et al. (2010; \emph{J. Evol. Biol.}). This is not a phylogenetic ANOVA, in which we want to compare the means of different sets of species on the tree. Instead, we are interested in the difference between two characters, or two measures of a character within a species, and we want to know if this difference is significantly different from zero controlling for the phylogenetic non-independence of species. 
}
\details{
	Likelihood optimization is performed using \code{\link{optim}} with \code{method="L-BFGS-B"} with box constraints on \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} (0,1).
}
\value{
	An object of class \code{"phyl.pairedttest"} with the following components:
	\item{dbar}{phylogenetic mean difference.}
	\item{se}{standard error of \code{dbar}.}
	\item{sig2}{estimated evolutionary variance (of the difference).}
	\item{lambda}{fitted (or fixed) value of \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}}.}
	\item{logL}{log-likelihood of the fitted model.}
	\item{t.dbar}{t-value (\code{(dbar-h0)/se} where \code{se} is computed from the Hessian).}
	\item{P.dbar}{P-value.}
	\item{df}{the degrees of freedom.}
	\item{h0}{the null hypothesis that was tested.}
}
\references{
	Lindenfors, P., L. J. Revell, and C. L. Nunn (2010) Sexual dimorphism in primate aerobic capacity: A phylogenetic test. \emph{J. Evol. Biol.}, \bold{23}, 1183-1194.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{statistics}
\keyword{least squares}
