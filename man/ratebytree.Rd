\name{ratebytree}
\alias{ratebytree}
\alias{posthoc.ratebytree}
\title{Likelihood test for rate variation among trees}
\usage{
ratebytree(trees, x, ...)
\method{posthoc}{ratebytree}(x, ...)
}
\arguments{
	\item{trees}{an object of class \code{"multiPhylo"}.}
	\item{x}{a list of trait vectors for a continuous trait in which the names of each vectors correspond to the tip labels of \code{trees}. This is not used if \code{type="diversification"}. In the case of \code{posthoc.ratebytree}, an object of class \code{"ratebytree"}.}
	\item{...}{optional arguments, including the argument \code{type} (\code{"continuous"}, \code{"discrete"}, or \code{"diversification"}), which, if not specified, the function will attempt to ascertain. See \code{Details} for more information.}
}
\description{
	This function essentially implements three different methods for comparing the rate or process of evolution between trees: one for continuously-valued traits, a second for discrete characters, and a third for the rate of diversification (speciation & extinction).
	
	In all cases, the function takes an object of class \code{"multiPhylo"} containing two or more phylogenies (\code{trees}), and, for the first two analyses, a list of trait vectors (\code{x}).

	For continuous traits, the function then proceeds to fit two models: one in which the rate (or regime, for models \code{"OU"} and \code{"EB"}) of trait evolution is equal among all trees; and a second in which the rates or regimes can differ between trees.
	
	The latter model corresponds to an extension the \emph{censored} approach of O'Meara et al. (2006; Revell et al. \emph{In review}) and should also be related to the method of Adams (2012) for comparing rates among traits. See \code{\link{brownie.lite}} for a different implementation of the \emph{noncensored} approach of O'Meara et al. (2006).
	
	For discrete traits, the function instead proceeds to fit two variants of the M<i>k</i> model (Lewis 2001): one in which the parameters values (transition rates) of the process are free to vary between trees, and a second in which they are fixed to be the same. 
	
	For diversification alone, the function fits two different diversification (speciation & extinction) models (Nee et al. 1994; Stadler 2012): one in which the birth (speciation) and death (extinction) rates are identical between the trees, and a second in which they are permitted to differ in various ways depending on the value of \code{"model"}.
	
	The method \code{posthoc} conducts a post-hoc comparison of parameter estimates between trees in the multi-rate or multi-process model. The parameter that is compared depends on the fitted model. For instance, in \code{model="BM"} posthoc comparison is made of \code{sig2}; if \code{model="OU"} fitted values of \code{alpha} are compared; and so on. The argument \code{p.adjust.method} can be used to specify a method for adjusting P-values for multiple tests following \code{p.adjust} (defaults to \code{p.adjust.method="none"}.
}
\details{
	At present it is not possible to specify different models to fit for the different trees - although if (for instance) character evolution on tree 1 proceeded by a strong \emph{OU} process while character evolution on tree 2 was by \emph{BM}, we would probably reject a constant-process model and tree 2 should show a very low value of \code{alpha}.
	
	To compute the standard errors for each fitted paramater value, the function computes the negative inverse of the Hessian matrix at the MLEs; however, if this matrix is computationally singular generalized inverse (\code{\link{ginv}}) will be used instead without warning.

	The function also conducts a likelihood-ratio test to compare the two models.
	
	For continuous character, optional arguments presently include the following: \code{model}, the model of continuous trait evolution (options are \code{"BM"}, the default, \code{"OU"}, and \code{"EB"}). \code{tol}, used as a minimum value for the fitting rates, to prevent problems in optimization. \code{trace}, a logical value indicating whether or not to report progress in the optimization. \code{test}, the method for hypothesis testing (options are \code{"chisq"} and \code{"simulation"}). \code{quiet}, a logical value indicating whether or not to run perfectly quietly. Finally, \code{se}, a list of vectors containing the standard errors for each value of \code{x}. 
	
	For \code{type="discrete"} the optional arguments are slightly different. The argument \code{model} can be used, but it must assume the values \code{"ER"}, \code{"SYM"}, \code{"ARD"}, or a numeric matrix following \code{\link{ace}}. 
	
	Finally, for \code{type= "diversification"} models are so far \code{model="birth-death"}, \code{"equal-extinction"}, and \code{"equal-specation"}, and \code{"Yule"}. It is also important to consider supplying the sampling fractions, \code{rho}, which is a vector of values between 0 and 1 of the same length as \code{trees}. If not provided the method will assume a sampling fraction of 1.0 for all trees - which is seldom true of empirical studies.
}
\value{
	An object of class \code{"ratebytree"} or an object of class \code{"posthoc.ratebytree"} in the case of the method \code{posthoc}.
}
\references{
	Adams, D. C. (2012) Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. \emph{Syst. Biol.}, \bold{62}, 181-192.
	
	Lewis, P. O. (2001) A likelihood approach to estimating phylogeny from discrete morphological character data. \emph{Systematic Biology}, \bold{50}, 913-925.
	
	Nee, S., May, R. M. and Harvey, P. H. (1994) The reconstructed evolutionary process. \emph{Philosophical Transactions of the Royal Society of London B}, \bold{344}, 305-311.

	O'Meara, B. C., C. Ane, M. J. Sanderson, and P. C. Wainwright. (2006) Testing for different rates of continuous trait evolution using likelihood. \emph{Evolution}, \bold{60}, 922-933.

	Stadler, T. (2012) How can we improve the accuracy of macroevolutionary rate estimates? \emph{Systematic Biology}, \bold{62}, 321-329.
	
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
	
	Revell, L. J. (2018) Comparing the rates of speciation and extinction between phylogenetic trees. emph{Ecology and Evolution}, \bold{8}, 5303-5312.
	
	Revell, L. J., Gonzalez-Valenzuela, L. E., Alfonso, A., Castellanos-Garcia, L. A., Guarnizo, C. E., and Crawford, A. J. (2018) Comparing evolutionary rates between trees, clades, & traits. \emph{Methods Ecol. Evol.}, \bold{9}, 994-1005.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{fitMk}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
