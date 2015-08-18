\name{phylANOVA}
\alias{phylANOVA}
\title{Phylogenetic ANOVA and post-hoc tests}
\usage{
phylANOVA(tree, x, y, nsim=1000, posthoc=TRUE, p.adj="holm")
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{x}{a vector containing the groups.}
	\item{y}{a vector containing the response variable (continuously valued).}
	\item{nsim}{an integer specifying the number of simulations (including the observed data).}
	\item{posthoc}{a logical value indicating whether or not to conduct posthoc tests to compare the mean among groups.}
	\item{p.adj}{method to adjust P-values for the posthoc tests to account for multiple testing. Options same as \code{\link{p.adjust}}.}
}
\description{
	This function performs the simulation-based phylogenetic ANOVA of Garland et al. (1993) and (optionally) conducts all posthoc comparisons of means among groups (also obtaining the P-values by phylogenetic simulation).
}
\details{
	Uses a little bit of code from \code{phy.anova} in the "geiger" package as well as \code{\link{pairwise.t.test}}.
}
\value{
	A list containing the following elements:
	\item{F}{F from observed data.}
	\item{Pf}{P-value for F from simulation.}
	\item{T}{matrix of t-values.}
	\item{Pt}{matrix of multiple test corrected P-values from posthoc t-tests.}
}
\references{
	Garland, T., Jr., A. W. Dickerman, C. M. Janis, & J. A. Jones. (1993) Phylogenetic analysis of covariance by computer simulation. \emph{Systematic Biology}, \bold{42}, 265-292.

	Harmon, L. J., J. T. Weir, C. D. Brock, R. E. Glor, W. Challenger. (2008) GEIGER: investigating evolutionary radiations. \emph{Bioinformatics}, 24, 129-131.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anova}}, \code{\link{pairwise.t.test}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{statistics}
\keyword{least squares}
\keyword{simulation}
