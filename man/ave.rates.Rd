\name{ave.rates}
\alias{ave.rates}
\title{Average the posterior rates}
\usage{
ave.rates(tree, shift, tips, sig1, sig2, ave.shift, showTree=TRUE)
}
\arguments{
	\item{tree}{a tree.}
	\item{shift}{the shift point for this sample.}
	\item{tips}{tip names tipward of \code{shift}.}
	\item{sig1}{\eqn{\sigma_1^2}.}
	\item{sig2}{\eqn{\sigma_2^2}.}
	\item{ave.shift}{average shift from all samples.}
	\item{showTree}{logical value indicating whether to plot the rate-stretched tree.}

}
\description{
	Internal function for \code{\link{posterior.evolrate}}.
}
\value{
	A list of the rates.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evol.rate.mcmc}}, \code{\link{minSplit}}, \code{\link{posterior.evolrate}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
\keyword{continuous characters}
\keyword{math}
