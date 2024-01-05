\name{density.multiSimmap}
\alias{density.multiSimmap}
\alias{plot.changesMap}
\title{Computes a posterior distribution for the number and types of changes on the tree}
\usage{
\method{density}{multiSimmap}(x, ...)
\method{plot}{changesMap}(x, ...)
}
\arguments{
	\item{x}{object of class \code{"multiSimmap"} (see \code{\link{make.simmap}}), or, in the case \code{plot.changesMap}, an object of class \code{"changesMap"} produced via a call to \code{density.multiSimmap}.}
	\item{...}{optional arguments. For \code{density.multiSimmap} these consist of \code{bw} (bandwidth) & \code{method} (\code{"changes"}, \code{"densityMap"}, or \code{"timings"}).}
}
\details{
	\code{density.multiSimmap} \code{method="changes"} (the default) results in a posterior distribution of the number & types of changes on the tree.
	
	If the package \pkg{coda} has been installed, then the function \code{HPD.interval} is used to compute a 95\% high probability density interval for the number of changes of each type on the tree. Otherwise, the central 95\% of the posterior sample is returned as an estimate of the 95\% HPD interval for each change type. 
	
	The method also computes the full posterior density for each change type using a bandwidth specified by the user. 
	
	\code{method="densityMap"} computes a standard \code{"\link{densityMap}"} object, and thus only permits binary characters.
	
	\code{plot.changesMap} plots the posterior density returned by \code{density.multiSimmap} for \code{method=} \code{"changes"}.
}
\description{
	S3 methods to summarize the results of one or more stochastic maps.
}
\value{
	For \code{method="changes"} \code{density.multiSimmap} returns an object of class \code{"changesMap"}.
	
	For \code{method="densityMap"} \code{density.multiSimmap} returns an object of class \code{"\link{densityMap}"}.
	
	\code{plot.changesMap} generates a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{ancestral states}
\keyword{comparative method}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{discrete character}
