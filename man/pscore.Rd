\name{pscore}
\alias{pscore}
\title{Compute the parsimony score}
\usage{
pscore(tree, x, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{vector (e.g., factor vector), matrix, or data frame. Should contain names or row names.}
	\item{...}{optional arguments.}
}
\description{
	Calculates the parsimony score using the Fitch algorithm.
}
\details{
	Mostly for diagnostic purposes. Users interested in using Maximum Parsimony for phylogeny inference or ancestral state reconstruction should refer to the \pkg{phangorn} package.
}
\value{
	A numerical value or vector of values.
}
\references{
	Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sinauer.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}, \code{\link{sim.Mk}}
}
\examples{
## load tree and data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)

## extract discrete character (feeding mode)
fmode<-setNames(sunfish.data$feeding.mode,
    rownames(sunfish.data))

## compute the parsimony score
pscore(sunfish.tree,fmode)
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
