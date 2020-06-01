\name{allFurcTrees}
\alias{allFurcTrees}
\title{Generate all bi- and multifurcating unrooted trees}
\usage{
allFurcTrees(n, tip.label=NULL, to.plot=TRUE)
}
\arguments{
	\item{n}{an integer giving the desired number of species.}
	\item{tip.label}{an optional vector of length n containing the tip names.}
	\item{to.plot}{an optional logical value indicating whether or not to plot the trees.}
}
\description{
	This function creates all possible unrooted bi- and multifurcating trees and returns a list of trees as an object of class \code{"multiPhylo"}.
}
\details{
	This function should be used with caution for \code{n} greater than about 8, as in this case the number of possible trees is extremely large.
}
\value{
	A list of trees as an object of class \code{"multiPhylo"}.
}
\references{
	Felsenstein, J. 2004. \emph{Inferring Phylogenies}. Sinauer.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{add.everywhere}}, \code{\link{exhaustiveMP}}
}
\examples{
## compute & plot all bi- and multifurcating
## trees for six taxa
trees<-allFurcTrees(n=6)
par(mfrow=c(16,15))
nulo<-sapply(trees,plot,type="unrooted",
    no.margin=TRUE)
par(mfrow=c(1,1))
}
\keyword{phylogenetics}
\keyword{phylogeny inference}
