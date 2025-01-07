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
	This function should be used with extreme caution for \code{n} greater than about 8, as in this case the number of possible trees is extremely large.
}
\value{
	A list of trees as an object of class \code{"multiPhylo"}.
}
\references{
	Felsenstein, J. 2004. \emph{Inferring Phylogenies}. Sinauer.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{add.everywhere}}, \code{\link{exhaustiveMP}}
}
\examples{
## compute & plot all bi- and multifurcating
## trees for six taxa
trees<-allFurcTrees(n=6)

## plot trees
par(mfrow=c(16,15))
nulo<-sapply(trees,plot,type="unrooted",
    no.margin=TRUE)

## reset par
par(mfrow=c(1,1))
}
\keyword{phylogenetics}
\keyword{phylogeny inference}
