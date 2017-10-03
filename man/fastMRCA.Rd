\name{fastMRCA}
\alias{fastMRCA}
\alias{fastHeight}
\alias{fastDist}
\title{Get the MRCA (or height above the root of the MRCA) of a pair of tip taxa}
\usage{
fastMRCA(tree, sp1, sp2)
fastHeight(tree, sp1, sp2)
fastDist(tree, sp1, sp2)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{sp1}{species name.}
	\item{sp2}{species name.}
}
\description{
	This function returns the most recent common ancestor (node number) for a pair of taxa; or, in the case of \code{fastHeight}, the height above the root of the MRCA of a pair of taxa; or, in the case of \code{fastDist}, the patristic distance between a pair of taxa. 
}
\details{
	This function is mostly redundant with \code{\link{findMRCA}} (or \code{findMRCA(...,type="height")} in the case of \code{fastHeight}) but for very large trees will be considerably faster. (Also see \code{\link{getMRCA}} in the ape package.)
}
\value{
	The node number of the MRCA or the height above the root (for \code{fastHeight}).
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{getMRCA}}, \code{\link{findMRCA}}, \code{\link{mrca}}
}
\examples{
tree<-pbtree(n=2000)
anc<-fastMRCA(tree,"t1","t15")
}
\keyword{phylogenetics}
\keyword{utilities}
