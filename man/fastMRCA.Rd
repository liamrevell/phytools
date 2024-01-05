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
	\item{tree}{an object of class \code{"phylo"}.}
	\item{sp1}{species one name.}
	\item{sp2}{species two name.}
}
\description{
	Finds the most recent common ancestor (MRCA) for a pair of tip taxa.
}
\details{
	Function (\code{fastMRCA}) returns the most recent common ancestor (node number) for a pair of taxa; or, in the case of \code{fastHeight}, the height above the root of the MRCA of a pair of taxa; or, in the case of \code{fastDist}, the patristic distance between a pair of taxa. 
	
	This function is mostly redundant with \code{\link{findMRCA}} (or \code{findMRCA(...,type="height")} in the case of \code{fastHeight}), but for very large trees will be considerably faster.
	
	(Also see \code{\link{getMRCA}} in the \pkg{ape} package.)
}
\value{
	The node number of the MRCA, the height above the root (for \code{fastHeight}), or the patristic distance between two taxa (for \code{fastDist}).
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
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
