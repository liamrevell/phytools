\name{findMRCA}
\alias{findMRCA}
\title{Get the MRCA of a set of taxa}
\usage{
findMRCA(tree, tips=NULL, type=c("node","height"))
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{tips}{a vector containing a set of tip labels.}
	\item{type}{either \code{"node"} to return the node of the MRCA; or \code{"height"} to return the height above the root of the MRCA of \code{tips}.}
}
\description{
	This function returns the most recent common ancestor (node number) for a set of taxa. If \code{tips=NULL} will be redundant with \code{\link{mrca}} (for \code{type="node"}) or \code{\link{vcv.phylo}}, but much slower (for \code{type="height"}).
}
\details{
	If \code{tips==NULL} will return the result of a normal function call to \code{\link{mrca}}. If \code{tips=NULL} will return a matrix equal to \code{\link{vcv.phylo}}.
}
\value{
  The node number of the MRCA, or a matrix of node numbers (if \code{tips==NULL}) - for \code{type="node"}; or the height of the MRCA, or a matrix of heights (if \code{tips==NULL}) - for \code{type="height"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{findMRCA}}, \code{\link{mrca}}
}
\examples{
tree<-pbtree(n=20)
anc<-findMRCA(tree,c("t1","t10","t15"))
}
\keyword{phylogenetics}
\keyword{utilities}
