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
	This function returns node number of the most recent common ancestor of a set of taxa. If \code{tips=NULL} the function is redundant with \code{\link{mrca}} (for \code{type="node"}) or \code{\link{vcv.phylo}}, but much slower (for \code{type="height"}).
}
\details{
	If \code{tips==NULL} and \code{type="node"} (the default) it will return the result of a normal function call to \code{\link{mrca}}. 
	
	If \code{tips=NULL} and \code{type="height"} it will return a matrix equal to that produced by \code{\link{vcv.phylo}}.
	
	From \code{phytools 0.5-66} forward \code{findMRCA} uses \code{\link{getMRCA}} in the \emph{ape} package internally, which results in a big speed-up. Even though the two functions are thus totally redundant I have left \code{findMRCA} in the package to ensure backward compatibility.
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
