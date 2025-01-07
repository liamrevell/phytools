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
	Finds the most recent common ancestor (MRCA) of a set of tips.
}
\details{
	This function returns node number of the most recent common ancestor of a set of taxa.
	
	If \code{tips==NULL} and \code{type="node"} (the default) it will return the result of a normal function call to \code{\link[ape]{mrca}}. 
	
	If \code{tips=NULL} and \code{type="height"} it will return a matrix equal to that produced by \code{\link[ape]{vcv.phylo}}.
	
	From \pkg{phytools} 0.5-66 forward, \code{findMRCA} uses \code{\link[ape]{getMRCA}} in the \pkg{ape} package internally, which results in a big speed-up. Even though the two functions are thus totally redundant I have left \code{findMRCA} in the package namespace to ensure backward compatibility.
}
\value{
  The node number of the MRCA, or a matrix of node numbers (if \code{tips==NULL}) - for \code{type="node"}; or the height of the MRCA, or a matrix of heights (if \code{tips==NULL}) - for \code{type="height"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{findMRCA}}, \code{\link[ape]{mrca}}
}
\examples{
## load tree
data(anoletree)

## find the common ancestor of three taxa
anc<-findMRCA(anoletree,c("cristatellus","cooki",
    "gundlachi"))

## graph the tree
plotTree(anoletree,type="fan",fsize=0.7,lwd=1)

## label the identified ancestor
nodelabels(node=anc,frame="circle",pch=21,cex=1.5,
    bg="blue")
legend("topleft",
    "most recent common ancestor\nof Puerto Rican TG anoles",
    pch=21,pt.cex=1.5,pt.bg="blue",cex=0.7,bty="n")
  
## reset margin to default
par(mar=c(5.1,4.1,4.1,2.1)) 
}
\keyword{phylogenetics}
\keyword{utilities}
