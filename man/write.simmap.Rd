\name{write.simmap}
\alias{write.simmap}
\title{Write a stochastic character mapped tree to file}
\usage{
write.simmap(tree, file=NULL, append=FALSE, map.order=NULL)
}
\arguments{
	\item{tree}{a phylogenetic tree as a modified object of class \code{"phylo"}. See \code{\link{make.simmap}} and \code{\link{read.simmap}}.}
	\item{file}{an optional filename.}
	\item{append}{a logical value indicating whether to append to file.}
	\item{map.order}{a optional value specifying whether to write the map in left-to-right or right-to-left order.  Acceptable values are "left-to-right" or "right-to-left" or some abbreviation of either.  If not provided, \code{write.simmap} will use \code{attr(tree,"map.order")} if available.}
}
\description{
	This function writes stochastic character mapped trees to file using the Newick style format of SIMMAP v1.0 (Bollback 2006). Note, can only write one tree at a time to file (hence the \code{append} option).
}
\value{
	A file or string (if \code{file=NULL}).
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\examples{
# simulate a tree & data
tree<-sim.history(pbtree(n=100,scale=1),Q=matrix(c(-1,1,1,-1),2,2))
# generate stochastic character maps
mtrees<-make.simmap(tree,tree$states,nsim=20)
# write them to file
\dontrun{
lapply(mmtrees,write.simmap,file="treefile.tre",append=TRUE) 
}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{input/output}
