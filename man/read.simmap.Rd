\name{read.simmap}
\alias{read.simmap}
\title{Read SIMMAP style trees from file}
\usage{
read.simmap(file="", text, format="nexus", rev.order=TRUE, version=1)
}
\arguments{
	\item{file}{name of text file with one or multiple SIMMAP v1.0 or v1.5 style trees.}
	\item{text}{character string containing the tree. If \code{version=1.5} this argument is ignored. (This format tree can only be read from file in the present version.)}
	\item{format}{format of the trees: either \code{"phylip"} or \code{"nexus"} - the latter is the default output from SIMMAP. If \code{version=1.5} this argument is ignored.}
	\item{rev.order}{a logical value indicating whether the states and times along each branch is given (from root to tip) in right-to-left order (if TRUE) or in left-to-right order. If \code{version=1.5} this argument is ignored.}
	\item{version}{version of SIMMAP for input tree. If the tree(s) was/were simulated in SIMMAP v1.0 or written to file by \code{link{make.simmap}} then \code{version=1.0}; if the tree(s) was/were simulated using SIMMAP v1.5 then \code{version=1.5}.}
}
\description{
	This reads one or multiple SIMMAP style trees from file.
}
\details{
	This function now accepts trees in both SIMMAP v1.0 and SIMMAP v1.5 format. In addition, it can read a more flexible format than is produced by SIMMAP (for instance, multi-character mapped states and more than 7 mapped states).  Uses some modified code from \code{\link{read.nexus}} from the "ape" package to read the NEXUS block created by SIMMAP.  Also creates the attribute \code{"map.order"} which indicates whether the stochastic map was read in from left to right or right to left.  This attribute is used by default by \code{\link{write.simmap}} to write the tree in the same order. 
}
\value{
	A modified object of class \code{"phylo"} (or list of class \code{"multiPhylo"}) with the following additional elements:
	\item{maps}{a list of named vectors containing the times spent in each state on each branch, in the order in which they occur.}
	\item{mapped.edge}{a matrix containing the total time spent in each state along each edge of the tree.}
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{evol.vcv}}, \code{\link{read.tree}}, \code{\link{read.nexus}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{input/output}
\keyword{discrete character}
