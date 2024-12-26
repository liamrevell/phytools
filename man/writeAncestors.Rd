\name{writeAncestors}
\alias{writeAncestors}
\title{Write a tree to file with ancestral states and (optionally) CIs at nodes}
\usage{
writeAncestors(tree, Anc=NULL, file="", digits=6, format=c("phylip","nexus"), 
   ...)
}
\arguments{
	\item{tree}{a phylogenetic tree or set of trees as an object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{Anc}{a vector of ancestral states, a list containing the ancestral states and 95\% confidence intervals (as from \code{\link{fastAnc}} or \code{\link[ape]{ace}}, or a list of such results.}
	\item{file}{an optional string with the filename for output.}
	\item{digits}{an integer indicating the number of digits to print for branch lengths and ancestral character values.}
	\item{format}{a string indicating whether to output the result in simple Newick (i.e., \code{"phylip"}) or Nexus format.}
	\item{...}{additional arguments including \code{x}: a vector of character values, in which case ancestral states are estimated internally using \code{fastAnc}; and \code{CI}: a logical value indicating whether or not to estimate 95\% confidence intervals.}
}
\description{
	This function writes a tree to file with ancestral character states and (optionally) 95\% confidence intervals stored as node value.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\value{
	A file, string, or vector of strings.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{ace}}, \code{\link{fastAnc}}, \code{\link[ape]{write.tree}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{input/output}
