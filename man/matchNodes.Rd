\name{matchNodes}
\alias{matchNodes}
\alias{matchLabels}
\title{Matches nodes between two trees}
\usage{
matchNodes(tr1, tr2, method=c("descendants","distances"), ...)
matchLabels(tr1, tr2)
}
\arguments{
	\item{tr1}{first tree.}
	\item{tr2}{second tree.}
	\item{method}{method to use to match nodes between trees. \code{"descendants"} uses the tip species descended from each node; \code{"distances"} uses the distances from the nodes to the tips. Any umambiguous shortening of \code{"descendants"} or \code{"distances"} is also permitted.}
	\item{...}{optional arguments which may or may not be used depending on the value of \code{method}. \code{tol} is a tolerance value for the difference from exact matching that is allowed for \code{method="distances"}. \code{corr}, which is \code{FALSE} by default, indicates whether to match nodes under \code{method="distances"} using the correlation (\code{corr=TRUE}) or the absolute similarity of distances.}
}
\description{
	This function returns a matrix in which the first column contains \emph{all} of the internal nodes of \code{tr1} and the second column contains the matching nodes from \code{tr2}, inasmuch as they can be identified. 
	
	For \code{method="descendants"}, pairs of matching nodes are defined by sharing all descendant leaves in common. 
	
	For \code{method="distances"}, nodes are considered to matched if the share the same set of distances (or proportional distances, for optional argument \code{corr=TRUE}) to all tips. 
	
	\code{matchLabels} is functionally equivalent but matches node (tip) indices based on identifying matching in the labels only.
}
\value{
	A matrix in which the first column contains the nodes of \code{tr1} with the second column containing matching nodes in \code{tr2}, with the criterion for matching defined by \code{method}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
