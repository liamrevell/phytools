\name{force.ultrametric}
\alias{force.ultrametric}
\title{Coerces a phylogenetic tree to be ultrametric}
\usage{
force.ultrametric(tree, method=c("nnls","extend"), ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{method}{the method to use to force the tree to be ultrametric. Options are \code{"nnls"} (which uses the \pkg{phangorn} function \code{\link[phangorn]{nnls.tree}} internally), or \code{"extend"}.}
	\item{...}{optional arguments: principally, \code{message}. This argument (if set to \code{FALSE}) can be used to suppress the default warning message that \code{force.ultrametric} \emph{should not} be used as a formal statistical method to ultrametricize a tree.}
}
\description{
	Coerces an object of class \code{"phylo"} to be ultrametric. 
}
\details{
	\code{force.ultrametric} coerces a non-ultrametric tree to be ultrametric.
	
	This is achieved either by using \code{\link[phangorn]{nnls.tree}} from the \pkg{phangorn} package to compute the set of edge lengths that result in a minimized sum-of-squares distance between the patristic distance of the output and input trees (\code{method="nnls"}); or by simply extending all the external edges of the tree to match the external edge with the greatest total height (\code{method="extend"}). 
	
	Note that neither of these should be treated as formal statistical methods for inferring an ultrametric tree. Rather, this method can be deployed when a genuinely ultrametric tree read from file fails \code{\link[ape]{is.ultrametric}} for reasons of numerical precision. 
}
\value{
	An ultrametric tree in an object of class \code{"phylo"}.
}
\seealso{
	\code{\link[ape]{is.ultrametric}}, \code{\link[phangorn]{nnls.tree}}
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{statistics}
\keyword{utilities}
