\name{midpoint_root}
\alias{midpoint_root}
\alias{midpoint.root}
\title{Midpoint root a phylogeny}
\usage{
midpoint_root(tree)
\method{midpoint}{root}(tree, node.labels="support", ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{node.labels}{argument to be passed to \code{\link{midpoint}}.}
	\item{...}{optional arguments to be passed to \code{\link{midpoint}}.}
}
\description{
	These functions midpoint root a rooted or unrooted tree (Farris 1972).
}
\details{
	Midpoint rooting involves locating the midpoint of the longest path between any two tips and putting the root in that location. 
	
	The function \code{midpoint_root} performs the same operation as \code{\link{midpoint}} in the \pkg{phangorn} package, but uses no \pkg{phangorn} (Schliep, 2011) code internally.
	
	The function \code{midpoint.root} is a pseudo S3 method for the object class \code{"root"} that exists because when \code{\link{midpoint}} was created in \pkg{phangorn} it was not defined as a generic method. This function merely points to \code{\link{midpoint}} and is being deprecated out.
}
\value{
	An object of class \code{"phylo"} containing a rooted phylogenetic tree.
}
\references{
	Farris, J. (1972) Estimating phylogenetic trees from distance matrices. \emph{American Naturalist}, \bold{106}, 645-667.

	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
	
	Schliep K. P. (2011) phangorn: phylogenetic analysis in R. \emph{Bioinformatics}, \bold{27}, 592-593.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{midpoint}}, \code{\link{reroot}}, \code{\link{root}}
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{phylogeny inference}
