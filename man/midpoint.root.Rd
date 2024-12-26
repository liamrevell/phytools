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
	\item{node.labels}{argument to be passed to \code{\link[phangorn]{midpoint}}.}
	\item{...}{optional arguments to be passed to \code{\link[phangorn]{midpoint}}.}
}
\description{
	These functions midpoint root a rooted or unrooted tree (Farris 1972).
}
\details{
	Midpoint rooting involves locating the midpoint of the longest path between any two tips and putting the root in that location. 
	
	The function \code{midpoint_root} performs the same operation as \code{\link[phangorn]{midpoint}} in the \pkg{phangorn} package, but uses no \pkg{phangorn} (Schliep, 2011) code internally.
	
	The function \code{midpoint.root} is a pseudo S3 method for the object class \code{"root"} that exists because when \code{\link[phangorn]{midpoint}} was created in \pkg{phangorn} it was not defined as a generic method. This function merely points to \code{\link[phangorn]{midpoint}} and is being deprecated out.
}
\value{
	An object of class \code{"phylo"} containing a rooted phylogenetic tree.
}
\references{
	Farris, J. (1972) Estimating phylogenetic trees from distance matrices. \emph{American Naturalist}, \bold{106}, 645-667.

	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
	
	Schliep K. P. (2011) phangorn: phylogenetic analysis in R. \emph{Bioinformatics}, \bold{27}, 592-593.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[phangorn]{midpoint}}, \code{\link{reroot}}, \code{\link[ape]{root}}
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{phylogeny inference}
