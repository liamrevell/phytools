\name{vcvPhylo}
\alias{vcvPhylo}
\title{Calculates cophenetic (i.e., phylogenetic VCV) matrix}
\usage{
vcvPhylo(tree, anc.nodes=TRUE, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{anc.nodes}{logical value indicating whether or not to include ancestral nodes.}
	\item{...}{optional arguments including \code{internal} (synonym of \code{anc.nodes}) and \code{model} (can be \code{"BM"}, \code{"OU"}, or \code{"lambda"}.}
}
\description{
	Computes a phylogenetic variance-covariance matrix.
}
\details{
	This function returns a so-called \emph{phylogenetic variance covariance matrix} (e.g., see \code{\link{vcv.phylo}}), but (optionally) including ancestral nodes, as well as under multiple evolutionary models.
	
	\code{vcvPhylo} is designed primarily for internal use by other \emph{phytools} functions.
}
\value{
	A matrix.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{statistics}
\keyword{utilities}
