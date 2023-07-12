\name{exhaustiveMP}
\alias{exhaustiveMP}
\title{Exhaustive and branch & bound MP optimization}
\usage{
exhaustiveMP(data, tree=NULL, method="branch.and.bound")
}
\arguments{
	\item{data}{is a \code{\link{phyDat}} (Schliep 2011) object containing DNA or other data.}
	\item{tree}{an optional input tree (used only with \code{method="branch.and.bound"}).}
	\item{method}{an optional string indicating which method to use: \code{"branch.and.bound"}, implementing a branch-and-bound search (obviously), or \code{"exhaustive"}.}
}
\description{
	This function does exhaustive and branch & bound MP searches.
}
\details{
	Should probably not be used for more than about 8 species (and definitely not more than 10 species).  Performs parsimony calculations using \code{\link{parsimony}} in the \pkg{phangorn} package (Schliep, 2011).
}
\value{
	A \code{"phylo"} or \code{"multiPhylo"} object that is the MP tree or set of MP trees. It also returns the parsimony scores in \code{attr(trees,"pscore")} or \code{attr(trees[[i]],"pscore")} for the \emph{i}th tree.
}
\references{
	Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sinauer.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Schliep, K. P. (2011) phangorn: phylogenetic analysis in R. \emph{Bioinformatics}, \bold{27}, 592-593.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{mrp.supertree}}, \code{\link{optim.parsimony}}, \code{\link{pratchet}}
}
\keyword{phylogenetics}
\keyword{phylogeny inference}
\keyword{parsimony}
