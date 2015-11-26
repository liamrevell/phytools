\name{mrp.supertree}
\alias{mrp.supertree}
\alias{compute.mr}
\title{Matrix representation parsimony supertree estimation}
\usage{
mrp.supertree(trees, method=c("pratchet","optim.parsimony"), ...)
compute.mr(trees, type=c("phyDat","matrix"))
}
\arguments{
	\item{trees}{an object of class \code{"multiPhylo"} (i.e., a list of trees).}
	\item{method}{an argument specifying whether to optimize the tree using \code{\link{pratchet}} or \code{\link{optim.parsimony}}.}
	\item{type}{for \code{compute.mr}, the type of object to return.}
	\item{...}{optional arguments - mostly to be passed to \code{\link{pratchet}} or \code{\link{optim.parsimony}}.}
}
\description{
	This function estimates the MRP (matrix representation parsimony) supertree from a set of trees (Baum 1992; Ragan 1992).
}
\details{
	Function uses \code{\link{pratchet}} or \code{\link{optim.parsimony}} from the "phangorn" package (Schliep 2011) and \code{\link{prop.part}} from the "ape" package (Paradis et al. 2004). See \code{\link{pratchet}} or \code{\link{optim.parsimony}} for optional arguments, which vary slightly depending on the method. All optional arguments of these methods are available to the user with one exception. The argument \code{tree} in \code{\link{optim.parsimony}} is supplied instead as \code{start}. In addition to being an object of class \code{"phylo"}, \code{start} can also be assigned the string values of \code{"NJ"} or \code{"random"}, in which case either a neighbor-joining or random tree will be used as the starting tree for optimization.

	The function \code{compute.mr} computes the matrix-representation matrix of the input trees. It is used internally by \code{mrp.supertree}, but can also be used to export an object that can be written to file if desired.
}
\value{
	A \code{"phylo"} or \code{"multiPhylo"} object that is the MP or set of MP MRP trees.

	In the case of \code{compute.mr}, an object of class \code{"phyDat"} or a matrix.
}
\references{
	Baum, B. R., (1992) Combining trees as a way of combining data sets for phylogenetic inference, and the desirability of combining gene trees. \emph{Taxon}, \bold{41}, 3-10.

	Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sinauer.

	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Ragan, M. A. (1992) Phylogenetic inference based on matrix representation of trees. \emph{Molecular Phylogenetics and Evolution}, \emph{1}, 53-58.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Schliep, K. P. (2011) phangorn: phylogenetic analysis in R. \emph{Bioinformatics}, \bold{27}, 592-593.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{exhaustiveMP}}, \code{\link{optim.parsimony}}, \code{\link{pratchet}}
}
\keyword{phylogenetics}
\keyword{inference}
\keyword{parsimony}
\keyword{supertree}
