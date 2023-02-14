\name{optim.phylo.ls}
\alias{optim.phylo.ls}
\title{Phylogeny inference using the least squares method}
\usage{
optim.phylo.ls(D, stree=NULL, set.neg.to.zero=TRUE, fixed=FALSE,
   tol=1e-10, collapse=TRUE)
}
\arguments{
	\item{D}{a distance matrix.}
	\item{stree}{an optional starting tree for the optimization.}
	\item{set.neg.to.zero}{a logical value indicating whether to set negative branch lengths to zero (default \code{TRUE}).}
	\item{fixed}{a logical value indicating whether to estimate the topology - if \code{TRUE} only the branch lengths will be computed.}
	\item{tol}{a tolerance value used to assess whether the optimization has converged.}
	\item{collapse}{ a logical indicating whether to collapse branches with zero length.}
}
\description{
	Phylogenetic inference using the method of least-squares (Cavalli-Sforza & Edwards, 1967).
}
\details{
	Function uses \code{\link{nni}} from the \emph{phangorn} package (Schliep 2011) to conduct NNIs for topology estimation.
	
	Since topology optimization is performed using NNIs, convergence to the true least-squares topology is not guaranteed. It is consequently probably wise to start with a very good tree - such as a NJ tree.
}
\value{
	An object of class \code{"phylo"} that (may be) the least-squares tree with branch lengths; also returns the sum of squares in \code{attr(tree,"Q-score")}.
}
\references{
	Cavalli-Sforza, L. L., and A. W. F. Edwards. (1967) Phylogenetic analysis: Models and estimation procedures. \emph{American Journal of Human Genetics}, \bold{19}, 233-257.

	Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sinauer.

	Paradis, E., J. Claude, and K. Strimmer. (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-290.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Schliep, K. P. (2011) phangorn: phylogenetic analysis in R. \emph{Bioinformatics}, \bold{27}, 592-593.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{exhaustiveMP}}, \code{\link{nni}}
}
\keyword{phylogenetics}
\keyword{phylogeny inference}
\keyword{distance matrix}
\keyword{least squares}
