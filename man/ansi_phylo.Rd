\name{ansi_phylo}
\alias{ansi_phylo}
\title{Compute the parsimony score}
\usage{
ansi_phylo(tree, vertical=c("|","-"), ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{vertical}{character for vertical lines}
	\item{...}{optional arguments.}
}
\description{
	Plots a phylogeny in a silly, ANSI graphics style.
}
\value{
	A plotted tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotTree}}, \code{\link{splinePhylogram}}
}
\examples{
data(salamanders)
ansi_phylo(salamanders)
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
