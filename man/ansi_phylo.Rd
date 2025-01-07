\name{ansi_phylo}
\alias{ansi_phylo}
\title{Graph a phylogeny using only ANSI characters}
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
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotTree}}, \code{\link{splinePhylogram}}
}
\examples{
## load trees
data(salamanders)
data(sunfish.tree)

## plot using defaults
ansi_phylo(salamanders)

## plot using special characters
ansi_phylo(sunfish.tree,vertical="//",horizontal="*",
    xlim=c(0,0.3),fsize=0.9)

## reset par to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
