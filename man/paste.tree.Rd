\name{paste.tree}
\alias{paste.tree}
\title{Paste two trees together}
\usage{
paste.tree(tr1, tr2)
}
\arguments{
	\item{tr1}{receptor tree.}
	\item{tr2}{donor clade.}
}
\description{
	Primarily internal function for \code{\link{posterior.evolrate}}; can be used to graft a clade into a receptor tree, at the "sticky tip" labeled with \code{"NA"}.
}
\details{
	The donor clade needs to have a root edge, even if it is zero length.
}
\value{
	A tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
tr1<-rtree(10)
tr2<-rtree(10)
tr1$tip.label[1]<-"NA"
tr2$root.edge<-0
tr3<-paste.tree(tr1,tr2)
}
\keyword{phylogenetics}
\keyword{utilities}
