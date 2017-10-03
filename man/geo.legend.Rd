\name{geo.legend}
\alias{geo.legend}
\title{Adds a geological (or other temporal) legend to a plotted tree}
\usage{
geo.legend(leg=NULL, colors=NULL, alpha=0.2, ...)
}
\arguments{
	\item{leg}{a matrix with the starting & ending point of each plotted era in rows, & names of the time periods as rownames.}
	\item{colors}{a vector of colors for the time periods of the rows in \code{leg}.}
	\item{alpha}{transparency level to apply to \code{colors}.}
	\item{...}{optional arguments.}
}
\description{
	This function adds a geological (or other temporal) legend to a plotted tree.
}
\value{
	Functions adds a visual element to a plotted tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
