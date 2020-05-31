\name{as.Qmatrix}
\alias{as.Qmatrix}
\alias{as.Qmatrix.fitMk}
\alias{plot.Qmatrix}
\alias{print.Qmatrix}
\title{Convert a fitted M\emph{k} model to a Q-matrix}
\usage{
as.Qmatrix(x, ...)
\method{as.Qmatrix}{fitMk}(x, ...)
\method{plot}{Qmatrix}(x, ...)
\method{print}{Qmatrix}(x, ...)
}
\arguments{
	\item{x}{fitted M\emph{k} model. (For instance, an object of class \code{"fitMk"}.) In the case of \code{print.Qmatrix}, an object of class \code{"Qmatrix"}.}
	\item{...}{optional arguments.}
}
\description{
	This function extracts a Q-matrix (in the form of an object of class \code{"Qmatrix"} from a fitted M\emph{k} model.
}
\value{
	An object of class \code{"Qmatrix"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}, Joan Maspons}
\seealso{
	\code{\link{fitMk}}
}
\keyword{phylogenetics}
\keyword{maximum likelihood}
\keyword{comparative method}
\keyword{discrete character}
\keyword{utilities}
