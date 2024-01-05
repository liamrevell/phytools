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
	Extracts a \bold{Q}-matrix from a fitted M\emph{k} model.
}
\details{
	This function extracts a \bold{Q}-matrix (in the form of an object of class \code{"Qmatrix"}) from a fitted M\emph{k} model.
}
\value{
	An object of class \code{"Qmatrix"}.
	
	\code{plot.Qmatrix} invisibly returns the coordinates of vertices of the plotted \bold{Q}-matrix.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}
}
\keyword{phylogenetics}
\keyword{maximum likelihood}
\keyword{comparative method}
\keyword{discrete character}
\keyword{utilities}
