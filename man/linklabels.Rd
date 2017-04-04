\name{linklabels}
\alias{linklabels}
\title{Function to add tip labels to a plotted tree with linking lines}
\usage{
linklabels(text,tips,link.type=c("bent","curved","straight"),
   ...)
}
\arguments{
	\item{text}{text string or vector to be used as labels.}
	\item{tips}{node numbers (indices) for the tips to be labeled.}
	\item{link.type}{manner in which to draw the linking lines.}
	\item{...}{optional arguments, including \code{cex}, \code{lty}, \code{lwd}, and \code{col}.}
}
\description{
	Function adds tip labels to a plotted tree by drawing curved, bent, or straight linking lines.
}
\details{
	The idea underlying this function is that the user should first plot the tree without tip labels, but set the area of the plotting device to be sufficient to accommodate the tip labels once they have been added. The function then can be called to add tip labels connected by linking lines to the tips of the plotted tree.
}
\value{
	This function annotates a plot.
}
\seealso{
	\code{\link{cladelabels}}, \code{\link{nodelabels}}, \code{\link{tiplabels}}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{plotting}
\keyword{utilities}
