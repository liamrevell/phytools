\name{labelnodes}
\alias{labelnodes}
\title{Function to interactively label nodes of a plotted tree}
\usage{
labelnodes(text, node=NULL, interactive=TRUE, shape=c("circle","ellipse",
   "rect"), ...)
}
\arguments{
	\item{text}{text string or vector to be used as labels.}
	\item{node}{node numbers (indices) for the labels.}
	\item{interactive}{logical value indicating whether or not nodes should be supplied interactively. (I.e., by clicking on the nodes.)}
	\item{shape}{shape to plot around the plotted node label(s).}
	\item{...}{optional arguments.}
}
\description{
	Adds node labels to a plotted object of class \code{"phylo"}.
}
\details{
	The nodes to be labels can be selected interactively by the user (i.e., by clicking on the corresponding nodes of the plotted tree).
}
\value{
	Invisibly returns a vector of the node indices for the labeled nodes.
}
\seealso{
	\code{\link{cladelabels}}, \code{\link{nodelabels}}
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{plotting}
\keyword{utilities}
