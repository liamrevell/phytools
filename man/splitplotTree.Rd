\name{splitplotTree}
\alias{splitplotTree}
\alias{plotTree.splits}
\title{Plots a phylogeny in two columns}
\usage{
splitplotTree(tree, fsize=1.0, ftype="reg", lwd=2, split=NULL, new.window=FALSE)
plotTree.splits(tree, splits=NULL, file=NULL, fn=NULL, ...)
}
\arguments{
 	\item{tree}{an object of class "phylo".}
 	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{lwd}{line width for plotting.}
	\item{split}{relative vertical position for splitting the tree (between 0 & 1).}
	\item{new.window}{whether or not to plot the split tree in a new window. If \code{FALSE} then the tree will be plotted in two columns within the same plotting window.}
	\item{splits}{for \code{plotTree.splits} relative positions (from 0 to 1) to split the tree across pages or devices.}
	\item{file}{filename if saving to a PDF file is desired. Otherwise will plot to the default plotting device.}
	\item{fn}{function to be executed on each plotted page. For instance, might be: \code{function()} \code{cladelabels()} if clade labels are desired.}
	\item{...}{other arguments to be passed to \code{\link{plotTree}}.}
}
\description{
	Function plots a tree in two columns or windows.
}
\value{
	Plots a tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotTree}}, \code{\link{plotSimmap}} 
}
\keyword{phylogenetics}
\keyword{plotting}
