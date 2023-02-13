\name{get.treepos}
\alias{get.treepos}
\alias{getnode}
\title{Get position or node of a plotted tree interactively}
\usage{
get.treepos(message=TRUE, ...)
getnode(...)
}
\arguments{
	\item{message}{for \code{get.treepos}, a logical value indicating whether or not to print an instructional message.}
	\item{...}{optional arguments.}
}
\description{
	Return the phylogenetic position of a mouse click on a plotted tree.
}
\details{
	Both \code{get.treepos} and \code{getnode} are primarily meant for internal use in other \emph{phytools} functions.
	
	\code{get.treepos} returns the index of the node at the end of the selected edge, along with the branch distance to that node.
	
	\code{getnode} simply returns the closest node to the user mouse click.
}
\value{
	A list for \code{get.treepos} and a node number for \code{getnode}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{utilities}
\keyword{plotting}
