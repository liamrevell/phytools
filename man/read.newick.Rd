\name{read.newick}
\alias{read.newick}
\title{Robust Newick style tree reader}
\usage{
read.newick(file="", text, ...)
}
\arguments{
	\item{file}{name of text file with single Newick style tree or multiple trees, one per line.}
	\item{text}{character string containing tree.}
	\item{...}{optional arguments to be passed to \code{\link{scan}}. Note that if the arguments \code{sep} or \code{what} are supplied this could generate an error. Useful optional arguments might include \code{skip} (number of lines to skip) and \code{nlines} (number of lines to read).}
}
\description{
	This function reads a Newick style tree from file.
}
\details{
	This function is almost completely redundant with \code{\link{read.tree}}; however it is 'robust' in that it does not fail if the tree contains so-called 'singles' (nodes with only one descendant).
}
\value{
	An object of class \code{"phylo"}, possibly containing singles (see \code{\link{collapse.singles}}).
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{read.tree}}, \code{\link{read.nexus}}
}
\examples{
tree<-"((Human,Chimp),Gorilla),Monkey);"
phy<-read.newick(text=tree)
}
\keyword{phylogenetics}
\keyword{input/output}
