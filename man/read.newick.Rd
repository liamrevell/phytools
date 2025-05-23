\name{read.newick}
\alias{read.newick}
\alias{readNexus}
\title{Newick or Nexus style tree reader}
\usage{
read.newick(file="", text, ...)
readNexus(file="", format=c("standard","raxml"))
}
\arguments{
	\item{file}{name of text file with single Newick style tree or multiple trees, one per line. For \code{readNexus} this should be a Nexus format tree.}
	\item{text}{character string containing tree.}
	\item{format}{file format (source) for \code{readNexus}. In the case of \code{format="standard"}, \code{\link[ape]{read.nexus}} from \pkg{ape} will be used internally. For \code{format="raxml"}, the parser assumes that bootstrap values have been stored as node labels in the format \code{[&label=bootstrap]}.}
	\item{...}{optional arguments to be passed to \code{\link{scan}}. Note that if the arguments \code{sep} or \code{what} are supplied this could generate an error. Useful optional arguments might include \code{skip} (number of lines to skip) and \code{nlines} (number of lines to read).}
}
\description{
	Reads a phylogenetic tree from file.
}
\details{
	The function \code{read.newick} reads a simple Newick style tree from file. This function is now almost completely redundant with \code{\link[ape]{read.tree}}. At the time of development, however, it was more 'robust' than \code{read.tree} in that it didn't fail if the tree contained so-called 'singles' (nodes with only one descendant); however, \code{read.tree} can now handle singleton nodes without difficulty.
	
	The function \code{readNexus} reads a Nexus formatted tree, optionally with bootstrap values as node labels. This function can read a simple Nexus formatted tree from file (like \code{\link[ape]{read.nexus}}); however, it can also parse the node labels as bootstrap values. This is the output format from the software \emph{RAxML}. For Nexus tree files with complex node labels (e.g., from the software \emph{MrBayes}) it will probably fail to parse node labels correctly, if at all.
}
\value{
	An object of class \code{"phylo"}, possibly containing singletons (see \code{\link[ape]{collapse.singles}}); or an object of class \code{"multiPhylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{read.nexus}}, \code{\link[ape]{read.tree}}
}
\examples{
## create text Newick string
tree<-"(((Human:7,Chimp:7):3,Gorilla:10):15,Monkey:25);"

## parse string
phy<-read.newick(text=tree)

## create plot
plotTree(phy,fsize=1.2,lwd=5,direction="downwards",
    mar=c(2.1,5.1,2.1,1.1),ylim=c(25,-3.8))
axis(2)
title(ylab="time before present (ma)")

}
\keyword{phylogenetics}
\keyword{input/output}
