\name{write.simmap}
\alias{write.simmap}
\title{Write a stochastic character mapped tree to file}
\usage{
write.simmap(tree, file=NULL, append=FALSE, map.order=NULL, quiet=FALSE, 
    format="phylip", version=1.0)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"} or \code{"multiSimmap"}.}
	\item{file}{an optional filename.}
	\item{append}{a logical value indicating whether to append to file.}
	\item{map.order}{a optional value specifying whether to write the map in left-to-right or right-to-left order.  Acceptable values are "left-to-right" or "right-to-left" or some abbreviation of either.  If not provided, \code{write.simmap} will use \code{attr(tree,"map.order")} if available.}
	\item{quiet}{logical value indicating whether or not to print a warning message when \code{map.order} is neither specified by a function argument or as an attribute of \code{tree}.}
	\item{format}{file format for output.}
	\item{version}{version of SIMMAP. Note that the options are \code{1.0} and \code{1.5}. \code{version=1.5} is generally recommended because in this format the tree can also be parsed by typical tree readers, but absent its mapped trait.}
}
\description{
	Writes one or multiple stochastic character mapped trees to file in several formats (following Bollback, 2006).
}
\value{
	A file.
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{input/output}
