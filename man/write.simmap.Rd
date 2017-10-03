\name{write.simmap}
\alias{write.simmap}
\title{Write a stochastic character mapped tree to file}
\usage{
write.simmap(tree, file=NULL, append=FALSE, map.order=NULL, quiet=FALSE)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"} or \code{"multiSimmap"}.}
	\item{file}{an optional filename.}
	\item{append}{a logical value indicating whether to append to file.}
	\item{map.order}{a optional value specifying whether to write the map in left-to-right or right-to-left order.  Acceptable values are "left-to-right" or "right-to-left" or some abbreviation of either.  If not provided, \code{write.simmap} will use \code{attr(tree,"map.order")} if available.}
	\item{quiet}{logical value indicating whether or not to print a warning message when \code{map.order} is neither specified by a function argument or as an attribute of \code{tree}.}
}
\description{
	This function writes stochastic character mapped trees to file using the Newick style format of SIMMAP v1.0 (Bollback 2006).
}
\value{
	A file or string (if \code{file=NULL}).
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}, \code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{input/output}
