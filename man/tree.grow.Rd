\name{tree.grow}
\alias{tree.grow}
\title{Creates an animation of a tree growing from left-to-right or upwards}
\usage{
tree.grow(..., res=200, direction="rightwards", ladderize=TRUE)
}
\arguments{
	\item{...}{arguments to pass to \code{\link{pbtree}}.}
	\item{res}{number of steps (the resolution of the animation). This also corresponds to the number of frames that will be created if the animation is to be converted to a .gif file.}
	\item{direction}{the direction to plot the tree. Only \code{direction="rightwards"} (the default) and \code{direction="upwards"} are supported.}
	\item{ladderize}{logical value indicating whether or not to 'ladderize' the plotted tree. (Defaults to \code{TRUE}.)}
}
\description{
	Animates a birth-death tree simulation.
}
\details{
	This function simulates a birth-death tree under user-defined conditions and then creates an animation of that tree growing from left-to-right in the plotting device, or upwards.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{pbtree}}
}
\examples{
\dontrun{
## to create a .gif with ImageMagick installed
png(file="pbtree-%03d.png",width=1000,height=600,res=144)
tree.grow(b=0.06,d=0.02,t=100)
dev.off()
system("ImageMagick convert -delay 5 -loop 0 *.png pbtree-anim.gif")}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{simulation}
\keyword{diversification}
