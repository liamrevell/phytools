\name{markChanges}
\alias{markChanges}
\title{Add marked changes to a plotted tree with mapped discrete character}
\usage{
markChanges(tree, colors=NULL, cex=1, lwd=2, plot=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}.}
	\item{colors}{a named vector of colors used to plot the stochastically mapped character on the tree.}
	\item{cex}{expansion factor for line height.}
	\item{lwd}{line width.}
	\item{plot}{logical value indicating whether the changes should be plotted or not.}
}
\description{
	Adds the reconstructed changes to a plotted tree with a stochastically mapped discrete character.
}
\value{
	This function returns (invisibly) a matrix containing the x & y coordinates of the marked changes on the plotted tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}
}
\examples{
## load tree and data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)
## extract discrete character (feeding mode)
fmode<-setNames(sunfish.data$feeding.mode,
  rownames(sunfish.data))
## fit model
er_model<-fitMk(sunfish.tree,fmode,model="ER",
  pi="fitzjohn")
## generate single stochastic map
sunfish_smap<-simmap(er_model,nsim=1)
## plot stochastic map & mark changes
cols<-setNames(c("blue","red"),levels(fmode))
plot(sunfish_smap,cols,ftype="i")
markChanges(sunfish_smap,colors=cols,lwd=6)
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{ancestral states}
\keyword{bayesian}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{discrete character}
