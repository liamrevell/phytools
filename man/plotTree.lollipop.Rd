\name{plotTree.lollipop}
\alias{plotTree.lollipop}
\title{Plot a tree with one or more matched lollipop plots}
\usage{
plotTree.lollipop(tree, x, args.plotTree=list(), args.lollipop=list(), ...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a matrix, data frame, or vector of numeric values.}
	\item{args.plotTree}{list of arguments to be passed to \code{\link{plotTree}}.}
	\item{args.lollilop}{list of arguments for the lollipop plot(s).}
	\item{...}{optional arguments.}
}
\description{
	Plots an upward facing phylogenetic tree with one or more matched lollipop plots.
}
\value{
	Creates a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.heatmap}}, \code{\link{plotTree.barplot}}
}
\examples{
data(anoletree)
anole_tree<-as.phylo(anoletree)
data(anole.data)
anole_data<-cbind(phyl.resid(anole_tree,
  x=as.matrix(anole.data[,"SVL",drop=FALSE]),
  Y=as.matrix(anole.data[,c(6,4,2)]))$resid,
  exp(anole.data[,"SVL",drop=FALSE]))
plotTree.lollipop(anole_tree,anole_data,
  ylab=c("relative (TL)","relative (FLL)",
    "relative (HL)","SVL"))
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{discrete character}
