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
	\item{args.lollipop}{list of arguments for the lollipop plot(s).}
	\item{...}{optional arguments.}
}
\description{
	Plots an upward facing phylogenetic tree with one or more matched lollipop plots. Inspired by a plotting style used by Title et al. (2024).
}
\value{
	Creates a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
	
	Title, P. O., Singhal, S., Grundler, M. C., Costa, G. C., Pyron, R. A., Colston, T. J., Grundler, M. R., Prates, I., Stepanova, N., Jones, M. E. H., Cavalcanti, L. B. Q., Colli, G. R., Di-Poi, N., Donnellan, S. C., Moritz, C., Mesquita, D. O., Pianka, E. R., Smith, S. A., Vitt, L. J., and Rabosky, D. L. (2024) The macroevolutionary singularity of snakes. \emph{Science}, \bold{383}, 918-923.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phylo.heatmap}}, \code{\link{plotTree.barplot}}
}
\examples{
## load tree
data(anoletree)

## convert to simple "phylo" object
anole_tree<-as.phylo(anoletree)

## load trait data
data(anole.data)

## create new data frame with size &
## residuals of other traits
anole_data<-cbind(phyl.resid(anole_tree,
    x=as.matrix(anole.data[,"SVL",drop=FALSE]),
    Y=as.matrix(anole.data[,c(6,4,2)]))$resid,
    exp(anole.data[,"SVL",drop=FALSE]))

## create plot
plotTree.lollipop(anole_tree,anole_data,
    ylab=c("relative (TL)","relative (FLL)",
        "relative (HL)","SVL"))

## reset margins to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{continuous character}
