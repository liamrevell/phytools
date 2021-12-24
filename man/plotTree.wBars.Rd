\name{plotTree.wBars}
\alias{plotTree.wBars}
\alias{plotTree.barplot}
\alias{plotTree.boxplot}
\title{Plot a tree with bars at the tips}
\usage{
plotTree.wBars(tree, x, scale=NULL, width=NULL, type="phylogram", 
   method="plotTree", tip.labels=FALSE, col="grey", border=NULL, 
   ...)
plotTree.barplot(tree, x, args.plotTree=list(), args.barplot=list(), 
	...)
plotTree.boxplot(tree, x, args.plotTree=list(), args.boxplot=list(),
	...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a named vector or matrix of trait values. For \code{plotTree.boxplot}, the names should repeat for multiple observations per species. For \code{plotTree.boxplot} \code{x} can also be supplied as a formula, though in that case the factor levels need to be provided in a valid cladewise order of the tips in \code{tree}. This order doesn't need to correspond with the current order of the tip labels. For \code{plotTree.barplot} \code{x} can be a matrix (or a data frame) in which columns are the values of multiple traits to be simultaneously plotted on the tree.}
 	\item{scale}{scaling factor for the tip bars (relative to the total tree height). If left as \code{NULL} a reasonable scaling factor is computed automatically.}
	\item{width}{width of the tip bars.}
	\item{type}{plot type. Can be \code{"phylogram"} or \code{"fan"}.}
	\item{method}{plotting method to use. Can be \code{"plotTree"} (for \code{\link{plotTree}}) or \code{"plotSimmap"} (for \code{\link{plotSimmap}}).}
	\item{tip.labels}{argument indicating whether or not tip labels should be plotted. Defaults to \code{tip.labels=FALSE}.}
	\item{col}{colors of the plotted bars. Can be a single value or a vector with length equal to the number of tips in the tree.}
	\item{border}{single value specifying the color of the border for the plotted bars. Defaults to \code{border=NULL}, which means that black borders will be plotted.}
	\item{args.plotTree}{in \code{plotTree.barplot}, arguments to be passed to \code{\link{plotTree}}.}
	\item{args.barplot}{in \code{plotTree.barplot}, arguments to be passed to \code{\link{barplot}}.}
	\item{args.boxplot}{in \code{plotTree.boxplot}, arguments to be passed to \code{\link{boxplot}}.}
	\item{...}{optional arguments to be passed to \code{\link{plotTree}} or \code{\link{plotSimmap}} in the case of \code{plotTree.wBars}. For \code{plotTree.barplot}, the only optional arguments are \code{add} and \code{ylim}. Generally \code{add} should not be used; however it can be employed to tell the function to draw the tree & barplot, respectively, in the next two open plotting devices - rather than creating a table of figures in the current plotting device. \code{ylim} (which is also an optional argument for \code{plotTree.boxplot} should be supplied here rather than using \code{args.plotTree} or \code{args.boxplot}\\code{args.barplot} because \emph{y} axis limits must match exactly between the two plots.}
}
\description{
	\code{plotTree.wbars} plots a phylogeny in phylogram or fan style with bars at the tips representing the values for a phenotypic trait.
	
	\code{plotTree.barplot} creates a split plot in which a right-facing phylogram is on the left, and a bar plot is shown on the right.
	
	\code{plotTree.boxplot} creates a split plot in which a right-facing phylogram is on the left, and a box plot is shown on the right.
}
\value{
	Plots a tree with an associated bar plot for a continuously valued character at the tips.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{barplot}}, \code{\link{dotTree}}, \code{\link{plotSimmap}}, \code{\link{plotTree}}
}
\examples{
## load data from Mahler et al. (2010)
data(anoletree)
data(anole.data)
## extract overall body size (SVL)
svl<-setNames(anole.data$SVL,rownames(anole.data))
## plotTree.wBars
plotTree.wBars(anoletree,svl,type="fan",scal=0.5)
par(mar=c(5.1,4.1,4.1,2.1))
## plotTree.barplot
plotTree.barplot(anoletree,exp(svl),
    args.plotTree=list(fsize=0.5),
    args.barplot=list(xlab="SVL (mm)"))
	
## load vertebrate tree and data
data(vertebrate.tree)
data(vertebrate.data)
## plotTree.barplot
options(scipen=3) ## change sci-notation
plotTree.barplot(vertebrate.tree,
	setNames(vertebrate.data$Mass,
	rownames(vertebrate.data)),
	args.barplot=list(log="x",
	xlab="mass (kg)",
	xlim=c(0.01,500000),
	col=palette()[4]))
options(scipen=0)

## reset par to defaults
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{continuous character}
