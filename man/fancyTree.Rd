\name{fancyTree}
\alias{fancyTree}
\alias{phyloScattergram}
\alias{phenogram95}
\title{Plots special types of phylogenetic trees}
\usage{
fancyTree(tree, type=c("extinction","traitgram3d","droptip","densitymap",
   "contmap","phenogram95","scattergram"), ..., control=list())
phyloScattergram(tree, X=NULL, ...)
phenogram95(tree, x=NULL, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{type}{the type of special plot to create. See Description.}
	\item{...}{arguments to be passed to different methods. See Description.}
	\item{control}{a list of control parameters, depending on \code{type}.}
	\item{X}{in \code{phyloScattergram}, a matrix of continuous trait values. Row names in the matrix should correspond to species names in the tree.}
	\item{x}{in \code{phenogram95}, a named vector with values for a continuously distributed trait.}
}
\description{
	Plots phylogenies (or phylogenetic trees and comparative data) in a variety of different styles.
}
\details{
	This function plots a phylogeny or phylogenetic tree and comparative data in a variety of different styles, depending on the value of \code{type}. In some instances, \code{fancyTree} is now just a wrappe for other \pkg{phytools} functions, such as \code{\link{contMap}} and \code{\link{densityMap}}.
	
	If \code{type="extinction"} (or any unambiguous abbreviation) the function will plot a tree in which branches preceding the MRCA of all extant taxa and branches leading only to extnct lineages are plotted with dashed red lines. 
	
	If \code{type="traitgram3d"} the function will plot a three dimensional traitgram (that is, a projection of the tree into three dimensional morphospace where two dimensions are the phenotypic trait and the third axis is time since the root). In this case, the additional argument \code{X}, a matrix containing the tip values of all species (with species IDs as row names) should be supplied. Optionally, the user can also supply the matrix \code{A}, which contains the ancestral states in the tree with rows labeled by node number. 
	
	If \code{type="droptip"} the function will create a two panel figure in which the first panel is the tree with lineages to be pruned highlighted; and the second panel is the pruned tree. In this case, the additional argument \code{tip}, the tip name or vector of tip names to be dropped, must be supplied. 
	
	If \code{type="densitymap"}, a posterior probability density "heat-map" is created based on a set of trees in a \code{"multiSimmap"} object containing a binary [0,1] mapped character. (See \code{\link{densityMap}} for additional optional arguments if \code{type="densitymap"}.)
	
	If \code{type="contmap"}, reconstructed continuous trait evolution is mapped on the tree. Again, see \code{\link{contMap}} for additional arguments if \code{type="contmap"}.

	If \code{type="phenogram95"} a 95-percent phenogram is plotted using transparency to visualize uncertainty at ancestral nodes and along branches. Most of the options of \code{\link{phenogram}} are available. 
	
	Finally, if \code{type="scattergram"} a phylogenetic scatter plot matrix containing \code{\link{contMap}} style trees on the diagonal and \code{\link{phylomorphospace}} plots in non-diagonal panels is produced. For this type a trait matrix \code{X} must also be supplied. The only additional arguments available for this type are \code{ftype}, \code{fsize}, \code{colors}, and \code{label}. (See \code{\link{phylomorphospace}} for details on how these arguments should be used.) This function calls \code{\link{phyloScattergram}} (which is also now exported to the name space) internally. In addition to creating a plot, \code{phyloScattergram} also returns an object of class \code{"phyloScattergram"} which can be replotted using different options if desired.

	Presently only \code{type="traitgram3d"} uses the list \code{control} which can be supplied the same set of control parameters as \code{\link{phylomorphospace3d}}, as well as the control parameter \code{maxit} which will be passed to \code{\link{anc.ML}}. 
	
	Finally, the optional argument \code{hold} will be passed to multiple methods if supplied. It is a logical value that indicates whether or not the output to the graphical device should be held using \code{\link{dev.hold}} before plotting (defaults to \code{hold=TRUE}).
}
\value{
	This function plots different types of phylogenetic trees. For \code{type="droptip"} the function also returns the pruned tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{contMap}}, \code{\link{densityMap}}, \code{\link{drop.tip}}, \code{\link{phenogram}}, \code{\link{phylomorphospace3d}}, \code{\link{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
# plot tree with extinction
set.seed(10)
tree<-pbtree(b=1,d=0.4,t=4)
fancyTree(tree,type="extinction")

\dontrun{
# plot 3D traitgram
tree<-pbtree(n=50,scale=10)
Y<-sim.corrs(tree,vcv=matrix(c(1,0.75,0.75,1),2,2))
fancyTree(tree,type="traitgram3d",X=Y,
    control=list(spin=FALSE))
	
# plot with internal nodes from simulation
Y<-sim.corrs(tree,vcv=matrix(c(1,0.75,0.75,1),
    2,2),internal=TRUE)
B<-Y[length(tree$tip)+1:tree$Nnode,]
Y<-Y[1:length(tree$tip),]
fancyTree(tree,type="traitgram3d",X=Y,A=B,
    control=list(simple.axes=TRUE,spin=FALSE))}
	
# plot with dropped tips
tree<-pbtree(n=30)
tips<-sample(tree$tip.label)[1:10]
pruned<-fancyTree(tree,type="droptip",tip=tips)

\dontrun{
# plot 95-percent CI phenogram
tree<-pbtree(n=30)
x<-fastBM(tree)
fancyTree(tree,type="phenogram95",x=x)}

## reset par to defaults
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
