\name{fitThresh}
\alias{fitThresh}
\title{Fits multi-state threshold model using the discrete approximation}
\usage{
fitThresh(tree, x, sequence=NULL, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{a named vector.}
	\item{sequence}{assumed ordering of the discrete character state. If not supplied and then alphanumeric order is assumed.}
	\item{...}{optional arguments. Optional arguments include \code{levs} (the number of levels for the discrete approximation, defaults to \code{levs=200}), \code{root} (a root prior assumption, defaults to \code{root="fitzjohn"} which involves treating the root like a nuisance parameter), and \code{rand_start} (defaults to \code{rand_start=TRUE}).}
}
\description{
	This function uses the discrete approximation of Boucher and Demery (2016) to fit a multi-state discrete character threshold model following Felsenstein (2012; Revell 2014).
}
\value{
	An object of class \code{"fitThresh"}.
}
\references{
	Boucher, F. C., and V. Demery (2016) Inferring bounded evolution in phenotypic characters from phylogenetic comparative data. \emph{Systematic Biology}, \bold{65}, 651-661.

	Felsenstein, J. (2012) A comparative method for both discrete and continuous characters using the threshold model. \emph{American Naturalist}, \bold{179}, 145-156.

	Revell, L. J. (2014) Ancestral character estimation under the threshold model from quantitative genetics. \emph{Evolution}, \bold{68}, 743-759.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ancThresh}}, \code{\link{bounded_bm}}, \code{\link{fitmultiBM}}, \code{\link{threshBayes}}
}
\examples{
\dontrun{
## load data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)

## extract character of interest
fmode<-setNames(sunfish.data$feeding.mode,
    rownames(sunfish.data))

## fit model
sunfish_thresh<-fitThresh(sunfish.tree,fmode,
    levs=100)
sunfish_thresh

## fit ARD model for comparison
sunfish_ard<-fitMk(sunfish.tree,fmode,model="ARD")
sunfish_ard

## compare results
anova(sunfish_thresh,sunfish_ard)

## estimate ancestral states under threshold model
anc_sunfish<-ancr(sunfish_thresh)

## plot results
plot(anc_sunfish,mar=c(0.1,0.1,4.1,0.1),
    args.plotTree=list(lwd=2,type="cladogram",
        nodes="inner",offset=1),
    args.nodelabels=list(piecol=hcl.colors(n=2),
        cex=0.6),
    args.tiplabels=list(cex=0.5),legend=FALSE)
legend("topleft",c("non-piscivorous","piscivorous"),
    pch=16,col=hcl.colors(n=2),pt.cex=2,bty="n")
	
## reset par to defaults
par(mar=c(5.1,4.1,4.1,2.1),cex.axis=1)}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{discrete character}
