\name{ancr}
\alias{ancr}
\alias{hide.hidden}
\title{Compute marginal ancestral states}
\usage{
ancr(object, ...)
hide.hidden(object, ...)
}
\arguments{
	\item{object}{fitted M\emph{k} model (presently object of class \code{"fitMk"}, \code{"fitpolyMk"}, and \code{"fitHRM"}). For \code{hide.hidden}, object of class \code{"ancr"}.}
	\item{...}{optional arguments.}
}
\description{
	Computes marginal ancestral states, also known as empirical Bayes posterior probabilities, conditioning on the fitted (or set) model of \code{object}.
	
	\code{hide.hidden} merges hidden states (if any).
}
\details{
	If the optional argument \code{tips=TRUE}, then the matrix returned contains empirical Bayes posterior probabilities (marginal scaled likelihoods) for both tips \emph{and} internal nodes. Otherwise (the default) only node marginal states are returned.
}
\value{
	An object of class \code{"ancr"}, which includes a matrix of marginal (empirical Bayes) probabilities and a likelihood.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitMk}}
}
\examples{
## load tree and data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)
## extract discrete character (feeding mode)
fmode<-setNames(sunfish.data$feeding.mode,
    rownames(sunfish.data))
## fit model
ard_fmode<-fitMk(sunfish.tree,fmode,model="ARD",
    pi="fitzjohn")
## compute ancestral states
anc_fmode<-ancr(ard_fmode)
## plot the results
cols<-setNames(c("blue","red"),levels(fmode))
plotTree(sunfish.tree,ftype="i",offset=0.5)
nodelabels(pie=anc_fmode$ace,piecol=cols,cex=0.6)
tiplabels(pie=to.matrix(fmode[sunfish.tree$tip.label],
    levels(fmode)),piecol=cols,cex=0.4)
legend("topleft",levels(fmode),pch=21,pt.cex=2,
    pt.bg=cols,bty="n")
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
