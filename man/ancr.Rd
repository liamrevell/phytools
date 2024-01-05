\name{ancr}
\alias{ancr}
\alias{hide.hidden}
\alias{plot.ancr}
\title{Compute marginal or joint ancestral state estimates}
\usage{
ancr(object, ...)
hide.hidden(object, ...)
\method{plot}{ancr}(x, args.plotTree=list(...), args.nodelabels=list(...), ...)
}
\arguments{
	\item{object}{fitted M\emph{k} model (presently object of class \code{"fitMk"}, \code{"fitpolyMk"}, and \code{"fitHRM"}), or a set of models in the form of a table from \code{anova} comparison. For \code{hide.hidden}, object of class \code{"ancr"}.}
	\item{...}{optional arguments.}
	\item{x}{in the case of \code{plot.ancr}, an object of class \code{"ancr"}.}
	\item{args.plotTree}{arguments to be passed to \code{\link{plotTree}}, in a list.}
	\item{args.nodelabels}{arguments to be passed to \code{\link{nodelabels}}, in a list.}
}
\description{
	By default, \code{ancr} computes marginal ancestral states, also known as empirical Bayes posterior probabilities, conditioning on the fitted (or set) model of \code{object}.
	
	Can also perform \emph{joint} ancestral state estimation, if the optional argument \code{type} is set to \code{type="joint"}.
	
	\code{hide.hidden} merges hidden states (if any).
}
\details{
	If the optional argument \code{tips=TRUE}, then the matrix returned contains empirical Bayes posterior probabilities (marginal scaled likelihoods) for both tips \emph{and} internal nodes. Otherwise (the default) only node marginal states are returned.
	
	If the input object is a set of models (in the form of an \code{anova} table), then \code{ancr} will compute model-averaged marginal ancestral states (for \code{type="marginal"}, unless the optional argument \code{weighted=FALSE}, in which case only the best-supported model is used.
}
\value{
	An object of class \code{"ancr"}. In the case of \code{type="marginal"}, this object consists of a matrix of marginal (empirical Bayes) probabilities and a likelihood. In the case of \code{type="joint"}, the object contains a set of most-likely internal node states stored in a data frame.
}
\references{
	Pagel, M. (1999) The Maximum Likelihood approach to reconstructing ancestral character states of discrete characters on phylogenies. \emph{Systematic Biology}, \bold{3}, 612-622.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
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
## fit ARD model
ard_fmode<-fitMk(sunfish.tree,fmode,model="ARD",
  pi="fitzjohn")
## compute ancestral states
anc_fmode<-ancr(ard_fmode)
## plot the results
par(mfrow=c(2,1))
cols<-setNames(c("blue","red"),levels(fmode))
plot(anc_fmode,
  args.plotTree=list(lwd=2,direction="upwards",
    mar=c(0.1,1.1,2.1,1.1),fsize=0.8),
  args.nodelabels=list(piecol=cols),
  args.tiplabels=list(cex=0.3),
  legend="bottomright")
mtext("a) marginal states under ARD model",adj=0)
## fit ER model
er_fmode<-fitMk(sunfish.tree,fmode,model="ER",
  pi="fitzjohn")
## compare models
aov_fmode<-anova(er_fmode,ard_fmode)
## compute model-averaged ancestral states
anc_fmode_model.averaged<-ancr(aov_fmode)
plot(anc_fmode_model.averaged,
  args.plotTree=list(lwd=2,direction="upwards",
    mar=c(0.1,1.1,2.1,1.1),fsize=0.8),
  args.nodelabels=list(piecol=cols),
  args.tiplabels=list(cex=0.3),
  legend="bottomright")
mtext("b) marginal states model-averaging ER & ARD models",
  adj=0)
par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
