\name{branching.diffusion}
\alias{branching.diffusion}
\title{Animation of branching random diffusion}
\usage{
branching.diffusion(sig2=1, b=0.0023, time.stop=1000, ylim=NULL, 
   smooth=TRUE, pause=0.02, record=NULL, path=NULL)
}
\arguments{
	\item{sig2}{variance of BM process, \eqn{\sigma^2}.}
	\item{b}{birthrate for branching process.}
	\item{time.stop}{number of generations to run.}
	\item{ylim}{\emph{y} limits (for plotting).}
	\item{smooth}{no longer used.}
	\item{pause}{pause (in s) between generations.}
	\item{record}{filename for video file output (no video if \code{NULL}).}
	\item{path}{full path to file for video rendering. (By default \code{branching.diffusion} will look for the executable \code{ffmpeg.exe} in the directory \code{C:/Program Files/ffmpeg/bin}, even though this will not make sense on non-Windows machines.)}
}
\description{
	This function creates an animation of branching random diffusion (i.e., BM with speciation).
}
\value{
	An animated plot and (optionally) a recorded video file. For animation to be recorded to file, the function requires the package \emph{animation} as well as a video renderer.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{bmPlot}}, \code{\link{fastBM}}
}
\keyword{phylogenetics}
\keyword{animation}
\keyword{plotting}
\keyword{simulation}
\keyword{continuous character}
