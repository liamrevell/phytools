\name{posthoc}
\alias{posthoc}
\title{Generic post-hoc test}
\usage{
posthoc(x, ...)
}
\arguments{
	\item{x}{an object on which to conduct a post-hoc test.}
	\item{...}{optional arguments to be passed to method.}
}
\description{
	Conducts posthoc test.
}
\details{
	So far is only implemented for object class \code{"ratebytree"}.
}
\value{
	An object of the appropriate class containing the results of a posthoc test.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ratebytree}}
}
\keyword{statistics}
