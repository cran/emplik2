\name{emplik2.internal}
\alias{myWCY}
\alias{myWKM}
\alias{myWdataclean2}
\alias{myWdataclean3}
\title{Internal emplik2 functions
}
\description{
These are internal functions called by \code{el2.test.wts} and \code{el2.test.wtm}.  They are not intended
 to be called by the user.
}
\usage{
myWCY(x, d, zc = rep(1, length(d)), wt = rep(1, length(d)), maxit = 25, error = 1e-09)
myWKM(x, d, zc = rep(1, length(d)), w = rep(1, length(d)))
myWdataclean2(z, d, wt = rep(1, length(z)))
myWdataclean3(z, d, zc = rep(1, length(z)), wt = rep(1, length(z)))
}
\details{
\code{WCY} calculates the weighted Chang-Yang self-consistent estimator for doubly-censored data.

\code{WKM} calculates the weighted Kaplan-Meier estimator for right-censored data.

\code{myWdataclean2} sorts the data, collapses the true ties, and puts the number of tied values as the weights.

\code{myWdataclean3} sorts the data, collapses the true ties, and puts the number of tied values as the weights.
}
\keyword{internal}
