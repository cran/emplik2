\name{el2.test.wts}
\alias{el2.test.wts}
\title{
Computes maximium-likelihood probability jumps for a single mean-type hypothesis, based on two independent uncensored samples}
\description{
This function computes the maximum-likelihood probability jumps for a single mean-type hypothesis, based on two samples
 that are independent, uncensored, and weighted.  The target function for the maximization is the constrained
log(empirical likelihood) which can be expressed as,
\deqn{\sum_{dx_i=1} wx_i \log{\mu_i} + \sum_{dy_j=1} wy_j \log{\nu_j} - \eta ( 1 - \sum_{dx_i=1} \mu_i ) - \delta
( 1 -\sum_{dy_j=1} \nu_j ) - \lambda \sum_{dx_i=1} \sum_{dy_j=1} ( g(x_i,y_j)- mean ) \mu_i \nu_j }

where the variables are defined as follows:

\eqn{x} is a vector of data for the first sample

\eqn{y} is a vector of data for the second sample

\eqn{wx} is a vector of estimated weights for the first sample

\eqn{wy} is a vector of estimated weights for the second sample

\eqn{\mu} is a vector of estimated probability jumps for the first sample

\eqn{\nu} is a vector of estimated probability jumps for the second sample
}
\usage{
el2.test.wts(u,v,wu,wv,mu0,nu0,indicmat,mean,lamOld=0)
}
\arguments{
\item{u}{a vector of uncensored data for the first sample}
\item{v}{a vector of uncensored data for the second sample}
\item{wu}{a vector of estimated weights for \code{u}}
\item{wv}{a vector of estimated weights for \code{v}}
\item{mu0}{a vector of estimated probability jumps for \code{u}}
\item{nu0}{a vector of estimated probability jumps for \code{v}}
\item{indicmat}{a matrix \eqn{[g(u_i,v_j)-mean]} where \eqn{g(u, v)} is a user-chosen function}
\item{mean}{a hypothesized value of \eqn{E(g(u,v))}, where \eqn{E} indicates ``expected value.''}
\item{lamOld}{The previous solution of lambda, used as the starting point to search for new solution of lambda.}
}
\details{
This function is called by \code{el2.cen.EMs}.  It is listed here because the user may find it useful elsewhere.

The value of \eqn{mean} should be chosen between the maximum and minimum values of
 \eqn{(u_i,v_j)}; otherwise there may be no distributions for \eqn{u} and \eqn{v} that
 will satisfy the the mean-type hypothesis.  If \eqn{mean} is inside this interval, but the convergence is
 still not satisfactory, then the value of \eqn{mean} should be moved closer to the NPMLE for \eqn{E(g(u,v))}.
(The NPMLE itself should always  be a feasible value for \eqn{mean}.)  The calculations for this function
 are derived in Owen (2001).
}
\value{
\code{el2.test.wts} returns a list of values as follows:

\item{u}{the vector of uncensored data for the first sample}
\item{wu}{the vector of weights for \eqn{u}}
\item{jumpu}{the vector of probability jumps for \eqn{u} that maximize the weighted empirical likelihood}
\item{v}{the vector of uncensored data for the second sample}
\item{wv}{the vector of weights for \eqn{v}}
\item{jumpv}{the vector of probability jumps for \eqn{v} that maximize the weighted empirical likelihood}
\item{lam}{the value of the Lagrangian multipler found by the calculations}
}
\references{
Owen, A.B. (2001). \code{Empirical Likelihood}. Chapman and Hall/CRC, Boca Raton, pp.223-227.
}
\author{
William H. Barton <bbarton@lexmark.com> and modified by Mai Zhou.
}
\examples{ 
u<-c(10, 209, 273, 279, 324, 391, 566, 785)
v<-c(21, 38, 39, 51, 77, 185, 240, 289, 524)
wu<-c(2.442931, 1.122365, 1.113239, 1.113239, 1.104113, 1.104113, 1.000000, 1.000000)
wv<-c( 1, 1, 1, 1, 1, 1, 1, 1, 1)
mu0<-c(0.3774461, 0.1042739, 0.09649724, 0.09649724, 0.08872055, 0.08872055, 0.0739222, 0.0739222)
nu0<-c(0.1013718, 0.1013718, 0.1013718, 0.1013718, 0.1013718, 0.1013718, 0.1095413, 0.1287447,
 0.1534831)
mean<-0.5

#let fun=function(x,y){x>=y}
indicmat<-matrix(nrow=8,ncol=9,c(
-0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5, 
-0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
-0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
-0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
-0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
-0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
-0.5, -0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,
-0.5, -0.5, -0.5, -0.5,  0.5,  0.5,  0.5,  0.5,
-0.5, -0.5, -0.5, -0.5, -0.5, -0.5,  0.5,  0.5))
el2.test.wts(u,v,wu,wv,mu0,nu0,indicmat,mean)

# jumpu
# [1] 0.3774461, 0.1042739, 0.09649724, 0.09649724, 0.08872055, 0.08872055, 0.0739222, 0.0739222

# jumpv
# [1] 0.1013718, 0.1013718, 0.1013718, 0.1013718, 0.1013718, 0.1013718, 0.1095413, 0.1287447,
# [9] 0.1534831

# lam
# [1] 7.055471
}
\keyword{nonparametric}
