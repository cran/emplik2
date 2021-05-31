\name{el2.cen.EMs}
\alias{el2.cen.EMs}
\title{
Computes p-value for a single mean-type hypothesis, based on two independent samples that may contain censored data.
}
\description{
This function uses the EM algorithm to calculate a maximized empirical likelihood ratio for the hypothesis
\deqn{
H_o:  E(g(x,y)-mean)=0 
}
where \eqn{E} indicates expected value; \eqn{g(x,y)} is a user-defined function of \eqn{x} and \eqn{y}; and
 \eqn{mean} is the hypothesized value of \eqn{E(g(x,y))}. The samples \eqn{x} and \eqn{y} are assumed independent.
 They may be uncensored, right-censored, left-censored, or left-and-right (``doubly'') censored. A p-value for 
 \eqn{H_o} is also calculated, based on the assumption that -2*log(empirical likelihood ratio) is approximately
 distributed as chisq(1).  
}
\usage{
el2.cen.EMs(x,dx,y,dy,fun=function(x,y){x>=y}, mean=0.5, maxit=25)
}
\arguments{
\item{x}{a vector of the data for the first sample}
\item{dx}{a vector of the censoring indicators for \code{x}: 0=right-censored, 1=uncensored, 2=left-censored}
\item{y}{a vector of the data for the second sample}
\item{dy}{a vector of the censoring indicators for \code{y}: 0=right-censored, 1=uncensored, 2=left-censored}
\item{fun}{a user-defined, continuous-weight-function \eqn{g(x,y)} used to define the mean in the hypothesis \eqn{H_o}.
  The default is \code{fun=function(x,y){x>=y}}.}
\item{mean}{the hypothesized value of \eqn{E(g(x,y))}; default is 0.5}
\item{maxit}{a positive integer used to set the number of iterations of the EM algorithm; default is 25}
}
\details{
The value of \eqn{mean} should be chosen between the maximum and minimum values of \eqn{g(x_i,y_j)}; otherwise
 there may be no distributions for \eqn{x} and \eqn{y} that will satisfy \eqn{H_o}.  If \eqn{mean} is inside 
 this interval, but the convergence is still not satisfactory, then the value of \eqn{mean} should be moved 
 closer to the NPMLE for \eqn{E(g(x,y))}.  (The NPMLE itself should always be a feasible value for \eqn{mean}.)
}
\value{
\code{el2.cen.EMs} returns a list of values as follows:

\item{xd1}{a vector of the unique, uncensored \eqn{x}-values in ascending order}

\item{yd1}{a vector of the unique, uncensored \eqn{y}-values in ascending order}

\item{temp3}{a list of values returned by the \code{el2.test.wts} function (which is called by \code{el2.cen.EMs})}

\item{mean}{the hypothesized value of \eqn{E(g(x,y))}}

\item{funNPMLE}{the non-parametric-maximum-likelihood-estimator of \eqn{E(g(x,y))}}

\item{logel00}{the log of the unconstrained empirical likelihood}

\item{logel}{the log of the constrained empirical likelihood}

\item{"-2LLR"}{\code{-2*(logel-logel00)}}

\item{Pval}{the estimated p-value for \eqn{H_o}, computed as \code{1 - pchisq(-2LLR, df = 1)}}

\item{logvec}{the vector of successive values of \code{logel} computed by the EM algorithm (should
 converge toward a fixed value)}

\item{sum_muvec}{sum of the probability jumps for the uncensored \eqn{x}-values, should be 1}

\item{sum_nuvec}{sum of the probability jumps for the uncensored \eqn{y}-values, should be 1}

\item{constraint}{the realized value of \eqn{\sum_{i=1}^n \sum_{j=1}^m (g(x_i,y_j) - mean) \mu_i \nu_j}, 
where \eqn{mu_i} and \eqn{nu_j} are the probability jumps at \eqn{x_i} and \eqn{y_j}, respectively, 
that maximize the empirical likelihood ratio. The value of \code{constraint} should be close to 0.} 
}
\references{

Barton, W. (2010). Comparison of two samples by a nonparametric likelihood-ratio test.  
PhD dissertation at University of Kentucky.

Chang, M. and Yang, G. (1987). ``Strong Consistency of a Nonparametric Estimator of the Survival Function
 with Doubly Censored Data.''  \code{Ann. Stat.},15, pp. 1536-1547.

Dempster, A., Laird, N., and Rubin, D. (1977). ``Maximum Likelihood from Incomplete Data via the EM Algorithm.'' \code{J. Roy. Statist. Soc.}, Series B, 39, pp.1-38.

Gomez, G., Julia, O., and Utzet, F. (1992). ``Survival Analysis for Left-Censored Data.'' In Klein, J. and Goel, P. (ed.),
 \code{Survival Analysis: State of the Art.} Kluwer Academic Publishers, Boston, pp. 269-288.

Li, G. (1995). ``Nonparametric Likelihood Ratio Estimation of Probabilities for Truncated Data.''
   \code{J. Amer. Statist. Assoc.}, 90, pp. 997-1003.

Owen, A.B. (2001). \code{Empirical Likelihood}. Chapman and Hall/CRC, Boca Raton, pp.223-227.

Turnbull, B. (1976). ``The Empirical Distribution Function with Arbitrarily Grouped, Censored and Truncated Data.''
   \code{J. Roy. Statist. Soc.}, Series B, 38, pp. 290-295.

Zhou, M. (2005). ``Empirical likelihood ratio with arbitrarily censored/truncated data by EM algorithm.''
   \code{J. Comput. Graph. Stat.}, 14, pp. 643-656.

Zhou, M. (2009) \code{emplik} package on CRAN website.  
  The \code{el2.cen.EMs} function extends \code{el.cen.EM}
  function from one-sample to two-samples.
}
\author{
William H. Barton <bbarton@lexmark.com>
}
\examples{ 
x<-c(10,80,209,273,279,324,391,415,566,785,852,881,895,954,1101,
1133,1337,1393,1408,1444,1513,1585,1669,1823,1941)
dx<-c(1,2,1,1,1,1,1,2,1,1,1,1,1,1,1,0,0,1,0,0,0,0,1,1,0)
y<-c(21,38,39,51,77,185,240,289,524,610,612,677,798,881,899,946,
1010,1074,1147,1154,1199,1269,1329,1484,1493,1559,1602,1684,1900,1952)
dy<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0)

# Ho1:  X is stochastically equal to Y
el2.cen.EMs(x, dx, y, dy, fun=function(x,y){x>=y}, mean=0.5, maxit=25)
# Result: Pval = 0.7090658, so we cannot with 95 percent confidence reject Ho1

# Ho2: mean of X equals mean of Y
el2.cen.EMs(x, dx, y, dy, fun=function(x,y){x-y}, mean=0.5, maxit=25)
# Result: Pval = 0.9695593, so we cannot with 95 percent confidence reject Ho2
}
\keyword{nonparametric}
