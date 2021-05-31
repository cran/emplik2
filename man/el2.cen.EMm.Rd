\name{el2.cen.EMm}
\alias{el2.cen.EMm}
\title{
Computes p-value for multiple mean-type hypotheses, based on two independent samples that may contain censored data.
}
\description{
This function uses the EM algorithm to calculate a maximized empirical likelihood ratio for a set of \eqn{p} hypotheses
as follows:
\deqn{
H_o:  E(g(x,y)-mean)=0 
}
where \eqn{E} indicates expected value; \eqn{g(x,y)} is a vector of user-defined functions \eqn{g_1(x,y), \ldots,
g_p(x,y)}; and \eqn{mean} is a vector of \eqn{p} hypothesized values of \eqn{E(g(x,y))}. The two samples \eqn{x} and \eqn{y}
are assumed independent. They may be uncensored, right-censored, left-censored, or left-and-right (``doubly'') 
censored. A p-value for \eqn{H_o} is also calculated, based on the assumption that -2*log(empirical likelihood ratio)
is asymptotically distributed as chisq(p). 
}
\usage{
el2.cen.EMm(x, dx, wx, y, dy, wy, p, H, xc=1:length(x), yc=1:length(y),
  mean, maxit=15)
}
\arguments{
\item{x}{a vector of the data for the first sample}
\item{dx}{a vector of the censoring indicators for x: 0=right-censored, 
1=uncensored, 2=left-censored}
\item{wx}{a vector of data case weight for x}
\item{y}{a vector of the data for the second sample}
\item{dy}{a vector of the censoring indicators for y: 0=right-censored, 
1=uncensored, 2=left-censored}
\item{wy}{a vector of data case weight for y}
\item{p}{the number of hypotheses}
\item{H}{a matrix defined as \eqn{H = [H_1, H_2, \ldots, H_p]}, where 
\eqn{H_k = [g_k(x_i,y_j)-mu_k], k=1, \ldots, p}}
\item{xc}{a vector containing the indices of the \code{x} datapoints, controls if tied x collapse or not}
\item{yc}{a vector containing the indices of the \code{y} datapoints, ditto}
\item{mean}{the hypothesized value of \eqn{E(g(x,y)})}
\item{maxit}{a positive integer used to control the maximum number of iterations of the EM algorithm; default is 15}
}
\details{
The value of \eqn{mean_k} should be chosen between the maximum and minimum values of \eqn{g_k(x_i,y_j)}; otherwise
 there may be no distributions for \eqn{x} and \eqn{y} that will satisfy \eqn{H_o}.  If \eqn{mean_k} is inside 
 this interval, but the convergence is still not satisfactory, then the value of \eqn{mean_k} should be moved 
 closer to the NPMLE for \eqn{E(g_k(x,y))}.  (The NPMLE itself should always be a feasible value for \eqn{mean_k}.)
}
\value{
\code{el2.cen.EMm} returns a list of values as follows:

 
\item{xd1}{a vector of unique, uncensored \eqn{x}-values in ascending order}

\item{yd1}{a vector of unique, uncensored \eqn{y}-values in ascending order}

\item{temp3}{a list of values returned by the \code{el2.test.wtm} function (which is called by \code{el2.cen.EMm})}

\item{mean}{the hypothesized value of \eqn{E(g(x,y))}}
\item{NPMLE}{a non-parametric-maximum-likelihood-estimator vector of \eqn{E(g(x,y))}}

\item{logel00}{the log of the unconstrained empirical likelihood}

\item{logel}{the log of the constrained empirical likelihood}

\item{"-2LLR"}{-2*(log-likelihood-ratio) for the \code{p} simultaneous hypotheses}

\item{Pval}{the p-value for the \code{p} simultaneous hypotheses, equal to \code{1 - pchisq(-2LLR, df = p)}}

\item{logvec}{the vector of successive values of \code{logel} computed by the EM algorithm
 (should converge toward a fixed value)}

\item{sum_muvec}{sum of the probability jumps for the uncensored \eqn{x}-values, should be 1}

\item{sum_nuvec}{sum of the probability jumps for the uncensored \eqn{y}-values, should be 1}
}
\references{

Barton, W. (2010). Comparison of two samples by a nonparametric likelihood-ratio test.  
PhD dissertation at University of Kentucky.

Chang, M. and Yang, G. (1987). ``Strong Consistency of a Nonparametric Estimator of the Survival Function with Doubly Censored Data.''  \code{Ann. Stat.},15, pp. 1536-1547.

Dempster, A., Laird, N., and Rubin, D. (1977). ``Maximum Likelihood from Incomplete Data via the EM Algorithm.'' \code{J. Roy. Statist. Soc.}, Series B, 39, pp.1-38.

Gomez, G., Julia, O., and Utzet, F. (1992). ``Survival Analysis for Left-Censored Data.'' In Klein, J. and Goel, P. (ed.),
 \code{Survival Analysis: State of the Art.} Kluwer Academic Publishers, Boston, pp. 269-288.

Li, G. (1995). ``Nonparametric Likelihood Ratio Estimation of Probabilities for Truncated Data.''
   \code{J. Amer. Statist. Assoc.}, 90, pp. 997-1003.

Owen, A.B. (2001). \code{Empirical Likelihood}. Chapman and Hall/CRC, Boca Raton, pp. 223-227.

Turnbull, B. (1976). ``The Empirical Distribution Function with Arbitrarily Grouped, Censored and Truncated Data.''
   \code{J. Roy. Statist. Soc.}, Series B, 38, pp. 290-295.

Zhou, M. (2005). ``Empirical likelihood ratio with arbitrarily censored/truncated data by EM algorithm.''
   \code{J. Comput. Graph. Stat.}, 14, pp. 643-656.

Zhou, M. (2009) \code{emplik} package on CRAN website. 
     The function \code{el2.cen.EMm} here extends \code{el.cen.EM2} 
	 inside emplik from one-sample to two-samples.
}
\author{
William H. Barton <bbarton@lexmark.com>
}
\examples{ 
x<-c(10, 80, 209, 273, 279, 324, 391, 415, 566, 85, 852, 881, 895, 954, 1101, 1133,
1337, 1393, 1408, 1444, 1513, 1585, 1669, 1823, 1941)
dx<-c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0)
y<-c(21, 38, 39, 51, 77, 185, 240, 289, 524, 610, 612, 677, 798, 881, 899, 946, 1010,
1074, 1147, 1154, 1199, 1269, 1329, 1484, 1493, 1559, 1602, 1684, 1900, 1952)
dy<-c(1,1,1,1,1,1,2,2,1,1,1,1,1,2,1,1,1,1,1,1,0,0,1,1,0,0,1,0,0,0)
nx<-length(x)
ny<-length(y)
xc<-1:nx
yc<-1:ny
wx<-rep(1,nx)
wy<-rep(1,ny)
mu=c(0.5,0.5)
p <- 2
H1<-matrix(NA,nrow=nx,ncol=ny)
H2<-matrix(NA,nrow=nx,ncol=ny)
for (i in 1:nx) {
  for (j in 1:ny) {
   H1[i,j]<-(x[i]>y[j])
   H2[i,j]<-(x[i]>1060) } }
H=matrix(c(H1,H2),nrow=nx,ncol=p*ny)

# Ho1: X is stochastically equal to Y
# Ho2: mean of X equals mean of Y

el2.cen.EMm(x=x, dx=dx, y=y, dy=dy, p=2, H=H, mean=mu, maxit=10)

# Result: Pval is 0.6310234, so we cannot with 95 percent confidence reject the two
# simultaneous hypotheses Ho1 and Ho2
}
\keyword{nonparametric}
