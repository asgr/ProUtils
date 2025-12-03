\name{quan_wt}
\alias{quan_wt}
\alias{quan_wt_mat_col}
\alias{quan_wt_mat_row}
\title{Compute Weighted Quantiles for a Numeric Vector of Matrix}
\description{
Computes one or more quantiles for a numeric vector or matrix, optionally using weights.
}
\usage{
quan_wt(x, q = as.numeric(c(0.5)), w = NULL)
quan_wt_mat_col(mat, q = as.numeric(c(0.5)), w = NULL)
quan_wt_mat_row(mat, q = as.numeric(c(0.5)), w = NULL)
}
\arguments{
  \item{x}{Numeric vector. Missing values (\code{NA}) are ignored).}
  \item{mat}{Numeric matrix. Missing values (\code{NA}) are ignored).}
  \item{q}{Numeric vector of probabilities in [0, 1]. Defaults to \code{0.5} (median).}
  \item{w}{An optional numeric vector of weights of the same length/dims as \option{x}/\option{mat}. If \code{NULL}, equal weights are assumed.}
}
\details{
For \code{quan_wt}, the function computes weighted quantiles using the midpoint rule for cumulative probabilities (see Examples). This rule behaves differently to no-weight \code{\link{quantile}} with the default \option{type} = 7, which is not to say it is wrong, there are just different ways to consider quantile treatment. The method used here is basically akin to \option{5} when not using weights. If weights are not provided, all weights (\option{w}) are assumed to be the same and internally valued 1. If all values are \code{NA} or weights sum to zero, \code{NA} is returned.

For \code{quan_wt_mat_col} and \code{quan_wt_mat_row} matrix functions, when \code{length(q) == 1} the functions return a numeric vector of length equal to the number of rows/columns (as relevant). If multiple quantiles are requested, the function returns a numeric matrix with:
\itemize{
  \item Rows corresponding to rows/columns of \code{mat} (as relevant).
  \item Columns corresponding to quantiles in \code{q}.
}
}
\value{
For \code{quan_wt}: A numeric vector of length equal to \code{length(q)} containing the requested quantiles. \cr

For \code{quan_wt_mat_col} and \code{quan_wt_mat_row}: A numeric vector (if one quantile requested) or a numeric matrix (if multiple quantiles requested).
}
\examples{
#vector examples
x = c(5, 1, 4, 7, 3)
quan_wt(x)                      # Median
quan_wt(x, c(0.25, 0.5, 0.75))  # Multiple quantiles

#We are the same as default quantile with option=5:
tempx = runif(100)
plot(quantile(tempx, seq(0,1,len=101), type = 5), quan_wt(tempx, seq(0,1,len=101)))
abline(0,1,col='red')

w = c(3, 5, 2, 1, 4)
out = quan_wt(x, c(0.25, 0.5, 0.75), w) # Weighted quantiles
print(out)

#to see explictly how we use weights:
x_ord = order(x)
w_mod = w[x_ord]
w_tot = w_mod[1]/2 + sum(w_mod[1:4] + w_mod[2:5])/2 + w_mod[5]/2
w_mod = c(w_mod[1]/2, (w_mod[2:5] + w_mod[1:4])/2)
plot(x[x_ord], cumsum(w_mod)/w_tot, type='l')
rect(0, 0, out[1], 0.25, lty=3)
rect(0, 0, out[2], 0.5, lty=3)
rect(0, 0, out[3], 0.75, lty=3)

#matrix examples
mat = matrix(1:12, nrow = 4)

#quantiles by column:
quan_wt_mat_col(mat)                      # Column medians
quan_wt_mat_col(mat, c(0.25, 0.5, 0.75))  # Multiple column quantiles
w = matrix(runif(12), nrow = 4)
quan_wt_mat_col(mat, c(0.25, 0.5, 0.75), w) # Weighted column quantiles

#quantiles by row:
quan_wt_mat_row(mat)                      # Row medians
quan_wt_mat_row(mat, c(0.25, 0.5, 0.75))  # Multiple row quantiles
quan_wt_mat_row(mat, c(0.25, 0.5, 0.75), w) # Weighted row quantiles
}
\seealso{
\code{\link{quantile}}, \code{\link{weighted.mean}}
}
