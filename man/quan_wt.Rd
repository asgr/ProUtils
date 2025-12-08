\name{quan_wt}
\alias{quan_wt}
\alias{quan_wt_mat_col}
\alias{quan_wt_mat_row}
\title{Compute Weighted Quantiles for a Numeric Vector of Matrix}
\description{
Computes one or more quantiles for a numeric vector or matrix, optionally using weights.
}
\usage{
quan_wt(x, probs = as.numeric(c(0.5)), wt = NULL, type = 7L)
quan_wt_mat_col(mat, probs = as.numeric(c(0.5)), wt = NULL, type = 7L)
quan_wt_mat_row(mat, probs = as.numeric(c(0.5)), wt = NULL, type = 7L)
}
\arguments{
  \item{x}{Numeric vector. Missing values (\code{NA}) are ignored).}
  \item{mat}{Numeric matrix. Missing values (\code{NA}) are ignored).}
  \item{probs}{Numeric vector of probabilities in [0, 1]. Defaults to \code{0.5} (median).}
  \item{wt}{An optional numeric vector of weights of the same length/dims as \option{x}/\option{mat}. If \code{NULL}, equal weights are assumed. If using errors, then \option{wt} should be the inverse variance.}
  \item{type}{Integer scalar between 4 and 9 selecting one of the continuous quantile algorithms detailed below to be used (see \code{\link{quantile}}). Default is 7, the same as base \code{\link{quantile}}.}
}
\details{
For \code{quan_wt}, the function computes weighted quantiles using different rules that capture the classic continuous quantile \option{type} (see Examples).

For \code{quan_wt_mat_col} and \code{quan_wt_mat_row} matrix functions, when \code{length(probs) == 1} the functions return a numeric vector of length equal to the number of columns/rows (as relevant). If multiple quantiles are requested, the function returns a numeric matrix with:
\itemize{
  \item Rows corresponding to rows/columns of \code{mat} (as relevant).
  \item Columns corresponding to quantiles in \code{probs}.
}

For clarity, if the matrix is [N,M] dimensions, then \code{quan_wt_mat_row} with defaults (so median) returns a vector length [N] and \code{quan_wt_mat_col} with defaults (so median) returns a vector length [M].
}
\value{
For \code{quan_wt}: A numeric vector of length equal to \code{length(probs)} containing the requested quantiles. If all values are \code{NA} or weights sum to zero, \code{NA} is returned. \cr

For \code{quan_wt_mat_col} and \code{quan_wt_mat_row}: A numeric vector (if one quantile requested) or a numeric matrix (if multiple quantiles requested). In the matrix case, the column order corresponds to the requested \option{probs} and rows are the outputs of \code{quan_wt}.
}
\examples{
#vector examples
x = c(5, 1, 4, 7, 3)
quan_wt(x)                      # Median
quan_wt(x, c(0.25, 0.5, 0.75))  # Multiple quantiles

#check how we compare to base quantile (only type 5:9 supported):
tempx = runif(100)
for(type in 4:9){
  plot(quantile(tempx, seq(0,1,len=101), type=type),
    quan_wt(tempx, seq(0,1,len=101), type=type), xlab='base', ylab='quan_wt')
  abline(0,1)
  legend('topleft', legend=type)
}

#now weighted quantiles:
#The weight free version should return 5:
quan_wt(0:10)

#But with weights 0:10 we should expect close to 10*sqrt(2)/2 (why is hwk for user!):
quan_wt(0:10, wt=0:10)

w = c(3, 5, 2, 1, 4)
quan_wt(x, c(0.25, 0.5, 0.75), w) # Weighted quantiles

#to see explictly how we use weights, this is how type = 4-9 are computed:
for(type in 4:9){
  out = quan_wt(x, c(0.25, 0.5, 0.75), w, type=type)
  start_mod = switch(type - 3, 1, 0.5, 1, 0, 0.6666667, 0.625)
  end_mod = switch(type - 3, 0, 0.5, 1, 0, 0.6666667, 0.625)
  x_ord = order(x)
  w_mod = w[x_ord]
  w_tot = w_mod[1]*start_mod + sum(w_mod[1:4] + w_mod[2:5])/2 + w_mod[5]*end_mod
  w_mod = c(w_mod[1]*start_mod, (w_mod[1:4] + w_mod[2:5])/2)

  plot(x[x_ord], cumsum(w_mod)/w_tot, type='l', xlim=c(0.5,7.5), ylim=c(0,1))
  rect(-1, -1, out[1], 0.25, lty=3)
  rect(-1, -1, out[2], 0.5, lty=3)
  rect(-1, -1, out[3], 0.75, lty=3)
  legend('topleft', legend=paste('type:', type))
}

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

#these should be the same:
set.seed(666); temp = matrix(rnorm(1e6),1e5,10)
set.seed(666); temp2 = matrix(rnorm(1e6),10,1e5,byrow=T)
all(quan_wt_mat_row(temp) == quan_wt_mat_col(temp2))
}
\seealso{
\code{\link{quantile}}, \code{\link{weighted.mean}}
}
