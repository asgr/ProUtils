\name{median_stream}
\alias{addVector}
\alias{getMedian}
\alias{clearHeaps}
\alias{getMinHeap}
\alias{getMaxHeap}
\title{
Maintain a Running Median Using Heaps
}
\description{
These functions allow you to maintain a running median of numeric values using two heaps (a max-heap and a min-heap). You can add values incrementally, retrieve the current median, inspect the heap contents, and clear the state.
}
\usage{
addVector(vec)
getMedian()
clearHeaps()
getMinHeap()
getMaxHeap()
}
\arguments{
  \item{vec}{
A numeric vector of values to add to the running median calculation.
}
}
\details{
The algorithm uses two heaps:
\itemize{
  \item A max-heap for the lower half of the data.
  \item A min-heap for the upper half of the data.
}
This allows efficient computation of the median after each insertion without sorting the entire dataset.

\code{addVector()} adds all elements of \code{vec} to the heaps.
\code{getMedian()} returns the current median based on all values added so far.
\code{clearHeaps()} resets both heaps, removing all stored values.
\code{getMinHeap()} and \code{getMaxHeap()} return the contents of the respective heaps as numeric vectors (sorted by heap order, not overall order).
}
\value{
\itemize{
  \item \code{addVector()} returns \code{NULL} (invisible).
  \item \code{getMedian()} returns a single numeric value (the current median) or \code{NA} if no values have been added.
  \item \code{clearHeaps()} returns \code{NULL} (invisible).
  \item \code{getMinHeap()} and \code{getMaxHeap()} return numeric vectors of heap contents.
}
}
\author{
Your Name
}
\note{
These functions maintain global state in C++ heaps. If you need multiple independent median trackers, consider using an Rcpp module or an R6 wrapper instead.
}
\seealso{
\code{\link{median}}, \code{\link{quantile}}
}
\examples{
# Add values to the heaps
addVector(c(1, 3, 5, 7))
getMedian()  # Returns 4

# Add more values
addVector(c(10, 2))
getMedian()  # Returns 4 (balanced heaps)

# Inspect heaps
getMaxHeap()
getMinHeap()

# Clear heaps
clearHeaps()
getMedian()  # Returns NA
}
\keyword{statistics}
\keyword{median}
