% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calkmeans.R
\name{calkmeans}
\alias{calkmeans}
\title{Calibrating kmeans}
\usage{
calkmeans(
  img,
  k = NULL,
  iter.max = 10,
  algo = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
  iter = 30,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{img}{RasterStack or RasterBrick.}

\item{k}{This argument is \code{NULL} when the objective is to obtain the bes \code{k} value. If the objective
is to select the best algorithm embedded in kmeans, please especify a \code{k} value.}

\item{iter.max}{The maximum number of iterations allowed. See \link[stats]{kmeans}.}

\item{algo}{It can be "Hartigan-Wong", "Lloyd", "Forgy" or "MacQueen". See \link[stats]{kmeans}.}

\item{iter}{Iterations number to obtain the best k value. \code{iter} must be greater than the number of classes
expected to be obtained in the classification. Default is 30.}

\item{verbose}{This paramater is Logical. It Prints progress messages during execution.}

\item{...}{Options to be passed to the function. See 'Details'.}
}
\description{
This function allows to calibrate the kmeans algorithm. It is possible to obtain
the best k value and the best embedded algorithm in kmeans.
}
\section{References}{

Tarazona, Y., Maria, Miyasiro-Lopez. (2020). Monitoring tropical forest degradation using
remote sensing. Challenges and opportunities in the Madre de Dios region, Peru. Remote
Sensing Applications: Society and Environment, 19, 100337.

Gareth James, Daniela Witten, Trevor Hastie, Robert Tibshirani. (2013).
An introduction to statistical learning : with applications in R. New York: Springer.
}

\section{Details}{
 If we want to find the optimal value of \code{k} (clusters or classes),
so we must put \code{k = NULL} as an argument of the function. Here, we are finding k for
which the intra-class inertia is stabilized. If we know the \code{k} value and the idea
is to find the best algorithm embedded in kmeans, that maximizes inter-class distances,
we must put \code{k = n}, where \code{n} is a specific class number.
}

\examples{
library(ForesToolboxRS)

# Load the dataset
data(FTdata)

# Selecting the best k value
best_k <- calkmeans(img = image[[1:2]], k = NULL, iter.max = 10,
           algo = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), iter = 30)
# Jambu Elbow
plot(best_k)


# Selecting the best embedded algorithm in kmeans
best_algo <- calkmeans(img = image[[1:2]], k = 4, iter.max = 10,
           algo = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), iter = 30)

# Choose the algorithm with the highest value
best_algo

}
\author{
Yonatan Tarazona
}
