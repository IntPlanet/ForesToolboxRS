% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calmla.R
\name{calmla}
\alias{calmla}
\title{Calibrating Supervised classification in Remote Sensing}
\usage{
calmla(
  img,
  endm,
  model = c("svm", "randomForest", "naiveBayes", "LMT", "nnet", "knn"),
  training_split = 50,
  approach = "Set-Approach",
  k = 5,
  iter = 10,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{img}{RasterStack or RasterBrick.}

\item{endm}{SpatialPointsDataFrame or SpatialPolygonsDataFrame (typically shapefile)
containing the training data.}

\item{model}{Model to use. It can be Support Vector Machine (\link[e1071]{svm}) like
\code{model = 'svm'}, Random Forest (\link[randomForest]{randomForest})
like \code{model = 'randomForest'}, Naive Bayes (\link[e1071]{naiveBayes})
like \code{model = 'naiveBayes'}, Decision Tree (\link[caret]{train})
like \code{model = 'LMT'}, Neural Networks (\link[nnet]{nnet})
like \code{model = 'nnet'}, K-nearest Neighbors (\link[caret]{knn3}) like \code{model = 'knn'}.}

\item{training_split}{For splitting samples into two subsets, i.e. training data and
for testing data.}

\item{approach}{Calibration method. There are for options: Simple training and testing
(Set-Approach) like \code{approach = 'Set-Approach'}, Leave One Out Cross-Validation (LOOCV) like
\code{approach = 'LOOCV'}, Cross-Validation (K-fold) like \code{approach = 'K-fold'} and
Monte Carlo Cross-Validation (MCCV) like \code{approach = 'MCCV'}.}

\item{k}{Number of groups for splitting samples. It must be used only with the
Cross-Validation (k-fold) approach.}

\item{iter}{Number of iterations, i.e number of times the analysis is executed.}

\item{verbose}{This paramater is Logical. It Prints progress messages during execution.}
}
\description{
This function allows to calibrate supervised classification in satellite images
through various algorithms and using approches such as Set-Approach,
Leave-One-Out Cross-Validation (LOOCV), Cross-Validation (k-fold) and
Monte Carlo Cross-Validation (MCCV).
}
\details{
If the "Set-Approach" method is being used, it is not necessary to use parameter \code{k}.
 \code{k} only can be used when the Cross-Validation (k-fold) method is used. On the other hand,
 to create groups in Cross-Validation, the \code{createFolds} function of the "confusionMatrix" is used.
 See \link[confusionMatrix]{createFolds} for more details. In addition, to generate random splits
 in Monte Carlos Cross-Validation the \code{generate.split} function of the "WilcoxCV" package was used.
 Please see \link[WilcoxCV]{generate.split} for more details.
}
\section{References}{

Gareth James, Daniela Witten, Trevor Hastie, Robert Tibshirani. (2013).
An introduction to statistical learning : with applications in R. New York: Springer.

Thomas G. Dietterich. (2006).Approximate Statistical Tests for Comparing Supervised
Classification Learning Algorithms. The MIT Press Journal, 10 (7).

Mountrakis, G., Im, J., Ogole, C. (2011). Support vector machines in remote sensing:
A review. ISPRS Journal of Photogrammetry and Remote Sensing, 66, 247-259.

Belgiu, M., Dragut., L. (2016). Random Forest in Remote Sensing: A Review of Applications
and Future Directions. ISPRS Journal of Photogrammetry and Remote Sensing, 114, 24-31.

Maxwell, A.E., Warner, T.A., Fang, F. (2018). Implementation of machine-learning
classification in remote sensing: an applied review. International Journal of Remote
Sensing, 29(9), 2784-2817.

Pradhan, R., Ghose, M.K., Jeyaram, A. (2010). Land Cover Classification of Remotely
Sensed Satellite Data using Bayesian and Hybrid classifier. International Journal
of Computer Applications, 7(11).

Holloway, J., Mengersen, K. (2018). Statistical Machine Learning Methods and Remote
Sensing for Sustainable Development Goals: A Review. Remote Sensing, 10(9), 1365.
}

\section{Note}{
 At the moment, only one calibration approach can be used.
}

\examples{
library(ForesToolboxRS)
library(raster)
library(caret)

# Load the dataset
data(FTdata)


# Support Vector Machine and Random Forest Classifiers
# Calibrating using "Set-Approach"
knn <- calmla(img = image, endm = endm, model = c("svm", "randomForest"), training_split = 80,
           approach = "Set-Approach", iter = 10)

}
\author{
Yonatan Tarazona
}
