

#' Samples from a mixture model
#'
#' A matrix of 2-dimensional Gaussian mixture samples, the samples came from 4 different Gaussian distributions
#'
#' @docType data
#'
#' @usage data(mmData)
#'
#' @format A matrix of Gaussian samples, each row is one sample.
"mmData"

#' Samples from a hierarchical mixture model
#'
#' A list of 2-dimensional Gaussian mixture samples and the corresponding group labels. The samples are generated by 4 different Gaussian distributions. The samples are separated into 30 groups.
#'
#' @docType data
#'
#' @usage data(mmhData)
#'
#' @format A list of two elements:
#' \describe{
#'   \item{x}{matrix, 2-dimensional Gaussian samples, each row is a sample}
#'   \item{groupLabel}{integer, the group label of each Gaussian sample}
#' }
"mmhData"

#' Samples from a hierarchical mixture model with two layers of hierarchies
#'
#' A list of 2-dimensional Gaussian mixture samples and the corresponding group labels. The samples are generated by 4 different Gaussian distributions. Each sample is assigned to a group and a subgroup. There are 2 groups, group 1 has 10 subgroups, group 2 has 20 subgroups.
#'
#' @docType data
#'
#' @usage data(mmhhData)
#'
#' @format A list of two elements:
#' \describe{
#'   \item{x}{matrix, 2-dimensional Gaussian samples, each row is a sample}
#'   \item{groupLabel}{integer, the group label of each Gaussian sample}
#'   \item{subGroupLabel}{integer, the subgroup label of each Gaussian sample}
#' }
"mmhhData"

#' Samples from a simple linear model
#'
#' A list of two elements, a matrix "X" and a numeric vector "x". They came from a linear model: x = X*0.3 + epsilon, where epsilon is Gaussian distributed with mean 0 and variance 25.
#'
#' @docType data
#'
#' @usage data(lrData)
#'
#' @format A list of two elements:
#' \describe{
#'   \item{x}{numeric, linear samples}
#'   \item{X}{matrix, the "locations" of the linear samples}
#' }
"lrData"

#' @title Samples from a hierarchical linear model
#' @description
#' The data was part of the 2002 Educational Longitudinal Study (ELS), a survey of students from a large sample of schools across the United States. This dataset includes a population of schools as well as a population of students within each school.
#' @docType data
#'
#' @usage data(lrData)
#'
#' @format A list of three elements:
#' \describe{
#'   \item{mathScore}{numeric, the mathScore of each student}
#'   \item{socioeconomicStatus}{numeric, the socioeconomic status score of each student}
#'   \item{schoolID}{integer, the school ID of each student}
#' }
#' @references Hoff, Peter D. A first course in Bayesian statistical methods. Vol. 580. New York: Springer, 2009.
"hlrData"


#' farm ads data
#'
#' A subset of farm ads data from https://archive.ics.uci.edu/ml/datasets/Farm+Ads
#'
#' @docType data
#'
#' @usage data(farmadsData)
#'
#' @format A list of two elements:
#' \describe{
#'   \item{word}{character, the words}
#'   \item{document}{integer, document id of each word}
#' }
#' @source \href{https://archive.ics.uci.edu/ml/datasets/Farm+Ads}{Farm-Ads}
"farmadsData"

#' Cancer mortality of 20 cities
#'
#' Cancer mortality data from Ordinal Data Modeling(1999) p24. This dataset is a list of 20 character vectors, representing the individuals in 20 cities' cancer mortality information. There can be only two possible values "death" and "no death" in each character vector. For example "death" occurs 3 times the 10th character vector, while "no death" occurs 582 times, this mean there are 3 death among the 585 cancer patients in city 10.
#'
#' @docType data
#'
#' @usage data(cancerData)
#'
#' @format A list of 20 character vectors.
#' @references Johnson, Valen E., and James H. Albert. Ordinal data modeling. Springer Science & Business Media, 2006.
"cancerData"



