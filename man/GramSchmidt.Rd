% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_GramSchmidt.R
\name{GramSchmidt}
\alias{GramSchmidt}
\title{Gram-Schmidt orthogonalization of vectors}
\usage{
GramSchmidt(x, otherVecs, weights = rep(1, length(x)))
}
\arguments{
\item{x}{The vector that is to be orthogonalized}

\item{otherVecs}{a matrix; x is orthogonalized with respect to its rows}

\item{weights}{The weights used in the orthogonalization}
}
\value{
The orthogonalized vector
}
\description{
Gram-Schmidt orthogonalization of vectors
}
