% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{unique_row}
\alias{unique_row}
\title{Reduce matrix to unique rows}
\usage{
unique_row(mat)
}
\arguments{
\item{mat}{A matrix}
}
\value{
A list containing \code{mat} and \code{count}, the reduced matrix and a vector of the number of each occurrence.
}
\description{
The function \code{unique_row} is used to reduce the input matrix
to a matrix of unique rows.
}
\details{
It uses a radix sort algorithm to efficiently reorder the rows by size
and then linearly reduces each rows to a single occurrence, and returns
each row's number of occurrences.
Adapted from https://stackoverflow.com/a/29829228
}
