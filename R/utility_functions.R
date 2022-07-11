#' Calculate expit (inverse logit) of \code{phi}
#'
#' @param phi
#'
expit <- function(phi) {
  return(exp(phi) / (1 + exp(phi)))
}

#' @describeIn expit Nexpit of \code{phi} for inverse Hessian
nexpit <- function(phi) {
  return(exp(phi) / (1 + exp(phi)) ^ 2)
}

#' Indexing matrix for converting parameter matrix to vector.
#'
#' The function \code{indexing} is used to convert the \code{p} by \code{p}
#' parameter matrix \code{sigma} to a \code{p * (p + 1) / 2} vector \code{eta}
#' during optimization.
#'
#' The function \code{indexing} is used by
#' \code{\link{compute_standard_deviation}},
#' \code{\link{optimize_pseudoposterior}},
#' and ....
#'
#' @param p The number of rows and columns of the parameter matrix \code{sigma}.
#'
#' @return A \code{p * (p + 1) / 2} by \code{3} matrix of integers. The first
#'   column contains the index of the vector \code{eta}, the second and third
#'   column the corresponding row and column in \code{sigma}.
#'
#' @export
indexing <- function(p) {
  index <- matrix(0, nrow = p * (p - 1) / 2, ncol = 3)
  counter <- 0
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      counter <- counter + 1
      index[counter, 1] <- counter
      index[counter, 2] <- s
      index[counter, 3] <- t
    }
  }
  return(index)
}

#' Matching the spike and slab intersection.
#'
#' The function \code{xi_delta_matching} is used to match the intersection
#' of the spike and slab prior distributions to a predefined number
#' \code{delta}, as a function of a penalty parameter \code{xi}.
#'
#' The root of \code{xi_delta_matching} is searched by
#' \code{\link{set_spike_and_slab}}. The root value \code{xi} for a given
#' \code{delta} and \code{n} is used to set the spike variance.
#'
#' @param xi A positive penalty parameter.
#'
#' @param delta A positive number. Taken to be the inverse of the intended
#'   precision, e.g., \code{delta = qnorm(.975, lower.tail = TRUE)}.
#'
#' @param n A positive integer. Number of observations.
#'
#' @return A numerical value.
xi_delta_matching <- function(xi, delta, n) {
  n * log(n / xi) / (n / xi - 1) - delta ^ 2
}

#' Specification of the spike and slab variance parameters.
#'
#' The function \code{set_spike_and_slab} is used to specify the variance
#' parameters of the spike and slab prior distribution that is stipulated
#' on the association parameters of the Ising model. It is used by
#' \code{\link{edge_screening}} and \code{\link{structure_selection}}.
#'
#' The function \code{set_spike_and_slab} runs
#' \code{\link{optimize_pseudoposterior}} to produce \code{slab_var},
#' a \code{p} by \code{p} matrix of unit information prior variances. It
#' furthermore aims to find the root of \code{\link{xi_delta_matching}} to
#' match the spike and slab intersections to a predefined precision. If
#' succesfull, \code{spike_var} contains the corresponding variances of the
#' spike components.
#'
#' @param x An \code{n} by \code{p} matrix containing binary coded
#'   responses (i.e., coded \code{0,1}) for \code{n} independent observations
#'   on \code{p} variables in the network or graph.
#'
#' @param precision A number between zero and one. The prior precision that is
#' desired for edge selection. Equal to one minus the desired type-1 error.
#'
#' @importFrom stats qnorm uniroot
#'
#' @return A list containing \code{spike_var} and \code{slab_var}, which are
#'   \code{p} by \code{p} matrices of variances used for the spike and slab
#'   prior components, \code{sigma_ml}, a \code{p} by \code{p} matrix of
#'   maximum pseudolikelihood estimates, with association estimates on the off-
#'   diagonal and threshold estimates on the diagonal.
set_spike_and_slab <- function (x, precision) {
  p <- ncol(x)
  n <- nrow(x)

  # determine xi penalty ------------------------------------------------------
  delta <- qnorm(precision, lower.tail = TRUE)
  xi <- uniroot (f = xi_delta_matching,
                 interval = c(.Machine$double.eps,
                              n - sqrt(.Machine$double.eps)),
                 delta = delta,
                 n = n)$root

  # determine spike and slab variance -----------------------------------------
  sigma_ml <- try(optimize_pseudoposterior(x = x,
                                           prior_var = Inf)$sigma,
                  silent = TRUE)
  if(class(sigma_ml)[1] != "try-error") {
    index <- indexing(p)

    inv_hessian <- try(invert_hessian(sigma = sigma_ml, index = index, x = x,
                                      prior_var = Inf), silent = TRUE)
    if(class(inv_hessian)[1] == "try-error") {
      stop("Could not compute unit information matrix.")
    }
  } else {
    stop("Could not estimate maximum pseudolikelihood estimates.")
  }

  slab_var <- matrix(0, nrow = p, ncol = p)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- p + index[index[, 2] == s & index[, 3] == t, 1]
      slab_var[s, t] <- - n * inv_hessian[row, row]
      slab_var[t, s] <- slab_var[s, t]
    }
  }
  spike_var <- slab_var * xi / n

  output <- list(spike_var = spike_var, slab_var = slab_var,
                 sigma_ml = sigma_ml)
  return(output)
}

#' Reduce matrix to unique rows
#'
#' The function \code{unique_row} is used to reduce the input matrix
#' to a matrix of unique rows.
#'
#' It uses a radix sort algorithm to efficiently reorder the rows by size
#' and then linearly reduces each rows to a single occurrence, and returns
#' each row's number of occurrences.
#' Adapted from https://stackoverflow.com/a/29829228
#'
#' @param mat A matrix
#'
#' @return A list containing \code{mat} and \code{count}, the reduced matrix and a vector of the number of each occurrence.
#'
#' @export
unique_row = function(mat){
    n = nrow(mat)
    p = ncol(mat)

    # radix sort indices of rows
    idx = 1:n
    for (i in p:1){
        c = mat[idx, i]
        ix = sort(c, method = "radix", index.return = TRUE)
        idx = idx[ix$ix]
    }
    mat = mat[idx, ]  # reorder all rows

    # count occurrences of each row
    occ = integer(length = n)
    lastrow = integer(length = p); lastrow[1] = 1L
    id = 1

    for (j in 1:n){
        row = mat[j, ]
        if (!all(row == lastrow)){
            id = j
        }
        occ[id] = occ[id] + 1
        lastrow = row
    }

    ind = as.logical(occ)  # row inclusion vector

    return(list(mat = mat[ind, ],
                count = occ[ind]))
}
