#' Invert the log joint pseudoposterior Hessian matrix.
#'
#' The function \code{invert_hessian_emvs} computes the inverse Hessian of the
#' log of the joint pseudoposterior distribution.
#'
#' The function \code{invert_hessian_emvs} is used by
#' \code{multidimensional_update_emvs} for the multidimensional Newton Raphson
#' step.
#'
#' @param x An \code{n} by \code{p} matrix containing binary coded
#'   responses (i.e., coded \code{0,1}) for \code{n} independent observations
#'   on \code{p} variables in the network or graph.
#'
#' @param sigma A \code{p} by \code{p} numeric matrix with pairwise association
#'   estimates on the off-diagonal elements and threshold estimates on the
#'   diagonal elements.
#'
#' @param gamma A \code{p} by \code{p} numeric matrix with expected inclusion
#'   probabilities for the pairwise association parameters (i.e., edges).
#'
#' @param index An indexing matrix used to convert the matrix \code{sigma} to a
#' vector \code{eta} during optimization. See \code{\link{indexing}}.
#'
#' @param prior_var_intercepts The variance of the prior distribution
#'   stipulated on the Ising model's threshold parameters. Currently a normal
#'   distribution is used for all threshold parameters, with a mean equal to
#'   zero and a variance equal to \code{prior_var_intercepts}. Defaults to
#'   \code{1}.
#'
#' @param spike_var,slab_var The \code{p} by \code{p} matrices of variances
#'   that are used in the specification of the spike and slab prior
#'   distributions that are stipulated on the association parameters.
#'
#' @return A \code{p(p-1)/2} by \code{p(p-1)/2} numeric matrix.
#'
#' @export
invert_hessian_emvs <- function(x, sigma, gamma, index,
                                prior_var_intercepts = 1,
                                spike_var, slab_var) {
  # the Hessian matrix is build up as
  #           (main_hessian    , cross_hessian)
  #           (cross_hessian^T , int_hessian)
  # main_hessian contains partial derivatives w.r.t the main effects.
  # int_hessian contains partial derivatives w.r.t the interaction effects.

  p <- ncol(x)
  n <- nrow(x)

  if(!is.matrix(spike_var)) {
    spike_var <- matrix(spike_var[1], nrow = p, ncol = p)
  }
  if(!is.matrix(slab_var)) {
    slab_var <- matrix(slab_var[1], nrow = p, ncol = p)
  }

  # precomputation ------------------------------------------------------------
  pq <- matrix(0, nrow = n, ncol = p)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    pq[, s] <- nexpit(phi)
  }

  # compute int_hessian -------------------------------------------------------
  int_hessian <- matrix(0, nrow = p * (p - 1) / 2, ncol = p * (p-1) / 2)
  for(row in 1:(p * (p - 1) / 2)) {
    s <- index[row, 2]
    t <- index[row, 3]

    post.expectation <- gamma[s, t] / slab_var[s, t]
    post.expectation <- post.expectation + (1 - gamma[s, t]) / spike_var[s, t]

    # partial derivative w.r.t. sigma[s, t]; s < t ----------------------------
    int_hessian[row, row] <- -x[, t] %*% pq[, s] - x[, s] %*% pq[, t]
    int_hessian[row, row] <- int_hessian[row, row] - post.expectation

    # partial derivative w.r.t. sigma[s, t] & sigma [s, r]; s < t -------------
    I <- which(index[, 2] == s & index[, 3] != t)
    for(i in I) {
      r <- index[i, 3]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, t] * x[, r]) %*% pq[, s]
    }

    # partial derivative w.r.t. sigma[s, t] & sigma [r, s]; s < t -------------
    I <- which(index[, 3] == s)
    for(i in I) {
      r <- index[i, 2]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, t] * x[, r]) %*% pq[, s]
    }

    # partial derivative w.r.t. sigma[s, t] & sigma [k, t]; s < t -------------
    I <- which(index[, 3] == t & index[, 2] != s)
    for(i in I) {
      k <- index[i, 2]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, s] * x[, k]) %*% pq[, t]
    }

    # partial derivative w.r.t. sigma[s, t] & sigma [t, k]; s < t -------------
    I <- which(index[, 2] == t)
    for(i in I) {
      k <- index[i, 3]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, s] * x[, k]) %*% pq[, t]
    }
  }

  cross_hessian <- matrix(0, nrow = p, ncol = p * (p - 1) / 2)
  for(s in 1:p) {
    I <- which(index[, 2] == s)
    for(i in I) {
      column <- index[i, 1]
      t <- index[i, 3]
      cross_hessian[s, column] <- - x[, t] %*% pq[, s]
    }

    I <- which(index[, 3] == s)
    for(i in I) {
      column <- index[i, 1]
      t <- index[i, 2]
      cross_hessian[s, column] <- - x[, t] %*% pq[, s]
    }
  }

  # compute inv_hessian -------------------------------------------------------
  inv_hessian <- matrix(data = 0,
                        nrow = p * (p + 1) / 2,
                        ncol = p * (p + 1) / 2)

  # indices for main effects and interactions
  index_main <- 1:p
  index_int <- (1:(p * (p + 1) / 2))[-index_main]

  # inverse of main for Schur complement --------------------------------------
  inv_main <- matrix(0, nrow = p, ncol = p)
  for(s in 1:p)
    inv_main[s, s] <- 1 / (-sum(pq[, s]) - 1 / prior_var_intercepts)

  # use Schur complement for inverting ----------------------------------------
  inv_hessian[index_int, index_int] <- int_hessian -
    t(cross_hessian) %*% inv_main %*% cross_hessian
  inv_hessian[index_int, index_int] <- solve(inv_hessian[index_int, index_int])
  inv_hessian[index_main, index_int] <-
    -inv_main %*% cross_hessian %*% inv_hessian[index_int, index_int]
  inv_hessian[index_int, index_main] <- t(inv_hessian[index_main, index_int])
  inv_hessian[index_main, index_main] <- inv_main -
    inv_hessian[index_main, index_int] %*% t(cross_hessian) %*% inv_main

  return(inv_hessian)
}

#' Multidimensional Newton Raphson joint pseudoposterior updates for EMVS.
#'
#' The function \code{multidimensional_update_emvs} updates all pseudoposterior
#' values for the Ising model parameters at once.
#'
#' The function \code{multidimensional_update_emvs} is used by
#' \code{edge_screening_emvs} to optimize the Ising model parameters.
#'
#' @inheritParams invert_hessian_emvs
#'
#' @param suff_stat A \code{p} by \code{p} matrix of sufficient statistics.
#'
#' @return A \code{p} by \code{p} numeric matrix with updates pairwise
#'   association estimates on the off-diagonal elements and updated threshold
#'   estimates on the diagonal elements.
#'
#' @export
multidimensional_update_emvs <- function(x, sigma, gamma, index, suff_stat,
                                         prior_var_intercepts = 1,
                                         spike_var, slab_var) {
  p <- ncol(x)
  n <- nrow(x)

  # compute vector of first-order derivatives (main effects) ------------------
  prob <- matrix(0, nrow = n, ncol = p)
  derivatives <- vector(length = p * (p + 1) / 2)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    prob[, s] <- expit(phi)
    derivatives[s] <- suff_stat[s, s] - sum(prob[, s])
    derivatives[s] <- derivatives[s] - sigma[s, s] / prior_var_intercepts
  }

  # compute vector of first-order derivatives (interaction effects) -----------
  for(s in 1:(p - 1)) {
    for(t in (s + 1): p) {
      row <- p + index[index[, 2] == s & index[, 3] == t, 1]
      post_expectation <- gamma[s, t] / slab_var[s, t]
      post_expectation <- post_expectation + 1 / spike_var[s, t]
      post_expectation <- post_expectation - gamma[s, t] / spike_var[s, t]

      derivatives[row] <- 2 * suff_stat[s, t] - x[, t] %*% prob[, s]
      derivatives[row] <- derivatives[row] - x[, s] %*% prob[, t]
      derivatives[row] <- derivatives[row] - sigma[s, t] * post_expectation
    }
  }

  #Compute inverse Hessian ----------------------------------------------------
  inv_hessian <- invert_hessian_emvs(sigma = sigma, gamma = gamma,
                                     index = index, x = x,
                                     prior_var_intercepts =
                                       prior_var_intercepts,
                                     spike_var = spike_var,
                                     slab_var = slab_var)

  # convert to eta values -----------------------------------------------------
  eta <- vector(length = p * (p + 1) / 2)
  eta[1:p] <- diag(sigma)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- index[index[,2] == s & index[, 3] == t, 1] + p
      eta[row] <- sigma[s, t]
    }
  }

  # Newton - Raphson update ---------------------------------------------------
  eta <- eta - inv_hessian %*% derivatives

  # revert to sigma values ----------------------------------------------------
  diag(sigma) <- eta[1:p]
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- p + index[index[,2] == s & index[, 3] == t, 1]
      sigma[s, t] <- eta[row]
      sigma[t, s] <- eta[row]
    }
  }
  return(sigma)
}

#' Compute the proportional log joint pseudoposterior density value.
#'
#' The function \code{\link{compute_log_pseudoposterior_emvs}} computes the
#' proportional log joint pseudoposterior density.
#'
#' Used by \code{\link{edge_screening_emvs}} to check for convergence of the
#' Newton Raphson procedure.
#'
#' @inheritParams invert_hessian_emvs
#'
#' @param theta The prior inclusion probability. The value \code{theta = 0.5},
#'   in combination with \code{hierarchical = FALSE} stipulates a uniform prior
#'   on the space of network structures.
#'
#' @param hierarchical Logical. If TRUE, a beta prior distribution is
#'  imposed on the prior inclusion probability \code{theta} with
#'  hyperparameters \code{alpha} and \code{beta}. A uniform prior on the
#'  inclusion probability, a beta with \code{alpha = beta = 1}, stipulates a
#'  uniform prior on network structure complexity.
#'
#' @param alpha,beta The hyperparameters of the beta prior distribution
#'   stipulated on the prior inclusion probability \code{theta} if
#'   \code{hierarchical = TRUE}. Default to \code{1}.
#'
#' @return A single numeric value.
#'
#' @importFrom stats dbeta dnorm
#'
#' @export
compute_log_pseudoposterior_emvs <- function(x, sigma, theta = 0.5,
                                             prior_var_intercepts = 1,
                                             spike_var, slab_var,
                                             hierarchical = FALSE,
                                             alpha = 1, beta = 1) {
  p <- ncol(x)

  log_pseudoposterior <- 0
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    log_pseudoposterior <- log_pseudoposterior + x[, s] %*% phi
    log_pseudoposterior <- log_pseudoposterior - sum(log(1 + exp(phi)))
    log_pseudoposterior <- log_pseudoposterior +
      dnorm(x = sigma[s, s], sd = sqrt(prior_var_intercepts), log = TRUE)
  }

  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      tmp1 <- theta * dnorm(x = sigma[s, t], sd = sqrt(slab_var[s, t]))
      tmp2 <- (1 - theta) * dnorm(x = sigma[s, t], sd = sqrt(spike_var[s, t]))
      log_pseudoposterior <- log_pseudoposterior + log(tmp1 + tmp2)
    }
  }
  if(hierarchical == TRUE)
    log_pseudoposterior <- log_pseudoposterior + dbeta(x = theta,
                                                       shape1 = alpha,
                                                       shape2 = beta,
                                                       log = TRUE)

  return(log_pseudoposterior)
}

#' EM edge inclusion procedure.
#'
#' The function \code{edge_screening_emvs} estimates the
#' (local) posterior probability that an edge should be included using the EM
#' algorithm.
#'
#' The function \code{\link{screen_edges}} is a wrapper for
#' \code{edge_screening_emvs}.
#'
#' @inheritParams compute_log_pseudoposterior_emvs
#'
#' @param sigma A \code{p} by \code{p} numeric matrix with pairwise association
#'   estimates on the off-diagonal elements and threshold estimates on the
#'   diagonal elements. Optional. Can be used to specify starting values.
#'
#' @param iteration_max The maximum number of Newton-Raphson iterations used in
#'   optimization. Defaults to \code{1e2}. A warning is issued if procedure has
#'   not converged in \code{iteration_max} iterations.
#'
#' @return A list containing the \code{p} by \code{p} matrices \code{sigma} and
#'   \code{gamma}. The matrix \code{sigma} is a numeric matrix with pairwise
#'   association estimates on the off-diagonal elements and threshold
#'   estimates on the diagonal elements. The matrix \code{gamma} contains the
#'   expected values of edge inclusion variables (i.e., the local posterior
#'   probability of edge inclusion. If \code{hierarchical = TRUE}, the modal
#'   estimate of the prior inclusion probability \code{theta} is also provided.
#'
#' @export
edge_screening_emvs <- function(x, sigma, theta = 0.5, prior_var_intercepts = 1,
                                spike_var, slab_var, hierarchical = FALSE,
                                alpha = 1, beta = 1,
                                iteration_max = 1e2) {
  p <- ncol(x)
  n <- nrow(x)

  # parameter matrix and sufficient statistics --------------------------------
  if(!hasArg("sigma")) {
    sigma <- matrix(0, nrow = p, ncol = p)
  }
  gamma <- matrix(1, nrow = p, ncol = p)
  diag(gamma) <- 0
  suff_stat <- t(x) %*% x

  # index allows us to convert sigma (matrix) to eta (vector) and revert back -
  index <- indexing (p)

  if(hierarchical == TRUE) {
    log_pseudoposterior_emvs <- compute_log_pseudoposterior_emvs(sigma = sigma,
                                               theta = theta, x = x,
                                               prior_var_intercepts =
                                                 prior_var_intercepts,
                                               spike_var = spike_var,
                                               slab_var = slab_var,
                                               hierarchical = hierarchical,
                                               alpha = alpha, beta = beta)
  } else {
    log_pseudoposterior_emvs <- compute_log_pseudoposterior_emvs(sigma = sigma,
                                              theta = theta, x = x,
                                              prior_var_intercepts =
                                                prior_var_intercepts,
                                              spike_var = spike_var,
                                              slab_var = slab_var,
                                              hierarchical = FALSE)
  }

  for(iteration in 1:iteration_max) {
    # E-step - update selection variables -------------------------------------
    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        tmp1 <- theta * dnorm(x = sigma[s, t], mean = 0,
                              sd = sqrt(slab_var[s, t]))
        tmp2 <- (1 - theta) * dnorm(x = sigma[s, t],
                                    mean = 0,
                                    sd = sqrt(spike_var[s, t]))
        gamma[s, t] <- gamma[t, s] <- tmp1 / (tmp1 + tmp2)
      }
    }

    # M-step - update prior inclusion probability -----------------------------
    if(hierarchical == TRUE) {
      tmp <- sum(gamma[lower.tri(gamma)])
      theta <- (tmp + alpha - 1) / (alpha + beta - 2 + p * (p - 1) / 2)
    }

    # M-step - update model parameters ----------------------------------------
    sigma <- multidimensional_update_emvs(sigma = sigma, gamma = gamma,
                                          index = index, x = x,
                                          suff_stat = suff_stat,
                                          prior_var_intercepts =
                                            prior_var_intercepts,
                                          spike_var = spike_var,
                                          slab_var = slab_var)

    # recompute log-pseudoposterior -------------------------------------------
    old <- log_pseudoposterior_emvs
    if(hierarchical == TRUE) {
      new <- compute_log_pseudoposterior_emvs(sigma = sigma,
                                              theta = theta, x = x,
                                              prior_var_intercepts =
                                                prior_var_intercepts,
                                              spike_var = spike_var,
                                              slab_var = slab_var,
                                              hierarchical = hierarchical,
                                              alpha = alpha, beta = beta)
    } else {
      new <- compute_log_pseudoposterior_emvs(sigma = sigma, theta = theta,
                                              x = x,
                                              prior_var_intercepts =
                                                prior_var_intercepts,
                                              spike_var = spike_var,
                                              slab_var = slab_var,
                                              hierarchical = hierarchical)
    }
    if(abs(new - old) < sqrt(.Machine$double.eps))
      break

    log_pseudoposterior_emvs <- new

    if(iteration == iteration_max)
      warning(paste("The optimization procedure did not convergence in", iteration_max, "iterations.",
                    sep = " "), call. = FALSE)
  }

  if(hierarchical == TRUE) {
    return(list(sigma = sigma, gamma = gamma, theta = theta))
  } else {
    return(list(sigma = sigma, gamma = gamma))
  }
}

#' Bayesian edge screening for the Ising model.
#'
#' The function \code{screen_edges} screens promising edges for the Ising
#' model, using the joint pseudolikelihood and a spike and slab prior
#' distribution stipulated on the Ising model's association parameters.
#'
#' The function \code{\link{screen_edges}} is a wrapper for
#' \code{edge_screening_emvs}.
#'
#' @param x An \code{n} by \code{p} matrix containing binary coded
#'   responses (i.e., coded \code{0,1}) for \code{n} independent observations
#'   on \code{p} variables in the network or graph.
#'
#' @param theta The prior inclusion probability. The value \code{theta = 0.5},
#'   in combination with \code{hierarchical = FALSE} stipulates a uniform prior
#'   on the space of network structures.
#'
#' @param precision A number between zero and one. The prior precision that is
#' desired for edge selection. Equal to one minus the desired type-1 error.
#' Defaults to \code{.975}.
#'
#' @param prior_var_intercepts The variance of the prior distribution
#'   stipulated on the Ising model's threshold parameters. Currently a normal
#'   distribution is used for all threshold parameters, with a mean equal to
#'   zero and a variance equal to \code{prior_var_intercepts}. Defaults to
#'   \code{1}.
#'
#' @param hierarchical Logical. If TRUE, a beta prior distribution is
#'  imposed on the prior inclusion probability \code{theta} with
#'  hyperparameters \code{alpha} and \code{beta}. A uniform prior on the
#'  inclusion probability, a beta with \code{alpha = beta = 1}, stipulates a
#'  uniform prior on network structure complexity.
#'
#' @param alpha,beta The hyperparameters of the beta prior distribution
#'   stipulated on the prior inclusion probability \code{theta} if
#'   \code{hierarchical = TRUE}. Default to \code{1}.
#'
#' @param iteration_max The maximum number of Newton-Raphson iterations used in
#'   optimization. Defaults to \code{1e2}. A warning is issued if procedure has
#'   not converged in \code{iteration_max} iterations.
#'
#' @return A list containing the \code{p} by \code{p} matrices \code{sigma} and
#'   \code{gamma}. The matrix \code{sigma} is a numeric matrix with pairwise
#'   association estimates on the off-diagonal elements and threshold
#'   estimates on the diagonal elements. The matrix \code{gamma} contains the
#'   expected values of edge inclusion variables (i.e., the local posterior
#'   probability of edge inclusion. If \code{hierarchical = TRUE}, the modal
#'   estimate of the prior inclusion probability \code{theta} is also provided.
#'
#' @examples
#' library("IsingSampler")
#' ### Simulate dataset ###
#' # Input:
#' p <- 6 # Number of nodes
#' n <- 1000 # Number of samples
#' # Ising parameters:
#' Graph <- matrix(data = 0, nrow = p, ncol = p)
#' Graph[lower.tri(Graph)] <- rbinom(n = p * (p - 1) / 2,
#'                                   size = 1,
#'                                   prob = 0.2)
#' Graph <- Graph * runif(n = p ^ 2, min = 0.5, max = 2)
#' Graph <- Graph + t(Graph)
#' Thresholds <- -rowSums(Graph) / 2
#' # Simulate:
#' Data <- IsingSampler(n = n, graph = Graph, thresholds = Thresholds)
#' ### Fit using fit_pseudoposterior ###
#' screened <- screen_edges(x = Data)
#' graph <- screened$sigma
#' graph[screened$gamma < 0.5] <- 0
#' # Plot results:
#' library("qgraph")
#' layout(t(1:2))
#' qgraph(graph,fade = FALSE)
#' title("Network of screened edges")
#' qgraph(Graph,fade = FALSE)
#' title("Original network")
#'
#' @export
screen_edges <- function(x = x, theta = 0.5, precision = .975, prior_var_intercepts = 1,
                         hierarchical = FALSE, alpha = 1,
                         beta = 1, iteration_max = 1e2) {

  spike_and_slab <- set_spike_and_slab (x, precision)
  spike_var <- spike_and_slab$spike_var
  slab_var <- spike_and_slab$slab_var
  sigma_ml <- spike_and_slab$sigma_ml

  # edge screening using emvs -------------------------------------------------
  screened_edges <- try(edge_screening_emvs(x = x, sigma = sigma_ml, theta = theta,
                             prior_var_intercepts = prior_var_intercepts,
                             spike_var = spike_var, slab_var = slab_var,
                             hierarchical = hierarchical, alpha = alpha,
                             beta = beta, iteration_max = iteration_max),
                        silent = TRUE)
  if(class(screened_edges) != "try-error") {
    return(screened_edges)
  } else {
    stop("Edge screening (EMVS) failed.")
  }
}
