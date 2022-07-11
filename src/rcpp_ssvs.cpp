// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include "rcpp_pgdraw.h"
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;

// Sample latent Polya-Gamma variates & main effects
void sample_omega_mu(IntegerMatrix x, NumericMatrix sigma, NumericMatrix omega,
                     double prior_var_intercepts, NumericVector sufC){
    // constants
    int p = x.ncol();
    int n = x.nrow();
    //NumericMatrix omega(n, p);
    int i, j, v;

    double mu, sum_v, omega_i, chi, post_var, post_mean;

    // loop over nodes
    for (i = 0; i < p; i++){
        mu = sigma(i, i);
        sum_v = 0;  // sum over \sigma_{ij}x_{jv}\omega_{iv}
        omega_i = 0;  // \omega_{i+}
        // loop over participants
        for (v = 0; v < n; v++){
            chi = 0;
            for (j = 0; j < p; j++){
                if ((j != i) & (x(v, j) == 1)){
                    chi += sigma(i, j);
                }
            }
            // sample omegas from PG
            omega(v, i) = samplepg(mu + chi);
            sum_v += omega(v, i) * chi;
            omega_i += omega(v, i);
        }
        // sample main effects
        post_var = 1.0 + omega_i * prior_var_intercepts;
        post_var = prior_var_intercepts / post_var;

        post_mean = sufC[i] - 0.5 * n - sum_v;
        post_mean = post_var * post_mean;

        sigma(i, i) = R::rnorm(post_mean, std::sqrt(post_var));
    }
}

// Sample inclusion variables & interaction effects
void sample_sigma_gamma(IntegerMatrix x, NumericMatrix sigma, NumericMatrix omega, IntegerMatrix gamma,
                        IntegerMatrix include, NumericMatrix slab_var, NumericMatrix spike_var,
                        double prior_var_intercepts, double theta, double r, double df, NumericVector sufC,
                        CharacterVector components){
    // constants
    int p = x.ncol();
    int n = x.nrow();

    // placeholder initiates
    double tmp1, tmp2, prob, post_mean, post_var, x_sum, sum_i, sum_j, q_i, q_j, lambda, rr, phi;
    int i, j, v, q;

    double theta_odds = theta / (1 - theta);

    for (i = 0; i < p - 1; i++){
        for (j = i + 1; j < p; j++){
            if (include(i, j) == 1){
                // sample inclusion variable
                tmp2 = (1 / slab_var(i, j) - 1 / spike_var(i, j)) * sigma(i, j) * sigma(i, j) * 0.5;
                prob = theta_odds / (theta_odds + std::sqrt(slab_var(i, j) / spike_var(i, j)) * std::exp(tmp2));
                //tmp1 = theta * R::dnorm(sigma(i, j), 0, std::sqrt(slab_var(i, j)), false);
                //tmp2 = (1 - theta) * R::dnorm(sigma(i, j), 0, std::sqrt(spike_var(i, j)), false);
                //prob = tmp1 / (tmp1 + tmp2);
                gamma(i, j) = R::rbinom(1, prob);
                gamma(j, i) = gamma(i, j);

                // sample interaction effect
                tmp1 = slab_var(i, j) * gamma(i, j) + spike_var(i, j) * (1 - gamma(i, j));
                tmp2 = 0, x_sum = 0, sum_i = 0, sum_j = 0;
                for (v = 0; v < n; v++){
                    // complicated sums
                    x_sum += x(v, i) * x(v, j);
                    q_i = 0, q_j = 0;
                    for (q = 0; q < p; q++){
                        if ((q != i) & (q != j) & (x(v, q) == 1)){
                            q_i += sigma(i, q);
                            q_j += sigma(j, q);
                        }
                    }
                    // tmp2 += omega(v, i) * x(v, j) + omega(v, j) * x(v, i);
                    if (x(v, j) == 1){
                        sum_i += omega(v, i) * (sigma(i, i) + q_i);
                        tmp2 += omega(v, i);
                    }
                    if (x(v, i) == 1){
                        sum_j += omega(v, j) * (sigma(j, j) + q_j);
                        tmp2 += omega(v, j);
                    }
                }
                post_var = tmp1 / (1 + tmp1 * tmp2);

                post_mean = 2 * x_sum - 0.5 * sufC[i] - 0.5 * sufC[j] - sum_i - sum_j;
                post_mean = post_var * post_mean;

                sigma(i, j) = R::rnorm(post_mean, std::sqrt(post_var));
                sigma(j, i) = sigma(i, j);

                // sample hierarchical variance
                if(components[0] != "normal"){
                    if (components[0] == "t"){
                        rr = gamma(i, j) + (1 - gamma(i, j)) * r;
                        phi = 1.0 / R::rgamma(df + 0.5,
                                              slab_var(j, i) * (df - 1) + sigma(i, j) * sigma(i, j) / (2.0 * rr));
                    }
                    if (components[0] == "laplace"){
                        lambda = std::sqrt(slab_var(j, i) * 0.5);
                        rr = gamma(i, j) + (1 - gamma(i, j)) * r;
                        phi = 1.0 / randinvg(std::sqrt(rr) / (lambda * std::abs(sigma(i, j))),
                                             1.0 / (lambda * lambda));
                    }
                    //update subdiagonal elements of spike and slab variance matrices --
                    spike_var(i, j) = r * phi;
                    slab_var(i, j) = phi;
                }
            }
        }
    }
}

// [[Rcpp::export]]
List gibbs_hierarchical_samples(int number_burnin_iterations, int number_iterations,
          IntegerMatrix x, NumericMatrix sigma, IntegerMatrix gamma, //NumericMatrix omega,
          NumericMatrix spike_var, NumericMatrix slab_var, IntegerMatrix include, NumericVector sufC,
          int alpha, int beta, int number_edges,
          double theta, double prior_var_intercepts, double df,
          CharacterVector components,
          bool display_progress){

    // Setup
    int idx_sigma, idx_gamma, iteration, row, col;

    int p = x.ncol();
    int n = x.nrow();
    NumericMatrix omega(n, p);

    double r = spike_var(0, 1) / slab_var(0, 1);

    // output
    NumericMatrix sigma_samples(number_iterations, p + R::choose(p, 2));
    IntegerMatrix structures(number_iterations, R::choose(p, 2));
    NumericVector theta_out(number_iterations);

    // Progress bar
    Progress bar(number_burnin_iterations + number_iterations, display_progress);

    for (iteration = -number_burnin_iterations; iteration < number_iterations; iteration++){
        bar.increment();

        //sample PG & main effects
        sample_omega_mu(x, sigma, omega, prior_var_intercepts, sufC);

        // sample inclusion variables & interaction effects
        sample_sigma_gamma(x, sigma, omega, gamma,
                           include, slab_var, spike_var,
                           prior_var_intercepts, theta, r, df, sufC,
                           components);

        // sample prior inclusion probability
        theta = R::rbeta(alpha + 0.5 * sum(gamma),
                         beta + number_edges - 0.5 * sum(gamma));

        // collect output
        if (iteration >= 0){
            idx_sigma = 0, idx_gamma = 0;
            theta_out[iteration] = theta;
            for (col = 0; col < p; col++){
                for (row = col; row < p; row++){
                    sigma_samples(iteration, idx_sigma) = sigma(row, col);
                    idx_sigma++;
                    if (row != col){
                        structures(iteration, idx_gamma) = gamma(row, col);
                        idx_gamma++;
                    }
                }
            }
        }
    }
    return List::create(structures, sigma_samples, theta_out);
}

// [[Rcpp::export]]
List gibbs_hierarchical_eap(int number_burnin_iterations, int number_iterations,
           IntegerMatrix x, NumericMatrix sigma, IntegerMatrix gamma, //NumericMatrix omega,
           NumericMatrix spike_var, NumericMatrix slab_var, IntegerMatrix include, NumericVector sufC,
           int alpha, int beta, int number_edges,
           double theta, double prior_var_intercepts, double df,
           CharacterVector components,
           bool display_progress){

    // Setup
    int p = x.ncol();
    int n = x.nrow();
    NumericMatrix omega(n, p);
    int idx_gamma, iteration, row, col;
    double r = spike_var(0, 1) / slab_var(0, 1);

    // output
    NumericMatrix sigma_out(p, p);
    IntegerMatrix structures(number_iterations, R::choose(p, 2));
    double theta_out = 0;

    Progress bar(number_burnin_iterations + number_iterations, display_progress);

    for (iteration = -number_burnin_iterations; iteration < number_iterations; iteration++){
        bar.increment();

        //sample PG & main effects
        sample_omega_mu(x, sigma, omega, prior_var_intercepts, sufC);

        // sample inclusion variables & interaction effects
        sample_sigma_gamma(x, sigma, omega, gamma,
                           include, slab_var, spike_var,
                           prior_var_intercepts, theta, r, df, sufC,
                           components);

        // sample prior inclusion probability
        theta = R::rbeta(alpha + 0.5 * sum(gamma),
                         beta + number_edges - 0.5 * sum(gamma));

        // collect output
        if (iteration >= 0){
            sigma_out += sigma;
            idx_gamma = 0;
            theta_out += theta;
            for (col = 0; col < p; col++){
                for (row = col + 1; row < p; row++){
                    structures(iteration, idx_gamma) = gamma(row, col);
                    idx_gamma++;
                }
            }
        }
    }
    sigma_out = sigma_out / number_iterations;
    theta_out /= number_iterations;

    return List::create(structures, sigma_out, theta_out);
}

// [[Rcpp::export]]
List gibbs_uniform_samples(int number_burnin_iterations, int number_iterations,
           IntegerMatrix x, NumericMatrix sigma, IntegerMatrix gamma, //NumericMatrix omega,
           NumericMatrix spike_var, NumericMatrix slab_var, IntegerMatrix include, NumericVector sufC,
           int alpha, int beta, int number_edges,
           double theta, double prior_var_intercepts, double df,
           CharacterVector components,
           bool display_progress){

    // Setup
    int p = x.ncol();
    int n = x.nrow();
    NumericMatrix omega(n, p);
    int idx_sigma, idx_gamma, iteration, row, col;

    double r = spike_var(0, 1) / slab_var(0, 1);

    // output
    NumericMatrix sigma_samples(number_iterations, p + R::choose(p, 2));
    IntegerMatrix structures(number_iterations, R::choose(p, 2));

    Progress bar(number_burnin_iterations + number_iterations, display_progress);

    for (iteration = -number_burnin_iterations; iteration < number_iterations; iteration++){
        bar.increment();

        //sample PG & main effects
        sample_omega_mu(x, sigma, omega, prior_var_intercepts, sufC);

        // sample inclusion variables & interaction effects
        sample_sigma_gamma(x, sigma, omega, gamma,
                           include, slab_var, spike_var,
                           prior_var_intercepts, theta, r, df, sufC,
                           components);

        // collect output
        if (iteration >= 0){
            idx_sigma = 0, idx_gamma = 0;
            for (col = 0; col < p; col++){
                for (row = col; row < p; row++){
                    sigma_samples(iteration, idx_sigma) = sigma(row, col);
                    idx_sigma++;
                    if (row != col){
                        structures(iteration, idx_gamma) = gamma(row, col);
                        idx_gamma++;
                    }
                }
            }
        }
    }
    return List::create(structures, sigma_samples);
}
// [[Rcpp::export]]
List gibbs_uniform_eap(int number_burnin_iterations, int number_iterations,
           IntegerMatrix x, NumericMatrix sigma, IntegerMatrix gamma, //NumericMatrix omega,
           NumericMatrix spike_var, NumericMatrix slab_var, IntegerMatrix include, NumericVector sufC,
           int alpha, int beta, int number_edges,
           double theta, double prior_var_intercepts, double df,
           CharacterVector components,
           bool display_progress){

    // Setup
    int p = x.ncol();
    int n = x.nrow();
    NumericMatrix omega(n, p);
    int idx_gamma, iteration, row, col;
    double r = spike_var(0, 1) / slab_var(0, 1);

    // output
    NumericMatrix sigma_out(p, p);
    IntegerMatrix structures(number_iterations, R::choose(p, 2));

    Progress bar(number_burnin_iterations + number_iterations, display_progress);

    for (iteration = -number_burnin_iterations; iteration < number_iterations; iteration++){
        bar.increment();

        //sample PG & main effects
        sample_omega_mu(x, sigma, omega, prior_var_intercepts, sufC);

        // sample inclusion variables & interaction effects
        sample_sigma_gamma(x, sigma, omega, gamma,
                           include, slab_var, spike_var,
                           prior_var_intercepts, theta, r, df, sufC,
                           components);

        // collect output
        if (iteration >= 0){
            sigma_out += sigma;
            idx_gamma = 0;
            for (col = 0; col < p; col++){
                for (row = col + 1; row < p; row++){
                    structures(iteration, idx_gamma) = gamma(row, col);
                    idx_gamma++;
                }
            }
        }
    }
    sigma_out = sigma_out / number_iterations;
    return List::create(structures, sigma_out);
}
