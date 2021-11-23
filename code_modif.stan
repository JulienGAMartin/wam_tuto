// generated with brms 2.15.0
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_bwt;  // number of observations
  vector[N_bwt] Y_bwt;  // response variable
  int<lower=1> N_tarsus;  // number of observations
  vector[N_tarsus] Y_tarsus;  // response variable
  int<lower=1> nresp;  // number of responses
  int nrescor;  // number of residual correlations
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1_bwt[N_bwt];  // grouping indicator per observation
  int<lower=1> J_1_tarsus[N_tarsus];  // grouping indicator per observation
  matrix[N_1, N_1] Lcov_1;  // cholesky factor of known covariance matrix
  // group-level predictor values
  vector[N_bwt] Z_1_bwt_1;
  vector[N_tarsus] Z_1_tarsus_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  vector[nresp] Y[N];  // response array
  for (n in 1:N) {
    Y[n] = transpose([Y_bwt[n], Y_tarsus[n]]);
  }
}
parameters {
  real Intercept_bwt;  // temporary intercept for centered predictors
  real<lower=0> sigma_bwt;  // residual SD
  real Intercept_tarsus;  // temporary intercept for centered predictors
  real<lower=0> sigma_tarsus;  // residual SD
  cholesky_factor_corr[nresp] Lrescor;  // parameters for multivariate linear models
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[N_1, M_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_bwt_1;
  vector[N_1] r_1_tarsus_2;
  // compute actual group-level effects
  r_1 = (Lcov_1 * z_1) * diag_pre_multiply(sd_1, L_1)';
  r_1_bwt_1 = r_1[, 1];
  r_1_tarsus_2 = r_1[, 2];
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N_bwt] mu_bwt = Intercept_bwt + rep_vector(0.0, N_bwt);
    // initialize linear predictor term
    vector[N_tarsus] mu_tarsus = Intercept_tarsus + rep_vector(0.0, N_tarsus);
    // multivariate predictor array
    vector[nresp] Mu[N];
    vector[nresp] sigma = transpose([sigma_bwt, sigma_tarsus]);
    // cholesky factor of residual covariance matrix
    matrix[nresp, nresp] LSigma = diag_pre_multiply(sigma, Lrescor);
    for (n in 1:N_bwt) {
      // add more terms to the linear predictor
      mu_bwt[n] += r_1_bwt_1[J_1_bwt[n]] * Z_1_bwt_1[n];
    }
    for (n in 1:N_tarsus) {
      // add more terms to the linear predictor
      mu_tarsus[n] += r_1_tarsus_2[J_1_tarsus[n]] * Z_1_tarsus_2[n];
    }
    // combine univariate parameters
    for (n in 1:N) {
      Mu[n] = transpose([mu_bwt[n], mu_tarsus[n]]);
    }
    target += multi_normal_cholesky_lpdf(Y | Mu, LSigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept_bwt | 3, 7.3, 2.8);
  target += student_t_lpdf(sigma_bwt | 3, 0, 2.8)
    - 1 * student_t_lccdf(0 | 3, 0, 2.8);
  target += student_t_lpdf(Intercept_tarsus | 3, 20.4, 5.2);
  target += student_t_lpdf(sigma_tarsus | 3, 0, 5.2)
    - 1 * student_t_lccdf(0 | 3, 0, 5.2);
  target += lkj_corr_cholesky_lpdf(Lrescor | 1);
  target += student_t_lpdf(sd_1[1] | 3, 0, 2.8)
    - 1 * student_t_lccdf(0 | 3, 0, 2.8);
  target += student_t_lpdf(sd_1[2] | 3, 0, 5.2)
    - 1 * student_t_lccdf(0 | 3, 0, 5.2);
  target += std_normal_lpdf(to_vector(z_1));
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
}
generated quantities {
  // actual population-level intercept
  real b_bwt_Intercept = Intercept_bwt;
  // actual population-level intercept
  real b_tarsus_Intercept = Intercept_tarsus;
  // residual correlations
  corr_matrix[nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor);
  vector<lower=-1,upper=1>[nrescor] rescor;
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:nresp) {
    for (j in 1:(k - 1)) {
      rescor[choose(k - 1, 2) + j] = Rescor[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
