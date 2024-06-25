library(mvtnorm)
library(MCMCglmm)
library(tidyverse)
library(ggplot2)
library(gridExtra)

## see functions below

prior_G <- list(V = diag(1), nu = 1000, alpha.mu = 0, alpha.V = diag(1))

# getting simulated priors
pr <- data.frame(
  ag = sim_prior1(prior_G),
  yr = sim_prior1(prior_G),
  pe = sim_prior1(prior_G),
  col = sim_prior1(prior_G),
  res = 1
) %>%
  mutate(
    herit = ag / (rowSums(.))
  )
ggplot(dat, aes(x = pr_herit)) +
  geom_density(bounds = c(0, Inf))

## Writing the prior model for multi-response
dim <- 4
R <- list(V = diag(dim), nu = dim + 0.002, fix = 2)
G <- list(V = diag(dim) * 0.02, nu = dim + 1) # , alpha.mu = rep(0,dim), alpha.V = diag(dim)*1000)#*c(1000,10))

# R <- list(V = diag(c(1, 1, 1, 0.0001), 4, 4), nu = 3.002, fix = 2)
# G <- list(V = diag(4), nu = 4.002, fixed = 4)

spr_G <- sim_prior(G)
spr_R <- sim_prior(R)

spr_P <- spr_G + spr_R

## looking at prior for covaraince and correlation
ggplot(spr_P, aes(x = G4.3)) +
  geom_density(bounds = c(-Inf, Inf))
a <- cov2cor.spr(spr_P)
ggplot(a,aes(x=V15))  + geom_density(bounds = c(-1, 1))





## simulating the prior functions
# Prior for univariate only
sim_prior1 <- function(prior, n = 100000) {
  pr_var_eta <- rIW(n = n, V = prior$V, nu = prior$nu)
  if (!is.null(prior$alpha.mu)) {
    pr_alpha <- rnorm(n, prior$alpha.mu, sqrt(prior$alpha.V))
    pr_var <- (pr_alpha^2) * pr_var_eta
  } else {
    pr_var <- pr_var_eta
  }
  pr_var
}
# Prior for multivariate cases
sim_prior <- function(G, K = 100000) {
  dimG <- ncol(G$V)
  if (!is.null(G$fix)) {
    pr_var_eta <- rIW(n = K, V = G$V, nu = G$nu, fix = G$fix)
  } else {
    pr_var_eta <- rIW(n = K, V = G$V, nu = G$nu)
  }
  if (is.null(G$alpha.mu)) {
    pr_var_G <- lapply(1:K, function(i) {
      matrix(pr_var_eta[i, ], ncol = dimG)
    })
  } else {
    pr_alpha <- rmvnorm(K, G$alpha.mu, G$alpha.V)
    pr_tot <- cbind(pr_alpha, pr_var_eta)

    pr_var_G <- lapply(1:K, function(i) {
      vec <- pr_tot[i, ]
      alpha <- diag(dimG) * vec[1:dimG]
      ETA <- matrix(vec[-c(1:dimG)], ncol = dimG)
      t(alpha) %*% ETA %*% alpha
    })
  }
  ## Combining the priors in a dataframe
  pr_V_G <- as.data.frame(matrix(unlist(pr_var_G), ncol = dimG^2, byrow = TRUE))
  colnames(pr_V_G) <- paste0("G", rep(1:dimG, each = dimG), ".", rep(1:dimG, dimG))
  pr_V_G
}

cov2cor.spr <- function(spr) {
  as.data.frame(t(apply(spr, 1, function(x) {
    c(cov2cor(matrix(unlist(x), nrow = sqrt(length(x)), ncol = sqrt(length(x)), byrow = TRUE)))
  })))
}
