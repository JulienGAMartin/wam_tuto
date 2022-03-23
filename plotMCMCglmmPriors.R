library(mvtnorm)
library(MCMCglmm)
library(ggplot2)
library(gridExtra)

##Writing the prior model for monovariate
K <- 100000	#Number of iterations
nu <- 1000	#nu parameter
alpha.mu <- 0	#alpha.mu parameter
alpha.V <- 1	#alpha.V parameter

#Prior for extended parameter Va
pr_var_eta <- rIW(n = K,V = diag(1),nu = nu)[,1]
pr_alpha <- rnorm(K, alpha.mu, sqrt(alpha.V))
pr_var_a <- (pr_alpha*2) * pr_var_eta

#Prior for residual variance
pr_var_r <- rIW(n = K, V = diag(1), nu = 0.01)[,1]

#Resulting prior for heritability
pr_herit <- pr_var_a / (pr_var_a + 1)

#Plotting the density of the prior
qplot(x=pr_herit,geom="density")

##Writing the prior model for multi-response
K <- 100000		#Number of iterations
nu <- 2.002         	#nu parameter
dim <- 2		#dimensions of the multivariate model
V <- diag(dim)*5.002
alpha.mu <- c(0,0)	#alpha.mu parameter
alpha.V <- diag(dim)*c(1000)	#alpha.V parameter

dim <- 4
R <- list(V = diag(dim), nu = dim + 0.002)#, fix = 2)
G <- list(V = diag(dim) * 0.02, nu = dim + 1)#, alpha.mu = rep(0,dim), alpha.V = diag(dim)*1000)#*c(1000,10))

  R = list(V = diag(c(1, 1, 1, 0.0001), 4, 4), nu = 3.002, fix = 2)
  G = list(V = diag(4), nu = 4.002, fixed = 4)

spr <- sim.priors(R = R, G = G)
a<-cov2cor.spr(spr$G)
qplot(a[,15], geom="density")
plot.priors(spr$G,"cor", dim = dim)
plot.priors(spr,"cov")



sim.priors <- function(K = 100000, R, G){
#Prior for extended parameter Va
    dimG <- ncol(G$V)
    dimR <- ncol(R$V)
    if(!is.null(G$alpha.mu)){
        pr_var_eta <- rIW(n = K, V = G$V, nu = G$nu)
        pr_alpha <- rmvnorm(K, G$alpha.mu, G$alpha.V)
        pr_tot <- cbind(pr_alpha, pr_var_eta)

#Prior for G matrix
        pr_var_G <- lapply(1:K, function(i){
            vec <- pr_tot[i,]
            alpha <- diag(vec[1:dimG])
            ETA <- matrix(vec[-c(1:dimG)], ncol = dimG)
            t(alpha) %*% ETA %*% alpha
        })
    }
    else {
        if(!is.null(G$fix)) {
          pr_var_eta <- rIW(n=K, V = G$V, nu = G$nu, fix = G$fix)
        } else {
          pr_var_eta <- rIW(n = K, V = G$V, nu = G$nu)
        }
        pr_var_G <- lapply(1:K, function(i){
            matrix(pr_var_eta[i,], ncol = dimG)
        })
    }

## Prior for R matrix
    if(is.null(R$fix)) {
        tmp <- rIW(n = K, V = R$V, nu = R$nu)
    } else {
        tmp <- rIW(n = K, V = R$V, nu = R$nu, fix = R$fix)
    }
    pr_var_R <- lapply(1:K, function(i){
        matrix(tmp[i,], ncol = dimR)
    })

## Combining the priors in a dataframe
    pr_V_G <- as.data.frame(matrix(unlist(pr_var_G), ncol = dimG^2, byrow = TRUE))
    colnames(pr_V_G) <- paste0("G",rep(1:dimG, each = dimG), ".", rep(1:dimG, dimG))
    pr_V_R <- as.data.frame(matrix(unlist(pr_var_R), ncol = dimR^2, byrow = TRUE))
    colnames(pr_V_R) <- paste0("R",rep(1:dimR, each = dimR), ".", rep(1:dimR, dimR))
    pr <- list(G = pr_V_G, R = pr_V_R)
    pr
}



cov2cor.spr <- function(spr) {
  as.matrix(t(apply(spr, 1, function(x) {
    c(cov2cor(matrix(unlist(x), nrow = sqrt(length(x)), ncol = sqrt(length(x)), byrow = TRUE)))
  }
  )))
}



K <- 100000 # Number of iterations
nu <- 0.01 # nu parameter
alpha.mu <- 0 # alpha.mu parameter
alpha.V <- 1
# Prior for extended parameter Va
pr_var_a <- rIW(n = K, V = diag(1), nu = 1)[, 1]


# Prior for residual variance
pr_var_r <- rIW(n = K, V = diag(1), nu = 1)[, 1]

# Resulting prior for heritability
pr_herit <- pr_var_a / (pr_var_a + pr_var_r)

# Plotting the density of the prior
qplot(x = pr_herit, geom = "density") #+ xlim(0,1)
