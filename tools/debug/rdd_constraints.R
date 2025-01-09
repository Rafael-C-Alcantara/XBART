set.seed(0)
## DGP
n <- 20000
ate <- 0.5
k1 <- 4 ## variability in mu0.x
k2 <- 0.75 ## amplitude of tau0.w relative to ATE (set from 0 to 1)
k3 <- 2 ## sd of mu0.w relative to sd of tau
sig_error <- 0.5 ## relative to sd of tau
p <- 2 # Dim of w
rho <- 0.5
c <- 1
pts_in_window <- 50
x.center <- 0 ## In case we want to center x at a different location
### Functions
mu0.x <- function(x,k,ate)
{
  -0.03*x^5 + k*0.5*x^3 + 0.1*x^2 - 0.1*x + 2
}
mu0.w <- function(w,k) 
{
  k*cos(rowMeans(w))
}
tau0.x <- function(x)
{
  0.25/(3+x-min(x))
}
tau0.w <- function(w,k)
{
  k*sin(3*rowMeans(w))
}
mu0 <- function(x,w,k1,k3,ate)
{
  mu0.x(x,k1,ate) + mu0.w(w,k3)
}
tau0 <- function(x,w,ate,k2)
{
  tau0.x(x) + tau0.w(w,k2)
}
mu <- function(x,w,k1,k3,ate,mu.mean)
{
  mu0(x,w,k1,k3,ate) - mu.mean
}
tau <- function(x,w,ate,k2,tau.mean)
{
  tau0(x,w,ate,k2) - tau.mean + ate
}
h.grid <- function(x,c,grid)
{
  abs.x <- sort(abs(x-c))
  out <- rep(0,length(grid))
  names(out) <- grid
  for(total in grid)
  {
    i <- 1
    sum.right <- sum.left <- 0
    while(sum.right < total | sum.left < total) 
    {
      sum.left <- sum(c-abs.x[i] <= x & x < c)
      sum.right <- sum(c < x & x <= c+abs.x[i])
      if (sum.left == sum(x<c) & sum.right == sum(c<x)) break
      i <- i+1
    }
    out[as.character(total)] <- abs.x[i]
  }
  return(out)
}
### Set parameters and demean data
x0 <- rnorm(n,x.center,1)
w0 <- as.matrix(rnorm(n,x.center,1))
k2 <- k2*ate
k3 <- k3*sd(tau0(c,w0,ate,k2))/sd(mu0(c,w0,k1,1,ate))
mu.bar <- mean(mu0(c,w0,k1,k3,ate))
tau.bar <- mean(tau0(c,w0,ate,k2))
sig_error <- sig_error*max(abs(mean(tau(c,w0,ate,k2,tau.bar))),2*sd(tau(c,w0,ate,k2,tau.bar)))
### Samples
x <- rnorm(n,x.center,1)
z <- as.numeric(x>=c)
w <- matrix(rnorm(n*p,rep((x-mean(x)),p)*rho,sqrt(1-rho^2)),n,p)
prog <- mu(x,w,k1,k3,ate,mu.bar)
cate <- tau(x,w,ate,k2,tau.bar)
y <- prog + cate*z + rnorm(n,0,sig_error)
## Fit model
p_categorical <- 0
burnin <- 0
num_sweeps <- 100
Omin <- 1
Opct <- 0.3
ntrees_con <- 5
ntrees_mod <- 5
Owidth <- h.grid(x,c,pts_in_window)
sample <- 3*Owidth
train <- c-sample < x & x < c+sample
fit <- XBART::XBCF.rd(y[train], w[train,], x[train], c,
                      Owidth = Owidth, Omin = Omin, Opct = Opct,
                      num_cutpoints = n,
                      num_trees_con = ntrees_con, num_trees_mod = ntrees_mod,
                      num_sweeps = num_sweeps,
                      burnin = burnin,
                      Nmin = 1,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical,
                      tau_con = var(y[train])/ntrees_con, tau_mod = var(y[train])/ntrees_mod)
fit$cutoff_nodes_con
fit$invalid_nodes_1_con
fit$invalid_nodes_2_con
fit$cutoff_nodes_mod
fit$invalid_nodes_1_mod
fit$invalid_nodes_2_mod
colSums(fit$invalid_nodes_2_con)
colSums(fit$invalid_nodes_2_mod)