#######################################################################
# set parameters of XBART
get_XBART_params <- function(y) {
    # XBART_params = list(num_trees = 1, # number of trees 
    #                   num_sweeps = 1, # number of sweeps (samples of the forest)
    #                   n_min = 1, # minimal node size
    #                   alpha = 0.95, # BART prior parameter 
    #                   beta = 1.25, # BART prior parameter
    #                   mtry = 10, # number of variables sampled in each split
    #                   burnin = 0,
    #                   no_split_penality = "Auto"
    #                   ) # burnin of MCMC sample
  XBART_params = list(num_trees = 1, # number of trees 
                      num_sweeps = 1, # number of sweeps (samples of the forest)
                      n_min = 1, # minimal node size
                      alpha = 0.95, # BART prior parameter 
                      beta = 1.25, # BART prior parameter
                      mtry = 10, # number of variables sampled in each split
                      burnin = 0,
                      no_split_penality = "Auto"
                      ) # burnin of MCMC sample
  num_tress = XBART_params$num_trees
  XBART_params$max_depth = 250
  XBART_params$num_cutpoints = 50;
  # number of adaptive cutpoints
  XBART_params$tau = 0.01 # var(y) / num_tress # prior variance of mu (leaf parameter)
  return(XBART_params)
}


#######################################################################
library(XBART)

set.seed(100)
new_data = TRUE # generate new data
run_dbarts = FALSE # run dbarts
run_xgboost = FALSE # run xgboost
run_lightgbm = FALSE # run lightgbm
parl = TRUE # parallel computing


small_case = TRUE # run simulation on small data set
verbose = TRUE # print the progress on screen


if (small_case) {
  n = 200 # size of training set
  nt = 100 # size of testing set
  d = 2 # number of TOTAL variables
  dcat = 2 # number of categorical variables
  # must be d >= dcat
  # (X_continuous, X_categorical), 10 and 10 for each case, 20 in total
} else {
  n = 1000000
  nt = 10000
  d = 50
  dcat = 0
}






#######################################################################
# Data generating process

#######################################################################
# Have to put continuous variables first, then categorical variables  #
# X = (X_continuous, X_cateogrical)                                   #
#######################################################################
if (new_data) {
  if (d != dcat) {
    x = matrix(runif((d - dcat) * n, -2, 2), n, d - dcat)
    if (dcat > 0) {
      x = cbind(x, matrix(as.numeric(sample(-2:2, dcat * n, replace = TRUE)), n, dcat))
    }
  } else {
    # x = matrix(as.numeric(sample(-2:2, dcat * n, replace = TRUE)), n, dcat)
    x = matrix(as.numeric(sample(0:1, dcat * n, replace = TRUE)), n, dcat)
  }

  xtest = unique(x)
  # if (d != dcat) {
  #   xtest = matrix(runif((d - dcat) * nt, -2, 2), nt, d - dcat)
  #   if (dcat > 0) {
  #     xtest = cbind(xtest, matrix(as.numeric(sample(-2:2, dcat * nt, replace = TRUE)), nt, dcat))
  #   }
  # } else {
  #   # xtest = matrix(as.numeric(sample(-2:2, dcat * nt, replace = TRUE)), nt, dcat)
  #   xtest = matrix(as.numeric(sample(0:1, dcat * nt, replace = TRUE)), nt, dcat)
  # }

  f = function(x) {
    sin(x[, 2] ^ 2) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 1] + x[, 2] ^ 2) / (3 + x[, 1])
    # sin(x[, 3] ^ 2) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 1] + x[, 2] ^ 2) / (3 + x[, 3])
    # sin(rowSums(x[, 3:4] ^ 2)) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 1] + x[, 2] ^ 2) / (3 + x[, 3] + x[, 4] ^ 2)
    # sin(rowSums(x[, 3:4] ^ 2)) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 5] + x[, 6]) ^ 2 * (x[, 1] + x[, 2] ^ 2) / (3 + x[, 3] + x[, 4] ^ 2)
  }

  # to test if ties cause a crash in continuous variables
  x[, 1] = round(x[, 1], 4)
  #xtest[,1] = round(xtest[,1],2)
  ftrue = f(x)
  ftest = f(xtest)
  sigma = sd(ftrue)

  #y = ftrue + sigma*(rgamma(n,1,1)-1)/(3+x[,d])
  #y_test = ftest + sigma*(rgamma(nt,1,1)-1)/(3+xtest[,d])

  y = ftrue + sigma * rnorm(n)
  y_test = ftest + sigma * rnorm(nrow(xtest)) # nt
  # sample prior from y
  y_prior = y[sample(n, 10)]
}

#######################################################################
# XBART
categ <- function(z, j) {
  q = as.numeric(quantile(x[, j], seq(0, 1, length.out = 100)))
  output = findInterval(z, c(q, + Inf))
  return(output)
}


params = get_XBART_params(y)
time = proc.time()
fit = XBART.density(as.matrix(y), as.matrix(x), as.matrix(xtest), as.matrix(y_prior), p_categorical = dcat,
            params$num_trees, params$num_sweeps, params$max_depth,
            params$n_min, alpha = params$alpha, beta = params$beta, tau = params$tau, s = 1, kap = 1,
            mtry = params$mtry, verbose = verbose, burnin = params$burnin,
            num_cutpoints = params$num_cutpoints, parallel = parl, random_seed = 100, no_split_penality = params$no_split_penality)

################################
# two ways to predict on testing set

# 1. set xtest as input to main fitting function
# fhat.1 = apply(fit$yhats_test[, params$burnin:params$num_sweeps], 1, mean)
fhat.1 = apply(fit$yhats_test, 1, mean)
time = proc.time() - time
print(time[3])

# 2. a separate predict function
pred = predict(fit, xtest)
# pred = rowMeans(pred[, params$burnin:params$num_sweeps])
pred = rowMeans(pred)

time_XBART = round(time[3], 3)

pred2 = predict(fit, xtest)
# pred2 = rowMeans(pred2[, params$burnin:params$num_sweeps])
pred2 = rowMeans(pred2)
stopifnot(pred == pred2)

#######################################################################
# dbarts
if (run_dbarts) {
  library(dbarts)

  time = proc.time()
  fit = bart(x, y, xtest, verbose = FALSE, numcut = 100, ndpost = 1000, nskip = 500)
  time = proc.time() - time
  print(time[3])
  fhat.db = fit$yhat.test.mean
  time_dbarts = round(time[3], 3)
} else {
  fhat.db = fhat.1
  time_dbarts = time_XBART
}


#######################################################################
# XGBoost
if (run_xgboost) {
  library(xgboost)
}



#######################################################################
# LightGBM
if (run_lightgbm) {
  library(xgboost)
}


#######################################################################
# print
xbart_rmse = sqrt(mean((fhat.1 - ftest) ^ 2))
print(paste("rmse of fit xbart: ", round(xbart_rmse, digits = 4)))
# print(paste("rmse of fit dbart: ", round(sqrt(mean((fhat.db - ftest) ^ 2)), digits = 4)))

# print(paste("running time, dbarts", time_dbarts))
print(paste("running time, XBART", time_XBART))


# plot(ftest, fhat.db, pch = 20, col = 'orange')
# points(ftest, fhat.1, pch = 20, col = 'slategray')
# legend("topleft", c("dbarts", "XBART"), col = c("orange", "slategray"), pch = c(20, 20))


# For Travis
# stopifnot(xbart_rmse < 1)
# stopifnot(time_XBART < 5)

# distribution of specific categories
cat_match = function(x, cat){
  if (length(x) != length(cat)){cat('dimension not match')}
  return(all(x==cat))
}

# xtest[1:4, ] = xtest[c(1, 4, 2, 3),]
color = c("lightblue1","darkolivegreen1","lightpink1", "darkseagreen1")
color_dark = c("blue", "green", "red", "darkseagreen4")
h <- hist(y, plot=FALSE, breaks = 15)
h$counts=h$counts/sum(h$counts)
plot(h, main = "", xlim = c(min(y), max(y)), ylim = c(0, 0.5))
color_ind=1
for (test_ind in 1:nrow(xtest)){
  cat1 = xtest[test_ind, ]
  ind = apply(x, 1, cat_match, cat=cat1)
  indt = apply(xtest, 1, cat_match, cat=cat1)
  
  # plot.new()
  h <- hist(y[ind], plot=FALSE, breaks = 15)
  h$counts=h$counts/sum(h$counts)
  plot(h, main = "", xlim = c(min(y), max(y)), col=color[color_ind], add=TRUE)
  # abline(v = ftest[indt], col = 'red')
  # title(toString(cat1))
  color_ind = color_ind+1
}
for (test_ind in 1:nrow(xtest))
{
  lines(x = seq(-2, 4, length.out=100), y = exp(fit$yhats_test[test_ind, 1, ]), type = 'l', col=color_dark[test_ind])
}


# library(PBDE)
# tau=0.01
# y_test = as.matrix(seq(min(y), max(y), length.out=n))
# density_x2_1_x1_0 = density_x2_1_x1_1 = density_x2_0 = rep(n, 0)
# for (i in 1:length(y_test)){
#   density_x2_1_x1_0[i] = exp(p_n(as.matrix(y_test[i]), as.matrix(y[x[,1]==0&x[,2]==1]),tau,TRUE)[[1]])
#   density_x2_1_x1_1[i] = exp(p_n(as.matrix(y_test[i]), as.matrix(y[x[,2]==1&x[,1]==1]),tau,TRUE)[[1]])
#   density_x2_0[i] = exp(p_n(as.matrix(y_test[i]), as.matrix(y[x[,2]==0]), tau,TRUE)[[1]])
# }
# h <- hist(y, plot=FALSE)
# h$counts=h$counts/sum(h$counts)
# plot(h, main = toString(tau), ylim = c(0, 0.5))
# h <- hist(y[x[,2]==0], plot=FALSE)
# h$counts=h$counts/sum(h$counts)
# plot(h, main = toString(tau), col = "lightpink1", add=T)
# h <- hist(y[x[,1]==0&x[,2]==1], plot=FALSE)
# h$counts=h$counts/sum(h$counts)
# plot(h, main = toString(tau), col = "darkolivegreen1", add=T)
# h <- hist(y[x[,1]==1&x[,2]==1], plot=FALSE)
# h$counts=h$counts/sum(h$counts)
# plot(h, main = toString(tau), col = "lightskyblue1", add=T)
# lines(y_test, density_x2_1_x1_0, type = 'l', col = "green")
# lines(y_test, density_x2_1_x1_1, type = 'l', col = 'blue')
# lines(y_test, density_x2_0, type = 'l', col = 'red')




