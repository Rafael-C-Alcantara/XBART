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
require(ggplot2)
require(tidyr)

set.seed(90)
new_data = TRUE # generate new data
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

  f = function(x) {
    sin(x[, 2] ^ 2) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 1] + x[, 2] ^ 2) / (3 + x[, 1])
    # sin(x[, 3] ^ 2) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 1] + x[, 2] ^ 2) / (3 + x[, 3])
    # sin(rowSums(x[, 3:4] ^ 2)) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 1] + x[, 2] ^ 2) / (3 + x[, 3] + x[, 4] ^ 2)
    # sin(rowSums(x[, 3:4] ^ 2)) + sin(rowSums(x[, 1:2] ^ 2)) + (x[, 5] + x[, 6]) ^ 2 * (x[, 1] + x[, 2] ^ 2) / (3 + x[, 3] + x[, 4] ^ 2)
  }

  # to test if ties cause a crash in continuous variables
  x[, 1] = round(x[, 1], 4)
  ftrue = f(x)
  sigma = sd(ftrue)

  #y = ftrue + sigma*(rgamma(n,1,1)-1)/(3+x[,d])

  y = ftrue + sigma * rnorm(n)
  # sample prior from y
  y_prior = y[sample(n, 10)]
  y_range = c(min(y), max(y))
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
fit = XBART.density(as.matrix(y), as.matrix(x), as.matrix(y_prior), as.matrix(y_range), p_categorical = dcat,
            params$num_trees, params$num_sweeps, params$max_depth,
            params$n_min, alpha = params$alpha, beta = params$beta, tau = params$tau, s = 1, kap = 1,
            mtry = params$mtry, verbose = verbose, burnin = params$burnin,
            num_cutpoints = params$num_cutpoints, parallel = parl, random_seed = 100, no_split_penality = params$no_split_penality)

################################
# two ways to predict on testing set

# 1. set xtest as input to main fitting function

time = proc.time() - time
print(time[3])

#######################################################################
# # print
# xbart_rmse = sqrt(mean((fhat.1 - ftest) ^ 2))
# print(paste("rmse of fit xbart: ", round(xbart_rmse, digits = 4)))
# # print(paste("rmse of fit dbart: ", round(sqrt(mean((fhat.db - ftest) ^ 2)), digits = 4)))

# # print(paste("running time, dbarts", time_dbarts))
# print(paste("running time, XBART", time_XBART))


# # plot(ftest, fhat.db, pch = 20, col = 'orange')
# # points(ftest, fhat.1, pch = 20, col = 'slategray')
# # legend("topleft", c("dbarts", "XBART"), col = c("orange", "slategray"), pch = c(20, 20))


# # For Travis
# # stopifnot(xbart_rmse < 1)
# # stopifnot(time_XBART < 5)

# distribution of specific categories

getcutpoints = function(cutpoints){
  cutpoints = matrix(cutpoints, nrow = length(cutpoints)/3, byrow = TRUE)
  output = ''
  for(i in 1:nrow(cutpoints)){
    if (cutpoints[i, 3]==1){
      output <- paste(output, 'X', toString(cutpoints[i, 1]),'<=', toString(cutpoints[i, 2]), sep="")
    } 
    else{output = paste(output, 'X', toString(cutpoints[i, 1]),'>', toString(cutpoints[i, 2]), sep="")}
    output = paste(output, "; ", sep="")
  }
  return(output)
}

n_bot = length(fit$density_info$cutpoints[[1]][[1]])
cutpoints = rep('', n_bot)
for(i in 1:n_bot) {cutpoints[i] = getcutpoints(fit$density_info$cutpoints[[1]][[1]][[i]])}
density = as.data.frame(fit$density_info$density[[1]][[1]], col.names = 1:n_bot)
density = gather(density, key = "group", value = "density")
p <- ggplot(data=as.data.frame(y), aes(y)) + 
  stat_bin(aes(y=..density..), fill = "grey69") +
  geom_line(data=density, aes(x = rep(seq(y_range[1], y_range[2], length.out=100), n_bot),
                      y = density, group = group, colour=group)) +
  scale_color_brewer(palette="Dark2")
p

