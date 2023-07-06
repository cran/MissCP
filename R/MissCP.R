## usethis namespace: start
#' @useDynLib MissCP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end

#' @title MCAR
#' @description  function to do the missing assuming the missing completely at random
#' @param data data before the missing case
#' @param alpha the percentage of missing compared to whole data
#' @return the data matrix with missing values
#' @export
MCAR <- function(data, alpha){
  if (!is.matrix(data)){
    data = as.matrix(data)
  }
  if (alpha > 1 || alpha < 0 ){
    stop("alpha has to between 0 and 1")
  }
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      if(rbinom(1,1,alpha) == 1)(
        data[i,j] = NA
      )
    }
  }
  return(data)
}

#' @title Heter_missing
#' @description function to do the missing assuming the missing completely at random
#' @param data data before the missing case
#' @param alpha the list of percentage of missing compared to whole data
#' @return the data matrix with missing values
#' @export
Heter_missing <- function(data, alpha){
  if (!is.matrix(data)){
    data = as.matrix(data)
  }
  if (ncol(data) != length(alpha)){
    stop("length of alpha should equal to the dimension")
  }
  for(i in 1:nrow(data)){
    for(j in 1:ncol(data)){
      if (alpha[j] > 1 || alpha[j] < 0 ){
        stop("alpha has to between 0 and 1")
      }
      if(rbinom(1,1,alpha[j]) == 1)(
        data[i,j] = NA
      )
    }
  }
  return(data)
}

#' @title imputation
#' @description function to do the imputation based on block size
#' @param data data before the imputation
#' @param block.size the block size that are used to impute the missing
#' @return the data matrix without missing values after imputation
#'
imputation <- function(data, block.size){
  if (!is.matrix(data)){
    data = as.matrix(data)
  }
  TT  <- length(data[, 1])
  # blocks <- seq(1, TT + 1, block.size);
  b_n_bound = 2*block.size  #block size for boundary
  blocks <- c(seq(1, b_n_bound*2+1 , b_n_bound),
              seq(b_n_bound*2+block.size+1, TT+1-2*b_n_bound, block.size),
              seq(TT+1-b_n_bound, TT+1,  b_n_bound))

  if(blocks[length(blocks)] < TT+1){
    blocks <- c(blocks[-length(blocks)], TT+1)
  }
  n.new <- length(blocks) - 1
  for(i in 1 : n.new){
    sum = matrix(0, nrow = 1, ncol = ncol(data))
    count = matrix(0, nrow = 1, ncol = ncol(data))
    mean = matrix(0, nrow = 1, ncol = ncol(data))
    for(j in blocks[i]:(blocks[i+1]-1)){
      for(k in 1: ncol(data)){
        if(!is.na(data[j,k])){
          sum[1,k] = sum[1,k] + data[j,k]
          count[1,k] = count[1,k] + 1
        }
      }
    }
    for(k in 1 : ncol(data)){
      mean[1,k] = sum[1,k]/count[1,k]
      if(is.na(mean[1,k])){
        mean[1,k] = 0
      }
    }
    for(j in blocks[i]:(blocks[i+1]-1)){
      for(k in 1: ncol(data)){
        if(is.na(data[j,k])){
          data[j,k] = mean[1,k]
        }
      }
    }
  }
  return(data)
}

#' @title constant_generation
#' @description function to generate constant given jump size and break points
#' @param n the sample size
#' @param p the data dimension
#' @param d the number of nonzero coeddficients
#' @param vns the jump size. It can be a vector or a single value. If single value, it is same for all break points
#' @param brk the break points' locations
#' @return the parameter matrix used to generate data
#' @import mvtnorm
#' @export
constant_generation <- function(n, p, d, vns, brk){
  p.y = p
  m = length(brk)
  constant.full <- matrix(0, p.y, m)
  constant.full[sample(1:p.y, d, replace = FALSE), 1] <- runif(d, 0.5, 1);
  if(length(vns) == 1 && m > 2){
    vns = rep(vns, (m-1))
  }
  vns = c(0, vns)
  if(m == 1){
    return(constant.full)
  }
  else{
    for(temp in 2:m){
      constant.full[, temp] <- constant.full[, (temp - 1)];
      a = sample(1:p.y, d, replace = FALSE)
      diff = sqrt((vns[temp]/16)**2 / d/2)
      if(temp%%2 == 0){
        constant.full[a, temp] <- constant.full[a, (temp - 1)] + diff;
      }
      else{
        constant.full[a, temp] <- constant.full[a, (temp - 1)] - diff;
      }
    }
  }
  return(constant.full)
}


#' @title imputation2
#' @description function to do the imputation based on change point candidate
#' @param data data before the imputation
#' @param cp.candidate the change point candidate that are used to impute the missing
#' @return the data matrix without missing values after imputation
#'
imputation2 <- function(data, cp.candidate){
  n = length(cp.candidate)/2
  for(i in 1:n){
    sum = matrix(0, nrow = 1, ncol = ncol(data))
    count = matrix(0, nrow = 1, ncol = ncol(data))
    mean = matrix(0, nrow = 1, ncol = ncol(data))
    for(j in cp.candidate[2*i-1] : (cp.candidate[2*i]-1)){
      for(k in 1: ncol(data)){
        if(!is.na(data[j,k])){
          sum[1,k] = sum[1,k] + data[j,k]
          count[1,k] = count[1,k] + 1
        }
      }
    }
    for(k in 1 : ncol(data)){
      mean[1,k] = sum[1,k]/count[1,k]
    }
    for(j in cp.candidate[2*i-1] : (cp.candidate[2*i]-1)){
      for(k in 1: ncol(data)){
        if(is.na(data[j,k])){
          data[j,k] <- mean[1,k]
        }
      }
    }
  }
  return(data)
}

#' @title pred
#' @description function to do the prediction
#' @param X data for prediction
#' @param phi parameter matrix
#' @param j the start time point for prediction
#' @param p.x the dimension of data X
#' @param p.y the dimension of data Y
#' @param h the length of observation to predict
#' @return prediction matrix
#'
pred <- function(X, phi, j, p.x, p.y, h = 1){
  concat.X <- matrix(0, p.x, 1);
  concat.Y <- matrix(0, p.y, 1);
  concat.X[,1] <- as.matrix(X[, j]);
  temp <- matrix(0, p.y, 1);
  if(p.y == 1){
    temp <- temp +  t(as.matrix(phi[, 1:p.x]))%*%concat.X[, 1];
  }else{
    temp <- temp +  as.matrix(phi[, 1:p.x])%*%concat.X[, 1];

  }
  concat.Y[, 1] <- temp;
  return(as.matrix(concat.Y[, 1]))
}

#' @title pred.block
#' @description Prediction function (block)
#' @param X data for prediction
#' @param phi parameter matrix
#' @param j the start time point for prediction
#' @param p.x the dimension of data X
#' @param p.y the dimension of data Y
#' @param h the length of observation to predict
#' @return prediction matrix
#'
pred.block <- function(X, phi, j, p.x, p.y, h){
  concat.X <- matrix(0, p.x, h);
  concat.Y <- matrix(0, p.y, h);
  concat.X[, 1:h] <- as.matrix(X[, (j):(j+h-1)]);
  for ( i in 1:h){
    temp <- matrix(0, p.y, 1);
    if(p.y == 1){
      temp <- temp +  t(as.matrix(phi[, (1):(p.x)]))%*%concat.X[, i];
    }else{
      temp <- temp +  as.matrix(phi[, (1):(p.x)])%*%concat.X[, i];
    }

    concat.Y[, i] <- temp;
  }
  return(as.matrix(concat.Y[, 1:h]))
}

#' @title BIC
#' @description BIC  and HBIC function
#' @param residual residual matrix
#' @param phi estimated coefficient matrix of the model
#' @return A list object, which contains the followings
#' \describe{
#'   \item{BIC}{BIC value}
#'   \item{HBIC}{HBIC value}
#' }
BIC <- function(residual, phi){
  p.y <- length(phi[, 1]);
  p.x <- length(phi[1, ]);
  T.new <- length(residual[1, ]);
  count <- 0;
  for (i in 1:p.y){
    for (j in 1:p.x){
      if(phi[i,j] != 0){
        count <- count + 1;
      }
    }
  }
  sigma.hat <- 0*diag(p.y);
  for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[, i]%*%t(residual[, i]);  }
  sigma.hat <- (1/(T.new))*sigma.hat;
  ee.temp <- min(eigen(sigma.hat)$values);
  if(ee.temp <= 10^(-8)){
    sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp) + 10^(-3))*diag(p.y);
  }
  log.det <- log(det(sigma.hat));
  return(list(BIC = log.det + log(T.new)*count/T.new , HBIC = log.det + 2*log(p.x*p.y)*count/T.new))
}

#' @title BIC_threshold
#' @description BIC threshold for final parameter estimation
#' @param beta.final estimated parameter coefficient matrices
#' @param k dimensions of parameter coefficient matrices
#' @param m.hat number of estimated change points
#' @param brk vector of estimated change points
#' @param data_y input data matrix (response), with each column representing the time series component
#' @param data_x input data matrix (predictor), with each column 1
#' @param b_n the block size
#' @param nlam number of hyperparameters for grid search
#' @return lambda.val.best, the tuning parameter lambda selected by BIC.
#'
BIC_threshold <- function(beta.final, k, m.hat, brk, data_y, data_x = NULL,
                          b_n = 2, nlam = 20){
  brk.full <- c(1, brk)
  jj = 1; flag = 0
  if(!is.null(beta.final)){
    lambda.val.best <- c()
    flag = 1
    for(i in 1:m.hat){
      temp <- unlist(beta.final[, ((i-1)*k+1):(i*k)] )
      lambda.max <-  max(abs(temp))
      if(lambda.max > 0){
        lambda.min <-  min(abs(temp[temp!=0]))
        if(lambda.max/lambda.min >= 10^4){
          nlam <- 50
        }
        if(lambda.max/lambda.min >= 10^8){
          lambda.min <-  lambda.max*10^(-4)
        }
        delata.lam <- (log(lambda.max)-log(lambda.min))/(nlam -1)
        lambda.val.full <-  sapply(1:(nlam), function(jjj) lambda.min*exp(delata.lam*(nlam-jjj)))
        mse.res <- c()
        BIC.res <- c()
        for(j in 1:nlam){
          lambda.val = lambda.val.full[j]
          beta.temp <- beta.final[((i-1)*k+1):(i*k)]
          beta.temp[abs(beta.temp) < lambda.val] <- 0
          data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
          data_x.temp <- matrix(1, nrow = brk.full[i+1] - brk.full[i])
          data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, nrow = 1)
          residual.temp <- data_y.temp - data_y.est
          BIC.res <- c(BIC.res, BIC(t(residual.temp[b_n : (nrow(residual.temp)-b_n), ]), phi = matrix(beta.temp, ncol = 1) )$BIC )
        }
      }else{
        lambda.min <- 0
        lambda.val.full <- c(0)
        BIC.res <- c(0)
        flag = 0
      }
      lambda.val.best <- c(lambda.val.best, lambda.val.full[which.min(BIC.res)])
    }
  }
  return(lambda.val.best)
}


#' @title data_generation
#' @description The function to generate mean shift data
#' @param n the number of data points
#' @param mu the matrix of mean parameter
#' @param sigma covariance matrix of the white noise
#' @param brk vector of change points
#' @return data_y matrix of generated mean shift data
#' @import mvtnorm
#' @export
data_generation <- function (n, mu, sigma, brk = n+1) {
  if (!is.matrix(sigma))
    sigma = as.matrix(sigma)
  p.y <-  nrow(sigma)
  m <- length(brk)
  # error term
  data_e = rmvnorm(n, rep(0, p.y), sigma)
  # data y
  data_y = matrix(0, n, p.y)
  if (m == 1){
    for (i in 1:n) {
      tmp = matrix(data_e[i, ], 1, p.y)
      data_y[i, ] = mu + tmp
    }
  }

  if (m > 1){
    for (i in 1:(brk[1]-1)) {
      tmp = matrix(data_e[i, ], 1, p.y)
      mu.temp = mu[, 1]
      data_y[i, ] = tmp + mu.temp
    }
    for (mm in 1:(m-1)){
      for (i in (brk[mm]):(brk[mm+1]-1) ) {
        tmp = matrix(data_e[i, ], 1, p.y)
        mu.temp = mu[, mm + 1]
        data_y[i, ] = tmp + mu.temp
      }
    }
  }
  data_y = data_y[1:n, ]
  return(data_y)
}




#' @title first.step
#' @description Perform the block fused lasso with thresholding to detect candidate break points.
#' @param data_y input data matrix Y, with each column representing the time series component
#' @param data_x input data matrix X
#' @param lambda1 tuning parmaeter lambda_1 for fused lasso
#' @param lambda2 tuning parmaeter lambda_2 for fused lasso
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param blocks the blocks
#' @param cv.index the index of time points for cross-validation
#' @param fixed_index index for linear regression model with only partial compoenents change.
#' @param nonfixed_index index for linear regression model with only partial compoenents change.
#' @return A list object, which contains the followings
#' \describe{
#'   \item{jump.l2}{estimated jump size in L2 norm}
#'   \item{jump.l1}{estimated jump size in L1 norm}
#'   \item{pts.list}{estimated change points in the first step}
#'   \item{beta.full}{estimated parameters in the first step}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
#' @import factoextra
#'
first.step <- function(data_y, data_x, lambda1, lambda2, max.iteration = max.iteration, tol = tol,
                                  blocks, cv.index, fixed_index = NULL, nonfixed_index = NULL){
  lambda.full <- expand.grid(lambda1, lambda2)
  kk <- length(lambda.full[, 1]);

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  cv.l <- length(cv.index);
  cv <- rep(0, kk);
  phi.final <- vector("list", kk);
  phi.final.2 <- vector("list", kk);

  data.y.temp <- data_y
  data.x.temp <- data_x
  TT <- length(data.y.temp[,1]);
  p.y <- length(data.y.temp[1,]); p.x.all <- length(data.x.temp[1, ]);
  p.x <- p.x.all - length(fixed_index);

  flag.full <- rep(0, kk);

  for (i in 1:kk) {
    if(!is.null(fixed_index)){
      if ( i == 1){
        test <- lm_partial_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1],lambda.full[i,2],
                                           max_iteration = max.iteration, tol = tol,
                                           initial_phi =  0.0+matrix(0.0,p.y,p.x*n.new),
                                           initial_phi_2 =  0.0+matrix(0.0,p.y,(p.x.all-p.x)),
                                           blocks = blocks, cv.index, fixed_index, nonfixed_index)
        flag.full[i] <- test$flag;
      }
      if ( i > 1 ){
        initial.phi <- phi.final[[i-1]]
        initial.phi.2 <- phi.final.2[[i-1]]
        if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
        test <- lm_partial_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1], lambda.full[i,2],
                                           max_iteration = max.iteration, tol = tol,
                                           initial_phi =  initial.phi, initial_phi_2 =  initial.phi.2,
                                           blocks = blocks, cv.index, fixed_index, nonfixed_index)
        flag.full[i] <- test$flag;
      }
      phi.hat.full <- test$phi.hat;
      phi.final[[i]] <- phi.hat.full;
      phi.hat.full.2 <- test$phi.hat.2;
      phi.final.2[[i]] <- phi.hat.full.2;
      phi.full.all <- vector("list",n.new);
      forecast <- matrix(0,p.y,TT);
      forecast.new <- matrix(0,p.y,cv.l);
      phi.hat <- phi.hat.full;
      phi.full.all[[1]] <- matrix(phi.hat[,(1):(p.x)],ncol = p.x);
      for(i.1 in 2:n.new){
        phi.full.all[[i.1]] <- matrix(phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        forecast[,(blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x[,nonfixed_index]),
                                                                 phi.full.all[[i.1-1]], blocks[i.1],
                                                                 p.x, p.y, blocks[i.1+1]-blocks[i.1]);
      }

      forecast.new <- matrix(0,p.y,cv.l);
      for(j in (1):cv.l){
        forecast.new[,j] <- pred(t(data_x[,nonfixed_index]),phi.full.all[[(cv.index[j])]],
                                 blocks[cv.index[j]+1]-1,p.x, p.y)
        forecast.new[,j] <- forecast.new[,j] + phi.hat.full.2%*%as.matrix(data_x[blocks[cv.index[j]+1]-1,
                                                                          fixed_index])
      }

      temp.index <- rep(0,cv.l);
      for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1]-1;}
      cv[i] <- (1/(p.y*cv.l))*sum( (forecast.new - t(data_y[temp.index,])  )^2 );


    }
    if(is.null(fixed_index)){
      if ( i == 1){
        test <- lm_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1], lambda.full[i,2],
                                   max_iteration = max.iteration, tol = tol,
                                   initial_phi =  0.0+matrix(0.0,p.y,p.x*n.new), blocks = blocks, cv.index)
        flag.full[i] <- test$flag;
      }
      if ( i > 1 ){
        initial.phi <- phi.final[[i-1]]
        if(max(abs(phi.final[[(i-1)]])) > 10^3 | flag.full[i-1] == 1  ){
          initial.phi <- 0*phi.final[[(i-1)]];
        }
        test <- lm_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1], lambda.full[i,2],
                                   max_iteration = max.iteration, tol = tol, initial_phi =  initial.phi,
                                   blocks = blocks, cv.index)
        flag.full[i] <- test$flag;
      }
      phi.hat.full <- test$phi.hat;
      phi.final[[i]] <- phi.hat.full;


      #forecast the time series based on the estimated matrix Phi (beta for the linear regression model)
      #and compute the forecast error
      phi.full.all <- vector("list", n.new);
      forecast <- matrix(0, p.y, TT);
      forecast.new <- matrix(0, p.y, cv.l);

      phi.hat <- phi.hat.full;
      phi.full.all[[1]] <- matrix(phi.hat[,(1):(p.x)], ncol = p.x);
      forecast[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data_x), phi.full.all[[1]], blocks[1],
                                                          p.x, p.y, blocks[2]-blocks[1]);
      for(i.1 in 2:n.new){
        phi.full.all[[i.1]] <- matrix(phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        forecast[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x), phi.full.all[[i.1]],
                                                                  blocks[i.1], p.x, p.y,
                                                                  blocks[i.1+1]-blocks[i.1]);
      }
      forecast.new <- matrix(0, p.y, cv.l);
      for(j in (1):cv.l){
        forecast.new[, j] <- pred(t(data_x), phi.full.all[[(cv.index[j])]], blocks[cv.index[j]+1]-1, p.x, p.y)
      }
      temp.index <- rep(0, cv.l);
      for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1]-1;}
      cv[i] <- (1/(p.y*cv.l))*sum( (forecast.new - t(data_y[temp.index, ])  )^2 );
    }
  }
  lll <- min(which(cv==min(cv)));


  phi.hat.full <- phi.final[[lll]];
  beta.fixed.full <- phi.final.2[[lll]];

  jumps.l2 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l2[i] <- (sum((phi.hat.full[,((i-1)*p.x+1):(i*p.x)] )^2 ));
  }
  jumps.l1 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l1[i] <- (sum(abs(phi.hat.full[,((i-1)*p.x+1):(i*p.x)] ) ));
  }

  jumps <- jumps.l2
  ignore_num <- min(5 , round(200/mean(blocks.size)))
  jumps[1:ignore_num]  <- 0
  jumps[(length(jumps)-(ignore_num-1)):length(jumps)]  <- 0



  BIC.diff <- 1
  BIC.old <- 10^8
  pts.sel <- c()
  loc.block.full <- c()
  while(BIC.diff > 0 & length(unique(jumps)) > 1 ){
    pts.sel.old <- pts.sel
    #use kmeans to hard threshold the jumps
    if( length(unique(jumps)) > 2 ){
      clus.2 <- kmeans(jumps, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss;
      if(fit.2 < 0.20){
        pts.sel <- c(pts.sel);
      }
      if( fit.2 >= 0.20 ){
        loc <- clus.2$cluster;
        if( clus.2$centers[1] > clus.2$centers[2]  ){
          loc.block <- which(loc==1);
        }
        if( clus.2$centers[1] < clus.2$centers[2]  ){
          loc.block <- which(loc==2);
        }
        pts.sel <- c(pts.sel, blocks[loc.block]);
        loc.block.full <- c(loc.block.full, loc.block)
      }
    }
    if( length(unique(jumps)) <= 2 ){
      if(length(unique(jumps)) == 2){
        loc.block <- which.max(jumps)
        pts.sel <- c(pts.sel, blocks[loc.block])
        loc.block.full <- c(loc.block.full, loc.block)
      }else{
        pts.sel <- c(pts.sel);
      }
    }

    phi.hat.full.new <- phi.hat.full
    for(i in 2:n.new){
      if(!(i %in% loc.block.full)){
        phi.hat.full.new[, ((i-1)*p.x+1):(i*p.x)] <- matrix(0, p.y, p.x)
      }
    }


    phi.full.all.new <- vector("list", n.new);
    phi.full.all.new.temp <- vector("list", n.new);
    forecast.all.new <- matrix(0, p.y, TT);

    phi.hat.new <- phi.hat.full.new;
    phi.full.all.new[[1]] <- matrix(phi.hat.new[, (1):(p.x)], ncol = p.x);
    phi.full.all.new.temp[[1]] <- matrix(phi.hat.new[, (1):(p.x)], ncol = p.x);

    if(!is.null(fixed_index)){
      forecast.all.new[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data_x[, nonfixed_index]), phi.full.all.new[[1]], blocks[1], p.x, p.y, blocks[2] - blocks[1]);
      forecast.all.new[, (blocks[1]):(blocks[2]-1)] <- forecast.all.new[, (blocks[1]):(blocks[2]-1)] + beta.fixed.full %*%t(as.matrix(data_x[(blocks[1]):(blocks[2]-1), fixed_index]))

      for(i.1 in 2:n.new){
        phi.full.all.new[[i.1]] <- matrix(phi.full.all.new[[i.1-1]] + phi.hat.new[,((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x[, nonfixed_index]), phi.full.all.new[[i.1]], blocks[i.1], p.x, p.y, blocks[i.1+1] - blocks[i.1]);
        forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] + beta.fixed.full %*%t(as.matrix(data_x[(blocks[i.1]):(blocks[i.1+1]-1), fixed_index]))
      }
    }
    if(is.null(fixed_index)){

      forecast.all.new[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data_x), phi.full.all.new[[1]],
                                                                  blocks[1], p.x, p.y, blocks[2] - blocks[1]);
      for(i.1 in 2:n.new){
        #phi.full.all.new.temp keeps adding
        phi.full.all.new.temp[[i.1]] <- matrix(phi.full.all.new.temp[[i.1-1]] + phi.hat.full[, ((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        if((i.1 %in% loc.block.full)){
          phi.full.all.new[[i.1]] <- phi.full.all.new.temp[[i.1]]
        }
        if(!(i.1 %in% loc.block.full)){
          phi.full.all.new[[i.1]] <- phi.full.all.new[[i.1-1]]
        }
        forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x), phi.full.all.new[[i.1]],
                                                                          blocks[i.1], p.x, p.y,
                                                                          blocks[i.1+1] - blocks[i.1]);
      }
    }
    residual <- t(data.y.temp[(1:TT), ]) - forecast.all.new;


    BIC.new <- BIC(residual, phi = phi.hat.full.new )$BIC
    BIC.diff <- BIC.old - BIC.new
    if(is.na(BIC.diff)){
      BIC.diff = 0
    }
    BIC.old <- BIC.new
    if(BIC.diff <= 0){
      pts.sel <- sort(pts.sel.old)
      break
    }
    jumps[loc.block] <- 0
  }


  if ( length(pts.sel) == 0){cp.final <- c(); pts.list <-  vector("list", 0);}
  if ( length(pts.sel) > 0){
    cp.final <- pts.sel;
    cp.final <- cp.final[which(cp.final > sum(blocks.size[1:3]))];
    cp.final <- cp.final[which(cp.final < (TT-sum(blocks.size[(length(blocks.size)-2):length(blocks.size)])))];
    cp.final <- sort(cp.final);

    if(length(cp.final) >=2){
      gap.temp <- sapply(2:length(cp.final), function(jjj) cp.final[jjj]-cp.final[jjj-1])
    }



    if(length(cp.final) > 5){
      if(length(unique(gap.temp)) > 1  ){

        cl <- fviz_nbclust(matrix(cp.final, length(cp.final), 1), kmeans, nstart = 25,  method = "gap_stat",
                           k.max = min(50, length(cp.final)-1), nboot = 100)+
          labs(subtitle = "Gap statistic method")
        cl.data <- cl$data;
        gap <- cl.data$gap;
        cl.number <- which.max(gap)
        gap.order <- order(gap, decreasing = TRUE)

        if(median(blocks.size) <= sqrt(TT)/4){
          cnst <- 2*4+1
        }else if(median(blocks.size) <= 2*sqrt(TT)/4){
          cnst <- 2*3+1
        }else{
          cnst <- 2*2+1
        }

        flag = TRUE
        idx <- 1
        while(flag & idx <= length(gap.order)){
          cl.number <- gap.order[idx]
          if(cl.number > 1){
            cluster.pos <- sort(order(gap.temp, decreasing = TRUE)[1:(cl.number-1)])
            cluster.pos <- c(0, cluster.pos, length(cp.final))
          }else{
            cluster.pos <- c(0, length(cp.final))
          }

          wide <- 0
          for (i in c(1:cl.number)) {
            pts.i <- cp.final[(cluster.pos[i]+1): cluster.pos[i+1]]
            block_avg <- mean(blocks.size[match(pts.i, blocks)])
            if ( max(pts.i) - min(pts.i) > cnst*block_avg ){
              wide <- 1
              break
            }
          }
          if (wide){
            idx <- idx + 1
          }else{
            idx <- idx
            flag <- FALSE
          }
        }
        if(flag){
          cl.number <- length(gap.temp) + 1
        }

        if(median(blocks.size) <= sqrt(TT)/4){
          cnst2 <- 7
        }else if(median(blocks.size) <= sqrt(TT)/2){
          cnst2 <- 5
        }else{
          cnst2 <- 3
        }

        if(cl.number > 1){
          cl.number <- sum(sort(gap.temp, decreasing = TRUE)[1:(cl.number - 1)] > cnst2*median(blocks.size)) + 1
        }
      }else if(unique(gap.temp) == median(blocks.size) ){
        cl.number <- 1
      }else{
        cl.number <- length(gap.temp) + 1
      }

      if(cl.number > 1){
        cluster.pos <- sort(order(gap.temp, decreasing = TRUE)[1:(cl.number-1)])
        cluster.pos <- c(0, cluster.pos, length(cp.final))
      }else{
        cluster.pos <- c(0, length(cp.final))
      }

      pts.list <-  vector("list", cl.number);
      for (i in c(1:cl.number)) {
        pts.i <- cp.final[(cluster.pos[i]+1): cluster.pos[i+1]]
        pts.list[[i]] <- pts.i
      }
    }
    if(length(cp.final) <= 5 & length(cp.final) > 1 ){
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))

      for (i in 2:length(cp.final)){
        if (cp.final[i]-cp.final[i-1]<= max(3*mean(blocks.size))){
          cl.number <-cl.number-1
          loc.new[i] <-loc.new[i-1]
        }else{
          loc.new[i] <- i
        }
      }

      pts.list <-  vector("list",cl.number);
      loc.new.unique <- unique(loc.new)
      for (i in 1:length(loc.new.unique)) {
        pts.i <- cp.final[which(loc.new==loc.new.unique[i])]
        pts.list[[i]] <- pts.i
      }
    }
    if(length(cp.final) == 1 ){
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))
      pts.list <-  vector("list",cl.number);
      for (i in unique(loc.new)) {
        pts.i <- cp.final[which(loc.new==i)]
        pts.list[[i]] <- pts.i
      }
    }
    if(length(cp.final) == 0 ){
      pts.list <-  vector("list", 0);
    }
  }
  phi.par.sum <- vector("list",n.new);
  phi.par.sum[[1]] <- phi.hat.full[,1:(p.x)];
  for(i in 2:n.new){
    phi.par.sum[[i]] <- phi.par.sum[[i-1]] + phi.hat.full[,((i-1)*p.x+1):(i*p.x)];
  }
  return(list(jumps.l2 = jumps.l2, jumps.l1 = jumps.l1, pts.list = pts.list, beta.full = phi.par.sum))

}


#' @title second.step
#' @description Reimputate the missing values and perform the exhaustive search to "thin out" redundant break points.
#'
#' @param data_y input data matrix, with each column representing the time series component
#' @param data_x input data matrix
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cp.first the selected break points after the first step
#' @param beta.est the estiamted parameters by block fused lasso
#' @param blocks the blocks
#' @param data_y_miss the data y matrix before the first imputation
#' @return A list object, which contains the followings
#' \describe{
#'   \item{cp.final}{a set of selected break point after the exhaustive search step}
#'   \item{beta.hat.list}{the estimated coefficient matrix for each segmentation}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
second.step <- function(data_y, data_x, max.iteration = max.iteration, tol = tol,
         cp.first, beta.est, blocks, data_y_miss){

  TT <- length(data_y[,1])
  p.y <- length(data_y[1,])
  p.x <- length(data_x[1,])

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );

  cp.search <- cp.first;
  cl.number <- length(cp.first)

  cp.list <- vector("list", cl.number + 2);
  cp.list[[1]] <- c(1);
  cp.list[[cl.number+2]] <- c(TT+1);

  cp.index.list <- vector("list", cl.number + 2);
  cp.index.list[[1]] <- c(1);
  cp.index.list[[cl.number+2]] <- c(n.new+1);

  for (i in 1:cl.number) {
    cp.list[[i+1]] <- cp.first[[i]]
    cp.index.list[[i+1]] <- match(cp.first[[i]], blocks)
  }

  cp.search <- rep(0, cl.number);
  cp.list.full <- cp.list

  cp.imputation <- rep(0, 2*cl.number+2)
  cp.imputation[1] <- 1
  cp.imputation[2*cl.number+2] <- TT+1
  for(i in 1 : (cl.number)){
    cp.imputation[2*i] = cp.list[[i+1]][1]
    cp.imputation[2*i + 1] = cp.list[[i+1]][length(cp.list[[i+1]])]-1
  }
  data_y2 = imputation2(data_y_miss, cp.imputation)
  for(i in 1:p.y){
    for(j in 1:TT){
      if(is.na(data_y2[j,i])){
        data_y2[j,i] <- data_y[j,i]
      }
    }
  }

  beta.hat.list <- vector("list", cl.number+1)
  idx1.lb <- 1
  idx1.up <- min(cp.index.list[[2]]) - 2
  temp.beta <- vector("list", idx1.up - idx1.lb+1)
  for(i.idx in idx1.lb:idx1.up){
    temp.beta[[i.idx - idx1.lb + 1]] <- beta.est[[i.idx]]
  }
  beta.hat.list[[1]] <- Reduce("+", temp.beta)/length(temp.beta)


  idx2.lb <- max(cp.index.list[[cl.number+1]]) + 1
  idx2.up <- length(blocks) - 1
  temp.beta2 <- vector("list", idx2.up - idx2.lb+1)
  for(i.idx in idx2.lb:idx2.up){
    temp.beta2[[i.idx - idx2.lb + 1]] <- beta.est[[i.idx]]
  }
  beta.hat.list[[cl.number+1]] <- Reduce("+", temp.beta2)/length(temp.beta2)


  if(cl.number>1){
    for(i in 2:(cl.number)){
      idx1.lb <- max(cp.index.list[[i]]) + 1
      idx1.up <- min(cp.index.list[[i+1]]) - 2
      temp.beta <- vector("list", idx1.up - idx1.lb+1)
      for(i.idx in idx1.lb:idx1.up){
        temp.beta[[i.idx - idx1.lb + 1]] <- beta.est[[i.idx]]
      }
      beta.hat.list[[i]] <- Reduce("+", temp.beta)/length(temp.beta)
    }
  }

  for(i in 1:(cl.number)){



    if(length(cp.list[[i+1]]) >1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] + 1)  :  (cp.list[[i+1]][length(cp.list[[i+1]])]-1  ) )
    }
    if(length(cp.list[[i+1]]) == 1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] -  (blocks.size[cp.index.list[[i+1]][1] ]) + 1) :  (cp.list[[i+1]][length(cp.list[[i+1]])] +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1   ) )
    }


    #compare the SSE of first num and last num
    num  <- cp.list.full[[i+1]][1]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    #idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2)
    #beta.hat <- beta.est[[idx.1]]
    beta.hat <- beta.hat.list[[i]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.1 == 1){
      temp.1 <- sum( (data_y2[lb.1:ub.1,]-forecast)^2)
    }else{
      temp.1 <- sum( (t(data_y2[lb.1:ub.1,])-forecast)^2 )
    }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    #idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    #beta.hat <- beta.est[[idx.2]]
    beta.hat <- beta.hat.list[[i+1]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.2 == 1){
      temp.2 <- sum( ( data_y2[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data_y2[lb.2:ub.2,])-forecast)^2 );
    }
    sse1 <- temp.1 + temp.2;


    num  <- cp.list.full[[i+1]][length(cp.list.full[[i+1]])]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    #idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    #beta.hat <- beta.est[[idx.1]]
    beta.hat <- beta.hat.list[[i]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.1 == 1){
      temp.1 <- sum( (data_y2[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data_y2[lb.1:ub.1,])-forecast)^2 );
    }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    #idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    #beta.hat <- beta.est[[idx.2]]
    beta.hat <- beta.hat.list[[i+1]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )

    if(len.2 == 1){
      temp.2 <- sum( ( data_y2[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data_y2[lb.2:ub.2,])-forecast)^2 );
    }
    sse2 <- temp.1 + temp.2;



    if(sse1 <= sse2){
      sse.full <- 0;
      ii <- 0
      for(num in cp.list.full[[i+1]]  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        beta.hat <- beta.hat.list[[i]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
        if(len.1 == 1){
          temp.1 <- sum( (data_y2[lb.1:ub.1,]-forecast)^2 );

        }else{
          temp.1 <- sum( (t(data_y2[lb.1:ub.1,])-forecast)^2 );
        }


        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        #idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        #beta.hat <- beta.est[[idx.1]]
        beta.hat <- beta.hat.list[[i+1]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x),matrix(beta.hat,ncol = p.x),jjj, p.x, p.y)  )
        if(len.2 == 1){
          temp.2 <- sum( (data_y2[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data_y2[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full,0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][min(which(sse.full == min(sse.full)))];

    }
    if(sse1 > sse2){
      sse.full <- 0;
      ii <- 0
      for(num in rev(cp.list.full[[i+1]])  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        #idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        #beta.hat <- beta.est[[idx.1]]
        beta.hat <- beta.hat.list[[i]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
        if(len.1 == 1){
          temp.1 <- sum( (data_y2[lb.1:ub.1,]-forecast)^2 );

        }else{
          temp.1 <- sum( (t(data_y2[lb.1:ub.1,])-forecast)^2 );
        }


        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        #idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        #beta.hat <- beta.est[[idx.1]]
        beta.hat <- beta.hat.list[[i+1]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x),matrix(beta.hat,ncol = p.x),jjj, p.x, p.y)  )
        if(len.2 == 1){
          temp.2 <- sum( (data_y2[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data_y2[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full,0.30)){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][length(cp.list.full[[i+1]]) + 1 - min(which(sse.full == min(sse.full)))];

    }

  }
  idx <- floor((min(cp.index.list[[cl.number+1+1]]) + max(cp.index.list[[cl.number+1]]))/2);
  beta.hat.list[[cl.number+1]] <- beta.est[[idx]]
  return(list(cp.final = cp.search, beta.hat.list = beta.hat.list))
}


#' @title BTIE
#' @description Perform the BTIE algorithm to detect the structural breaks
#' in large scale high-dimensional mean shift models.
#'
#' @param data_y input data matrix (response), with each column representing the time series component
#' @param lambda.1.cv tuning parmaeter lambda_1 for fused lasso
#' @param lambda.2.cv tuning parmaeter lambda_2 for fused lasso
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param block.size the block size
#' @param refit logical; if TRUE, refit the model, if FALSE, use BIC to find a thresholding value and then output the parameter estimates without refitting. Default is FALSE.
#' @param optimal.block logical; if TRUE, grid search to find optimal block size, if FALSE, directly use the default block size. Default is TRUE.
#' @param optimal.gamma.val hyperparameter for optimal block size, if optimal.blocks == TRUE. Default is 1.5.
#' @param block.range the search domain for optimal block size.
#' @importFrom Rcpp sourceCpp
#' @importFrom glmnet cv.glmnet
#' @importFrom sparsevar fitVAR
#' @export
#' @examples
#' set.seed(1)
#' n <- 1000;
#' p <- 50;
#' brk <-  c(333, 666, n+1)
#' m <- length(brk)
#' d <- 5
#' constant.full <- constant_generation(n, p, d, 50, brk)
#' e.sigma <- as.matrix(1*diag(p))
#' data_y <- data_generation(n = n, mu = constant.full, sigma = e.sigma, brk = brk)
#' data_y <- as.matrix(data_y, ncol = p.y)
#' data_y_miss <- MCAR(data_y, 0.3)
#' temp <- BTIE(data_y_miss, optimal.block = FALSE, block.size = 30)
#' temp$cp.final
#' @return A list object, which contains the followings
BTIE <- function(data_y, lambda.1.cv = NULL, lambda.2.cv = NULL,
                max.iteration = 100, tol = 10^(-2), block.size = NULL,
                refit = FALSE, optimal.block = TRUE, optimal.gamma.val = 1.5,
                block.range = NULL){
  TT  <- length(data_y[, 1]);
  p.y <- length(data_y[1, ]);
  data_x <- matrix(1, nrow = TT)
  p.x <- 1
  data_y_miss = data_y

  if( !is.null(block.size) && optimal.block == TRUE){
    stop("Set optimal.block to be FALSE if entered block")
  }

  if( is.null(block.size) && optimal.block == FALSE){
    stop("Set the block size if optimal.block is FALSE")
  }

  if(optimal.block == TRUE){
    if(sqrt(TT) > p.x*p.y){
      b_n.max <- ceiling(min(sqrt(TT), TT/20))
    }else{
      b_n.max <- ceiling(min(sqrt(TT)*log(p.x*p.y), TT/20))
    }
    b_n.min <- floor(min(log(TT*p.y), TT/20 ))
    b_n.range <- round(seq(b_n.min, b_n.max, length.out = 5))
    b_n.range <- unique(b_n.range)
    b_n.range <- b_n.range[b_n.range > 1]

    n.method <- length(b_n.range) #number of methods
    temp.full <- vector("list", n.method);
    pts.full <- vector("list", n.method);
    for(j.2 in 1:n.method){
      block.size = b_n.range[j.2]
      data_y <- imputation(data_y, block.size)
      b_n_bound = 2*block.size  #block size for boundary
      blocks <- c(seq(1, b_n_bound*2+1 , b_n_bound),
                  seq(b_n_bound*2+block.size+1, TT+1-2*b_n_bound, block.size),
                  seq(TT+1-b_n_bound, TT+1,  b_n_bound));
      # blocks <- seq(1, TT + 1, block.size);
      if(blocks[length(blocks)] < TT+1){
        blocks <- c(blocks[-length(blocks)], TT + 1)
      }
      n.new <- length(blocks) - 1;
      blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj + 1] - blocks[jjj]  );

      bbb <- floor(n.new/4);
      aaa <- 4;
      cv.index <- seq(aaa,n.new,floor(n.new/bbb));
      cv.l <- length(cv.index);


      ############# Tuning parameter ################
      if(is.null(lambda.1.cv)){
        lambda.1.max <- lambda_warm_up(data_y, data_x, blocks, cv.index)$lambda_1_max
        if(blocks[2] <= (p.x + p.y) ){
          epsilon <-  10^(-3)
        }
        if(blocks[2] >= (p.x + p.y) ){
          epsilon <-  10^(-4)
        }
        nlam <- 10
        lambda.1.min <-  lambda.1.max*epsilon
        delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
        lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
      }

      if(is.null(lambda.2.cv)){
        lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
      }
      ########## first step       #######################
      temp.first <- first.step(data_y, data_x, lambda.1.cv, lambda.2.cv,
                                          max.iteration = max.iteration, tol = tol, blocks,
                                          cv.index)
      cp.first <- temp.first$pts.list;
      cl.number <- length(cp.first);
      beta.est <- temp.first$beta.full
      temp.first.all<- temp.first
      jumps.l2 <- temp.first$jumps.l2
      ########## second step       #######################
      if(length(cp.first) > 0){
        temp.second<- second.step(data_y, data_x, max.iteration = max.iteration,
                                  tol = tol, cp.first, beta.est, blocks, data_y_miss = data_y_miss)
        temp.second.all<- temp.second
        cp.final<- temp.second$cp.final;
        beta.hat.list <- temp.second$beta.hat.list
      }else{
        cp.final <- c()
        beta.hat.list <- vector("list", 1)
        beta.hat.list[[1]] <- beta.est[[floor(n.new/2)]]
      }

      if(refit){
        cp.full <- c(1, cp.final, TT+1)
        m <- length(cp.final) + 1
        temp_beta <-  list(m)
        for(i in 1:m){
          # data_x_temp <- matrix(1, nrow = TT)
          data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
          temp_beta[[i]] <- apply(data_y_temp, 2, mean)
        }

        beta.final.temp <- matrix(c(do.call("cbind", temp_beta)), nrow = 1)
        lambda.val.best <- BIC_threshold(beta.final.temp, p.y,
                                         length(c(cp.final, TT+1)), c(cp.final, TT+1), data_y , b_n = block.size)

        for(j in 1:length(c(cp.final, TT+1))){
          temp_beta[[j]][abs(temp_beta[[j]]) < lambda.val.best[j]] <- 0
        }
        beta.final <- temp_beta

      }else{
        beta.final.temp <- matrix(c(do.call("cbind", beta.hat.list)), nrow= 1)
        lambda.val.best <- BIC_threshold(beta.final.temp, p.y, length(cp.final) + 1, c(cp.final, TT+1),
                                         data_y,  b_n = block.size, nlam = 20)
        temp <- beta.hat.list
        for(j in 1:(length(cp.final) + 1) ){
          temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
        }
        cp.full <- c(1, cp.final, TT+1)
        m = NROW(cp.final) + 1
        temp_beta <- list(m)
        for(i in 1:m){
          data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
          temp_beta[[i]] <- apply(data_y_temp, 2, mean)
        }
        beta.final <- temp_beta
      }
      temp.full[[j.2]] <- beta.final
      if(length(cp.final) > 0){
        pts.full[[j.2]] <- cp.final
      }


      temp.full[[j.2]] <- beta.final
      if(length(cp.final) > 0){
        pts.full[[j.2]] <- cp.final
      }
    }

    HBIC.full <- c()
    BIC.full <- c()
    for(j.2 in 1:n.method){
      brk.temp <- pts.full[[j.2]]
      brk.full <- c(1, brk.temp, TT+1)
      m.hat <- length(brk.full) - 1
      k <- p.y
      beta.final <- matrix(c(do.call("cbind", temp.full[[j.2]])), nrow = 1)
      HBIC.res <- c()
      residual.full <- c()
      for(i in 1:m.hat){
        beta.temp <- beta.final[((i-1)*k+1):(i*k)]
        data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
        data_x.temp <- matrix(1, nrow = brk.full[i+1] - brk.full[i])
        data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, nrow = 1)
        residual.temp <- data_y.temp - data_y.est
        residual.full <- rbind(residual.full, residual.temp)
      }

      phi = do.call("cbind", temp.full[[j.2]])
      residual <- t(residual.full)
      p.y <- length(phi[, 1]);
      p.x <- length(phi[1, ]);
      T.new <- length(residual[1, ]);
      count <- 0;
      for (i in 1:p.y){
        for (j in 1:p.x){
          if(phi[i,j] != 0){
            count <- count + 1;
          }
        }
      }

      sigma.hat <- 0*diag(p.y);
      for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[, i]%*%t(residual[, i]);  }
      sigma.hat <- (1/(T.new))*sigma.hat;
      ee.temp <- min(eigen(sigma.hat)$values);
      if(ee.temp <= 10^(-8)){
        sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp) + 10^(-3))*diag(p.y);
      }
      log.det <- log(det(sigma.hat));
      HBIC.res <- log.det*TT + 2*optimal.gamma.val*log(p.y)*count
      HBIC.full <- c(HBIC.full, sum(HBIC.res))
    }
    pts.res.HBIC<- pts.full[[which.min(HBIC.full)]]
    bn.res.HBIC <- b_n.range[which.min(HBIC.full)]

    return(list(cp.final = pts.res.HBIC, beta.final = temp.full[[which.min(HBIC.full)]], bn.optimal = bn.res.HBIC,
                bn.range = b_n.range, HBIC.full = HBIC.full, pts.full = pts.full))



  }else{
    b_n_bound = 2*block.size  #block size for boundary
    blocks <- c(seq(1, b_n_bound*2+1 , b_n_bound),
                seq(b_n_bound*2+block.size+1, TT+1-2*b_n_bound, block.size),
                seq(TT+1-b_n_bound, TT+1,  b_n_bound));

    if(blocks[length(blocks)] < TT+1){
      blocks <- c(blocks[-length(blocks)], TT + 1)
    }

    n.new <- length(blocks) - 1;
    blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj + 1] - blocks[jjj]  );
    if(is.null(block.size)){
      block.size <- median(blocks.size)
    }
    bbb <- floor(n.new/4);
    aaa <- 4;
    cv.index <- seq(aaa,n.new,floor(n.new/bbb));
    cv.l <- length(cv.index);

    p.y <- length(data_y[1, ]);
    data_x <- matrix(1, nrow = TT)
    data_y <- imputation(data_y, block.size)
    p.x <- 1
    if(is.null(lambda.1.cv)){
      lambda.1.max <- lambda_warm_up(data_y, data_x, blocks, cv.index)$lambda_1_max
      if(blocks[2] <= (p.x + p.y) ){
        epsilon <-  10^(-3)
      }
      if(blocks[2] >= (p.x + p.y) ){
        epsilon <-  10^(-4)
      }
      nlam <- 10
      lambda.1.min <-  lambda.1.max*epsilon
      delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
      lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
    }

    if(is.null(lambda.2.cv)){
      lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
    }

    temp.first <- first.step(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration,
                             tol = tol, blocks , cv.index)
    cp.first <- temp.first$pts.list;
    cl.number <- length(cp.first);
    beta.est <- temp.first$beta.full
    temp.first.all<- temp.first
    jumps.l2 <- temp.first$jumps.l2

    if(length(cp.first) > 0){
      temp.second<- second.step(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first,
                                beta.est, blocks, data_y_miss = data_y_miss)
      temp.second.all<- temp.second
      cp.final<- temp.second$cp.final;
      beta.hat.list <- temp.second$beta.hat.list
    }else{
      cp.final <- c()
      beta.hat.list <- vector("list", 1)
      beta.hat.list[[1]] <- beta.est[[floor(n.new/2)]]
    }

    if(refit){
      cp.full <- c(1, cp.final, TT+1)
      m <- length(cp.final) + 1
      temp_beta <-  list(m)
      for(i in 1:m){
        data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
        temp_beta[[i]] <- apply(data_y_temp, 2, mean)
      }

      beta.final.temp <- matrix(c(do.call("cbind", temp_beta)), nrow = 1)
      lambda.val.best <- BIC_threshold(beta.final.temp, p.y,
                                       length(c(cp.final, TT+1)), c(cp.final, TT+1), data_y , b_n = block.size)

      for(j in 1:length(c(cp.final, TT+1))){
        temp_beta[[j]][abs(temp_beta[[j]]) < lambda.val.best[j]] <- 0
      }
      beta.final <- temp_beta

    }else{
      beta.final.temp <- matrix(c(do.call("cbind", beta.hat.list)), nrow= 1)
      lambda.val.best <- BIC_threshold(beta.final.temp, p.y, length(cp.final) + 1, c(cp.final, TT+1),
                                       data_y,  b_n = block.size, nlam = 20)
      cp.full <- c(1, cp.final, TT+1)
      m = NROW(cp.final) + 1
      temp_beta <- list(m)
      for(i in 1:m){
        data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
        temp_beta[[i]] <- apply(data_y_temp, 2, mean)
      }
      beta.final <- temp_beta
    }
    return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
                beta.est = beta.est, beta.final = beta.final, jumps = jumps.l2))
  }
}

