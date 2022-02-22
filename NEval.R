insample_eval <- function(data, kernel_sigma, lm_Hh1, lm_Qh1, lm_Hh0, lm_Qh0, lm_Hq1, lm_Qq1, lm_Hq0, Nystroem = FALSE, appro_rate = 0.05){
  
  if (!Nystroem) {
    print("We do not use a low rank approximation.")
    n <- length(data$Y)
    ## Gaussian RBF kernel ############################################################
    h_arg <- cbind(data$W,data$X)
    q_arg <- cbind(data$Z,data$X)
    K_H <- rbfkernel(h_arg, sigma = kernel_sigma)
    K_Q <- rbfkernel(q_arg, sigma = kernel_sigma)
    ###################################################################################
    
    ## Optimization ###################################################################
    ident_n <- diag(rep(1,n))
    D1 <- diag(as.numeric(data$A),n,n)
    D0 <- diag(as.numeric(1-data$A),n,n)
    ## Optimization h ##
    Gamm <- 0.25*K_Q %*% ginv( K_Q/n + lm_Qh1/n^0.8*ident_n )
    alpha1h <- ginv(K_H %*% D1 %*% Gamm %*% D1 %*% K_H+n^2*lm_Hh1/n^0.8*K_H) %*%
      K_H %*% D1 %*% Gamm %*% (data$Y*data$A)
    Gamm <- 0.25*K_Q %*% ginv( K_Q/n + lm_Qh0/n^0.8*ident_n )
    alpha0h <- ginv(K_H %*% D0 %*% Gamm %*% D0 %*% K_H+n^2*lm_Hh0/n^0.8*K_H) %*%
      K_H %*% D0 %*% Gamm %*% (data$Y*(1-data$A))
    ## Optimization q ##
    Gamm <- 0.25*K_H %*% ginv( K_H/n + lm_Hq1/n^0.8*ident_n )
    alpha1q <- ginv(K_Q %*% D1 %*% Gamm %*% D1 %*% K_Q+n^2*lm_Qq1/n^0.8*K_Q) %*%
      K_Q %*% D1 %*% Gamm %*% matrix(1,n,1)
    Gamm <- 0.25*K_H %*% ginv( K_H/n + lm_Hq0/n^0.8*ident_n )
    alpha0q <- ginv(K_Q %*% D0 %*% Gamm %*% D0 %*% K_Q+n^2*lm_Qq0/n^0.8*K_Q) %*%
      K_Q %*% D0 %*% Gamm %*% matrix(1,n,1)
    ###################################################################################
    ## Estimate function h_hat ########################################################
    h_hat_1 <- function(w1,x1){
      W <- c(w1,x1)
      S <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W,1,length(W)))
      return(t(S) %*% alpha1h)
    }
    h_hat_0 <- function(w1,x1){
      W <- c(w1,x1)
      S <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W,1,length(W)))
      return(t(S) %*% alpha0h)
    }
    h_hat <- function(w1,a1,x1){
      return(a1*h_hat_1(w1,x1)+(1-a1)*h_hat_0(w1,x1))
    }
    ###################################################################################
    
    ## Estimate function q_hat ########################################################
    q_hat_1 <- function(z1,x1){
      Z <- c(z1,x1)
      S <- rbfkernel(q_arg, sigma = kernel_sigma,matrix(Z,1,length(Z)))
      return(t(S) %*% alpha1q)
    }
    q_hat_0 <- function(z1,x1){
      Z <- c(z1,x1)
      S <- rbfkernel(q_arg, sigma = kernel_sigma,matrix(Z,1,length(Z)))
      return(t(S) %*% alpha0q)
    }
    q_hat <- function(z1,a1,x1){
      return(a1*q_hat_1(z1,x1)+(1-a1)*q_hat_0(z1,x1))
    }
    ###################################################################################
    
  } else {
    print("We do use a low rank approximation.")
    
    n <- length(data$Y)
    r <- ceiling(n * appro_rate)
    S <- diag(1, n)[, sort(sample(1:n, r))]
    ## Gaussian RBF kernel ############################################################
    h_arg <- cbind(data$W,data$X)
    q_arg <- cbind(data$Z,data$X)
    K_H <- rbfkernel(h_arg, sigma = kernel_sigma)
    K_Q <- rbfkernel(q_arg, sigma = kernel_sigma)
    ###################################################################################
    
    ## Optimization ###################################################################
    ident_n <- diag(1, n)
    ident_r <- diag(1, r)
    D1 <- diag(as.numeric(data$A),n,n)
    D0 <- diag(as.numeric(1-data$A),n,n)
    ## Optimization h ##
    sqrtmtSK_QS <- sqrtm(ginv(t(S) %*% K_Q %*% S))
    D <- K_Q %*% S %*% sqrtmtSK_QS
    Gamm <- 0.25*D %*% ginv( t(D) %*% D/n + lm_Qh1/n^0.8*ident_r ) %*% t(D)
    sqrtmtSK_HS <- sqrtm(ginv(t(S) %*% K_H %*% S))
    V <- K_H %*% S %*% sqrtmtSK_HS
    alpha1h <- ginv(t(V) %*% D1 %*% Gamm %*% D1 %*% V+n^2*lm_Hh1/n^0.8*ident_r) %*%
      t(V) %*% D1 %*% Gamm %*% (data$Y*data$A)
    Gamm <- 0.25*D %*% ginv( t(D) %*% D/n + lm_Qh0/n^0.8*ident_r ) %*% t(D)
    alpha0h <- ginv(t(V) %*% D0 %*% Gamm %*% D0 %*% V+n^2*lm_Hh0/n^0.8*ident_r) %*%
      t(V) %*% D0 %*% Gamm %*% (data$Y*(1-data$A))
    ## Optimization q ##
    Gamm <- 0.25*V %*% ginv( t(V) %*% V/n + lm_Hq1/n^0.8*ident_r ) %*% t(V)
    alpha1q <- ginv(t(D) %*% D1 %*% Gamm %*% D1 %*% D+n^2*lm_Qq1/n^0.8*ident_r) %*%
      t(D) %*% D1 %*% Gamm %*% matrix(1,n,1)
    Gamm <- 0.25*V %*% ginv( t(V) %*% V/n + lm_Hq0/n^0.8*ident_r ) %*% t(V)
    alpha0q <- ginv(t(D) %*% D0 %*% Gamm %*% D0 %*% D+n^2*lm_Qq0/n^0.8*ident_r) %*%
      t(D) %*% D0 %*% Gamm %*% matrix(1,n,1)
    ###################################################################################
    ## Estimate function h_hat ########################################################
    h_hat_1 <- function(w1,x1){
      W <- c(w1,x1)
      K <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W,1,length(W)))
      return(t(K) %*% S %*% sqrtmtSK_HS %*% alpha1h)
    }
    h_hat_0 <- function(w1,x1){
      W <- c(w1,x1)
      K <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W,1,length(W)))
      return(t(K) %*% S %*% sqrtmtSK_HS %*% alpha0h)
    }
    h_hat <- function(w1,a1,x1){
      return(a1*h_hat_1(w1,x1)+(1-a1)*h_hat_0(w1,x1))
    }
    ###################################################################################
    
    ## Estimate function q_hat ########################################################
    q_hat_1 <- function(z1,x1){
      Z <- c(z1,x1)
      K <- rbfkernel(q_arg, sigma = kernel_sigma, matrix(Z,1,length(Z)))
      return(t(K) %*% S %*% sqrtmtSK_QS %*% alpha1q)
    }
    q_hat_0 <- function(z1,x1){
      Z <- c(z1,x1)
      K <- rbfkernel(q_arg, sigma = kernel_sigma, matrix(Z,1,length(Z)))
      return(t(K) %*% S %*% sqrtmtSK_QS %*% alpha0q)
    }
    q_hat <- function(z1,a1,x1){
      return(a1*q_hat_1(z1,x1)+(1-a1)*q_hat_0(z1,x1))
    }
    ###################################################################################
    
  }
  
  
  ## Evaluation #####################################################################
  gt <- data$truth
  ## Evaluation OR ##################################################################
  est <- 0
  for (j in seq(n)){
    est <- est+(1/n)*(h_hat(data$W[j,],1,data$X[j,])-h_hat(data$W[j,],0,data$X[j,]))
  }
  OR_est <- est-gt
  
  ## Evaluation IPW #################################################################
  est <- 0
  for (j in seq(n)){
    est <- est+(1/n)*(-1)^(1-data$A[j,])*q_hat(data$Z[j,],data$A[j,],data$X[j,])*
      data$Y[j,]
  }
  IPW_est <- est-gt
  
  ## Evaluation DR ##################################################################
  est <- 0
  for (j in seq(n)){
    est <- est+(1/n)*(-1)^(1-data$A[j,])*q_hat(data$Z[j,],data$A[j,],data$X[j,])*
      (data$Y[j,]-h_hat(data$W[j,],data$A[j,],data$X[j,]))+
      (1/n)*(h_hat(data$W[j,],1,data$X[j,])-h_hat(data$W[j,],0,data$X[j,]))
  }
  DR_est <- est-gt
  ###################################################################################
  
  bias <- c(OR_est, IPW_est, DR_est)
  return(bias)
   
}


CV_eval <- function(data, kernel_sigma_list, theta = 0, lm_Hh1_list, lm_Qh1_list, lm_Hh0_list, lm_Qh0_list, 
                    lm_Hq1_list, lm_Qq1_list, lm_Hq0_list, lm_Qq0_list, CV_K_fold, Nystroem = FALSE, appro_rate = 0.05){
  
  n_sample <- length(data$Y)
  working_data <- with(data, list(X = X, Z = Z, W = W, A = A, Y = Y))
  risk_sum <- c()
  hyperparameters_mat <- c()
  kernel_sigma_num <- length(kernel_sigma_list)
  for (kernel_rep in 1:kernel_sigma_num) {
    random_sample <- createFolds(1:n_sample, k = CV_K_fold)
    riskh1_result <- matrix(0, nrow = length(lm_Hh1_list), ncol = length(lm_Qh1_list))
    riskh0_result <- matrix(0, nrow = length(lm_Hh0_list), ncol = length(lm_Qh0_list))
    riskq1_result <- matrix(0, nrow = length(lm_Hq1_list), ncol = length(lm_Qq1_list))
    riskq0_result <- matrix(0, nrow = length(lm_Hq0_list), ncol = length(lm_Qq0_list))
    kernel_sigma <- kernel_sigma_list[kernel_rep]
    
    if (!Nystroem) {
      for (cv_rep in 1:CV_K_fold) {
        tr_data <- sapply(working_data, function(x) matrix(x[-random_sample[[cv_rep]], ], ncol = ncol(x)))
        te_data <- sapply(working_data, function(x) matrix(x[random_sample[[cv_rep]], ], ncol = ncol(x)))
        n_tr <- n_sample - length(random_sample[[cv_rep]])
        n_te <- length(random_sample[[cv_rep]])
        
        ## Gaussian RBF kernel ############################################################
        h_arg_tr <- cbind(tr_data$W, tr_data$X)
        q_arg_tr <- cbind(tr_data$Z, tr_data$X)
        K_H_tr <- rbfkernel(h_arg_tr, sigma = kernel_sigma)
        K_Q_tr <- rbfkernel(q_arg_tr, sigma = kernel_sigma)
        ###################################################################################
        
        ## Optimization ###################################################################
        ident_n_tr <- diag(rep(1, n_tr))
        D1_tr <- diag(as.numeric(tr_data$A))
        D0_tr <- diag(as.numeric(1 - tr_data$A))
    
        ###################################################################################
        
    
        ## Gaussian RBF kernel ############################################################
        h_arg_te <- cbind(te_data$W, te_data$X)
        q_arg_te <- cbind(te_data$Z, te_data$X)
        K_H_tr_te <- rbfkernel(X = h_arg_te, sigma = kernel_sigma, Y = h_arg_tr)
        K_Q_tr_te <- rbfkernel(X = q_arg_te, sigma = kernel_sigma, Y = q_arg_tr)
        K_H_te <- rbfkernel(h_arg_te, sigma = kernel_sigma)
        K_Q_te <- rbfkernel(q_arg_te, sigma = kernel_sigma)
        
        ident_n_te <- diag(rep(1, n_te))
        D1_te <- diag(as.numeric(te_data$A), n_te, n_te)
        D0_te <- diag(as.numeric(1 - te_data$A), n_te, n_te)
        
        Theta <- 0.25 * K_Q_te %*% ginv(K_Q_te / n_te + theta * ident_n_te)
        ## find the best combinations to minimize the inside maximization of h1
        for (lm_Hh1_index in 1:length(lm_Hh1_list)) {
          for (lm_Qh1_index in 1:length(lm_Qh1_list)) {
            lm_Hh1 <- lm_Hh1_list[lm_Hh1_index]
            lm_Qh1 <- lm_Qh1_list[lm_Qh1_index]
            Gamm <- 0.25 * K_Q_tr %*% ginv(K_Q_tr / n_tr + lm_Qh1 / n_tr^0.7 * ident_n_tr)
            alpha1h <- ginv(K_H_tr %*% D1_tr %*% Gamm %*% D1_tr %*% K_H_tr + n_tr^2 * lm_Hh1 / n_tr^0.7 * K_H_tr) %*%
              K_H_tr %*% D1_tr %*% Gamm %*% (tr_data$Y * tr_data$A)
            xi_n <- (D1_te %*% K_H_tr_te %*% alpha1h + te_data$Y * te_data$A) / n_te
            riskh1 <- t(xi_n) %*% Theta %*% xi_n
            riskh1_result[lm_Hh1_index, lm_Qh1_index] <- riskh1_result[lm_Hh1_index, lm_Qh1_index] + riskh1
          }
        }
        
        ## find the best combinations to minimize the inside maximization of h0
        for (lm_Hh0_index in 1:length(lm_Hh0_list)) {
          for (lm_Qh0_index in 1:length(lm_Qh0_list)) {
            lm_Hh0 <- lm_Hh0_list[lm_Hh0_index]
            lm_Qh0 <- lm_Qh0_list[lm_Qh0_index]
            Gamm <- 0.25 * K_Q_tr %*% ginv(K_Q_tr / n_tr + lm_Qh0 / n_tr^0.7 * ident_n_tr)
            alpha0h <- ginv(K_H_tr %*% D0_tr %*% Gamm %*% D0_tr %*% K_H_tr + n_tr^2 * lm_Hh0 / n_tr^0.7 * K_H_tr) %*%
              K_H_tr %*% D0_tr %*% Gamm %*% (tr_data$Y * (1 - tr_data$A))
            xi_n <- (D0_te %*% K_H_tr_te %*% alpha0h + te_data$Y * (1 - te_data$A)) / n_te
            riskh0 <- t(xi_n) %*% Theta %*% xi_n
            riskh0_result[lm_Hh0_index, lm_Qh0_index] <- riskh0_result[lm_Hh0_index, lm_Qh0_index] + riskh0
          }
        }
    
        Theta <- 0.25 * K_H_te %*% ginv(K_H_te / n_te + theta * ident_n_te)
        ## find the best combinations to minimize the inside maximization of q1
        for (lm_Hq1_index in 1:length(lm_Hq1_list)) {
          for (lm_Qq1_index in 1:length(lm_Qq1_list)) {
            lm_Hq1 <- lm_Hq1_list[lm_Hq1_index]
            lm_Qq1 <- lm_Qq1_list[lm_Qq1_index]
            Gamm <- 0.25 * K_H_tr %*% ginv(K_H_tr / n_tr + lm_Hq1 / n_tr^0.7 * ident_n_tr)
            alpha1q <- ginv(K_Q_tr %*% D1_tr %*% Gamm %*% D1_tr %*% K_Q_tr + n_tr^2 * lm_Qq1 / n_tr^0.7 * K_Q_tr) %*%
              K_Q_tr %*% D1_tr %*% Gamm %*% matrix(1, n_tr, 1)
            xi_n <- (D1_te %*% K_Q_tr_te %*% alpha1q + 1) / n_te
            riskq1 <- t(xi_n) %*% Theta %*% xi_n
            riskq1_result[lm_Hq1_index, lm_Qq1_index] <- riskq1_result[lm_Hq1_index, lm_Qq1_index] + riskq1
            
          }
        }
        
        ## find the best combinations to minimize the inside maximization of q0
        for (lm_Hq0_index in 1:length(lm_Hq0_list)) {
          for (lm_Qq0_index in 1:length(lm_Qq0_list)) {
            lm_Hq0 <- lm_Hq0_list[lm_Hq0_index]
            lm_Qq0 <- lm_Qq0_list[lm_Qq0_index]
            Gamm <- 0.25 * K_H_tr %*% ginv(K_H_tr / n_tr + lm_Hq0 * ident_n_tr)
            alpha0q <- ginv(K_Q_tr %*% D0_tr %*% Gamm %*% D0_tr %*% K_Q_tr + n_tr^2 * lm_Qq0 / n_tr^0.7 * K_Q_tr) %*%
              K_Q_tr %*% D0_tr %*% Gamm %*% matrix(1, n_tr, 1)
            xi_n <- (D0_te %*% K_Q_tr_te %*% alpha0q + 1) / n_te
            riskq0 <- t(xi_n) %*% Theta %*% xi_n
            riskq0_result[lm_Hq0_index, lm_Qq0_index] <- riskq0_result[lm_Hq0_index, lm_Qq0_index] + riskq0
            
          }
        }
      }
    } else {
      # print("We do use a low rank approximation.")
      
      for (cv_rep in 1:CV_K_fold) {
        tr_data <- sapply(working_data, function(x) matrix(x[-random_sample[[cv_rep]], ], ncol = ncol(x)))
        te_data <- sapply(working_data, function(x) matrix(x[random_sample[[cv_rep]], ], ncol = ncol(x)))
        n_tr <- n_sample - length(random_sample[[cv_rep]])
        n_te <- length(random_sample[[cv_rep]])
        r_tr <- max(10, ceiling(n_tr * appro_rate))
        r_te <- max(10, ceiling(n_te * appro_rate))
        S_tr <- diag(1, n_tr)[, sort(sample(1:n_tr, r_tr))]
        S_te <- diag(1, n_te)[, sort(sample(1:n_te, r_te))]
        
        ## Gaussian RBF kernel ############################################################
        h_arg_tr <- cbind(tr_data$W, tr_data$X)
        q_arg_tr <- cbind(tr_data$Z, tr_data$X)
        K_H_tr <- rbfkernel(h_arg_tr, sigma = kernel_sigma)
        K_Q_tr <- rbfkernel(q_arg_tr, sigma = kernel_sigma)
        ###################################################################################
        
        ## Optimization ###################################################################
        ident_n_tr <- diag(1, n_tr)
        ident_r_tr <- diag(1, r_tr)
        D1_tr <- diag(as.numeric(tr_data$A))
        D0_tr <- diag(as.numeric(1 - tr_data$A))
        
        ###################################################################################
        
        
        ## Gaussian RBF kernel ############################################################
        h_arg_te <- cbind(te_data$W, te_data$X)
        q_arg_te <- cbind(te_data$Z, te_data$X)
        K_H_tr_te <- rbfkernel(X = h_arg_te, sigma = kernel_sigma, Y = h_arg_tr)
        K_Q_tr_te <- rbfkernel(X = q_arg_te, sigma = kernel_sigma, Y = q_arg_tr)
        K_H_te <- rbfkernel(h_arg_te, sigma = kernel_sigma)
        K_Q_te <- rbfkernel(q_arg_te, sigma = kernel_sigma)
        
        ident_n_te <- diag(1, n_te)
        ident_r_te <- diag(1, r_te)
        
        D1_te <- diag(as.numeric(te_data$A))
        D0_te <- diag(as.numeric(1 - te_data$A))
        
        sqrtmtSK_QS_tr <- sqrtm(ginv(t(S_tr) %*% K_Q_tr %*% S_tr))
        D_tr <- K_Q_tr %*% S_tr %*% sqrtmtSK_QS_tr
        
        sqrtmtSK_HS_tr <- sqrtm(ginv(t(S_tr) %*% K_H_tr %*% S_tr))
        V_tr <- K_H_tr %*% S_tr %*% sqrtmtSK_HS_tr
        
        sqrtmtSK_QS_te <- sqrtm(ginv(t(S_te) %*% K_Q_te %*% S_te))
        D_te <- K_Q_te %*% S_te %*% sqrtmtSK_QS_te
        
        sqrtmtSK_HS_te <- sqrtm(ginv(t(S_te) %*% K_H_te %*% S_te))
        V_te <- K_H_te %*% S_te %*% sqrtmtSK_HS_te
        
        
        
        
        Theta <- 0.25 * D_te %*% ginv(t(D_te) %*% D_te / n_te + theta * ident_r_te) %*% t(D_te)
        ## find the best combinations to minimize the inside maximization of h1
        for (lm_Hh1_index in 1:length(lm_Hh1_list)) {
          for (lm_Qh1_index in 1:length(lm_Qh1_list)) {
            lm_Hh1 <- lm_Hh1_list[lm_Hh1_index]
            lm_Qh1 <- lm_Qh1_list[lm_Qh1_index]
            Gamm <- 0.25 * D_tr %*% ginv( t(D_tr) %*% D_tr/n_tr + lm_Qh1/n_tr^0.7*ident_r_tr ) %*% t(D_tr)
            alpha1h <- ginv(t(V_tr) %*% D1_tr %*% Gamm %*% D1_tr %*% V_tr + n_tr^2*lm_Hh1/n_tr^0.7*ident_r_tr) %*%
              t(V_tr) %*% D1_tr %*% Gamm %*% (tr_data$Y*tr_data$A)
            xi_n <- (D1_te %*% K_H_tr_te %*% S_tr %*% alpha1h + te_data$Y * te_data$A) / n_te
            riskh1 <- Re(t(xi_n) %*% Theta %*% xi_n)
            riskh1_result[lm_Hh1_index, lm_Qh1_index] <- riskh1_result[lm_Hh1_index, lm_Qh1_index] + riskh1
          }
        }
        
        
        ## find the best combinations to minimize the inside maximization of h0
        for (lm_Hh0_index in 1:length(lm_Hh0_list)) {
          for (lm_Qh0_index in 1:length(lm_Qh0_list)) {
            lm_Hh0 <- lm_Hh0_list[lm_Hh0_index]
            lm_Qh0 <- lm_Qh0_list[lm_Qh0_index]
            Gamm <- 0.25 * D_tr %*% ginv( t(D_tr) %*% D_tr/n_tr + lm_Qh0/n_tr^0.7*ident_r_tr ) %*% t(D_tr)
            alpha0h <- ginv(t(V_tr) %*% D0_tr %*% Gamm %*% D0_tr %*% V_tr + n_tr^2*lm_Hh0/n_tr^0.7*ident_r_tr) %*%
              t(V_tr) %*% D0_tr %*% Gamm %*% (tr_data$Y*(1 - tr_data$A))
            xi_n <- (D0_te %*% K_H_tr_te %*% S_tr %*% alpha0h + te_data$Y * (1 - te_data$A)) / n_te
            riskh0 <- Re(t(xi_n) %*% Theta %*% xi_n)
            riskh0_result[lm_Hh0_index, lm_Qh0_index] <- riskh0_result[lm_Hh0_index, lm_Qh0_index] + riskh0
          }
        }
        
  
        Theta <- 0.25 * V_te %*% ginv(t(V_te) %*% V_te / n_te + theta * ident_r_te) %*% t(V_te)
        ## find the best combinations to minimize the inside maximization of q1
        for (lm_Hq1_index in 1:length(lm_Hq1_list)) {
          for (lm_Qq1_index in 1:length(lm_Qq1_list)) {
            lm_Hq1 <- lm_Hq1_list[lm_Hq1_index]
            lm_Qq1 <- lm_Qq1_list[lm_Qq1_index]
            Gamm <- 0.25 * V_tr %*% ginv( t(V_tr) %*% V_tr/n_tr + lm_Hq1/n_tr^0.7*ident_r_tr) %*% t(V_tr)
            alpha1q <- ginv(t(D_tr) %*% D1_tr %*% Gamm %*% D1_tr %*% D_tr + n_tr^2*lm_Qq1/n_tr^0.7*ident_r_tr) %*%
              t(D_tr) %*% D1_tr %*% Gamm %*% matrix(1,n_tr,1)
            xi_n <- (D1_te %*% K_Q_tr_te %*% S_tr %*% alpha1q + 1) / n_te
            riskq1 <- Re(t(xi_n) %*% Theta %*% xi_n)
            riskq1_result[lm_Hq1_index, lm_Qq1_index] <- riskq1_result[lm_Hq1_index, lm_Qq1_index] + riskq1
          }
        }
        ## find the best combinations to minimize the inside maximization of q0
        for (lm_Hq0_index in 1:length(lm_Hq0_list)) {
          for (lm_Qq0_index in 1:length(lm_Qq0_list)) {
            lm_Hq0 <- lm_Hq0_list[lm_Hq0_index]
            lm_Qq0 <- lm_Qq0_list[lm_Qq0_index]
            Gamm <- 0.25 * V_tr %*% ginv( t(V_tr) %*% V_tr/n_tr + lm_Hq0/n_tr^0.7*ident_r_tr) %*% t(V_tr)
            alpha0q <- ginv(t(D_tr) %*% D0_tr %*% Gamm %*% D0_tr %*% D_tr + n_tr^2*lm_Qq0/n_tr^0.7*ident_r_tr) %*%
              t(D_tr) %*% D0_tr %*% Gamm %*% matrix(1,n_tr,1)
            xi_n <- (D0_te %*% K_Q_tr_te %*% S_tr %*% alpha0q + 1) / n_te
            riskq0 <- Re(t(xi_n) %*% Theta %*% xi_n)
            riskq0_result[lm_Hq0_index, lm_Qq0_index] <- riskq0_result[lm_Hq0_index, lm_Qq0_index] + riskq0
          }
        }
      }
    }
    result_lm_Hh1 <- lm_Hh1_list[which(riskh1_result == min(riskh1_result), arr.ind = TRUE)[1, 1]]
    result_lm_Qh1 <- lm_Qh1_list[which(riskh1_result == min(riskh1_result), arr.ind = TRUE)[1, 2]]
    result_lm_Hh0 <- lm_Hh0_list[which(riskh0_result == min(riskh0_result), arr.ind = TRUE)[1, 1]]
    result_lm_Qh0 <- lm_Qh0_list[which(riskh0_result == min(riskh0_result), arr.ind = TRUE)[1, 2]]
    result_lm_Hq1 <- lm_Hq1_list[which(riskq1_result == min(riskq1_result), arr.ind = TRUE)[1, 1]]
    result_lm_Qq1 <- lm_Qq1_list[which(riskq1_result == min(riskq1_result), arr.ind = TRUE)[1, 2]]
    result_lm_Hq0 <- lm_Hq0_list[which(riskq0_result == min(riskq0_result), arr.ind = TRUE)[1, 1]]
    result_lm_Qq0 <- lm_Qq0_list[which(riskq0_result == min(riskq0_result), arr.ind = TRUE)[1, 2]]
    # hyperparameters <- list(result_lm_Hh1 = result_lm_Hh1, 
    #                         result_lm_Qh1 = result_lm_Qh1, result_lm_Hh0 = result_lm_Hh0, 
    #                         result_lm_Qh0 = result_lm_Qh0,
    #                      result_lm_Hq1 = result_lm_Hq1,
    #                      result_lm_Qq1 = result_lm_Qq1,
    #                      result_lm_Hq0 = result_lm_Hq0,
    #                      result_lm_Qq0 = result_lm_Qq0)
    hyperparameters <- c(result_lm_Hh1, result_lm_Qh1, result_lm_Hh0, result_lm_Qh0,
                       result_lm_Hq1, result_lm_Qq1, result_lm_Hq0, result_lm_Qq0)
    risk_sum <- c(risk_sum, min(riskh1_result) + min(riskh0_result) + min(riskq1_result) + min(riskq0_result))
    hyperparameters_mat <- rbind(hyperparameters_mat, hyperparameters)
  }
  result <- list(kernel_sigma = kernel_sigma_list[which.min(risk_sum)], 
                 hyperparameters = hyperparameters_mat[which.min(risk_sum), ])
  return(result)
}


# cross fitting evaluation
CF_eval <- function(data, kernel_sigma, lm_Hh1, lm_Qh1, lm_Hh0, lm_Qh0, lm_Hq1, lm_Qq1, lm_Hq0, lm_Qq0, CF_K_fold, Nystroem = FALSE, appro_rate = 0.05){
  
  
  n_sample <- length(data$Y)
  working_data <- with(data, list(X = X, Z = Z, W = W, A = A, Y = Y))
  random_sample <- createFolds(1:n_sample, k = CF_K_fold)
  OR_IF <- c()
  IPW_IF <- c()
  DR_IF <- c()
  
  if (!Nystroem) {
    print("We do not use a low rank approximation.")
    
    for (cf_rep in 1:CF_K_fold) {
      tr_data <- sapply(working_data, function(x) matrix(x[-random_sample[[cf_rep]], ], ncol = ncol(x)))
      eval_data <- sapply(working_data, function(x) matrix(x[random_sample[[cf_rep]], ], ncol = ncol(x)))
      n_tr <- n_sample - length(random_sample[[cf_rep]])
      n_eval <- length(random_sample[[cf_rep]])
      
      ## Gaussian RBF kernel ############################################################
      h_arg <- cbind(tr_data$W, tr_data$X)
      q_arg <- cbind(tr_data$Z, tr_data$X)
      K_H <- rbfkernel(h_arg, sigma = kernel_sigma)
      K_Q <- rbfkernel(q_arg, sigma = kernel_sigma)
      ###################################################################################
      
      ## Optimization ###################################################################
      ident_n <- diag(rep(1, n_tr))
      D1 <- diag(as.numeric(tr_data$A), n_tr, n_tr)
      D0 <- diag(as.numeric(1 - tr_data$A), n_tr, n_tr)
      ## Optimization h ##
      Gamm <- 0.25 * K_Q %*% ginv(K_Q / n_tr + lm_Qh1 / n_tr^0.7 * ident_n)
      alpha1h <- ginv(K_H %*% D1 %*% Gamm %*% D1 %*% K_H + n_tr^2 * lm_Hh1 / n_tr^0.7 * K_H) %*%
        K_H %*% D1 %*% Gamm %*% (tr_data$Y * tr_data$A)
      Gamm <- 0.25 * K_Q %*% ginv(K_Q / n_tr + lm_Qh0 / n_tr^0.7 * ident_n)
      
      alpha0h <- ginv(K_H %*% D0 %*% Gamm %*% D0 %*% K_H + n_tr^2 * lm_Hh0 / n_tr^0.7 * K_H) %*%
        K_H %*% D0 %*% Gamm %*% (tr_data$Y * (1 - tr_data$A))
      ## Optimization q ##
      Gamm <- 0.25 * K_H %*% ginv( K_H / n_tr + lm_Hq1 / n_tr^0.7 * ident_n )
      alpha1q <- ginv(K_Q %*% D1 %*% Gamm %*% D1 %*% K_Q + n_tr^2 * lm_Qq1 / n_tr^0.7 * K_Q) %*%
        K_Q %*% D1 %*% Gamm %*% matrix(1, n_tr, 1)
      Gamm <- 0.25 * K_H %*% ginv( K_H / n_tr + lm_Hq0 / n_tr^0.7 * ident_n )
      alpha0q <- ginv(K_Q %*% D0 %*% Gamm %*% D0 %*% K_Q + n_tr^2 * lm_Qq0 / n_tr^0.7 * K_Q) %*%
        K_Q %*% D0 %*% Gamm %*% matrix(1, n_tr, 1)
      ###################################################################################
      
      ## Estimate function h_hat ########################################################
      h_hat_1 <- function(w1, x1){
        T_mat <- c(w1,x1)
        S <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(T_mat,1,length(T_mat)))
        return(t(S) %*% alpha1h)
      }
      h_hat_0 <- function(w1, x1){
        T_mat <- c(w1,x1)
        S <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(T_mat,1,length(T_mat)))
        return(t(S) %*% alpha0h)
      }
      h_hat <- function(w1, a1, x1){
        return(a1*h_hat_1(w1,x1)+(1-a1)*h_hat_0(w1,x1))
      }
      ###################################################################################
      
      ## Estimate function q_hat ########################################################
      q_hat_1 <- function(z1, x1){
        T_mat <- c(z1,x1)
        S <- rbfkernel(q_arg, sigma = kernel_sigma,matrix(T_mat, 1, length(T_mat)))
        return(t(S) %*% alpha1q)
      }
      q_hat_0 <- function(z1, x1){
        T_mat <- c(z1,x1)
        S <- rbfkernel(q_arg, sigma = kernel_sigma,matrix(T_mat, 1, length(T_mat)))
        return(t(S) %*% alpha0q)
      }
      q_hat <- function(z1, a1, x1){
        return(a1 * q_hat_1(z1, x1) + (1 - a1) * q_hat_0(z1, x1))
      }
      ###################################################################################
      
      ## Evaluation #####################################################################
      #gt <- data$truth
      ## Evaluation OR ##################################################################
      est <- 0
      for (j in seq(n_eval)){
        IF <- (1 / n_eval) * (h_hat(eval_data$W[j, ], 1, eval_data$X[j, ]) - h_hat(eval_data$W[j, ], 0, eval_data$X[j,]))
        OR_IF <- c(OR_IF, IF)
      }
  
      ## Evaluation IPW #################################################################
      est <- 0
      for (j in seq(n_eval)){
        IF <- (1 / n_eval) * (-1)^(1 - eval_data$A[j, ]) * q_hat(eval_data$Z[j, ], eval_data$A[j, ], eval_data$X[j, ]) *
          eval_data$Y[j, ]
        IPW_IF <- c(IPW_IF, IF)
      }
  
      ## Evaluation DR ##################################################################
      est <- 0
      var_est <- 0
      for (j in seq(n_eval)){
        IF <- (1 / n_eval) * (-1)^(1 - eval_data$A[j, ]) * q_hat(eval_data$Z[j, ], eval_data$A[j, ], eval_data$X[j, ]) *
          (eval_data$Y[j, ] - h_hat(eval_data$W[j, ], eval_data$A[j, ], eval_data$X[j, ])) +
          (1 / n_eval) * (h_hat(eval_data$W[j, ], 1, eval_data$X[j, ]) - h_hat(eval_data$W[j, ], 0, eval_data$X[j,]))
        DR_IF <- c(DR_IF, IF)
      }
  
      ###################################################################################
    }
  } else {
    print("We do use a low rank approximation.")
    
    for (cf_rep in 1:CF_K_fold) {
      tr_data <- sapply(working_data, function(x) matrix(x[-random_sample[[cf_rep]], ], ncol = ncol(x)))
      eval_data <- sapply(working_data, function(x) matrix(x[random_sample[[cf_rep]], ], ncol = ncol(x)))
      n_tr <- n_sample - length(random_sample[[cf_rep]])
      n_eval <- length(random_sample[[cf_rep]])
      r <- max(10, ceiling(n_tr * appro_rate))      
      S <- diag(1, n_tr)[, sort(sample(1:n_tr, r))]
      
      ## Gaussian RBF kernel ############################################################
      h_arg <- cbind(tr_data$W, tr_data$X)
      q_arg <- cbind(tr_data$Z, tr_data$X)
      K_H <- rbfkernel(h_arg, sigma = kernel_sigma)
      K_Q <- rbfkernel(q_arg, sigma = kernel_sigma)
      ###################################################################################
      
      ## Optimization ###################################################################
      ident_n <- diag(1, n_tr)
      ident_r <- diag(1, r)
      D1 <- diag(as.numeric(tr_data$A))
      D0 <- diag(as.numeric(1-tr_data$A))
      ## Optimization h ##
      sqrtmtSK_QS <- sqrtm(ginv(t(S) %*% K_Q %*% S))
      D <- K_Q %*% S %*% sqrtmtSK_QS
      Gamm <- 0.25*D %*% ginv( t(D) %*% D/n_tr + lm_Qh1/n_tr^0.7*ident_r ) %*% t(D)
      sqrtmtSK_HS <- sqrtm(ginv(t(S) %*% K_H %*% S))
      V <- K_H %*% S %*% sqrtmtSK_HS
      alpha1h <- ginv(t(V) %*% D1 %*% Gamm %*% D1 %*% V+n_tr^2*lm_Hh1/n_tr^0.7*ident_r) %*%
        t(V) %*% D1 %*% Gamm %*% (tr_data$Y*tr_data$A)
      Gamm <- 0.25*D %*% ginv( t(D) %*% D/n_tr + lm_Qh0/n_tr^0.7*ident_r ) %*% t(D)
      alpha0h <- ginv(t(V) %*% D0 %*% Gamm %*% D0 %*% V+n_tr^2*lm_Hh0/n_tr^0.7*ident_r) %*%
        t(V) %*% D0 %*% Gamm %*% (tr_data$Y*(1-tr_data$A))
      ## Optimization q ##
      Gamm <- 0.25*V %*% ginv( t(V) %*% V/n_tr + lm_Hq1/n_tr^0.7*ident_r ) %*% t(V)
      alpha1q <- ginv(t(D) %*% D1 %*% Gamm %*% D1 %*% D+n_tr^2*lm_Qq1/n_tr^0.7*ident_r) %*%
        t(D) %*% D1 %*% Gamm %*% matrix(1,n_tr,1)
      Gamm <- 0.25*V %*% ginv( t(V) %*% V/n_tr + lm_Hq0/n_tr^0.7*ident_r ) %*% t(V)
      alpha0q <- ginv(t(D) %*% D0 %*% Gamm %*% D0 %*% D+n_tr^2*lm_Qq0/n_tr^0.7*ident_r) %*%
        t(D) %*% D0 %*% Gamm %*% matrix(1,n_tr,1)
      ###################################################################################
      ## Estimate function h_hat ########################################################
      h_hat_1 <- function(w1,x1){
        W <- c(w1,x1)
        K <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W,1,length(W)))
        return(t(K) %*% S %*% sqrtmtSK_HS %*% alpha1h)
      }
      h_hat_0 <- function(w1,x1){
        W <- c(w1,x1)
        K <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W,1,length(W)))
        return(t(K) %*% S %*% sqrtmtSK_HS %*% alpha0h)
      }
      h_hat <- function(w1,a1,x1){
        return(a1*h_hat_1(w1,x1)+(1-a1)*h_hat_0(w1,x1))
      }
      ###################################################################################
      
      ## Estimate function q_hat ########################################################
      q_hat_1 <- function(z1,x1){
        Z <- c(z1,x1)
        K <- rbfkernel(q_arg, sigma = kernel_sigma, matrix(Z,1,length(Z)))
        return(t(K) %*% S %*% sqrtmtSK_QS %*% alpha1q)
      }
      q_hat_0 <- function(z1,x1){
        Z <- c(z1,x1)
        K <- rbfkernel(q_arg, sigma = kernel_sigma, matrix(Z,1,length(Z)))
        return(t(K) %*% S %*% sqrtmtSK_QS %*% alpha0q)
      }
      q_hat <- function(z1,a1,x1){
        return(a1*q_hat_1(z1,x1)+(1-a1)*q_hat_0(z1,x1))
      }
      ###################################################################################
      
      ## Evaluation #####################################################################
      #gt <- data$truth
      ## Evaluation OR ##################################################################
      est <- 0
      for (j in seq(n_eval)){
        IF <- (1 / n_eval) * (h_hat(eval_data$W[j, ], 1, eval_data$X[j, ]) - h_hat(eval_data$W[j, ], 0, eval_data$X[j,]))
        OR_IF <- c(OR_IF, IF)
      }
      
      ## Evaluation IPW #################################################################
      est <- 0
      for (j in seq(n_eval)){
        IF <- (1 / n_eval) * (-1)^(1 - eval_data$A[j, ]) * q_hat(eval_data$Z[j, ], eval_data$A[j, ], eval_data$X[j, ]) *
          eval_data$Y[j, ]
        IPW_IF <- c(IPW_IF, IF)
      }
      
      ## Evaluation DR ##################################################################
      est <- 0
      var_est <- 0
      for (j in seq(n_eval)){
        IF <- (1 / n_eval) * (-1)^(1 - eval_data$A[j, ]) * q_hat(eval_data$Z[j, ], eval_data$A[j, ], eval_data$X[j, ]) *
          (eval_data$Y[j, ] - h_hat(eval_data$W[j, ], eval_data$A[j, ], eval_data$X[j, ])) +
          (1 / n_eval) * (h_hat(eval_data$W[j, ], 1, eval_data$X[j, ]) - h_hat(eval_data$W[j, ], 0, eval_data$X[j,]))
        DR_IF <- c(DR_IF, IF)
      }
      
      ###################################################################################
    }
  }
  gt <- data$truth
  
  OR_IF <- OR_IF / CF_K_fold
  IPW_IF <- IPW_IF / CF_K_fold
  DR_IF <- DR_IF / CF_K_fold
  OR_est <- sum(OR_IF) - gt
  IPW_est <- sum(IPW_IF) - gt
  DR_est <- sum(DR_IF) - gt
  OR_IF <- OR_IF - mean(OR_IF)
  IPW_IF <- IPW_IF - mean(IPW_IF)
  DR_IF <- DR_IF - mean(DR_IF)
  
  result <- c(OR_est, IPW_est, DR_est, sqrt(sum(OR_IF^2)), sqrt(sum(IPW_IF^2)), sqrt(sum(DR_IF^2)))
  return(result)
  
}
