
DGP <- function(n, dim_x = 1, dim_z = 1, dim_w = 1, type){
  if (type == "parametric") {
    ## Parametric DGP by Yifan ############################################################
    p <- 2
    b0 <- 2
    bA <- 2
    bW <- 4
    bX <- c(0.25,0.25)
    t0 <- 0.25 
    # tA
    # tZ
    tX <- c(0.25,0.25)
    alpha0 <- 0.25
    alphaA <- 0.25
    alphaX <- c(0.25,0.25)
    mu0 <- 0.25
    # muA
    muX <- c(0.25,0.25)
    kappa0 <- 0.25
    kappaA <- 0.25
    kappaX <- c(0.25,0.25)
    
    sigmaZW <- 0.25
    sigmaZU <- 0.5
    sigmaWU <- 0.5
    sigmaZ <- 1
    sigmaW <- 1
    sigmaU <- 1
    sigmaY <- 0.25
    sigmaX <- 0.25
    omega <- 2
    
    theta0 <- alpha0-kappa0*sigmaZU/sigmaU^2
    thetaA <- alphaA-kappaA*sigmaZU/sigmaU^2
    thetaX <- alphaX-kappaX*sigmaZU/sigmaU^2
    thetaU <- sigmaZU/sigmaU^2
    # constraints
    muA <- sigmaWU*kappaA/sigmaU^2 
    tZ <- -(kappaA-sigmaWU*muA/sigmaW^2)/thetaU/(sigmaU^2-sigmaWU^2/sigmaW^2)
    tA <- -tZ^2*(sigmaZ^2-sigmaZU^2/sigmaU^2)-tZ*thetaA
    
    
    sigma1 <- diag(3)
    sigma1[1,1] <- sigmaZ^2
    sigma1[2,2] <- sigmaW^2
    sigma1[3,3] <- sigmaU^2
    sigma1[1,2] <- sigma1[2,1] <- sigmaZW
    sigma1[1,3] <- sigma1[3,1] <- sigmaZU 
    sigma1[2,3] <- sigma1[3,2] <- sigmaWU
    
    nsimu <- 1
    
    psi.truth <- bA
    ground_truth <- bA
    psi.out <- rep(NA,nsimu)
    psi.ipw <- rep(NA,nsimu)
    psi.ddr <- rep(NA,nsimu)
    psi.dr <- rep(NA,nsimu)
    
    
    ### X
    X <- matrix(rnorm(n*p,0,sigmaX),nrow = n,ncol = p)
    
    ### A
    temp1 <- t0 + tA + X%*%tX + tZ *(theta0+thetaA+X%*%thetaX) +1/2*(1-(sigmaZU^2/sigmaZ^2/sigmaU^2))*tZ^2*sigmaZ^2
    temp2 <- tZ*thetaU*(kappa0+kappaA+X%*%kappaX) + (thetaU*tZ*sigmaU)^2/2
    p.A <- 1/(1+exp(temp1)*exp(temp2))
    A <- rbinom(n,1,p.A)
    
    ### W Z U
    mu <- cbind(alpha0+alphaA*A+X%*%alphaX,mu0+muA*A+X%*%muX,kappa0+kappaA*A+X%*%kappaX)
    zwu <- t(sapply(1:n,function(x){mvrnorm(n = 1, mu = mu[x,], Sigma = sigma1)}))
    Z <- zwu[,1]
    W <- zwu[,2]
    U <- zwu[,3]
    
    # Y
    temp <- mu0+muA*A+X%*%muX+sigmaWU/sigmaU^2*(U-kappa0-kappaA*A-X%*%kappaX)
    Y <- rnorm(n,b0+bA*A+X%*%bX+(bW-omega)*temp+omega*W,sigmaY)
    
    A <- matrix(A,n,1)
    Z <- matrix(Z,n,1)
    W <- matrix(W,n,1)
    Y <- matrix(Y,n,1)
  } else if (type == "nonparametric") {
    ## Nonparametric DGP with the same dimension as in parametric ########################################
    p <- 2
    X <- matrix(rnorm(n * p, 0, 1), ncol = p)
    U <- rnorm(n, 0, 1)
    A <- rbinom(n, 1, 1/(1 + exp(-0.25 - X %*% c(0.2, 0.3) - X^3 %*% c(0.1, -0.05) - 0.25 * U + 0.1 * U^3)))
    Z <- rnorm(n, 0.5 + 0.5 * A + X %*% c(0.2, -0.2) + 0.75 * U, 1)
    W <- rnorm(n, 0.3 + X %*% c(0.35, 0.25) - 0.75 * U, 1)
    Y <- rnorm(n, -0.5 + A + X %*% c(0.25, -0.2) + A * X %*% c(-0.5, 0.3) - X^3 %*% c(0.025, -0.03) - 0.3 * U + 0.25 * A * U + 0.025 * U^3, 1)
    ground_truth <- 1
    A <- matrix(A,n,1)
    Z <- matrix(Z,n,1)
    W <- matrix(W,n,1)
    Y <- matrix(Y,n,1)
  } else if (type == "highd") {
    

  
    M_ux <- matrix(seq(0.4, 0.9, length.out = 5), 5, 1)
    sigma_x <- diag(rep(0.3, 5))
    M_az <- matrix(c(0.5, 0.6), 2, 1)
    M_xz <- matrix(seq(0.4, 0.9, length.out = 2 * 5), 2, 5)
    M_uz <- matrix(c(0.8, 0.9), 2, 1)
    sigma_z <- diag(rep(0.3, 2))
    M_xw <- matrix(seq(0.9, 0.4, length.out = 2 * 5), 2, 5)
    M_uw <- matrix(c(0.7, 0.8), 2, 1)
    sigma_w <- diag(rep(0.3, 2))
    M_ay <- 2#matrix(runif(1,0.3,1),1,1)
    M_wy <- matrix(seq(0.4, 0.9, length.out = 2), 1, 2)
    M_xy <- matrix(seq(0.4, 0.9, length.out = 5), 1, 5)
    M_uy <- matrix(seq(0.4, 0.9, length.out = 2), 1, 2)
    sigma_y <- 0.3
    t_z <- matrix(seq(0.9, 0.4, length.out = 2), 1, 2)
    t_x <- matrix(seq(0.9, 0.4, length.out = 5), 1, 5)
    t_a <- -t_z %*% M_az - t_z %*% sigma_z %*% t(t_z)
    
    ground_truth <- M_ay
    
    X <- matrix(nrow = n, ncol = 5)
    Z <- matrix(nrow = n, ncol = 2)
    W <- matrix(nrow = n, ncol = 2)
    U <- matrix(nrow = n, ncol = 1)
    A <- matrix(nrow = n, ncol = 1)
    Y <- matrix(nrow = n, ncol = 1)


    ground_truth <- M_ay

    for(i in seq(n)){
      U[i] <- rnorm(1, 0, 0.3)
      mean_x <- M_ux %*% U[i]
      X[i,] <- mvrnorm(1, mean_x, sigma_x)
      temp_a <- t_a + t_x %*% matrix(X[i,], 5, 1) + t_z %*% (M_az + M_xz %*% matrix(X[i, ], 5, 1) + M_uz %*% U[i]) +
        0.5 * t_z %*% sigma_z %*% t(t_z)
      p <- 1/(1 + exp(temp_a))
      A[i] <- rbinom(1, size = 1, prob = p)
      mean_z <- M_az %*% A[i] + M_xz %*% matrix(X[i, ], 5, 1) + M_uz %*% U[i]
      Z[i, ] <- mvrnorm(1, mean_z, sigma_z)
      mean_w <- M_xw %*% matrix(X[i, ], 5, 1) + M_uw %*% U[i]
      W[i, ] <- mvrnorm(1, mean_w, sigma_w)
      mean_y <- M_ay %*% A[i] + M_wy %*% matrix(W[i, ], 2, 1) + M_xy %*% matrix(X[i, ], 5, 1) + M_uy %*% mean_w
      Y[i, ] <- mvrnorm(1, mean_y, sigma_y)

  }
  
  } else if (type == "highd2") {
    
    mean_u <- c(0, 0)
    sigma_u <- diag(rep(0.2, 2))
    
    M_ux <- matrix(seq(0.4, 0.9, length.out = 10), 5, 2)
    sigma_x <- diag(rep(0.2, 5))
    M_az <- matrix(c(0.5, 0.6), 2, 1)
    M_xz <- matrix(seq(0.4, 0.9, length.out = 5), 2, 5)
    M_uz <- matrix(c(0.8, 0.9, 0.7, 0.6), 2, 2)
    sigma_z <- diag(rep(0.2, 2))
    M_xw <- matrix(seq(0.9, 0.4, length.out = 2 * 5), 2, 5)
    M_uw <- matrix(c(0.7, 0.8, 0.6, 0.7), 2, 2)
    sigma_w <- diag(rep(0.2, 2))
    M_ay <- 2#matrix(runif(1,0.3,1),1,1)
    M_wy <- matrix(seq(0.4, 0.9, length.out = 1 * 2), 1, 2)
    M_xy <- matrix(seq(0.4, 0.9, length.out = 1 * 5), 1, 5)
    M_uy <- matrix(seq(0.4, 0.9, length.out = 2), 1, 2)
    sigma_y <- 0.2
    t_z <- matrix(seq(0.9, 0.4, length.out = 1 * 2), 1, 2)
    t_x <- matrix(seq(0.9, 0.4, length.out = 5), 1, 5)
    t_a <- -t_z %*% M_az - t_z %*% sigma_z %*% t(t_z)
    
    ground_truth <- M_ay
    
    X <- matrix(nrow = n, ncol = 5)
    Z <- matrix(nrow = n, ncol = 2)
    W <- matrix(nrow = n, ncol = 2)
    U <- matrix(nrow = n, ncol = 2)
    A <- matrix(nrow = n, ncol = 1)
    Y <- matrix(nrow = n, ncol = 1)
    
    
    ground_truth <- M_ay
    
    for(i in seq(n)){
      U[i, ] <- mvrnorm(1, mean_u, sigma_u)
      mean_x <- M_ux %*% U[i, ]
      X[i,] <- mvrnorm(1, mean_x, sigma_x)
      temp_a <- t_a + t_x %*% matrix(X[i,], 5, 1) + t_z %*% (M_az + M_xz %*% matrix(X[i, ], 5, 1) + M_uz %*% U[i, ]) +
        0.5 * t_z %*% sigma_z %*% t(t_z)
      p <- 1/(1 + exp(temp_a))
      A[i] <- rbinom(1, size = 1, prob = p)
      mean_z <- M_az %*% A[i] + M_xz %*% matrix(X[i, ], 5, 1) + M_uz %*% U[i, ]
      Z[i, ] <- mvrnorm(1, mean_z, sigma_z)
      mean_w <- M_xw %*% matrix(X[i, ], 5, 1) + M_uw %*% U[i, ]
      W[i, ] <- mvrnorm(1, mean_w, sigma_w)
      mean_y <- M_ay %*% A[i] + M_wy %*% matrix(W[i, ], 2, 1) + M_xy %*% matrix(X[i, ], 5, 1) + M_uy %*% mean_w
      Y[i, ] <- mvrnorm(1, mean_y, sigma_y)
      
    }
    
  }
  data <- list(X=X, Z=Z, W=W, A=A, Y=Y, truth=ground_truth)
  
  
  return(data)
  
}