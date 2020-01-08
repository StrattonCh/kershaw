setwd("~/GitHub/Summer2019/Kershaw Model/hierarchical kershaw model/models")
set.seed(06112019)
packs <- c('tidyverse', 'lubridate', 'ggthemr', 'coda', 'arrayhelpers')
lapply(packs, require, character.only = T)
load("~/GitHub/Summer2019/Kershaw Model/hierarchical kershaw model/scraping and cleaning/kershaw_clean.Rdata")
my.prog <- function(print = .05*num.mcmc, begin, num.mcmc, i, model = 'model 1', chain = 'chain 1'){
  if(i %% print == 0){
    cat("\014")
    runtime <- (Sys.time() - begin)
    percent <- round(i/num.mcmc * 100, 2)
    message <- paste('\nIteration ', i, ' of ', num.mcmc, ' for ', model, ', ', chain, '; ', percent, '% done with this chain. Current runtime of ', round(runtime, 2), ' ', attr(runtime, 'units'), '.\n', sep = "")
    cat(message)
    txtProgressBar(min = 2, max = num.mcmc, initial = i, style = 3)
  }
}

kershaw.full <- kershaw
kershaw.full <- kershaw.full %>%
  dplyr::mutate(date = sapply(strsplit(gameday_link, '_'), FUN = function(x) paste(x[2:4], collapse = '-')) %>%
                  ymd(.)) %>%
  mutate(year = year(date),
         month = month(date),
         day = day(date)) %>%
  mutate(count = factor(count),
         bat_order = factor(bat_order))
kershaw.train <- sample_frac(kershaw.full, .7)
kershaw.test <- anti_join(kershaw.full, kershaw.train)

save(kershaw, kershaw.train, kershaw.test, file = 'data.Rdata')

nchains <- 3
num.mcmc <- 10000
for(chain in 1:nchains){
  ###############
  ### MODEL 1 ###
  ###############
  {
    # define model formula
    formula <- as.formula(pitch_type ~ count)
    
    # prep for hierarchical model
    df.list <- split(kershaw.train, f = kershaw.train$year)
    X.list <- lapply(df.list, FUN = function(x){model.matrix(formula, data = x)})
    y.list <- lapply(df.list, FUN = function(x){
      mf <- model.frame(formula = formula, data = x)
      resp.mat <- dummies::dummy(model.response(mf))
      colnames(resp.mat) <- levels(model.response(mf))
      resp.mat
    })
    
    # setup sampler
    M <- length(y.list)
    J <- length(unique(kershaw.train$pitch_type))
    p <- ncol(model.matrix(formula, kershaw.train))
    
    # storage
    theta.mcmc <- array(0, dim = c(num.mcmc, p, J))
    dimnames(theta.mcmc)[[2]] <- colnames(X.list[[1]])
    dimnames(theta.mcmc)[[3]] <- colnames(y.list[[1]])
    
    Sigma.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.mcmc)[[3]] <- colnames(y.list[[1]])
    Sigma.inv.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.inv.mcmc)[[3]] <- colnames(y.list[[1]])
    beta.mcmc <- list()
    for(m in 1:M){
      beta.mcmc[[m]] <- array(0, dim = c(num.mcmc, p, J))
      dimnames(beta.mcmc[[m]])[[2]] <- colnames(X.list[[1]])
      dimnames(beta.mcmc[[m]])[[3]] <- colnames(y.list[[1]])
    }
    names(beta.mcmc) <- names(X.list)
    
    # initialize 
    theta.mcmc[1,,] <- rnorm(p*J)
    Sigma.mcmc[1,,] <- c(diag(p))
    Sigma.inv.mcmc[1,,] <- c(diag(p))
    for(m in 1:M){
      beta.mcmc[[m]][1,,] <- rnorm(p*J)
    }
    
    # priors 
    mu0 <- matrix(rep(0, p), ncol = 1)
    Lambda0 <- diag(p); Lambda0.inv <- solve(Lambda0)
    prior.prod <- Lambda0.inv %*% mu0
    eta0 <- p
    S0 <- diag(p)
    
    # run sampler
    begin <- Sys.time()
    for(i in 2:num.mcmc){
      # update hyper-parameters - INEFFICIENT, NEED TO OPTIMIZE
      for(j in 1:(J-1)){
        Sigma.inv <- matrix(Sigma.inv.mcmc[i-1,,j], p, p)
        
        ## theta
        beta <- matrix(0, p, M)
        for(m in 1:M){
          beta[,m] <- beta.mcmc[[m]][i-1,,j]
        }
        beta.bar <- rowMeans(beta)
        Lambda.m <- solve(Lambda0.inv + M*Sigma.inv)
        mu.m <- Lambda.m %*% (prior.prod + M*Sigma.inv %*% beta.bar)
        theta <- matrix(mnormt::rmnorm(1, mu.m, Lambda.m), ncol = 1)
        
        ## Sigma
        tmp <- beta - c(theta)
        tmp2 <- array(0, dim = c(p, p, M))
        for(m in 1:M){
          tmp2[,,m] <- tmp[,m] %*% t(tmp[,m ])
        }
        S.theta <- apply(tmp2, c(1:2), sum)
        # Sn <- solve(S0 + S.theta)
        # Sigma <- round(MCMCpack::riwish(eta0 + M, Sn), 6)
        Sigma.inv <- round(rWishart(1, eta0 + M, S0 + S.theta), 6)
        
        # store hyper-parameters
        theta.mcmc[i,,j] <- c(theta)
        # Sigma.mcmc[i,,j] <- c(Sigma)
        Sigma.inv.mcmc[i,,j] <- c(Sigma.inv)
      }
      
      # update regression coefficients
      for(m in 1:M){
        # convenience
        N <- nrow(X.list[[m]])
        X <- X.list[[m]]
        y <- y.list[[m]]
        kappa <- y - 1/2
        
        for(j in 1:(J-1)){
          # convenience
          Sigma.inv <- matrix(Sigma.inv.mcmc[i,,j], p, p)
          theta <- matrix(theta.mcmc[i,,j], ncol = 1)
          prod <- Sigma.inv %*% theta
          
          # calculate matrix of linear predictors
          linpred <- X %*% beta.mcmc[[m]][i-1,,]
          
          # update latent omegas
          C <- log(rowSums(exp(linpred[,-j])))
          eta <- linpred[,j] - C
          omega <- BayesLogit::rpg(N, rep(1, N), eta)
          Omega <- diag(omega)
          
          # update beta
          V <- solve(t(X) %*% Omega %*% X + Sigma.inv)
          mean <- V %*% (t(X) %*% (kappa[,j] + Omega %*% C) + prod)
          beta.mcmc[[m]][i,,j] <- matrix(mvtnorm::rmvnorm(1, mean, V), ncol = 1)
        }
      }
      
      # progress
      if(T) my.prog(begin = begin, num.mcmc = num.mcmc, i = i, model = 'model 1', chain = paste0('chain ', chain))
    }
    mod.chain <- paste0('mod1_', chain)
    file.name <- paste0(mod.chain, '.Rdata')
    assign(mod.chain, list(beta = beta.mcmc, theta = theta.mcmc, Sigma = Sigma.mcmc))
    
    save(list = mod.chain, file = file.name)
  }
  
  ###############
  ### MODEL 2 ###
  ###############
  {
    # define model formula
    formula <- as.formula(pitch_type ~ count + prev_pitch)
    
    # prep for hierarchical model
    df.list <- split(kershaw.train, f = kershaw.train$year)
    X.list <- lapply(df.list, FUN = function(x){model.matrix(formula, data = x)})
    y.list <- lapply(df.list, FUN = function(x){
      mf <- model.frame(formula = formula, data = x)
      resp.mat <- dummies::dummy(model.response(mf))
      colnames(resp.mat) <- levels(model.response(mf))
      resp.mat
    })
    
    # setup sampler
    M <- length(y.list)
    J <- length(unique(kershaw.train$pitch_type))
    p <- ncol(model.matrix(formula, kershaw.train))
    
    # storage
    theta.mcmc <- array(0, dim = c(num.mcmc, p, J))
    dimnames(theta.mcmc)[[2]] <- colnames(X.list[[1]])
    dimnames(theta.mcmc)[[3]] <- colnames(y.list[[1]])
    
    Sigma.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.mcmc)[[3]] <- colnames(y.list[[1]])
    Sigma.inv.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.inv.mcmc)[[3]] <- colnames(y.list[[1]])
    beta.mcmc <- list()
    for(m in 1:M){
      beta.mcmc[[m]] <- array(0, dim = c(num.mcmc, p, J))
      dimnames(beta.mcmc[[m]])[[2]] <- colnames(X.list[[1]])
      dimnames(beta.mcmc[[m]])[[3]] <- colnames(y.list[[1]])
    }
    names(beta.mcmc) <- names(X.list)
    
    # initialize 
    theta.mcmc[1,,] <- rnorm(p*J)
    Sigma.mcmc[1,,] <- c(diag(p))
    Sigma.inv.mcmc[1,,] <- c(diag(p))
    for(m in 1:M){
      beta.mcmc[[m]][1,,] <- rnorm(p*J)
    }
    
    # priors 
    mu0 <- matrix(rep(0, p), ncol = 1)
    Lambda0 <- diag(p); Lambda0.inv <- solve(Lambda0)
    prior.prod <- Lambda0.inv %*% mu0
    eta0 <- p
    S0 <- diag(p)
    
    # run sampler
    begin <- Sys.time()
    for(i in 2:num.mcmc){
      # update hyper-parameters - INEFFICIENT, NEED TO OPTIMIZE
      for(j in 1:(J-1)){
        Sigma.inv <- matrix(Sigma.inv.mcmc[i-1,,j], p, p)
        
        ## theta
        beta <- matrix(0, p, M)
        for(m in 1:M){
          beta[,m] <- beta.mcmc[[m]][i-1,,j]
        }
        beta.bar <- rowMeans(beta)
        Lambda.m <- solve(Lambda0.inv + M*Sigma.inv)
        mu.m <- Lambda.m %*% (prior.prod + M*Sigma.inv %*% beta.bar)
        theta <- matrix(mnormt::rmnorm(1, mu.m, Lambda.m), ncol = 1)
        
        ## Sigma
        tmp <- beta - c(theta)
        tmp2 <- array(0, dim = c(p, p, M))
        for(m in 1:M){
          tmp2[,,m] <- tmp[,m] %*% t(tmp[,m ])
        }
        S.theta <- apply(tmp2, c(1:2), sum)
        # Sn <- solve(S0 + S.theta)
        # Sigma <- round(MCMCpack::riwish(eta0 + M, Sn), 6)
        Sigma.inv <- round(rWishart(1, eta0 + M, S0 + S.theta), 6)
        
        # store hyper-parameters
        theta.mcmc[i,,j] <- c(theta)
        # Sigma.mcmc[i,,j] <- c(Sigma)
        Sigma.inv.mcmc[i,,j] <- c(Sigma.inv)
      }
      
      # update regression coefficients
      for(m in 1:M){
        # convenience
        N <- nrow(X.list[[m]])
        X <- X.list[[m]]
        y <- y.list[[m]]
        kappa <- y - 1/2
        
        for(j in 1:(J-1)){
          # convenience
          Sigma.inv <- matrix(Sigma.inv.mcmc[i,,j], p, p)
          theta <- matrix(theta.mcmc[i,,j], ncol = 1)
          prod <- Sigma.inv %*% theta
          
          # calculate matrix of linear predictors
          linpred <- X %*% beta.mcmc[[m]][i-1,,]
          
          # update latent omegas
          C <- log(rowSums(exp(linpred[,-j])))
          eta <- linpred[,j] - C
          omega <- BayesLogit::rpg(N, rep(1, N), eta)
          Omega <- diag(omega)
          
          # update beta
          V <- solve(t(X) %*% Omega %*% X + Sigma.inv)
          mean <- V %*% (t(X) %*% (kappa[,j] + Omega %*% C) + prod)
          beta.mcmc[[m]][i,,j] <- matrix(mvtnorm::rmvnorm(1, mean, V), ncol = 1)
        }
      }
      
      # progress
      if(T) my.prog(begin = begin, num.mcmc = num.mcmc, i = i, model = 'model 2', chain = paste0('chain ', chain))
    }
    mod.chain <- paste0('mod2_', chain)
    file.name <- paste0(mod.chain, '.Rdata')
    assign(mod.chain, list(beta = beta.mcmc, theta = theta.mcmc, Sigma = Sigma.mcmc))
    
    save(list = mod.chain, file = file.name)
  }
  
  ###############
  ### MODEL 3 ###
  ###############
  {
    # define model formula
    formula <- as.formula(pitch_type ~ count + prev_pitch + pitch_count)
    
    # prep for hierarchical model
    df.list <- split(kershaw.train, f = kershaw.train$year)
    X.list <- lapply(df.list, FUN = function(x){model.matrix(formula, data = x)})
    y.list <- lapply(df.list, FUN = function(x){
      mf <- model.frame(formula = formula, data = x)
      resp.mat <- dummies::dummy(model.response(mf))
      colnames(resp.mat) <- levels(model.response(mf))
      resp.mat
    })
    
    # setup sampler
    M <- length(y.list)
    J <- length(unique(kershaw.train$pitch_type))
    p <- ncol(model.matrix(formula, kershaw.train))
    
    # storage
    theta.mcmc <- array(0, dim = c(num.mcmc, p, J))
    dimnames(theta.mcmc)[[2]] <- colnames(X.list[[1]])
    dimnames(theta.mcmc)[[3]] <- colnames(y.list[[1]])
    
    Sigma.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.mcmc)[[3]] <- colnames(y.list[[1]])
    Sigma.inv.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.inv.mcmc)[[3]] <- colnames(y.list[[1]])
    beta.mcmc <- list()
    for(m in 1:M){
      beta.mcmc[[m]] <- array(0, dim = c(num.mcmc, p, J))
      dimnames(beta.mcmc[[m]])[[2]] <- colnames(X.list[[1]])
      dimnames(beta.mcmc[[m]])[[3]] <- colnames(y.list[[1]])
    }
    names(beta.mcmc) <- names(X.list)
    
    # initialize 
    theta.mcmc[1,,] <- rnorm(p*J)
    Sigma.mcmc[1,,] <- c(diag(p))
    Sigma.inv.mcmc[1,,] <- c(diag(p))
    for(m in 1:M){
      beta.mcmc[[m]][1,,] <- rnorm(p*J)
    }
    
    # priors 
    mu0 <- matrix(rep(0, p), ncol = 1)
    Lambda0 <- diag(p); Lambda0.inv <- solve(Lambda0)
    prior.prod <- Lambda0.inv %*% mu0
    eta0 <- p
    S0 <- diag(p)
    
    # run sampler
    begin <- Sys.time()
    for(i in 2:num.mcmc){
      # update hyper-parameters - INEFFICIENT, NEED TO OPTIMIZE
      for(j in 1:(J-1)){
        Sigma.inv <- matrix(Sigma.inv.mcmc[i-1,,j], p, p)
        
        ## theta
        beta <- matrix(0, p, M)
        for(m in 1:M){
          beta[,m] <- beta.mcmc[[m]][i-1,,j]
        }
        beta.bar <- rowMeans(beta)
        Lambda.m <- solve(Lambda0.inv + M*Sigma.inv)
        mu.m <- Lambda.m %*% (prior.prod + M*Sigma.inv %*% beta.bar)
        theta <- matrix(mnormt::rmnorm(1, mu.m, Lambda.m), ncol = 1)
        
        ## Sigma
        tmp <- beta - c(theta)
        tmp2 <- array(0, dim = c(p, p, M))
        for(m in 1:M){
          tmp2[,,m] <- tmp[,m] %*% t(tmp[,m ])
        }
        S.theta <- apply(tmp2, c(1:2), sum)
        # Sn <- solve(S0 + S.theta)
        # Sigma <- round(MCMCpack::riwish(eta0 + M, Sn), 6)
        Sigma.inv <- round(rWishart(1, eta0 + M, S0 + S.theta), 6)
        
        # store hyper-parameters
        theta.mcmc[i,,j] <- c(theta)
        # Sigma.mcmc[i,,j] <- c(Sigma)
        Sigma.inv.mcmc[i,,j] <- c(Sigma.inv)
      }
      
      # update regression coefficients
      for(m in 1:M){
        # convenience
        N <- nrow(X.list[[m]])
        X <- X.list[[m]]
        y <- y.list[[m]]
        kappa <- y - 1/2
        
        for(j in 1:(J-1)){
          # convenience
          Sigma.inv <- matrix(Sigma.inv.mcmc[i,,j], p, p)
          theta <- matrix(theta.mcmc[i,,j], ncol = 1)
          prod <- Sigma.inv %*% theta
          
          # calculate matrix of linear predictors
          linpred <- X %*% beta.mcmc[[m]][i-1,,]
          
          # update latent omegas
          C <- log(rowSums(exp(linpred[,-j])))
          eta <- linpred[,j] - C
          omega <- BayesLogit::rpg(N, rep(1, N), eta)
          Omega <- diag(omega)
          
          # update beta
          V <- solve(t(X) %*% Omega %*% X + Sigma.inv)
          mean <- V %*% (t(X) %*% (kappa[,j] + Omega %*% C) + prod)
          beta.mcmc[[m]][i,,j] <- matrix(mvtnorm::rmvnorm(1, mean, V), ncol = 1)
        }
      }
      
      # progress
      if(T) my.prog(begin = begin, num.mcmc = num.mcmc, i = i, model = 'model 3', chain = paste0('chain ', chain))
    }
    mod.chain <- paste0('mod3_', chain)
    file.name <- paste0(mod.chain, '.Rdata')
    assign(mod.chain, list(beta = beta.mcmc, theta = theta.mcmc, Sigma = Sigma.mcmc))
    
    save(list = mod.chain, file = file.name)
  }
  
  ###############
  ### MODEL 4 ###
  ###############
  {
    # define model formula
    formula <- as.formula(pitch_type ~ count + prev_pitch + bat_order)
    
    # prep for hierarchical model
    df.list <- split(kershaw.train, f = kershaw.train$year)
    X.list <- lapply(df.list, FUN = function(x){model.matrix(formula, data = x)})
    y.list <- lapply(df.list, FUN = function(x){
      mf <- model.frame(formula = formula, data = x)
      resp.mat <- dummies::dummy(model.response(mf))
      colnames(resp.mat) <- levels(model.response(mf))
      resp.mat
    })
    
    # setup sampler
    M <- length(y.list)
    J <- length(unique(kershaw.train$pitch_type))
    p <- ncol(model.matrix(formula, kershaw.train))
    
    # storage
    theta.mcmc <- array(0, dim = c(num.mcmc, p, J))
    dimnames(theta.mcmc)[[2]] <- colnames(X.list[[1]])
    dimnames(theta.mcmc)[[3]] <- colnames(y.list[[1]])
    
    Sigma.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.mcmc)[[3]] <- colnames(y.list[[1]])
    Sigma.inv.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.inv.mcmc)[[3]] <- colnames(y.list[[1]])
    beta.mcmc <- list()
    for(m in 1:M){
      beta.mcmc[[m]] <- array(0, dim = c(num.mcmc, p, J))
      dimnames(beta.mcmc[[m]])[[2]] <- colnames(X.list[[1]])
      dimnames(beta.mcmc[[m]])[[3]] <- colnames(y.list[[1]])
    }
    names(beta.mcmc) <- names(X.list)
    
    # initialize 
    theta.mcmc[1,,] <- rnorm(p*J)
    Sigma.mcmc[1,,] <- c(diag(p))
    Sigma.inv.mcmc[1,,] <- c(diag(p))
    for(m in 1:M){
      beta.mcmc[[m]][1,,] <- rnorm(p*J)
    }
    
    # priors 
    mu0 <- matrix(rep(0, p), ncol = 1)
    Lambda0 <- diag(p); Lambda0.inv <- solve(Lambda0)
    prior.prod <- Lambda0.inv %*% mu0
    eta0 <- p
    S0 <- diag(p)
    
    # run sampler
    begin <- Sys.time()
    for(i in 2:num.mcmc){
      # update hyper-parameters - INEFFICIENT, NEED TO OPTIMIZE
      for(j in 1:(J-1)){
        Sigma.inv <- matrix(Sigma.inv.mcmc[i-1,,j], p, p)
        
        ## theta
        beta <- matrix(0, p, M)
        for(m in 1:M){
          beta[,m] <- beta.mcmc[[m]][i-1,,j]
        }
        beta.bar <- rowMeans(beta)
        Lambda.m <- solve(Lambda0.inv + M*Sigma.inv)
        mu.m <- Lambda.m %*% (prior.prod + M*Sigma.inv %*% beta.bar)
        theta <- matrix(mnormt::rmnorm(1, mu.m, Lambda.m), ncol = 1)
        
        ## Sigma
        tmp <- beta - c(theta)
        tmp2 <- array(0, dim = c(p, p, M))
        for(m in 1:M){
          tmp2[,,m] <- tmp[,m] %*% t(tmp[,m ])
        }
        S.theta <- apply(tmp2, c(1:2), sum)
        # Sn <- solve(S0 + S.theta)
        # Sigma <- round(MCMCpack::riwish(eta0 + M, Sn), 6)
        Sigma.inv <- round(rWishart(1, eta0 + M, S0 + S.theta), 6)
        
        # store hyper-parameters
        theta.mcmc[i,,j] <- c(theta)
        # Sigma.mcmc[i,,j] <- c(Sigma)
        Sigma.inv.mcmc[i,,j] <- c(Sigma.inv)
      }
      
      # update regression coefficients
      for(m in 1:M){
        # convenience
        N <- nrow(X.list[[m]])
        X <- X.list[[m]]
        y <- y.list[[m]]
        kappa <- y - 1/2
        
        for(j in 1:(J-1)){
          # convenience
          Sigma.inv <- matrix(Sigma.inv.mcmc[i,,j], p, p)
          theta <- matrix(theta.mcmc[i,,j], ncol = 1)
          prod <- Sigma.inv %*% theta
          
          # calculate matrix of linear predictors
          linpred <- X %*% beta.mcmc[[m]][i-1,,]
          
          # update latent omegas
          C <- log(rowSums(exp(linpred[,-j])))
          eta <- linpred[,j] - C
          omega <- BayesLogit::rpg(N, rep(1, N), eta)
          Omega <- diag(omega)
          
          # update beta
          V <- solve(t(X) %*% Omega %*% X + Sigma.inv)
          mean <- V %*% (t(X) %*% (kappa[,j] + Omega %*% C) + prod)
          beta.mcmc[[m]][i,,j] <- matrix(mvtnorm::rmvnorm(1, mean, V), ncol = 1)
        }
      }
      
      # progress
      if(T) my.prog(begin = begin, num.mcmc = num.mcmc, i = i, model = 'model 4', chain = paste0('chain ', chain))
    }
    mod.chain <- paste0('mod4_', chain)
    file.name <- paste0(mod.chain, '.Rdata')
    assign(mod.chain, list(beta = beta.mcmc, theta = theta.mcmc, Sigma = Sigma.mcmc))
    
    save(list = mod.chain, file = file.name)
  }
  
  ###############
  ### MODEL 5 ###
  ###############
  {
    # define model formula
    formula <- as.formula(pitch_type ~ count + prev_pitch + stand)
    
    # prep for hierarchical model
    df.list <- split(kershaw.train, f = kershaw.train$year)
    X.list <- lapply(df.list, FUN = function(x){model.matrix(formula, data = x)})
    y.list <- lapply(df.list, FUN = function(x){
      mf <- model.frame(formula = formula, data = x)
      resp.mat <- dummies::dummy(model.response(mf))
      colnames(resp.mat) <- levels(model.response(mf))
      resp.mat
    })
    
    # setup sampler
    M <- length(y.list)
    J <- length(unique(kershaw.train$pitch_type))
    p <- ncol(model.matrix(formula, kershaw.train))
    
    # storage
    theta.mcmc <- array(0, dim = c(num.mcmc, p, J))
    dimnames(theta.mcmc)[[2]] <- colnames(X.list[[1]])
    dimnames(theta.mcmc)[[3]] <- colnames(y.list[[1]])
    
    Sigma.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.mcmc)[[3]] <- colnames(y.list[[1]])
    Sigma.inv.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.inv.mcmc)[[3]] <- colnames(y.list[[1]])
    beta.mcmc <- list()
    for(m in 1:M){
      beta.mcmc[[m]] <- array(0, dim = c(num.mcmc, p, J))
      dimnames(beta.mcmc[[m]])[[2]] <- colnames(X.list[[1]])
      dimnames(beta.mcmc[[m]])[[3]] <- colnames(y.list[[1]])
    }
    names(beta.mcmc) <- names(X.list)
    
    # initialize 
    theta.mcmc[1,,] <- rnorm(p*J)
    Sigma.mcmc[1,,] <- c(diag(p))
    Sigma.inv.mcmc[1,,] <- c(diag(p))
    for(m in 1:M){
      beta.mcmc[[m]][1,,] <- rnorm(p*J)
    }
    
    # priors 
    mu0 <- matrix(rep(0, p), ncol = 1)
    Lambda0 <- diag(p); Lambda0.inv <- solve(Lambda0)
    prior.prod <- Lambda0.inv %*% mu0
    eta0 <- p
    S0 <- diag(p)
    
    # run sampler
    begin <- Sys.time()
    for(i in 2:num.mcmc){
      # update hyper-parameters - INEFFICIENT, NEED TO OPTIMIZE
      for(j in 1:(J-1)){
        Sigma.inv <- matrix(Sigma.inv.mcmc[i-1,,j], p, p)
        
        ## theta
        beta <- matrix(0, p, M)
        for(m in 1:M){
          beta[,m] <- beta.mcmc[[m]][i-1,,j]
        }
        beta.bar <- rowMeans(beta)
        Lambda.m <- solve(Lambda0.inv + M*Sigma.inv)
        mu.m <- Lambda.m %*% (prior.prod + M*Sigma.inv %*% beta.bar)
        theta <- matrix(mnormt::rmnorm(1, mu.m, Lambda.m), ncol = 1)
        
        ## Sigma
        tmp <- beta - c(theta)
        tmp2 <- array(0, dim = c(p, p, M))
        for(m in 1:M){
          tmp2[,,m] <- tmp[,m] %*% t(tmp[,m ])
        }
        S.theta <- apply(tmp2, c(1:2), sum)
        # Sn <- solve(S0 + S.theta)
        # Sigma <- round(MCMCpack::riwish(eta0 + M, Sn), 6)
        Sigma.inv <- round(rWishart(1, eta0 + M, S0 + S.theta), 6)
        
        # store hyper-parameters
        theta.mcmc[i,,j] <- c(theta)
        # Sigma.mcmc[i,,j] <- c(Sigma)
        Sigma.inv.mcmc[i,,j] <- c(Sigma.inv)
      }
      
      # update regression coefficients
      for(m in 1:M){
        # convenience
        N <- nrow(X.list[[m]])
        X <- X.list[[m]]
        y <- y.list[[m]]
        kappa <- y - 1/2
        
        for(j in 1:(J-1)){
          # convenience
          Sigma.inv <- matrix(Sigma.inv.mcmc[i,,j], p, p)
          theta <- matrix(theta.mcmc[i,,j], ncol = 1)
          prod <- Sigma.inv %*% theta
          
          # calculate matrix of linear predictors
          linpred <- X %*% beta.mcmc[[m]][i-1,,]
          
          # update latent omegas
          C <- log(rowSums(exp(linpred[,-j])))
          eta <- linpred[,j] - C
          omega <- BayesLogit::rpg(N, rep(1, N), eta)
          Omega <- diag(omega)
          
          # update beta
          V <- solve(t(X) %*% Omega %*% X + Sigma.inv)
          mean <- V %*% (t(X) %*% (kappa[,j] + Omega %*% C) + prod)
          beta.mcmc[[m]][i,,j] <- matrix(mvtnorm::rmvnorm(1, mean, V), ncol = 1)
        }
      }
      
      # progress
      if(T) my.prog(begin = begin, num.mcmc = num.mcmc, i = i, model = 'model 5', chain = paste0('chain ', chain))
    }
    mod.chain <- paste0('mod5_', chain)
    file.name <- paste0(mod.chain, '.Rdata')
    assign(mod.chain, list(beta = beta.mcmc, theta = theta.mcmc, Sigma = Sigma.mcmc))
    
    save(list = mod.chain, file = file.name)
  }
  
  ###############
  ### MODEL 6 ###
  ###############
  {
    # define model formula
    formula <- as.formula(pitch_type ~ count + prev_pitch + + on.1b + on.2b + on.3b)
    
    # prep for hierarchical model
    df.list <- split(kershaw.train, f = kershaw.train$year)
    X.list <- lapply(df.list, FUN = function(x){model.matrix(formula, data = x)})
    y.list <- lapply(df.list, FUN = function(x){
      mf <- model.frame(formula = formula, data = x)
      resp.mat <- dummies::dummy(model.response(mf))
      colnames(resp.mat) <- levels(model.response(mf))
      resp.mat
    })
    
    # setup sampler
    M <- length(y.list)
    J <- length(unique(kershaw.train$pitch_type))
    p <- ncol(model.matrix(formula, kershaw.train))
    
    # storage
    theta.mcmc <- array(0, dim = c(num.mcmc, p, J))
    dimnames(theta.mcmc)[[2]] <- colnames(X.list[[1]])
    dimnames(theta.mcmc)[[3]] <- colnames(y.list[[1]])
    
    Sigma.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.mcmc)[[3]] <- colnames(y.list[[1]])
    Sigma.inv.mcmc <- array(0, dim = c(num.mcmc, p^2, J))
    dimnames(Sigma.inv.mcmc)[[3]] <- colnames(y.list[[1]])
    beta.mcmc <- list()
    for(m in 1:M){
      beta.mcmc[[m]] <- array(0, dim = c(num.mcmc, p, J))
      dimnames(beta.mcmc[[m]])[[2]] <- colnames(X.list[[1]])
      dimnames(beta.mcmc[[m]])[[3]] <- colnames(y.list[[1]])
    }
    names(beta.mcmc) <- names(X.list)
    
    # initialize 
    theta.mcmc[1,,] <- rnorm(p*J)
    Sigma.mcmc[1,,] <- c(diag(p))
    Sigma.inv.mcmc[1,,] <- c(diag(p))
    for(m in 1:M){
      beta.mcmc[[m]][1,,] <- rnorm(p*J)
    }
    
    # priors 
    mu0 <- matrix(rep(0, p), ncol = 1)
    Lambda0 <- diag(p); Lambda0.inv <- solve(Lambda0)
    prior.prod <- Lambda0.inv %*% mu0
    eta0 <- p
    S0 <- diag(p)
    
    # run sampler
    begin <- Sys.time()
    for(i in 2:num.mcmc){
      # update hyper-parameters - INEFFICIENT, NEED TO OPTIMIZE
      for(j in 1:(J-1)){
        Sigma.inv <- matrix(Sigma.inv.mcmc[i-1,,j], p, p)
        
        ## theta
        beta <- matrix(0, p, M)
        for(m in 1:M){
          beta[,m] <- beta.mcmc[[m]][i-1,,j]
        }
        beta.bar <- rowMeans(beta)
        Lambda.m <- solve(Lambda0.inv + M*Sigma.inv)
        mu.m <- Lambda.m %*% (prior.prod + M*Sigma.inv %*% beta.bar)
        theta <- matrix(mnormt::rmnorm(1, mu.m, Lambda.m), ncol = 1)
        
        ## Sigma
        tmp <- beta - c(theta)
        tmp2 <- array(0, dim = c(p, p, M))
        for(m in 1:M){
          tmp2[,,m] <- tmp[,m] %*% t(tmp[,m ])
        }
        S.theta <- apply(tmp2, c(1:2), sum)
        # Sn <- solve(S0 + S.theta)
        # Sigma <- round(MCMCpack::riwish(eta0 + M, Sn), 6)
        Sigma.inv <- round(rWishart(1, eta0 + M, S0 + S.theta), 6)
        
        # store hyper-parameters
        theta.mcmc[i,,j] <- c(theta)
        # Sigma.mcmc[i,,j] <- c(Sigma)
        Sigma.inv.mcmc[i,,j] <- c(Sigma.inv)
      }
      
      # update regression coefficients
      for(m in 1:M){
        # convenience
        N <- nrow(X.list[[m]])
        X <- X.list[[m]]
        y <- y.list[[m]]
        kappa <- y - 1/2
        
        for(j in 1:(J-1)){
          # convenience
          Sigma.inv <- matrix(Sigma.inv.mcmc[i,,j], p, p)
          theta <- matrix(theta.mcmc[i,,j], ncol = 1)
          prod <- Sigma.inv %*% theta
          
          # calculate matrix of linear predictors
          linpred <- X %*% beta.mcmc[[m]][i-1,,]
          
          # update latent omegas
          C <- log(rowSums(exp(linpred[,-j])))
          eta <- linpred[,j] - C
          omega <- BayesLogit::rpg(N, rep(1, N), eta)
          Omega <- diag(omega)
          
          # update beta
          V <- solve(t(X) %*% Omega %*% X + Sigma.inv)
          mean <- V %*% (t(X) %*% (kappa[,j] + Omega %*% C) + prod)
          beta.mcmc[[m]][i,,j] <- matrix(mvtnorm::rmvnorm(1, mean, V), ncol = 1)
        }
      }
      
      # progress
      if(T) my.prog(begin = begin, num.mcmc = num.mcmc, i = i, model = 'model 6', chain = paste0('chain ', chain))
    }
    mod.chain <- paste0('mod6_', chain)
    file.name <- paste0(mod.chain, '.Rdata')
    assign(mod.chain, list(beta = beta.mcmc, theta = theta.mcmc, Sigma = Sigma.mcmc))
    
    save(list = mod.chain, file = file.name)
  }
}

files <- paste0(rep(paste0('mod', 1:6), each = 3), rep(paste0('_', 1:3), 6), '.Rdata')
lapply(files, load, .GlobalEnv)

###############
### MODEL 1 ###
###############
mod1_beta1 <- mod1_1$beta
mod1_beta2 <- mod1_2$beta
mod1_beta3 <- mod1_3$beta

array2df(mod1_beta1$`2008`)






