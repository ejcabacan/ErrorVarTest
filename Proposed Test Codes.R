

##--- Nonparametric Test for Change in Error Variance of Multiple Time Series ---##

## Data Generating Process

# This function generates the kth time series, with specified parameters:

# k	    = kth of the N time series
# L	    = length of time series
# brkpt = break/change point
# rho 	= autoregressive parameter
# epsd	= vector of standard deviation of the error terms with length 2
# cv    = coefficient of variation of the random effects

gendata <- function(k,L,brkpt,rho,epsd,cv){
  
  # Create the data frames for a time series before and during the change point
  ts1 <- NULL
  ts2 <- NULL
  
  # Initial values for the series [y - current value; yl - previous value]
  y <- 0
  yl <- 0
  
  # Random effect for individual i
  rmean <- runif(1,10,100)
  rsd <- rmean*cv
  lambda <- rnorm(1, rmean, rsd)
  
  for(n in -500:(brkpt-1)){
    # Simulate an AR(1) data before the change point
    eps<-rnorm(1,0,epsd[1])
    ytf<-rho*y+eps
    
    # Add the random effect for individual i
    ylt <- ytf + lambda
    
    # Data frame
    ts<-data.frame(indv=k, time=n, yt=ylt, yt_1=yl)
    ts1<- rbind(ts1, ts)
    
    # store values
    y<-ytf
    yl <- ylt
  }
  
  for(n in brkpt:L){
    # Simulate an AR(1) data during the change point
    eps<-rnorm(1,0,epsd[2])
    ytf<-rho*y+eps
    
    # Add the random effect for individual i
    ylt <- ytf + lambda
    
    # Data frame
    ts<-data.frame(indv=k, time=n, yt=ylt, yt_1=yl)
    ts2<- rbind(ts2, ts)
    
    # store values
    y<-ytf
    yl <- ylt
  }
  
  # Combine series before and during the change point
  series <- rbind(ts1,ts2)
  
  # Delete the first 500 observations in the series
  series <- subset(series,time>0)
  series$indv <- as.factor(series$indv)
  
  return(series)
}




## Estimation Procedure

# This function estimates the diferent parameters of the proposed model before and after the change point

# data     = data frame that includes all Y(t)s in a column and all Y(t-1)s in another column
# N        = number of time series in the data
# L        = length of each series
# brkpt    = known break/changepoint
# max_iter = number of maximum allowable iteration before force stopping the algorithm 
# cut      = convergence criterion (absolute percent change)

estimate <- function(data,brkpt,max_iter,cut){
  
  # Obtain unique series IDs for iteration
  series.id <- unique(data$indv)
  
  # Set initial set-up for convergence
  conv <- FALSE
  a <- 0
  
  # Estimation proper
  while(!conv  && a < max_iter){
    a <- a + 1
    if(a == 1){
      
      # Step 1 : Estimate random effects per time series
      rand_effects <- lme(yt~1,random=~1|indv,control=lmeControl(opt='optim'), method="REML", data=data)
      
      # Save estimated random effects
      reff_p <- as.numeric(unlist(coef(rand_effects)))
      
      # Step 2 : Define residual resid1
      data$reff <- as.numeric(fitted(rand_effects))
      data$resid1 <- data$yt - data$reff
      
      # Step 3 : Estimate autoregressive parameter rho thru modified forward search algorithm
      ar_comp <- NULL
      
      for(n in 1:length(series.id)){
        sub <- subset(data,indv==series.id[n])
        
        ar_n1 <- NULL
        
        if(is.error(try(arima(sub$resid1,order=c(1,0,0)),silent=TRUE))==TRUE){
          ar_n1 <- 1
        } else{
          ar_mod <- arima(sub$resid1,order=c(1,0,0))
          ar_n1 <- as.numeric(coef(ar_mod)[1])
        }
        
        ar_comp <- c(ar_comp,ar_n1)
      }
      
      # Obtain bootstrap resample mean for rho
      ar_bs_comp <- NULL
      
      for(c in 1:200){
        samp_ind <- sample(length(ar_comp),length(ar_comp),replace=TRUE)
        ar_samp <- ar_comp[samp_ind]
        ar_bs_comp <- c(ar_bs_comp,mean(ar_samp,na.rm=TRUE))
      }
      
      # Save the estimated bootstrap mean
      ar_p <- mean(ar_bs_comp,na.rm=TRUE)
      
      # Step 4 : Define new residual resid2
      data$resid2 <- data$yt - ar_p*data$yt_1
      
    } else{
      
      # Step 5 : Estimate random effects per time series
      rand_effects <- lme(resid2~1,random=~1|indv,control=lmeControl(opt='optim'), method="REML",data=data)
      
      # Save estimated random effects
      reff_n <- as.numeric(unlist(coef(rand_effects)))
      
      # Step 6 : Define residual resid1
      data$reff <- as.numeric(fitted(rand_effects))
      data$resid1 <- data$yt - data$reff
      
      # Step 7 : Estimate autoregressive parameter rho thru modified forward search algorithm
      ar_comp <- NULL
      
      for(n in 1:length(series.id)){
        sub <- subset(data,indv==series.id[n])
        
        ar_n1 <- NULL
        
        if(is.error(try(arima(sub$resid1,order=c(1,0,0)),silent=TRUE))==TRUE){
          ar_n1 <- 1
        } else{
          ar_mod <- arima(sub$resid1,order=c(1,0,0))
          ar_n1 <- as.numeric(coef(ar_mod)[1])
        }
        ar_comp <- c(ar_comp,ar_n1)
      }
      
      # Obtain bootstrap resample mean for rho
      ar_bs_comp <- NULL
      for(c in 1:200){
        samp_ind <- sample(length(ar_comp),length(ar_comp),replace=TRUE)
        ar_samp <- ar_comp[samp_ind]
        ar_bs_comp <- c(ar_bs_comp,mean(ar_samp,na.rm=TRUE))
      }
      
      # Save the estimated bootstrap mean
      ar_n <- mean(ar_bs_comp,na.rm=TRUE)
      
      # Step 8 : Define new residual resid2
      data$resid2 <- data$yt - ar_n*data$yt_1
      
      # Check for convergence of the estimation procedure
      check_p <- abs((ar_n - ar_p)/ar_p)
      
      # When at least one random effect mean is estimated to be zero, set convergence for the random effects to be TRUE
      
      if(!0 %in% reff_p){
        check_reff <- mean(abs((reff_n - reff_p)/reff_p))
      }else{
        check_reff <- 0        
      }
      if(check_p < cut && check_reff < cut){
        conv = TRUE
      }
      
      # Setting each estimate as previous estimates for the next iteration
      ar_p <- ar_n
      reff_p <- reff_n
    }
  }
  # Variance of lambdas
  lambsd <- as.numeric(VarCorr(rand_effects)[1,1])
  
  # Compute for the fitted values
  data$fitted <- ar_n*data$yt + data$reff
  
  # Final residuals
  data$residual <- data$yt - data$fitted
  
  # Compute the mean residual per series
  group_mean <- mean(data$residual, na.rm=TRUE)
  
  # Compute for the adjusted residuals
  data$gresid <- data$residual - group_mean
  
  # Squared residuals for MSE
  data$residual_sq <- data$gresid^2
  
  # Overall MSE before and after breakpoint
  mse_bef <- mean(data$residual_sq[data$time<brkpt], na.rm=TRUE)
  mse_aft <- mean(data$residual_sq[data$time>=brkpt], na.rm=TRUE)
  
  # Return a list of values
  return(list(ar=ar_n,reff=reff_n,lsd=lambsd,mse_bef=mse_bef,mse_aft=mse_aft,num_iter=a))
}




## Bootstrap Method

# This function generates the kth time series to be used in AR Sieve Bootstrap

# k       = kth time series
# L.s     = length of time series
# brkpt.s = break/change point
# rho.s   = estimated autoregressive parameter
# rmean.s = estimated random effects
# epsd.s  = estimated vector of standard deviation of the error terms with length 2

sieve_data <- function(k,L.s,brkpt.s,rho.s,rmean.s,lsd.s,epsd.s){
  
  # Create the data frame for a time series before and after the change point
  ts1 <- NULL
  ts2 <- NULL
  
  # Initial values for the series
  y <- 0
  yl <- 0
  
  # Bounds for the AR Sieve error component
  uppr <- sqrt(3*epsd.s)
  lwr <- (-uppr)
  
  # Bounds for lambdas
  upp_lambda <- rmean.s[k] + sqrt(3*lsd.s)
  low_lambda <- rmean.s[k] - sqrt(3*lsd.s)
  
  # Random effects for individual k
  lambda <- runif(1, low_lambda, upp_lambda)
  
  for(n in -500:(brkpt.s-1)){
    
    # Simulate the dgp of an AR(1) model before the changepoint
    eps <- runif(1,lwr[1],uppr[1])
    ytf <- rho.s*y+eps
    
    # Add random effects for individual k
    ylt <- ytf + lambda
    
    # Data frame
    ts<-data.frame(indv=k, time=n, yt=ylt, yt_1=yl, eps=eps, lambda=lambda)
    ts1<- rbind(ts1, ts)
    
    # store values
    y<-ytf
    yl <- ylt
  }
  
  for(n in brkpt.s:L.s){
    
    # Simulate the dgp of an AR(1) model during the changepoint
    eps <- runif(1,lwr[2],uppr[2])
    ytf <- rho.s*y+eps
    
    # Add random effects for individual k
    ylt <- ytf + lambda
    
    # Data frame
    ts<-data.frame(indv=k, time=n, yt=ylt, yt_1=yl, eps=eps, lambda=lambda)
    ts2<- rbind(ts2, ts)
    
    # store values
    y<-ytf
    yl <- ylt
  }
  
  # Combine the two series before and during changepoint
  series <- rbind(ts1,ts2)
  
  # Delete the first 500 observations in the series
  series <- subset(series,time>0)
  series$indv <- as.factor(series$indv)
  
  return(series)
}



## Test for Change in Error Variance

# This function performs the proposed test for change in error variance given a multiple time series data and a known possible change point.

# yn      = data frame that includes all Y(t)s in a column and all Y(t-1)s in another column
# N       = number of time series in the data
# L       = length of each series
# brkpt   = known break/changepoint
# j       = seed number

ErrorVarTest <- function(yn,N,L,brkpt,j){
  
  # Estimate the model
  est_model <- estimate(yn,brkpt=brkpt,max_iter=1000,cut=0.001)
  
  # Estimated model parameters for the sieve bootstrap
  rho.s   <- est_model$ar
  rmean.s <- est_model$reff
  epsd.s  <- c(est_model$mse_bef,est_model$mse_aft)
  lsd.s <- est_model$lsd
  
  # Set up parallel computing parameters
  no_cores <- detectCores()
  no_cores
  cl <- makeCluster(no_cores-1)
  registerDoParallel(cl)
  getDoParWorkers()
  
  # Replicate the current scenario (v-u) times in parallel
  ratio_data <- foreach(m=1:200,.combine='rbind', .export=c('estimate','sieve_data'),  .packages=c('nlme','assertthat')) %dopar% {
    
    # Set the random number generator
    seed_no <- ((3^7)*j+m)
    set.seed(seed_no, kind="L'Ecuyer-CMRG")
    
    # Regenerate the N time series
    sdata <- do.call(rbind,lapply(1:N,sieve_data,L,brkpt,rho.s,rmean.s,lsd.s,epsd.s))
    
    ## Re-estimate the model
    sievetrial <- estimate(sdata,brkpt=brkpt,max_iter=1000,cut=0.001)
    
    # Compute for the Ratio Statistic
    comp_ratio <- sievetrial$mse_bef / sievetrial$mse_aft
    
    return(comp_ratio)
  }
  
  stopCluster(cl)
  
  # Determine the 95% bootstrap confidence interval
  p2.5 <- quantile(ratio_data,prob=0.025,na.rm=TRUE)
  p97.5 <- quantile(ratio_data,prob=0.975,na.rm=TRUE)
  
  # Determine if the confidence interval contains 1 
  if(p2.5<1 && p97.5>1){
    tag <- "Do not reject Ho"
  } else{
    tag <- "Reject Ho"
  }
  
  # Test Statistic
  MSE_Ratio  <- round(est_model$mse_bef/est_model$mse_aft,3)
  
  # Output
  TestResult <- data.frame(value=as.numeric(MSE_Ratio),lcl=as.numeric(p2.5),ucl=as.numeric(p97.5),Result=tag)
  
  return(list(TestResult))
  
}

## Packages Needed

library(nlme)
library(assertthat)
library(doParallel)


## Example 1
# Set the random number generator
set.seed(1)

# Simulate Data
y1<-do.call(rbind,lapply(1:5,gendata,L=150,brkpt=75,rho=0.6,epsd=c(2,2),cv=0.05))

# Test for Change in Error Variance
ErrorVarTest(y1,N=5,L=150,brkpt=75,j=123)



## Example 2
# Set the random number generator
set.seed(2)

# Simulate Data
y2<-do.call(rbind,lapply(1:5,gendata,L=150,brkpt=75,rho=0.6,epsd=c(2,2.5),cv=0.05))

# Test for Change in Error Variance
ErrorVarTest(y2,N=5,L=150,brkpt=75,j=124)



## Example 3
# Set the random number generator
set.seed(3)

# Simulate Data
y3<-do.call(rbind,lapply(1:5,gendata,L=150,brkpt=75,rho=0.6,epsd=c(2,5),cv=0.05))

# Test for Change in Error Variance
ErrorVarTest(y3,N=5,L=150,brkpt=75,j=125)



## Example 4
# Set the random number generator
set.seed(4)

# Simulate Data
y4<-do.call(rbind,lapply(1:5,gendata,L=150,brkpt=75,rho=0.6,epsd=c(2,10),cv=0.05))

# Test for Change in Error Variance
ErrorVarTest(y4,N=5,L=150,brkpt=75,j=126)





