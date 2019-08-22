
library(quantmod)
library(xts)
library(geoR)
library(mvtnorm)
library(dplyr)

# Define asset shortnames and current weights and other variables
stock_shortnames <- c("FB","AAPL")
current_weights <- c(0, 0)
n_days_garch <- 300


# Get asset data
getSymbols(stock_shortnames,src="yahoo")

# Merge time series
close_prices <- merge.xts(AAPL$AAPL.Close, FB$FB.Close, join = "inner")
dates <- index(close_prices)
n_days <- length(dates)
dt <- 1/255
close_prices <- coredata(close_prices)


# Log returns
log_return <- tail(log(close_prices[-1,]/close_prices[-n_days,]),n_days_garch)

# Negative MLE
garch_mle<-function(lambda ,s ,dt ){
  # Param
  mu=lambda[1]
  beta0=lambda[2]
  beta1=lambda[3]
  beta2=lambda[4]
  alpha1=lambda[5]
  alpha2=lambda[6]
  
  # Initialisation
  l=0
  v <- integer(length(s)+1)
  v[1]=(sd(s)^2)*dt
  
  # Iterate over dt
  for (i in 1:length(s)){
    l <- l - 0.5*log(v[i]) - 0.5*(s[i] - mu*dt)^2/(v[i]*dt)
    
    if(s[i]<0){
      v[i+1]<-beta0 + beta1*v[i] + beta2/dt*(s[i]-alpha1*dt)^2
    }else{
      v[i+1]<-beta0 + beta1*v[i] + beta2/dt*(s[i]-alpha2*dt)^2
    }
  }
  return(-l)
}


#Define constants
lambda_start <- c(0.1,0.001,0.9,0.05,0.1,0.1)
ui <- matrix(c(0,1,0,0,0,0, #beta0 LB
               0,0,1,0,0,0, #beta1 LB
               0,0,0,1,0,0, #beta2 LB
               0,0,-1,0,0,0, #beta1 UB
               0,0,0,-1,0,0, #beta2 UB
               0,0,-1,-1,0,0), #beta1+beta2<=1
             nrow = 6,
             byrow = TRUE
)
ci <- c(0,0,0,-1,-1,-1)
cov <- cov(log_return)

#Run Optimization and list the results
res <- constrOptim(lambda_start,garch_mle, grad = NULL, ui = ui, ci = ci, s = log_return[,1], dt = dt )
for (i in 2:ncol(log_return)) {
  res <- list(res,  constrOptim(lambda_start,garch_mle, grad = NULL, ui = ui, ci = ci, s = log_return[,i], dt = dt ))
}

garch_params <- lapply(res, '[[','par')
names(garch_params) <- stock_shortnames
garch_params <- as_tibble(garch_params)

## Test 
# Iterate over returns and calculate eps
eps <- matrix(,nrow = n_days_garch,ncol = length(stock_shortnames))
for (j in 1:2){
  v <- integer(n_days_garch+1)
  v[1]=(sd(log_return[,j])^2)*dt
  
  mu=garch_params[[j]][1]
  beta0=garch_params[[j]][2]
  beta1=garch_params[[j]][3]
  beta2=garch_params[[j]][4]
  alpha1=garch_params[[j]][5]
  alpha2=garch_params[[j]][6]
  
  for (i in 1:n_days_garch){
    (log_return[i,j]-mu*dt)/(sqrt(v[i]*dt))
    eps[i,j] <- (log_return[i,j]-mu*dt)/(sqrt(v[i]*dt))
    
    if(log_return[i,j]<0){
      v[i+1]<-beta0 + beta1*v[i] + beta2/dt*(log_return[i]-alpha1*dt)^2
    }else{
      v[i+1]<-beta0 + beta1*v[i] + beta2/dt*(log_return[i]-alpha2*dt)^2
    }
  }
  qqnorm(eps[,j], main = stock_shortnames[j])
  qqline(eps[,j], main = stock_shortnames[j])
}

# Covariance from eps
sigma = cov(eps)

stock_sim <- function(){
  
}
  
#simulate 


#nu=lambda[1]
#beta0=lambda[2]
#beta1=lambda[3]
#eta2=lambda[4]
#lpha1=lambda[5]
#alpha2=lambda[6]

#s <- log_return[,1]

