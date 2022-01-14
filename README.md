# StatisticalComputingSolutions
### Numeric Optimization
```
library(bootstrap)
library(boot)
library(GA)

set.seed(2020)
J<- 10
p <- 0.3
A <- emails <- matrix(NA, J, J)
lambdas <- c(2,3)
for(i in 1:(J-1)){
  for(j in (i+1):J){
    A[i,j] <- A[j, i] <- sample(c(1,0), 1, prob = c(p,1-p))
    emails[i,j] <- emails[j,i] <- rpois(n=1, lambda = lambdas[A[i,j] +1])
  }
}
emails[1,]
```
### Log-Likelihood & Maximum Likelihood Estimates for Emailing Service
```
emails_ll <- function(par =c(x, y, p), X) {
  return(sum(log(par[3] * (par[2]^X) * exp(-1 * par[2]) / factorial(X) + (1-par[3]) *
           (par[1]^X) * exp(-1 * par[1]) / factorial(X)), na.rm = TRUE) /2)
}

emails_ll(par = c(2, 3, 0.3), X = emails)
emails_ll(par = c(1, 2, 0.5), X = emails)

opt <- optim(c(1, 3, 0.5), emails_ll, X = emails, control = list(fnscale = -1))

```

### EM Algorithm for Maximum Liklihood Parameter Estimates 
```
n_iter <- 1e3
lambda_hat <- matrix(NA, n_iter, 3)
lambda_hat[1, ] <- c(1, 2, 0.5)
for(iter in 2:n_iter){
  X_vec <- emails[upper.tri(emails)]
  p1 <- lambda_hat[iter - 1, 3] * dpois(X_vec, lambda_hat[iter - 1, 2])
  p0 <- (1 - lambda_hat[iter - 1, 3]) * dpois(X_vec, lambda_hat[iter - 1, 1])
  gamma_hat <- p1 / (p0 + p1)
  lambda_hat[iter, 2] <- sum(gamma_hat * X_vec) / 
    sum(gamma_hat)
  lambda_hat[iter, 1] <- sum((1 - gamma_hat) * X_vec) / 
    sum(1 - gamma_hat)
  lambda_hat[iter, 3] <- sum(gamma_hat) / length(X_vec)
}
lambda_hat[1000,]
```
### Bootstrapping | Jackknife and Quantile Estimation
```
knitr::opts_chunk$set(cache = T, fig.width = 4, fig.align = 'center')
library(bootstrap)
library(boot)
load("/Users/michaelherrington/Desktop/water_qual.RData")
sample_corr <- cor(water_qual$median_cl2, water_qual$median_income)
sample_corr

#jackknife for bias estimation w/90th quantile for chlorine levels

n <- nrow(water_qual)
theta.hat <- cor(water_qual$median_cl2, water_qual$median_income)
print(theta.hat)

theta.jack <- numeric(n)
for(i in 1:n)
  theta.jack[i] <- cor(water_qual$median_cl2[-i], water_qual$median_income[-i])
bias <- (n-1) * (mean(theta.jack) - theta.hat)
print(bias)

ninety <- quantile(water_qual$median_cl2, .90)

#bootstrapped CI @ 95th percentile for teu 90th quantile: 1000 estimates

set.seed(1234)
B <- 1000
q95 <- numeric(B)
for (i in 1:B) {
  water_qual.b = sample(water_qual$median_cl2, length(water_qual$median_cl2), replace=T)
  q95[i] <- quantile(water_qual.b, 0.90)
}
2*ninety - quantile(q95, c(.975, .025))

ten <- quantile(water_qual$median_cl2, 0.10)
set.seed(1234)
B <- 1000
q95 <- numeric(B)
for (i in 1:B) {
  water_qual.b = sample(water_qual$median_cl2, length(water_qual$median_cl2), replace=T)
  q95[i] <- quantile(water_qual.b, 0.1)
}
2*ten - quantile(q95, c(.975, .025))
```

### Linear-Model for Median Chlorine Levels by $population, $median_income, $prop_children, $L0_health
```
fit <- lm(water_qual$median_cl2 ~ water_qual$population + water_qual$median_income + water_qual$prop_children + water_qual$LO_health, data = water_qual)
summary(fit)

get_regression_coefs <- function(data, ind){
  fit <- lm(median_cl2 ~., data = data[ind,])
  coef(fit)
}
get_regression_coefs(water_qual, 1:10)
```
