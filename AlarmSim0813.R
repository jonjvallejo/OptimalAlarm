# Bayesian prediction for x_1 and x_2 given x_{-1} and x_0

source('C:/Users/Ben/Desktop/Statistics/TimeSeries/AlarmFuncs.R')

set.seed(12345)

x <- arima.sim(n = 1002, list(ar = c(-0.85,0.1)), sd = 1)
z <- x[3:1002]
n <- length(z) # 10000

store <- array(NA, c(n, 2, 7), 
               dimnames = list(NULL, c('upper', 'lower'), 
                               c('50% quantiles', '80% quantiles',
                                 '90% quantiles', '95% quantiles', '97.5% quantiles',
                                 '99% quantiles', '99.5% quantiles')))

for (i in 3:(n-2)){ # predicts next 2 values, only need until n-2
  # ys are standardized!
  y_tplus1s <- rt(10000, df = (i+1))
  y_tplus2s <- rt(10000, df = (i+2))
  # maybe just need unstand_next here...
  x_tplus1s <- mapply(unstand_next, y_tplus1s, i, 
                      MoreArgs = list(data = x))
  x_tplus2s <- mapply(create_x_t0, y_tplus2s, x_tplus1s, i+1, 
                      MoreArgs = list(data = x))
  uppers <- quantile(x_tplus1s, 
                     probs = c(0.50, 0.80, 0.90, 0.95, 0.975, 0.99, 0.995),
                     na.rm = TRUE)
  lowers <- quantile(x_tplus2s, 
                     probs = c(0.50, 0.20, 0.10, 0.05, 0.025, 0.001, 0.0005),
                     na.rm = TRUE)
  store[i,'upper',] <- uppers
  store[i,'lower',] <- lowers
  print(i)
}

# Plot some sections

# ts.plot(z[2:41])
# lines(store[,'upper','50% quantiles'], col = 2)

uppers <- store[, 'upper', '90% quantiles']
lowers <- store[, 'lower', '90% quantiles']

ts.plot(z[-1])
lines(uppers, col = 2)
lines(lowers, col = 3)
plot(uppers, type = 'l', col = 2)
lines(lowers, col = 3)

z_temp <- z[4:999]
up_temp <- uppers[3:998]
low_temp <- lowers[3:998]
sum(z_temp > up_temp)

# Not saved in AlarmSim0316

events <- (z_temp[-1] > 2)
n_events <- sum(events)
n_alarms <- sum(alarms)
n_obs <- length(z_temp)
alarms <- (up_temp[-n_obs] < 2 & low_temp[-n_obs] > 2)

prob_df <- data.frame(Event = mean(z_temp > 2),
                      AlarmSize = mean(alarms),
                      DetectedEvent = sum(alarms & events) / n_alarms,
                      CorrectAlarm = sum(alarms & events) / n_events,
                      FalseAlarm = sum(alarms & !events) / (n_obs - n_events),
                      UndetectedEvent = sum(!alarms & events) / (n_obs - n_alarms))

probs_df <- gen_probs_df(10, 2, z)

mk_alarm <- function(alpha, t, df, tol = 0.01, maxiter = 5)
{
  require(dplyr)
  require(sp)
  
  t <- 10
  # Check t vs t-1
  x1_new <- rt(100000, t)
  x2_new <- rt(100000, t + 1)
  
  k <- 0.5
  k_low <- 0
  k_high <- 1
  diff <- 1
  
  while (diff > tol)
  {
    probs <- df$prob
    
    probs_trunc <- df %>%
      filter(prob > quantile(ecdf(probs), k)) %>%
      mutate(index = interaction(x1, x2))
    
    indexes <- probs_trunc$index
    
    # Initial boundary points
    probs_bdry <- probs_trunc %>%
      mutate(bdry = mapply(get_bdry, x1, x2, MoreArgs = list(indexes = indexes))) %>%
      filter(bdry)
    
    bdry_mat <- cbind(probs_bdry$x1, probs_bdry$x2)
    dists <- dist(bdry_mat)
    dist_mat <- as.matrix(dists)
    min_dists <- apply(dist_mat, 1, function(x){
      temp <- x[x > 0]
      sorted <- sort(temp)
      sorted[1:5]
    })
    max_dists <- apply(min_dists, 2, max)
    keep_pts <- which(max_dists < 5)
    
    # New boundary, with disjoint regions removed
    bdry_new <- probs_bdry %>%
      mutate(idx = 1:nrow(.)) %>%
      filter(idx %in% keep_pts)
    
    pts_in <- point.in.polygon(x1_new, x2_new, probs_bdry$x1, probs_bdry$x2)
    alpha_hat <- mean(pts_in)
    diff <- abs(alpha - alpha_hat)
    print(diff)
    
    k0 <- k
    
    if (alpha_hat > alpha & diff > tol){
      k <- (k0 + k_high) / 2
      k_high <- k_high
      k_low <- k0
    } else if (alpha_hat < alpha & diff > tol){
      k <- (k0 + k_low) / 2
      k_high <- k0
      k_low <- k_low
    }
  }

  list(bdry_df = probs_bdry,
       k = k)
}

alarm <- mk_alarm(0.15, 10, probs_df)

contour(x1, x2, probs_mat)
points(probs_bdry$x1, probs_bdry$x2, col = 2)

library(sp)


contour(x1, x2, probs_mat)
points(bdry_new$x1, bdry_new$x2, col = 2)

pts_in <- point.in.polygon(x1_new, x2_new, bdry_new$x1, bdry_new$x2)


plot(x1_new, x2_new, xlim = c(-20, 20), ylim = c(-20, 20))
points(probs_bdry$x1, probs_bdry$x2, col = 2)

prob(z[8], z[9], 10, 2, z)
