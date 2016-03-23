# m is the length of the t.s. sequence

M <- function(m, data){
  m_act <- m + 2 # lag must be 2 in this example
  if (m < 1) print('m must be greater than 0')
  else {
    data1 <- data[2:(m_act - 1)]
    data2 <- data[1:(m_act - 2)]
    matrix(c(sum(data1*data1), sum(data1*data2), 
             sum(data1*data2), sum(data2*data2)),
           ncol = 2)}
}

C <- function(m, data){
  if (m < 1) print('m must be greater than 0')
  else {
    m_act <- m + 2
    data0 <- data[3:m_act]
    data1 <- data[2:(m_act - 1)]
    data2 <- data[1:(m_act - 2)]
    matrix(c(sum(data0*data1), sum(data0*data2)), ncol = 1)
  }
}

E_fun <- function(t, data){
  t_act <- t + 2
  M_t1 <- M(t+1, data) 
  X <- c(data[t_act], data[(t_act - 1)])
  denom <- (1 - t(X) %*% solve(M_t1) %*% X)
  E <- (t(C(t, data)) %*% solve(M_t1) %*% X) / denom
  E
}

V_fun <- function(t, data){
  t_act <- t + 2
  temp <- data[3:t_act]
  M_t1 <- M(t+1, data)
  C_t <- C(t, data)
  X <- c(data[t_act], data[(t_act - 1)])
  denom <- (1 - t(X) %*% solve(M_t1) %*% X)
  V <- ((sum(temp^2) - t(C_t) %*% solve(M_t1) %*% C_t) /
          denom - ((t(C_t) %*% solve(M_t1) %*% X) / denom)^2) #* (1 / t)
  V / t
}

stand_next <- function(y_next, t, data){
  # taken from eqs (15) - (17)
  t_act <- t + 2
  temp <- data[3:t_act]
  M_t1 <- M(t+1, data) # M^(t+1)
  C_t <- C(t, data) # C^(t)
  X <- c(data[t_act], data[(t_act-1)])
  denom <- (1 - t(X) %*% solve(M_t1) %*% X)
  E <- (t(C_t) %*% solve(M_t1) %*% X) / denom
  V <- ((sum(temp^2) - t(C_t) %*% solve(M_t1) %*% C_t) /
          denom - ((t(C_t) %*% solve(M_t1) %*% X) / denom)^2) * (1 / t)
  V^(-1/2) * (y_next - E)
}

unstand_next <- function(y_next, t, data){
  # taken from eqs (15) - (17)
  t_act <- t + 2
  temp <- data[3:t_act]
  M_t1 <- M(t+1, data) # M^(t+1)
  C_t <- C(t, data) # C^(t)
  X <- c(data[t_act], data[(t_act-1)])
  denom <- (1 - t(X) %*% solve(M_t1) %*% X)
  E <- (t(C_t) %*% solve(M_t1) %*% X) / denom
  V <- ((sum(temp^2) - t(C_t) %*% solve(M_t1) %*% C_t) /
          denom - ((t(C_t) %*% solve(M_t1) %*% X) / denom)^2) * (1 / t)
  V^(1/2) * y_next + E
}

# Use E_fun to see what next value should be
prob <- function(x_1, x_0, t, u, data){
  
  t_act <- t + 2
  temp <- data[1:(t_act-2)]
  data_new0 <- c(temp, x_1, x_0)
  u1 <- stand_next(u, t, data_new0)
  
  integrand <- function(z_stand){
    # x_1 is x_{t-1}, x_0 is x_t
    # need to deal with time misalignment better
    u2 <- function(z_stand){
      z <- unstand_next(z_stand, t, data_new0)
      data_new1 <- c(data_new0, z)
      stand_next(u, t+1, data_new1)
    }
    dt(z_stand, df = t) * (1 - pt(u2(z_stand), df = t+1))
  }
  
  integrate(integrand, -Inf, u1)$value
}

create_x_t0 <- function(y_t0, x_t1, t, data){
  t_act <- t + 2
  temp <- data[1:t_act]
  temp1 <- c(temp[-t_act], x_t1)
  unstand_next(y_t0, t, temp1)
}

create_next <- function(y, x_t0, x_t1, t, data){
  t_act <- t + 2
  idxs <- c(t_act - 1, t_act)
  temp <- data[1:t_act]
  temp1 <- c(temp[-idxs], x_t1, x_t0)
  unstand_next(y, t, temp1)
}

get_bdry <- function(x1, x2, indexes)
{
  
  up <- interaction(x1, x2 + 1)
  down <- interaction(x1, x2 - 1)
  left <- interaction(x1 - 1, x2)
  right <- interaction(x1 + 1, x2)
  
  dirs <- unlist(list(up, down, left, right))
  
  if (all(dirs %in% indexes)){
    val <- FALSE
  } else {
    val <- TRUE
  }
  
  val
}

gen_probs_df <- function(t, u, data)
{
  x1 <- seq(-10, 10, length.out = 100)
  x2 <- seq(-10, 10, length.out = 100)
  len <- length(x1)
  mesh <- expand.grid(x1 = x1, x2 = x2)
  probs <- mapply(prob, mesh$x1, mesh$x2, t, u, MoreArgs = list(data = data))
  probs_df <- data.frame(x1 = mesh$x1,
                         x2 = mesh$x2,
                         prob = probs)
  probs_df
}
