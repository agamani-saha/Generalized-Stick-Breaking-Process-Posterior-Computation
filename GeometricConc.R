###########################################################################################
library(MCMCprecision)
library(markovchain)
set.seed(1)
alpha = 0.5
beta = 2
b0 = 0.1
N = 2
Burn.in = 1000
M = 2000
epsilon <- 0.05
sample.vec = 10^(2:6)
d <- 40
############################################################################################

fkTilde <- function(x, J, A, B) x^(A - 1)*exp(-B*x)/(gamma(x))^J

SamplingFromCoverFunction <- function(N, J, A, B){
  ###initialize (2N+2) length vector of knots
  m <- rep(0, 2*N + 2)
  
  ### initialize (2N+2) length vectors of intercepts and slopes of the tangent lines
  a <- rep(0, 2*N + 2)
  lambda <- rep(0, 2*N + 2)
  
  ### initialize (2N+2) length vector of points of intersection between the tangent lines
  ### concatenate a 0 in front later on to make q[0] = 0
  q <- rep(0, 2*N + 2)
  
  ### log of the k-th density function
  hkTilde <- function(x) {
    retVal <- -J*log(gamma(x)) + (A - 1)*log(x) - B*x
    return(retVal)
  }
  
  ### and its derivative
  hkPrimeTilde <- function(x) {
    retVal <- -J*psigamma(x, deriv = 0) + (A - 1)/x - B
    return(retVal)
  }
  
  ### central knot
  M <- 1.5*(B > 0) + exp(1 - B/J)*(B < 0)
  m[N + 1] <- optimise(hkTilde, interval = c(0, M), maximum = TRUE, tol = 0.000001)$maximum
  
  ### last knot 
  m[2*N + 2] <- (m[N + 1] + 1.5)*(B > 0) + exp(1 - B/J)*(B < 0)
  
  ### first knot
  m[1] <- 0.5*m[N + 1]
  
  ### last-but-one knot
  m[2*N + 1] <- mean(c(m[N + 1], m[2*N + 2]))
  
  ### rest of the notes filled up in an equi-spaced manner
  if (N > 1){
    m[2:N] <- seq(m[1], m[N + 1], length.out = N + 1)[-c(1, N + 1)]
    m[(N + 2):(2*N)] <- seq(m[N + 1], m[2*N + 1], length.out = N + 1)[-c(1, N + 1)]
  }
  
  for (i in 1:(2*N + 2)) {
    a[i] <- hkTilde(m[i])
    lambda[i] <- hkPrimeTilde(m[i])
  }
  
  for (i in 1:(2*N + 1)) q[i] <- (a[i + 1] - a[i] + m[i]*lambda[i] - m[i + 1]*lambda[i + 1])/(lambda[i]- lambda[i + 1])
  q[2*N + 2] <- Inf
  q <- c(0, q)
  
  ### integrating constants
  Cgki <- rep(0, 2*N + 2)
  Cgki[N + 1] <- exp(a[N+1])*(q[N + 2] - q[N + 1])
  
  for (i in setdiff(1:(2*N + 2), N + 1)) Cgki[i] <- exp(a[i] - m[i]*lambda[i])*
    (exp(q[i+1]*lambda[i]) - exp(q[i]*lambda[i]))/lambda[i]
  Cgki[!is.finite(Cgki)] <- 1
  Cgk <- sum(Cgki)
  
  ii <- sample(1:(2*N+2), size = 1, prob = Cgki/Cgk)
  ### now sample from gki
  if (ii == N + 1){
    sim <- runif(1, q[ii], q[ii + 1])
  } else {
    u <- runif(1)
    sim <- log(u*exp(lambda[ii]*q[ii + 1]) + (1 - u)*exp(lambda[ii]*q[ii]))/lambda[ii]
  }
  if(!is.finite(sim)){sim <- 1e2}
  for (i in 1:(2*N + 2)) if (q[i] <= sim && sim < q[i + 1]) 
    gkTildeSim <- exp(a[i] + (sim - m[i])*lambda[i])
  
  return(list(sim = sim, gkTildeSim = gkTildeSim))
}

##### samples from fkTilde using the cover function and rejection sampling method
RejectionSampling <- function(N, J, A, B){
  retVal <- NA
  while (is.na(retVal)){
    ss <- SamplingFromCoverFunction(N, J, A, B)
    gkTildeSim <- ss$gkTildeSim 
    u <- runif(1)
    num = fkTilde(ss$sim, J, A, B)
    den = gkTildeSim
    if(is.nan(num) == T||is.nan(den) == T){retVal <- ss$sim}
    else if (u <= num/den) {retVal <- ss$sim}
  }
  return(retVal)
}

#################### BLOCKED GIBBS SAMPLER #####################################
# Function to draw samples from a dirichlet distribution.
# alpha : concentration parameters

update_Pi = function(n, t){
  
  Pi =  matrix(NA, nrow = d, ncol = d)
  
  for(i in 1:d){
    
    # n_j = (n_j1, ..., n_jL), n_jk = \sum_i I(Z_ji = k)
    n_i = n[i, ]
    
    # Pi |... follows Dirichlet (n_j + t)
    Pi.draw = rdirichlet(1, n_i + t)
    
    # setting an lower bound of 10^(-10) for \pi_jk's to avoid numerical issues
    Pi.ind = which(Pi.draw < 1e-10)
    
    if(length(Pi.ind) > 0){
      
      Pi.draw[Pi.ind] = 1e-10
      
      excess.P = sum(Pi.draw) - 1
      
      ind.max = which.max(Pi.draw)
      Pi.draw[ind.max] = Pi.draw[ind.max] - excess.P
    }
    
    Pi[i, ] = Pi.draw
  }
  return(Pi)
}

update_u = function(d, sum_t){
  u = rgamma(d, shape = sum_t, scale = 1)
  return(u)
}

update_w = function(d, alpha, t){
  w = rep(NA, d)
  for(j in 1:d){
    sum_t_j = sum(t[j:(length(t))])
    w[j] = rgamma(1, shape = alpha, scale = sum_t_j)
  }
  return(w)
}

update_Beta = function(n, N, Pi, alpha, beta, B, u, w){
  
  d <- ncol(n)
  delta <- c(rep(alpha, d - 1), beta)
  
  log_u = sum(log(u))
  log.p = colSums(log(Pi))
  
  # log.p tends to -\infty for unoccupied clusters
  # setting a threshold of -10^10 for log.p to avoid numerical issues
  log.p = sapply(1:d, function(x) ifelse(var(Pi[,x]) == 0, -1e+10, log.p))
  
  # draw t_j from f_j(.), j = 1, 2,..., d.
  t = rep(NA, d)
  for(i in 1:d){
    A = delta[i]
    sum_w = sum(w[1:i])
    B = b0 + sum_w - log.p[i] - log_u
    t[i] = RejectionSampling(N = N, J = d, A = A, B = B)
  }
  sum_t = sum(t)
  Beta = t/sum_t
  
  return(list("Beta" = Beta, "sum_t" = sum_t, "t" = t))
}

blocked_gibbs = function(N, n, Burn.in, M, d, alpha, beta, b0){
  
  # set initial values for running the Gibbs sampler
  Pi = matrix(1/d, nrow = d, ncol = d)
  
  t = rep(1/d, d)
  u = rgamma(n = d, shape = 1, rate = 1)
  w = rgamma(n = d, shape = 1, rate = 1)
  
  # list to store the posterior samples
  Iterates = vector(mode = "list", length = M)
  
  for(m in 1:(M + Burn.in)){
    
    # time at the beginning
    T1 = Sys.time()
    
    # update Pi
    Pi = update_Pi(n = n, t = t)
    
    # update t and Beta 
    res = update_Beta(n = n, N = N, Pi = Pi, alpha = alpha, beta = beta, 
                      B = B, u = u, w = w)
    Beta = res$Beta
    sum_t = res$sum_t
    t = res$t
    
    # update u
    u = update_u(sum_t = sum_t, d = d)
    
    # update w
    w = update_w(d = d, alpha = alpha, t = t)
    
    # time at the end of all updates
    T2 = Sys.time()
    Tdiff =  difftime(T2, T1, units = "secs")
    
    
    # print every 200th iteration
    if(m %% 200 == 0){
      print(paste("iteration :", m))
    }
    
    # store samples after Burn in
    if(m > Burn.in){
      Iterates[[m-Burn.in]] = list("Pi" = Pi, "Beta" = Beta, 
                                   "t" = t, "w" = w,
                                   "u" = u, "time" = Tdiff)
    }
  }
  return(Iterates)
  
}

########################### True TPM ###########################################

trans.matrix.true = matrix(rep(NA, d*d), nrow = d)
for(i in 0:(d-1)){
  for(j in 0:(d-1)){
    p = 1/(log(i+10))
    trans.matrix.true[i+1,j+1] = ((1-p)^j)*p
  }
}
trans.matrix.true
#####################################################################################
sample.mat <- list()
sample.prob.mat <- list()
GSBP.mat <- list()


Pi.error.mat <- matrix(0, nrow = M, ncol = length(sample.vec))


MAE_MLE <- rep(NA, length(sample.vec))
RMSE_MLE <- rep(NA, length(sample.vec))
MAE_GSBP <- rep(NA, length(sample.vec))
RMSE_GSBP <- rep(NA, length(sample.vec))
for(l in 1:length(sample.vec)){
  sample.size = sample.vec[l]
  # sample.size = 10^6
  random.start = 0
  sim.data = rep(NA, sample.size)
  sim.data[1] = random.start
  for(i in 2:sample.size){
    sim.data[i] = rgeom(n = 1, prob = 1/(log(sim.data[i-1] + 10)))
  }
  trans.matrix <- function(X, prob = F)
  {
    tt <- table(c(X[-length(X)]), c(X[-1]))
    if(prob==T) tt <- tt / rowSums(tt)
    tt
  }
  trans.matrix.sample = trans.matrix(sim.data)
  fullmat_dim <- 1 + max(as.numeric(c(rownames(trans.matrix.sample), 
                                      colnames(trans.matrix.sample))))
  fullnames <- as.character(0:(fullmat_dim-1))
  fullmat <- matrix(0, nrow = fullmat_dim, ncol = fullmat_dim, 
                    dimnames = list(fullnames, fullnames))
  fullmat[rownames(trans.matrix.sample), colnames(trans.matrix.sample)] <- 
    trans.matrix.sample
  trans.matrix.sample <- fullmat
  sample.mat[[l]] <- trans.matrix.sample
  ########################## Put the values of the hyperparameters #########################
  #d = nrow(trans.matrix.sample)
  n <- matrix(rep(0, d*d), nrow = d)
  dim_sample <- nrow(trans.matrix.sample)
  n[0:dim_sample, 0:dim_sample] <- trans.matrix.sample
  
  sample.prob.matrix <- n/rowSums(n)
  sample.prob.matrix[is.nan(sample.prob.matrix)] <- 0
  sample.prob.mat[[l]] <- sample.prob.matrix

  bg <- blocked_gibbs(N = N, n = n, Burn.in = Burn.in, M = M, d = d, alpha = alpha, 
                      beta = beta, b0 = b0)
  
  ################################# calculate posterior mean ###############################
  Pi_sum =  matrix(rep(0, d*d), nrow = d)
  Pi_error = rep(0, d)
  for(i in 1:M){
    Pi_sum = Pi_sum + bg[[i]]$Pi
    Pi_error[i] = sqrt(mean((bg[[i]]$Pi - trans.matrix.true)^2))
  }
  Pi_mean = Pi_sum/M
  Pi.error.mat[, l] = Pi_error 
  GSBP.mat[[l]] <- Pi_mean
  ################################################################################
  TrueVsMLE = trans.matrix.true[1:d,1:d] - sample.prob.matrix
  SampVsEst = sample.prob.matrix - Pi_mean
  TrueVsEst = trans.matrix.true[1:d,1:d] - Pi_mean
  
  MAE_MLE[l] <- mean(abs(TrueVsMLE))*100
  MAE_GSBP[l] <- mean(abs(TrueVsEst))*100
  ################################################################################
}
matplot(cbind(MAE_MLE, MAE_GSBP), type = "l")
data.frame(sample.vec, MAE_MLE, MAE_GSBP)
MAE_MLE

saveRDS(sample.mat, file = "sample.RData")
saveRDS(sample.prob.mat, file = "sampleProb.RData")
saveRDS(GSBP.mat, file = "GSBP.RData")
sample.mat[[1]]
GSBP.mat

sampleList <- lapply(sample.mat, as.data.frame)
names(sampleList) <- paste0("Sheet", seq_along(sampleList))

GSBPList <- lapply(GSBP.mat, as.data.frame)
names(GSBPList) <- paste0("Sheet", seq_along(GSBPList))

normalized_list <- lapply(sample.mat, function(mat) {
  t(apply(mat, 1, function(row) row / sum(row)))
})
padded_list <- lapply(normalized_list, function(mat) {
  # Current size
  r <- nrow(mat)
  c <- ncol(mat)
  # Create d x d zero matrix
  new_mat <- matrix(0, nrow = d, ncol = d)
  # Fill with existing values
  new_mat[1:r, 1:c] <- mat
  return(new_mat)
})
cleaned_list <- lapply(padded_list, function(mat) {
  mat[is.nan(mat)] <- 0
  return(mat)
})

sample.prob.list <- cleaned_list

sampleProbList <- lapply(sample.prob.list, as.data.frame)
names(sampleProbList) <- paste0("Sheet", seq_along(sampleList))

############################### GRAPH of Posterior Mean #########################
frobenius_dist <- function(A, B) sqrt(mean((A - B)^2))
distances1 <- sapply(sample.prob.list, function(mat) 
frobenius_dist(mat, trans.matrix.true))
# distances2 <- sapply(GSBP.mat, function(mat)
# frobenius_dist(mat, trans.matrix.true)) #alpha = 1
distances3 <- sapply(GSBP.mat, function(mat) #alpha = 0.5
frobenius_dist(mat, trans.matrix.true))

############################### Rate of Convergence ################################
distances1 <- distancesData$Distances1_MLE
distances2 <- distancesData$Distances2_HSBP
distances3 <- distancesData$Distances3_GHSBP
lm <- lm(log(distances1)~log(sample.vec))
coef_MLE <- -coef(lm)[2]
lm <- lm(log(distances2)~log(sample.vec))
coef_HSBP <- -coef(lm)[2]
lm <- lm(log(distances3)~log(sample.vec))
coef_GSBP <- -coef(lm)[2]
distances1 <- distancesData$Distances1_MLE * 100
distances2 <- distancesData$Distances2_HSBP * 100
distances3 <- distancesData$Distances3_GHSBP * 100
# --- Create publication-quality PDF ---
pdf("Figure_RMSE_vs_SampleSize.pdf", width = 6, height = 4.5)

# --- Adjust text and spacing for AoAS-style figure ---
par(
  mar = c(4.5, 5, 2, 1),   # more left margin for y label
  mgp = c(2.8, 0.8, 0),    # control label and tick spacing
  las = 1,                 # horizontal y-axis labels
  cex.axis = 1.1,          # axis label size
  cex.lab = 1.3,           # axis title size
  cex.main = 1.4,          # main title size
  tcl = -0.3               # tick length
)


# --- Main Plot ---
plot(sample.vec, distances1, type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "Sample Size",
     ylab = "RMSE (%)",
     ylim = range(c(distances1, distances2, distances3)),
     main = "")
lines(sample.vec, distances2, type = "b", pch = 17, col = "red", lwd = 2)
lines(sample.vec, distances3, type = "b", pch = 15, col = "darkgreen", lwd = 2)

legend("topright",
       legend = c("MLE", "HSBP", "GHSBP"),
       col = c("blue", "red", "darkgreen"),
       pch = c(19, 17, 15),
       lwd = 2,
       bty = "n",
       cex = 1.1)

dev.off()
###################################################################################
