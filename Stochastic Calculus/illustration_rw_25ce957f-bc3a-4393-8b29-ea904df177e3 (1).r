rm(list=ls())
gc()


n<-500
N<-100                        

walks<-replicate(N, c(0, cumsum(sample(c(-1,1), n, replace = TRUE))))

matplot(0:n, walks, type = "l", lty = 1,
        xlab = "k", ylab = "M_k",
        main = "Symmetric Random Walk")
abline(h = 0, lty = 2)

# Non-symmetric random walk

n<-500
N<-10000      
p<-0.49         


walks <- replicate(
  N,
  c(0, cumsum(sample(c(-1, 1), n, replace = TRUE, prob = c(1 - p, p))))
)

matplot(0:n, walks, type = "l", lty = 1,
        xlab = "k", ylab = "M_k",
        main = paste("Non-symmetric random walks (p =", p, ")"))
abline(h = 0, lty = 2)




# Scaled Random Walk

set.seed(42)


n <- 2000
X <- sample(c(-1, 1), n, replace = TRUE)
M <- cumsum(X)

t <- (1:n) / n
Wn <- M / sqrt(n)

plot(t, Wn, type = "l", col = "darkgreen", lwd = 2,
     xlab = "t", ylab = expression(W^{(n)}(t)),
     main = "Scaled Random Walk")
abline(h = 0, lty = 2)

# Martingale


k<-100         # time point of observation
L<-300       # Horizon  
N<-1000       # Number of future paths

X_base<-sample(c(-1, 1), L, replace = TRUE)
M_base<-c(0, cumsum(X_base))
Mk<-M_base[k + 1]
future_steps<-L - k

paths_future<-replicate(N, {
  future_X<-sample(c(-1, 1), future_steps, replace = TRUE)
  Mk+cumsum(future_X)
})


endpoints<-paths_future[future_steps, ]
mean_endpoint<-mean(endpoints)

plot(0:k,M_base[1:(k+1)],type = "l",lwd =3,xlim=c(0,L),ylim =range(c(M_base, paths_future)), xlab = "t", ylab = "M_t", main = "Martingale Property")
for (i in 1:N) {
  lines(k:L, c(Mk, paths_future[, i]),col= rgb(0, 0, 1, 0.15))
}
abline(h = mean_endpoint, col = "red", lwd = 3, lty = 2)
abline(v = k, lty = 3)
legend("topleft",legend = c("Observed path up to k","Future continuations", "E[M_L | F_k]"), col = c("black", "blue", "red"),
      lwd = c(3, 2, 3), lty = c(1, 1, 2), bty = "n")



# Brownian motion

T<-100
n<-1000
dt<-T/n
t<-seq(0, T, length.out = n + 1)

N<-1         


B_mat <- matrix(0, nrow = n + 1, ncol = N)

for (i in 1:N) {
  dB <- rnorm(n, mean = 0, sd = sqrt(dt))
  B_mat[, i] <- c(0, cumsum(dB))
}

matplot(t, B_mat, type = "l", lty = 1,
        xlab = "t", ylab = expression(B[t]),
        main = "Brownian Motion")
abline(h = 0, lty = 2)

# Geometric Brownian motion

mu<-0.2
sigma<-0.4
S0<-1

S_mat<-matrix(0, nrow = n + 1, ncol = N)
S_mat[1,]<-S0

for (i in 1:N) {
  dB<-diff(B_mat[, i])   
  for(k in 1:n) {
    S_mat[k+1,i]<-S_mat[k, i]*exp((mu-0.5*sigma^2)*dt+sigma*dB[k])
  }
}

matplot(t, S_mat, type = "l", lty = 1, xlab = "t", ylab = expression(S[t]), main = "Geometric Brownian Motion")
