library( "RProtoBuf")
setwd("/users/daniel/desktop/ar_gibbs")


#read proto vector file
vector_proto <- RProtoBuf::readProtoFiles(files = "proto/vector_it.proto")

## read beta file
file_path_beta <- "cmake-build-debug/beta_serialized.bin"  
binary_data_beta <- readBin(file_path_beta, "raw", file.info(file_path_beta)$size)

msg_beta <- read(vector.full_iter_vec, binary_data_beta)

## writeLines(as.character(msg_beta$vec_t[1]))
list_beta <- as.list(msg_beta$vec_t)

betas <- sapply(list_beta, function(x){x$vec_value});

beta_true <- c(
  0.52218,
  -1.09668,
  -0.415183,
  -1.27539,
  -0.032302
)

## beta only sample  
dev.off()
par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
for (i in 1:5) {
  plot(betas[i,1000:5000], type='l' ,main = "betas", col = 'grey')
  lines(1:4000, rep(beta_true[i], 4000), col = "green", type = 'l')  
}

plot(betas[1,1000:5000], type='l' ,main = "betas", col = 'grey')
lines(1:4000, rep(beta_true[1], 4000), col = "green", type = 'l')
#plot(betas[2,1000:5000], type='l', col = 'grey')
#lines(1:4000, rep(2, 4000), col = "green", type = 'l')
#plot(betas[3,1000:5000], type='l', col = 'grey')
#lines(1:4000, rep(3, 4000), col = "green", type = 'l')
#plot(betas[4,1000:5000], type='l', col = 'grey')
#lines(1:4000, rep(4, 4000), col = "green", type = 'l')
#plot(betas[5,1000:5000], type='l', col = 'grey')
#lines(1:4000, rep(5, 4000), col = "green", type = 'l')

dev.off()
par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
plot(betas[1,], type='l')
lines(1:1000, rep(1, 1000), col = "green", type = 'l')
plot(betas[2,], type='l')
lines(1:1000, rep(2, 1000), col = "green", type = 'l')
plot(betas[3,], type='l')
lines(1:1000, rep(3, 1000), col = "green", type = 'l')
plot(betas[4,], type='l')
lines(1:1000, rep(4, 1000), col = "green", type = 'l')
plot(betas[5,], type='l')
lines(1:1000, rep(5, 1000), col = "green", type = 'l')

par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
plot(betas[1,50:3000], type='l', ylab = expression(beta[1]))
plot(betas[2,50:3000], type='l', ylab = expression(beta[2]))
plot(betas[3,50:3000], type='l', ylab = expression(beta[3]))
plot(betas[4,50:3000], type='l', ylab = expression(beta[4]))
plot(betas[5,50:3000], type='l', ylab = expression(beta[5]))


mean(betas[1,2000:3000])
mean(betas[2,2000:3000])
mean(betas[3,2000:3000])
mean(betas[4,2000:3000])
mean(betas[5,2000:3000])

dev.off()


###### rho
#read proto scalar file
scalar_proto <- RProtoBuf::readProtoFiles(files = "proto/scalar_it.proto")

## read rho file
file_path_rho <- "cmake-build-debug/rho_serialized.bin"  
binary_data_rho <- readBin(file_path_rho, "raw", file.info(file_path_rho)$size)

msg_rho <- read(scalar.full_scalar_it, binary_data_rho)

list_rho <- as.list(msg_rho$scalar)

rho <- sapply(list_rho, function(x){x$value});
dev.off()
plot(rho[1000:5000], type = 'l', col = 'grey')
lines(1:4000, rep(0.5, 4000), type = 'l', col = 'green')
mean(rho[2000:3000])



#mu_0

## read file
file_path_mu0 <- "cmake-build-debug/mu0_serialized.bin"  
binary_data_mu0 <- readBin(file_path_mu0, "raw", file.info(file_path_mu0)$size)

msg_mu0 <- read(vector.full_iter_vec, binary_data_mu0)

list_mu0 <- as.list(msg_mu0$vec_t)

mu0 <- sapply(list_mu0, function(x){x$vec_value});
dev.off()
par(mfrow = c(5, 2))
par(mar = c(2, 4, 2, 1))
mu_0_true <- c(
  1.44256,
  0.540784,
  -0.00554277,
  -1.18828,
  -0.340807,
  0.815028,
  -1.37403,
  -0.486838,
  -0.272475,
  1.639
)
for (i in 1:10) {
  plot(mu0[i,1000:5000], main = paste("mu ", i) ,  type = 'l', col = 'grey');
  lines(1:4000, rep(mu_0_true[i], 4000), type = 'l', col = "green");
}
mean(mu0[2,])


##### plot of variance components

## sig_eps
dev.off()
file_path_sig_eps <- "cmake-build-debug/sig_eps_serialized.bin"  
binary_data_sig_eps <- readBin(file_path_sig_eps, "raw", file.info(file_path_sig_eps)$size)

msg_sig_eps <- read(scalar.full_scalar_it, binary_data_sig_eps)

list_sig_eps <- as.list(msg_sig_eps$scalar)

sig_eps <- sapply(list_sig_eps, function(x){x$value});
dev.off()
plot(sig_eps[1000:5000], type = 'l', col = 'grey', main = "sig_eps true value : 1")
lines(1:4000, rep(1, 4000), type = 'l', col = "green")
var(sig_eps)
mean(sig_eps)


## sig_w
dev.off()
file_path_sig_w <- "cmake-build-debug/sig_w_serialized.bin"  
binary_data_sig_w <- readBin(file_path_sig_w, "raw", file.info(file_path_sig_w)$size)

msg_sig_w <- read(scalar.full_scalar_it, binary_data_sig_w)

list_sig_w <- as.list(msg_sig_w$scalar)

sig_w <- sapply(list_sig_w, function(x){x$value});
plot(sig_w[1000:5000], type = 'l', main = "true value = 1", col ='grey')
lines(1:4000, rep(1, 4000), type = 'l', col = "green")
var(sig_w)
mean(sig_w)

## sig_0
dev.off()
file_path_sig_0 <- "cmake-build-debug/sig_0_serialized.bin"  
binary_data_sig_0 <- readBin(file_path_sig_0, "raw", file.info(file_path_sig_0)$size)

msg_sig_0 <- read(scalar.full_scalar_it, binary_data_sig_0)


list_sig_0 <- as.list(msg_sig_0$scalar)

sig_0 <- sapply(list_sig_0, function(x){x$value});
dev.off()
plot(sig_0[1000:5000],col = 'grey' ,type = 'l', main = "true value = 1, est. mean = 1.0053", ylim =c(0,10))
lines(1:4000, rep(1, 4000), type = 'l', col = "green")
mean(sig_0)


x = seq(0,100,0.01);
dev.off()
plot(x, dinvgamma(x,1, 100), type = 'l', xlim = c(0,10))


### phi
file_path_phi <- "cmake-build-debug/phi_serialized.bin"  
binary_data_phi <- readBin(file_path_phi, "raw", file.info(file_path_phi)$size)

msg_phi <- read(scalar.full_scalar_it, binary_data_phi)

list_phi <- as.list(msg_phi$scalar)

phi <- sapply(list_phi, function(x){x$value});

plot(phi[1000:5000], type = 'l', col = 'grey')



#### o's
matrix_proto <- RProtoBuf::readProtoFiles(files = "proto/o.proto")
file_path_o <- "cmake-build-debug/o_serialized.bin"  
binary_data_o <- readBin(file_path_o, "raw", file.info(file_path_o)$size)
msg_o <- read(o_data.full_o_it, binary_data_o)
list_o <- as.list(msg_o$m)
o <- sapply(list_o, function(x){x$vec});
o <- sapply(o, function(x){x$vec_value});

o0 <- matrix(0, ncol = 5000, nrow = 10);
col = 1;
for(i in seq(11,999811 , 200)){
  o0[,col] = o[,i];
  col = col + 1;
}
dev.off()
plot(o0[1,], type = "l")  
mean(o0[1,])
plot(o0[2,], type = "l")
mean(o0[2,])
plot(o0[3,], type = "l")
mean(o0[3,])
plot(o0[4,], type = "l")
mean(o0[4,])
plot(o0[5,], type = "l")
mean(o0[5,])
plot(o0[6,], type = "l")
mean(o0[6,])
plot(o0[7,], type = "l")
mean(o0[7,])
plot(o0[8,], type = "l")
mean(o0[8,])
plot(o0[9,], type = "l")
mean(o0[9,])
plot(o0[10,], type = "l")
mean(o0[10,])

mean <- c(
  0.790282,
  2.44915,
  3.55486,
  5.77145,
  5.61658,
  5.07802,
  7.88179,
  8.32403,
  9.91118,
  10.2209
)
mean2 <- c(
  -7.14982,
  6.82703,
  -1.31385,
  -9.71581,
  9.54782,
  0.260374,
  -8.24077,
  -1.04074,
  3.25668,
  -16.2563
)
par(mfrow = c(5, 2))
par(mar = c(2, 4, 2, 1))
for (i in 1:10) {
  plot(o0[i,1000:5000], main = paste("o10 ", i) ,  type = 'l', col = 'grey');
  lines(1:4000, rep(mean2[i], 4000), type = 'l', col = "green");
}
