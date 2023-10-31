library("RProtoBuf")
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


## beta only sample  
dev.off()
par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
plot(betas[1,2000:3000], type='l' ,main = "betas", col = 'grey')
lines(1:1000, rep(1, 1000), col = "green", type = 'l')
plot(betas[2,2000:3000], type='l', col = 'grey')
lines(1:1000, rep(2, 1000), col = "green", type = 'l')
plot(betas[3,2000:3000], type='l', col = 'grey')
lines(1:1000, rep(3, 1000), col = "green", type = 'l')
plot(betas[4,2000:3000], type='l', col = 'grey')
lines(1:1000, rep(4, 1000), col = "green", type = 'l')
plot(betas[5,2000:3000], type='l', col = 'grey')
lines(1:1000, rep(5, 1000), col = "green", type = 'l')

dev.off()
par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
plot(betas[1,2000:3000], type='l', ylim = c(0.6, 1.1))
lines(1:1000, rep(1, 1000), col = "green", type = 'l')
plot(betas[2,2000:3000], type='l', ylim = c(1.86, 2.1))
lines(1:1000, rep(2, 1000), col = "green", type = 'l')
plot(betas[3,2000:3000], type='l', ylim = c(2.88, 3.1))
lines(1:1000, rep(3, 1000), col = "green", type = 'l')
plot(betas[4,2000:3000], type='l', ylim = c(3.95, 4.11))
lines(1:1000, rep(4, 1000), col = "green", type = 'l')
plot(betas[5,2000:3000], type='l',ylim = c(4.59, 5.01))
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
plot(rho[50:3000], type = 'l')
plot(rho[2000:3000], type = 'l', col = 'grey')
lines(1:1000, rep(0.8, 1000), type = 'l', col = 'green')
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
for (i in 1:10) {
  plot(mu0[i,1:3000], main = paste("mu ", i) ,  type = 'l', col = 'grey');
  lines(1:3000, rep(i, 3000), type = 'l', col = "green");
}


##### plot of variance components

## sig_eps
dev.off()
file_path_sig_eps <- "cmake-build-debug/sig_eps_serialized.bin"  
binary_data_sig_eps <- readBin(file_path_sig_eps, "raw", file.info(file_path_sig_eps)$size)

msg_sig_eps <- read(scalar.full_scalar_it, binary_data_sig_eps)

list_sig_eps <- as.list(msg_sig_eps$scalar)

sig_eps <- sapply(list_sig_eps, function(x){x$value});
plot(sig_eps[2000:3000], type = 'l', col = 'grey', main = "sig_eps true value : 2")
lines(1:1000, rep(2, 1000), type = 'l', col = "green")
mean(sig_eps)

plot(seq(0,10,0.01), dinvgamma(seq(0,10,0.01), shape = 0.1, rate = 1), type = 'l')

## sig_w
dev.off()
file_path_sig_w <- "cmake-build-debug/sig_w_serialized.bin"  
binary_data_sig_w <- readBin(file_path_sig_w, "raw", file.info(file_path_sig_w)$size)

msg_sig_w <- read(scalar.full_scalar_it, binary_data_sig_w)

list_sig_w <- as.list(msg_sig_w$scalar)

sig_w <- sapply(list_sig_w, function(x){x$value});
plot(sig_w[2000:3000], type = 'l',col = "grey", main = "true value = 5")
lines(1:1000, rep(5, 1000), type = 'l', col = "green")
mean(sig_w)

## sig_0
dev.off()
file_path_sig_0 <- "cmake-build-debug/sig_0_serialized.bin"  
binary_data_sig_0 <- readBin(file_path_sig_0, "raw", file.info(file_path_sig_0)$size)

msg_sig_0 <- read(scalar.full_scalar_it, binary_data_sig_0)

list_sig_0 <- as.list(msg_sig_0$scalar)

sig_0 <- sapply(list_sig_0, function(x){x$value});
plot(sig_0[2000:3000],col = 'grey' ,type = 'l', main = "true value = 1")
lines(1:1000, rep(1, 1000), type = 'l', col = "green")
mean(sig_0)


#### o's
matrix_proto <- RProtoBuf::readProtoFiles(files = "proto/o.proto")
file_path_o <- "cmake-build-debug/o_serialized.bin"  
binary_data_o <- readBin(file_path_o, "raw", file.info(file_path_o)$size)
msg_o <- read(o_data.full_o_it, binary_data_o)
list_o <- as.list(msg_o$m)
list_o <- as.list(msg_o$m)
o <- sapply(list_o, function(x){x$vec});
o <- sapply(o, function(x){x$vec_value});

o0 <- matrix(0, ncol = 3000, nrow = 10);
col = 1;
for(i in seq(1, 149951, 50)){
  o0[,col] = o[,i];
  col = col + 1;
}
plot(o0[1,], type = "l", col = "grey")  
plot(o0[2,2000:3000], type = "l")
plot(o0[3,2000:3000], type = "l")
plot(o0[4,2000:3000], type = "l")
plot(o0[5,2000:3000], type = "l")








