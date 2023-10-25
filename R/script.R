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


par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
plot(betas[1,10:1000], type='l')
plot(betas[2,10:1000], type='l')
plot(betas[3,10:1000], type='l')
plot(betas[4,10:1000], type='l')
plot(betas[5,10:1000], type='l')
dev.off()
## o's
## read o file
matrix_proto <- RProtoBuf::readProtoFiles(files = "proto/o.proto")
file_path_o <- "cmake-build-debug/o_serialized.bin"  
binary_data_o <- readBin(file_path_o, "raw", file.info(file_path_o)$size)

msg_o <- read(o_data.full_o_it, binary_data_o)
list_o <- as.list(msg_o$m)

list_o <- as.list(msg_o$m)

writeLines(as.character(msg_o[1]))
o <- sapply(list_o, function(x){x$vec})
o <- sapply(o, function(x){x$vec_value})



#read proto scalar file
scalar_proto <- RProtoBuf::readProtoFiles(files = "proto/scalar_it.proto")

## read beta file
file_path_rho <- "cmake-build-debug/rho_serialized.bin"  
binary_data_rho <- readBin(file_path_rho, "raw", file.info(file_path_rho)$size)

msg_rho <- read(scalar.full_scalar_it, binary_data_rho)

list_rho <- as.list(msg_rho$scalar)

rho <- sapply(list_rho, function(x){x$value});
plot(rho[20:1000], type = 'l')
