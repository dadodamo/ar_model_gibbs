library("RProtoBuf")
setwd("/users/daniel/desktop/ar_gibbs")


my_proto_file <- "/proto/ydata.proto"  
proto_file <- RProtoBuf::readProtoFiles(files = "proto/ydata.proto")

file_path <- "cmake-build-debug/serialized_data.bin"  

binary_data <- readBin(file_path, "raw", file.info(file_path)$size)

y_file <- system.file(file_path, proto_file, package = 'RProtoBuf')

# example of messages for testing 

msg <- read(y_data.full_y, binary_data)
writeLines(as.character(msg))
as.matrix(msg$vec_t)

#read proto vector file
vector_proto_file <- "/proto/vector_it.proto"
vector_proto <- RProtoBuf::readProtoFiles(files = "proto/vector_it.proto")
matrix_proto <- RProtoBuf::readProtoFiles(files = "proto/o.proto")

## read beta file
file_path <- "cmake-build-debug/beta_serialized.bin"  
binary_data_beta <- readBin(file_path, "raw", file.info(file_path)$size)

## read o file
file_path <- "cmake-build-debug/beta_serialized.bin"  
binary_data_o <- readBin(file_path, "raw", file.info(file_path)$size)

msg_beta <- read(vector.full_iter_vec, binary_data_beta)
msg_o <- read(o_data.full_o_it, binary_data_o)
## writeLines(as.character(msg_beta$vec_t[1]))
list_beta <- as.list(msg_beta$vec_t)
list_o <- as.list(msg_o$m)

writeLines(as.character(msg_o[[1]]))
list_o[1]

betas <- sapply(list_beta, function(x){x$vec_value});
o <- sapply(list_o, function(x){x$vec})

par(mfrow = c(5, 1))
par(mar = c(2, 4, 2, 1))
plot(betas[1,], type='l')
plot(betas[2,], type='l')
plot(betas[3,], type='l')
plot(betas[4,], type='l')
plot(betas[5,], type='l')


