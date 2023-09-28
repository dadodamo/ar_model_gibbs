library("RProtoBuf")
library("sys")
library("system2")
setwd("/users/daniel/desktop/ar_gibbs")


my_proto_file <- "/proto/ydata.proto"  
proto_file <- RProtoBuf::readProtoFiles(files = "proto/ydata.proto")

file_path <- "cmake-build-debug/serialized_data.bin"  

binary_data <- readBin(file_path, "raw", file.info(file_path)$size)

y_file <- system.file(file_path, proto_file, package = 'RProtoBuf')

# example of messages for testing 

msg <- read(y_data.full_y, binary_data)
writeLines(as.character(msg))
as.matrix(msg$vec_t[[1]])

#read proto vector file
vector_proto_file <- "/proto/vector_it.proto"
vector_proto <- RProtoBuf::readProtoFiles(files = "proto/vector_it.proto")

## read beta file
file_path <- "cmake-build-debug/beta_serialized.bin"  
binary_data <- readBin(file_path, "raw", file.info(file_path)$size)

msg <- read(vector.full_iter_vec, binary_data)
writeLines(as.character(msg))

function read_bin_to_matrix(binary_file_path, proto_file_path) {
  binary_data <- readBin(binary_file_path, "raw", file.info(file_path)$size)
}
