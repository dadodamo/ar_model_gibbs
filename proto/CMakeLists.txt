INCLUDE(FindProtobuf)
FIND_PACKAGE(Protobuf REQUIRED)
INCLUDE_DIRECTORIES(${PROTOBUF_INCLUDE_DIR})
PROTOBUF_GENERATE_CPP(PROTO_SRC PROTO_HEADER
        ydata.proto
        xdata.proto
        vector_it.proto
        scalar_it.proto
        o.proto)
ADD_LIBRARY(proto ${PROTO_HEADER} ${PROTO_SRC} ../protocpp/serialize.h ../protocpp/serialize.cpp)