//
// Created by Daniel Adamovic on 28/10/23.
//

#include<iostream>
#include "../cmake-build-debug/proto/ydata.pb.h"
#include "../cmake-build-debug/proto/paramdata.pb.h"
#include <fstream>


#ifndef AR_GIBBS_SERIALIZE_H
#define AR_GIBBS_SERIALIZE_H

namespace proto {
    void serialize(sampler_data::samples& sample_stream);
    void serialize_y(y_data::full_y& y_stream);
}

#endif //AR_GIBBS_SERIALIZE_H
