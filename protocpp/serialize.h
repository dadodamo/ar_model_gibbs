//
// Created by Daniel Adamovic on 28/10/23.
//

#include<iostream>
#include "../cmake-build-debug/proto/ydata.pb.h"
#include "../cmake-build-debug/proto/o.pb.h"
#include "../cmake-build-debug/proto/scalar_it.pb.h"
#include "../cmake-build-debug/proto/vector_it.pb.h"
#include <fstream>


#ifndef AR_GIBBS_SERIALIZE_H
#define AR_GIBBS_SERIALIZE_H

namespace proto {
    void serialize_o(o_data::full_o_it& o_stream);
    void serialize_beta(vector::full_iter_vec& beta_stream);
    void serialize_mu0(vector::full_iter_vec& mu_0_stream);
    void serialize_rho(scalar::full_scalar_it& rho_stream);
    void serialize_sig_0(scalar::full_scalar_it& sig_0_stream);
    void serialize_sig_w(scalar::full_scalar_it& sig_w_stream);
    void serialize_sig_eps(scalar::full_scalar_it& sig_eps_stream);

}

#endif //AR_GIBBS_SERIALIZE_H
