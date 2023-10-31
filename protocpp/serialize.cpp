//
// Created by Daniel Adamovic on 28/10/23.
//
#include "serialize.h"

void proto::serialize_o(o_data::full_o_it& o_stream){
    std::string serialized_o;
    if (!o_stream.SerializeToString(&serialized_o)) {
        std::cerr << "Failed to write random effects data." << std::endl;
    }
    std::ofstream outputFile_o("o_serialized.bin", std::ios::binary);
    if (outputFile_o.is_open()) {
        outputFile_o.write(serialized_o.c_str(), serialized_o.size());
        outputFile_o.close();
        std::cout << "Serialized data written to o_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
};
void proto::serialize_beta(vector::full_iter_vec& beta_stream){
    std::string serialized_beta;
    if (!beta_stream.SerializeToString(&serialized_beta)) {
        std::cerr << "Failed to write beta data." << std::endl;
    }
    std::ofstream outputFile_beta("beta_serialized.bin", std::ios::binary);
    if (outputFile_beta.is_open()) {
        outputFile_beta.write(serialized_beta.c_str(), serialized_beta.size());
        outputFile_beta.close();
        std::cout << "Serialized data written to beta_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
}
;

void proto::serialize_mu0(vector::full_iter_vec& mu_0_stream){
    std::string serialized_mu0;
    if (!mu_0_stream.SerializeToString(&serialized_mu0)) {
        std::cerr << "Failed to write mu_0 data." << std::endl;
    }
    std::ofstream outputFile_mu0("mu0_serialized.bin", std::ios::binary);
    if (outputFile_mu0.is_open()) {
        outputFile_mu0.write(serialized_mu0.c_str(), serialized_mu0.size());
        outputFile_mu0.close();
        std::cout << "Serialized data written to mu0_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

};
void proto::serialize_rho(scalar::full_scalar_it& rho_stream){
    std::string serialized_rho;
    if (!rho_stream.SerializeToString(&serialized_rho)) {
        std::cerr << "Failed to write rho data." << std::endl;
    }
    std::ofstream outputFile_rho("rho_serialized.bin", std::ios::binary);
    if (outputFile_rho.is_open()) {
        outputFile_rho.write(serialized_rho.c_str(), serialized_rho.size());
        outputFile_rho.close();
        std::cout << "Serialized data written to rho_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }


};
void proto::serialize_sig_0(scalar::full_scalar_it& sig_0_stream){
    std::string serialized_sig_0;
    if (!sig_0_stream.SerializeToString(&serialized_sig_0)) {
        std::cerr << "Failed to write sigma_0 data." << std::endl;
    }
    std::ofstream outputFile_sig_0("sig_0_serialized.bin", std::ios::binary);
    if (outputFile_sig_0.is_open()) {
        outputFile_sig_0.write(serialized_sig_0.c_str(), serialized_sig_0.size());
        outputFile_sig_0.close();
        std::cout << "Serialized data written to sig_0_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

};
void proto::serialize_sig_w(scalar::full_scalar_it& sig_w_stream){
    std::string serialized_sig_w;
    if (!sig_w_stream.SerializeToString(&serialized_sig_w)) {
        std::cerr << "Failed to write sigma_w data." << std::endl;
    }
    std::ofstream outputFile_sig_w("sig_w_serialized.bin", std::ios::binary);
    if (outputFile_sig_w.is_open()) {
        outputFile_sig_w.write(serialized_sig_w.c_str(), serialized_sig_w.size());
        outputFile_sig_w.close();
        std::cout << "Serialized data written to sig_w_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

};
void proto::serialize_sig_eps(scalar::full_scalar_it& sig_eps_stream){
    std::string serialized_sig_eps;
    if (!sig_eps_stream.SerializeToString(&serialized_sig_eps)) {
        std::cerr << "Failed to write sigma_eps data." << std::endl;
    }
    std::ofstream outputFile_sig_eps("sig_eps_serialized.bin", std::ios::binary);
    if (outputFile_sig_eps.is_open()) {
        outputFile_sig_eps.write(serialized_sig_eps.c_str(), serialized_sig_eps.size());
        outputFile_sig_eps.close();
        std::cout << "Serialized data written to sig_eps_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
};