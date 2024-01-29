#ifndef AR_GIBBS_AR_CLASS_H
#define AR_GIBBS_AR_CLASS_H

#include <iostream>
#include "Eigen/Dense"
#include "eigenmvn.h"
#include "../coordinates/coordinates.h"
#include "../matern.h"
#include "../cmake-build-debug/proto/ydata.pb.h"
#include "../cmake-build-debug/proto/paramdata.pb.h"
#include "../protocpp/serialize.h"
#include <fstream>
#include<string>
#include "../calc_posterior/posterior.h"
#include <random>
#include "../debug_functions/debug.h"

// DEBUG MODE


class ar_model {
private:
    //dimensions
    unsigned int T;
    unsigned int N;
    unsigned int p;

    // data passed to class
    std::vector<Eigen::VectorXd> y;
    std::vector<Eigen::MatrixXd> X;
    std::vector<coord> coordinates;
    double nu;

    //debug
    std::vector<Eigen::VectorXd> o_store_true;
    Eigen::VectorXd beta_true;
    double rho_true;
    Eigen::VectorXd mu_0_true;
    double phi_true;
    double sigma_eps_true;
    double sigma_w_true;
    double sigma_0_true;

    //sampled parameters
    std::vector<Eigen::VectorXd> o_store;
    Eigen::VectorXd beta;
    double rho;
    Eigen::VectorXd mu_0;
    double phi;
    double sigma_eps;
    double sigma_w;
    double sigma_0;


    // priors chosen
    double beta_sig_prior = 1;
    double rho_sig_prior = 1.;
    double mu0_sig_prior = 1.;

    //inverse gamma group
    std::pair<double, double> ab_eps_prior = {2,4};
    std::pair<double, double> ab_w_prior = {2,4};
    std::pair<double, double> ab_0_prior = {2,4};
    // phi prior and candidate variance
    std::pair<double, double> ab_phi_prior = {8,  4};
    double phi_cand_var = 0.1;

    // matrices
    Eigen::MatrixXd coord_mat;
    Eigen::MatrixXd matern_cov;
    Eigen::MatrixXd matern_inv;
    Eigen::MatrixXd w_full_cov_inv;

    // algo options
    bool use_cholesky = true;
    u_int64_t seed = 1;
    std::mt19937 generator = std::mt19937(1);

    Eigen::VectorXd std_mean_y;
    double std_var_y;
    std::vector<Eigen::VectorXd> std_mean_X;
    std::vector<double> std_var_X;
    double std_mean_coord;
    double std_var_coord;
    double phi_accept_rate = 0;
    double iter_count = 0;



    //samplers
    Eigen::EigenMultivariateNormal<double> o_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(1), Eigen::MatrixXd::Identity(1,1), false, 1);
    Eigen::EigenMultivariateNormal<double> beta_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(1),Eigen::MatrixXd::Identity(1,1), false, 1);
    Eigen::EigenMultivariateNormal<double> mu0_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(1), Eigen::MatrixXd::Identity(1,1), false, 1);
    std::normal_distribution<double> rho_sampler = std::normal_distribution<double>(0,1);
    std::gamma_distribution<double> sig_eps_sampler = std::gamma_distribution<double>(1, 1);
    std::gamma_distribution<double> sig_w_sampler = std::gamma_distribution<double>(1,1);
    std::gamma_distribution<double> sig_0_sampler = std::gamma_distribution<double>(1,1);
    std::normal_distribution<double> phi_sampler = std::normal_distribution<double>(0,1);
    std::uniform_real_distribution<double> unif = std::uniform_real_distribution<double>(0,1);

    sampler_data::samples sample_stream;

public:
    ar_model(
            std::vector<Eigen::VectorXd>& y_store,
            std::vector<Eigen::MatrixXd>& x_store,
            std::vector<coord>& coord_vec,
            std::vector<Eigen::VectorXd>& ot_store,
            Eigen::VectorXd& beta,
            Eigen::VectorXd& mu_0,
            double& rho,
            double& sigma_eps_true,
            double& sigma_w_true,
            double& sigma_0_true,
            double& phi_true,
            double& nu
    );

    void init();
    void standardize();
    void sample();
    void write_curr_state();

    void serialize();
    double get_acceptance_rate(){
        return phi_accept_rate/iter_count;
    }
};

#endif