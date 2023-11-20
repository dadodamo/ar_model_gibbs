#ifndef AR_GIBBS_AR_CLASS_H
#define AR_GIBBS_AR_CLASS_H

#include <iostream>
#include "Eigen/Dense"
#include "eigenmvn.h"
#include "../coordinates/coordinates.h"
#include "../matern.h"
#include "../cmake-build-debug/proto/ydata.pb.h"
#include "../cmake-build-debug/proto/o.pb.h"
#include "../cmake-build-debug/proto/scalar_it.pb.h"
#include "../cmake-build-debug/proto/vector_it.pb.h"
#include <fstream>
#include "../calc_posterior/posterior.h"
#include <random>
#include<string>
#include "../debug_functions/debug.h"
#include "../protocpp/serialize.h"

// DEBUG MODE


class ar_model {
    private:
        //iterations
        unsigned int n_iter;
        //dimensions
        unsigned int T;
        unsigned int N;
        unsigned int p;

        // data passed to class.
        std::vector<Eigen::VectorXd> y;
        std::vector<Eigen::MatrixXd> X;
        std::vector<coord> coordinates;

        //debug
        std::vector<Eigen::VectorXd> ot;
        Eigen::VectorXd beta_true;
        double rho_true;
        Eigen::VectorXd mu_0_true;

        Eigen::VectorXd beta;
        Eigen::VectorXd mu_0;

        double phi;
        double nu;

        double sigma_eps_true;
        double sigma_w_true;
        double sigma_0_true;

        // matern matrix
        Eigen::MatrixXd matern_inv;


        // priors chosen as in spTimer paper suggested
        //beta prior
        double beta_sig_prior = 1000000.;
        
        // rho
        double rho_sig_prior = 1000000.;

        //mu_0 prior
        double mu0_sig_prior = 100000.;

        //inverse gamma group 
        std::pair<double, double> ab_eps_prior = {2,1};
        std::pair<double, double> ab_w_prior = {2,1};
        std::pair<double, double> ab_0_prior = {2,1};

        // algo options
        bool use_cholesky;
        u_int64_t seed;



    public:
        ar_model(unsigned int n, unsigned int T, 
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
        double& phi,
        double& nu,
        bool use_cholesky = false,
        u_int64_t seed = 1);

        void sample() const;
};

#endif