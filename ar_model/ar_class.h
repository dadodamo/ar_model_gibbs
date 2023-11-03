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
        double rho;
        Eigen::VectorXd mu_0;

        double phi = 1.;
        double nu = 1.;

        double sigma_eps_true;
        double sigma_w_true;
        double sigma_0_true;

        // matern matrix
        Eigen::MatrixXd matern_inv;

        //beta prior
        double beta_sig_prior = 1;
        
        // rho 
        double rho_mean_prior = 0;
        double rho_sig_prior = 1;

        //mu_0 prior
        double mu0_sig_prior = 100;

        //inverse gamma group 
        std::pair<double, double> ab_eps_prior = {1,1};
        std::pair<double, double> ab_w_prior = {1,1};
        std::pair<double, double> ab_0_prior = {1,1};

        // algo options
        bool use_cholesky;
        u_int64_t seed;



    public:
        ar_model(unsigned int n, unsigned int T, 
        std::vector<Eigen::VectorXd>& y_store,
        std::vector<Eigen::MatrixXd>& x_store,
        std::vector<coord>& coord_vec,
        std::vector<Eigen::VectorXd> ot_store,
        Eigen::VectorXd beta,
        Eigen::VectorXd mu_0,
        double rho,
        double sigma_eps_true,
        double sigma_w_true,
        double sigma_0_true,
        bool use_cholesky = false,
        u_int64_t seed = 1):
        y(y_store), X(x_store),
        coordinates(coord_vec),
        ot(ot_store),
        beta_true(beta),
        mu_0_true(mu_0),
        rho_true(rho),
        sigma_eps_true(sigma_eps_true),
        sigma_w_true(sigma_w_true),
        sigma_0_true(sigma_0_true),
        n_iter(n), T(T),
        use_cholesky(use_cholesky), seed(seed)
        {
            N = (X)[0].rows();
            p = (X)[0].cols();

            Eigen::MatrixXd matern_cov = Eigen::MatrixXd::Zero(N,N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    double dist = eucl_dist(coordinates[i], coordinates[j]) ;
                    matern_cov(i, j) = matern(dist, phi, nu);
                 }
            }
            Eigen::MatrixXd id_N = Eigen::VectorXd::Ones(N).asDiagonal();
            matern_inv = matern_cov.fullPivLu().solve(id_N);
        };
        void sample() const;
};

#endif