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

        // matern matrix
        Eigen::MatrixXd matern_inv;

        //beta prior
        double beta_sig_prior = 1;
        
        // rho 
        double rho_mean_prior = 0;
        double rho_sig_prior = 1;

        //mu_0 prior
        double mu0_sig_prior = 1;

        //inverse gamma group 
        std::pair<double, double> ab_eps_prior = {1,1};
        std::pair<double, double> ab_w_prior = {1,1};
        std::pair<double, double> ab_0_prior = {1,1};


    public:
        ar_model(unsigned int n, unsigned int T, 
        std::vector<Eigen::VectorXd>& y_store,
        std::vector<Eigen::MatrixXd>& x_store,
        std::vector<coord>& coord_vec,
        std::vector<Eigen::VectorXd> ot_store,
        Eigen::VectorXd beta,
        Eigen::VectorXd mu_0,
        double rho):
        y(y_store), X(x_store),
        coordinates(coord_vec),
        ot(ot_store),
        beta_true(beta),
        mu_0_true(mu_0),
        rho_true(rho),
        n_iter(n), T(T)
        {
            N = (X)[0].rows();
            p = (X)[0].cols();

            Eigen::MatrixXd matern_cov = Eigen::MatrixXd::Zero(N,N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    double dist = eucl_dist(coordinates[i], coordinates[j]) ;//eucl_dist((*coordinates)[i], (*coordinates)[j]);
                    matern_cov(i, j) = matern(dist, phi, nu);
                 }
            }
            std::cout << matern_cov << std::endl;
            check_eigenvalues(matern_cov);

            matern_inv = matern_cov.inverse();
        };
        
        void sample() const;
};

#endif