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



class ar_model {
    private:
        //iterations
        unsigned int n_iter;
        //dimensions
        unsigned int T;
        unsigned int N;
        unsigned int p;

        // data passed to class.
        std::vector<Eigen::VectorXf> y;
        std::vector<Eigen::MatrixXf> X;
        std::vector<coord> coordinates;
        float phi = 1;
        float nu = 0.5;

        // matern matrix
        Eigen::MatrixXf matern_inv;

        //beta prior
        float beta_sig_prior = 1;
        
        // rho 
        float rho_mean_prior = 0;
        float rho_sig_prior = 1;

        //mu_0 prior
        float mu0_sig_prior = 1;

        //inverse gamma group 
        std::pair<float, float> ab_eps_prior = {1,1};
        std::pair<float, float> ab_w_prior = {1,1};
        std::pair<float, float> ab_0_prior = {1,1};


    public:
        ar_model(unsigned int n, unsigned int T, 
        std::vector<Eigen::VectorXf>& y_store,
        std::vector<Eigen::MatrixXf>& x_store,
        std::vector<coord>& coord_vec):
        y(y_store), X(x_store), coordinates(coord_vec),
        n_iter(n), T(T)
        {
            N = (X)[0].rows();
            p = (X)[0].cols();

            Eigen::MatrixXf matern_cov = Eigen::MatrixXf::Zero(N,N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    float dist = eucl_dist(coordinates[i], coordinates[j]) ;//eucl_dist((*coordinates)[i], (*coordinates)[j]);
                    matern_cov(i, j) = matern(dist, phi, nu);
                 }
            }
            matern_inv = matern_cov.inverse();
        };
        
        void sample() const;
};

#endif