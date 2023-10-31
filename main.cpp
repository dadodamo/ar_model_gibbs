#include <iostream>
#include "Eigen/Dense"
// needed maybe later
#include<boost/math/distributions.hpp>
#include<boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "boost/random.hpp"
//
#include "eigenmvn.h"
#include "matern.h"
#include "calc_posterior/posterior.h"
#include "chrono"
#include "coordinates/coordinates.h"
#include "ar_model/ar_class.h"
#include<random>
#include"debug_functions/debug.h"



// Source files

int main(int argc,char* argv[]) {

    // data generation of model using T = 10, N = 10, p = 5
    // can be rewritten with generic param, as class+methods or function
    const unsigned int T = 50;
    const unsigned int p = 5;
    const unsigned int N = 10;


    // all fixed parameters
    double phi = 1.0;
    double nu = 1;
    double rho = 0.8;

    //
    Eigen::VectorXd beta(p) ;
    beta << 1, 2, 3, 4, 5;

    Eigen::VectorXd mu_0(N);
    for (int i = 0; i <N ; ++i) {
        mu_0(i) = i + 2;
    }

    // model variance components
    double sigma_eps = 2.;
    double sigma_w = 5.;
    double sigma_0 = 1.0;

    //priors, to pass if needed
    double sigma_mu_0_prior = 1;

    double sig_eps_a_prior = 0.5;
    double sig_eps_b_prior = 0.5;

    double sig_w_a_prior = 0.5;
    double sig_w_b_prior = 0.5;

    double sig_0_a_prior = 0.5;
    double sig_0_b_prior = 0.5;



////// DATA GENERATION ///////
    //X_t generation

    std::vector<Eigen::MatrixXd> xt_store_vec(T);

    // sampler from external header file (OK)
    {
        Eigen::VectorXd mean_X(N);
        Eigen::MatrixXd covar_X(N, N);
        for (int i = 0; i < N; ++i) {
            mean_X(i) = 0;
            covar_X(i, i) = 1;
        }
        Eigen::EigenMultivariateNormal<double> normal_sampler(mean_X, covar_X);
        for (int t = 0; t <= T-1; ++t) {
            Eigen::MatrixXd X(N, p);
            X = normal_sampler.samples(p);
            xt_store_vec[t] = X;
        }
    }

    // w_t generation + generate coord points
    Eigen::MatrixXd covar_w(N, N);
    // calculation of covar matrix
    std::vector<coord> coord_store_vec(10);
    {
        // hard code 10 coords
        coord c1(10, 40);
        coord c2(16, 40);
        coord c3(17, 50);
        coord c4(5, 10);
        coord c5(5, 40);
        coord c6(9, 20);
        coord c7(15, 25);
        coord c8(40, 10);
        coord c9(35, 5);
        coord c10(8, 5);
        coord_store_vec = {c1, c2, c3, c4, c5, c6, c7, c8, c9, c10};

        //matern correlation matrix
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double dist = eucl_dist(coord_store_vec[i], coord_store_vec[j]);
                covar_w(i, j) = matern(dist, phi, nu);
            }
        }
    }
    // sample w's
    Eigen::VectorXd mean_w(N);
    std::vector<Eigen::VectorXd> wt_store_vec(T);
    {
        Eigen::EigenMultivariateNormal<double> normal_sampler(mean_w, sigma_w*covar_w);
        for (int t = 0; t <= T-1; ++t) {
            wt_store_vec[t] = normal_sampler.samples(1);
        }
    }

    //O_t calculation (arbitrary beta parameter)


    std::vector<Eigen::VectorXd> ot_store_vec(T+1);
    {
        // normal sampler for initial O_0; take same cov matrix as for w_t, but mean different from 0
        Eigen::EigenMultivariateNormal<double> normal_sampler(mu_0, covar_w*sigma_0);
        Eigen::VectorXd o_initial = normal_sampler.samples(1);
        ot_store_vec[0] = o_initial;
        for (int t = 1; t <= T ; t++) {
            ot_store_vec[t] = rho * ot_store_vec[t-1] + xt_store_vec[t-1]*beta + wt_store_vec[t];
        }

    }


    // Eps data generation
    std::vector<Eigen::VectorXd> epst_store_vec(T);
    {
        Eigen::VectorXd mean_eps(N); // init is (I think) always to 0
        Eigen::MatrixXd covar_eps(N, N);
        for (int i = 0; i < N; ++i) {
            covar_eps(i, i) = sigma_eps;
        }
        Eigen::EigenMultivariateNormal<double> normal_sampler(mean_eps, covar_eps);
        for (int t = 0; t <= T-1; ++t) {
            epst_store_vec[t] = normal_sampler.samples(1);
        }

    }

    // Y_t data calculation
    std::vector<Eigen::VectorXd> yt_store_vec(T);
    {
        for (int t = 0; t < T; t++) {
            yt_store_vec[t] = ot_store_vec[t+1] + epst_store_vec[t];
        }
    }
    ///////// DATA GENERATION END /////////

    unsigned int n_iter = 3000;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    ar_model a(n_iter, T,yt_store_vec, xt_store_vec, coord_store_vec, ot_store_vec, beta, mu_0, rho, sigma_eps, sigma_w, sigma_0);
    a.sample();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

}