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
#include "cmake-build-debug/proto/ydata.pb.h"
#include <fstream>
#include "coordinates/coordinates.h"
#include "ar_model/ar_class.h"
#include<random>

// Source files

int main(int argc,char* argv[]) {
    

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // data generation of model using T = 10, N = 10, p = 5
    // can be rewritten with generic param, as class+methods or function
    const unsigned int T = 50;
    const unsigned int p = 5;
    const unsigned int N = 10;


    // all fixed parameters (for testing functions etc.)
    float sigma_eps = 2;
    float phi = 1.0;
    float nu = 0.5;
    float sigma_w = 5;
    float rho = 0.8;
    Eigen::VectorXf beta(p) ;
    beta << 1, 2, 3, 4, 5;

    // prior param
    float sigma_0 = 1.0;
    Eigen::VectorXf mu_0(N);
    for (int i = 0; i <N ; ++i) {
        mu_0(i) = i + 2;
    }
    float sigma_mu_0_prior = 1;

    float sig_eps_a_prior = 0.5;
    float sig_eps_b_prior = 0.5;

    float sig_w_a_prior = 0.5;
    float sig_w_b_prior = 0.5;

    float sig_0_a_prior = 0.5;
    float sig_0_b_prior = 0.5;



////// DATA GENERATION ///////
    //X_t generation
    std::vector<Eigen::MatrixXf> xt_store_vec(T);

    // sampler from external header file (OK)
    {
        Eigen::VectorXf mean_X(N);
        Eigen::MatrixXf covar_X(N, N);
        for (int i = 0; i < N; ++i) {
            mean_X(i) = 0;
            covar_X(i, i) = 1;
        }
        Eigen::EigenMultivariateNormal<float> normal_sampler(mean_X, covar_X);
        for (int t = 0; t <= T-1; ++t) {
            Eigen::MatrixXf X(N, p);
            X = normal_sampler.samples(p).cwiseAbs(); // abs: to have at least meaningful matrix (wrt choice of beta param)
            xt_store_vec[t] = X;
        }
    }
    // w_t generation + generate coord points
    Eigen::MatrixXf covar_w(N, N);
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
                float dist = eucl_dist(coord_store_vec[i], coord_store_vec[j]);
                covar_w(i, j) = sigma_w*matern(dist, phi, nu);
            }
        }
    }
    // sample w's
    Eigen::VectorXf mean_w(N);
    std::array<Eigen::VectorXf, T> wt_store_vec;
    {
        Eigen::EigenMultivariateNormal<float> normal_sampler(mean_w, covar_w);
        for (int t = 0; t <= T-1; ++t) {
            wt_store_vec[t] = normal_sampler.samples(1);
        }
    }

    //O_t calculation (aribtrary beta parameter)

    std::array<Eigen::VectorXf, T> ot_store_vec;
    {
        // normal sampler for initial O_0; take same cov matrix as for w_t, but mean different from 0
        Eigen::EigenMultivariateNormal<float> normal_sampler(mu_0, covar_w*sigma_0 /sigma_w );
        Eigen::VectorXf o_initial = normal_sampler.samples(1);
        ot_store_vec[0] = o_initial;
        for (int t = 1; t <= T ; t++) {
            ot_store_vec[t] = rho * ot_store_vec[t-1] + xt_store_vec[t-1]*beta;
        }

    }

    // Eps data generation
    std::array<Eigen::VectorXf, T> epst_store_vec;
    {
        Eigen::VectorXf mean_eps(N); // init is (I think) always to 0
        Eigen::MatrixXf covar_eps(N, N);
        for (int i = 0; i < N; ++i) {
            covar_eps(i, i) = sigma_eps;
        }
        Eigen::EigenMultivariateNormal<float> normal_sampler(mean_eps, covar_eps);
        for (int t = 0; t <= T-1; ++t) {
            epst_store_vec[t] = normal_sampler.samples(1);
        }

    }
    y_data::full_y y;
    // Y_t data calculation
    std::vector<Eigen::VectorXf> yt_store_vec(T);
    {
        for (int t = 0; t <= T-1; t++) {
            yt_store_vec[t] = ot_store_vec[t+1] + epst_store_vec[t];
            y_data::vector* y_vector = y.add_vec_t();
            y_vector->set_t(t);
            y_vector->mutable_vec_value()-> Add(yt_store_vec[t].begin(), yt_store_vec[t].end());
            // copy data vector to protobuf
        }
    }
    ///////// DATA GENERATION END /////////


    // some matrix calculations
    Eigen::MatrixXf covar_inv = covar_w.inverse();
    Eigen::MatrixXf S_0 = covar_w/sigma_w;
    Eigen::MatrixXf S_0_inv = S_0.inverse();


    //serialization
  {
    std::string serializedData;
    if (!y.SerializeToString(&serializedData)) {
      std::cerr << "Failed to write data." << std::endl;
      return -1;
    }
    std::ofstream outputFile("serialized_data.bin", std::ios::binary);
    if (outputFile.is_open()) {
        outputFile.write(serializedData.c_str(), serializedData.size());
        outputFile.close();
        std::cout << "Serialized data written to serialized_data.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
  }



    unsigned int n_iter = 100;


    ar_model a(n_iter, T,yt_store_vec, xt_store_vec, coord_store_vec);
    a.sample();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

}