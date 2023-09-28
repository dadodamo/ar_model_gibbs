#include<cmath>
#include "Eigen/dense"
#include <iostream>

/// data type specified. maybe good in future to template the data structure class ?

// some calc purpose const


namespace post {
    // full conditional calculation functions for mean and covariance
    //initialization of dimensions (and other things maybe later?
// beta
    Eigen::VectorXf calc_mean_beta(const std::vector<Eigen::MatrixXf>& x_store_vec, Eigen::MatrixXf& covar_w_inv, std::vector<Eigen::VectorXf>& ot_store_vec,float& rho);
    Eigen::MatrixXf calc_cov_beta(const std::vector<Eigen::MatrixXf>& x_store_vec, Eigen::MatrixXf& covar_w_inv,const float& sigma_beta_prior);

// rho
    float calc_mean_rho(const std::vector<Eigen::MatrixXf>& x_store_vec, Eigen::MatrixXf& covar_w_inv, std::vector<Eigen::VectorXf>& ot_store_vec, Eigen::VectorXf& beta);
    float calc_var_rho(std::vector<Eigen::VectorXf>& ot_store_vec, Eigen::MatrixXf& covar_w_inv, const float& sigma_rho_prior);

// O_T
    Eigen::VectorXf calc_mean_eff_T(const Eigen::VectorXf& y_T,const Eigen::MatrixXf& x_T, Eigen::MatrixXf& covar_w_inv, Eigen::VectorXf& o_prev,
                                           Eigen::VectorXf& beta, float& sigma_eps, float& rho);
    Eigen::MatrixXf calc_cov_eff_T(float& sigma_eps, Eigen::MatrixXf& covar_w_inv);

// O_0
    Eigen::VectorXf calc_mean_eff_0( const Eigen::MatrixXf&  x_1, Eigen::MatrixXf& covar_w_inv, Eigen::VectorXf& o_1, Eigen::VectorXf& beta, const Eigen::MatrixXf& S_0_inv,
                                    Eigen::VectorXf& mu_0, float& rho, float& sigma_0);
    Eigen::MatrixXf calc_cov_eff_0(Eigen::MatrixXf& covar_w_inv, const Eigen::MatrixXf& S_0_inv, float& rho, float& sigma_0);

// O_t
    Eigen::VectorXf calc_mean_eff_t(const Eigen::VectorXf& y_curr, const Eigen::MatrixXf& x_curr,const Eigen::MatrixXf x_next, Eigen::MatrixXf& covar_w_inv, Eigen::VectorXf& o_prev, Eigen::VectorXf& o_next,
                                           Eigen::VectorXf& beta, float& rho, float& sigma_eps);
// note: It is the exact same as for T, hence I will probably just use one declaration of the function
    Eigen::MatrixXf calc_cov_eff_t(float& sigma_eps, Eigen::MatrixXf& covar_w_inv);

    std::pair<float, float> calc_a_b_sigma_eps(float& a_prior, float& b_prior,
                                               const unsigned int& n,const  unsigned int& T , std::vector<Eigen::VectorXf>& y_store_vec, std::vector<Eigen::VectorXf>& o_store_vec);
    std::pair<float, float> calc_a_b_sigma_w(float& a_prior, float& b_prior,
                                             const unsigned int& n,const unsigned  int& T , std::vector<Eigen::MatrixXf>& x_store_vec,Eigen::MatrixXf& covar_w_inv ,
                                             std::vector<Eigen::VectorXf>& o_store_vec, Eigen::VectorXf& beta, float& rho);
    std::pair<float, float> calc_a_b_sigma_0(float& a_prior, float& b_prior,
                                             const unsigned int& n, const unsigned int& T, Eigen::MatrixXf& S_0_inv ,Eigen::VectorXf& o_0, Eigen::VectorXf& mu_0_prior);

    Eigen::VectorXf calc_mean_mu_0(Eigen::MatrixXf& S_0_inv, Eigen::VectorXf& o_0, float& sigma_0);
    Eigen::MatrixXf calc_cov_mu_0(Eigen::MatrixXf& S_0_inv, Eigen::VectorXf& o_0, float& sigma_0, float& sigma_mu_prior);

}
