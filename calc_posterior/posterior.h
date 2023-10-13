#include<cmath>
#include "Eigen/dense"
#include <iostream>

/// data type specified. maybe good in future to template the data structure class ?

// some calc purpose const


namespace post {
    // full conditional calculation functions for mean and covariance
    //initialization of dimensions (and other things maybe later?
// beta
    Eigen::VectorXd calc_mean_beta(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, std::vector<Eigen::VectorXd>& ot_store_vec,double& rho);
    Eigen::MatrixXd calc_cov_beta(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv,const double& sigma_beta_prior);

// rho
    double calc_mean_rho(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, std::vector<Eigen::VectorXd>& ot_store_vec, Eigen::VectorXd& beta);
    double calc_var_rho(std::vector<Eigen::VectorXd>& ot_store_vec, Eigen::MatrixXd& covar_w_inv, const double& sigma_rho_prior);

// O_T
    Eigen::VectorXd calc_mean_eff_T(const Eigen::VectorXd& y_T,const Eigen::MatrixXd& x_T, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_prev,
                                           Eigen::VectorXd& beta, double& sigma_eps, double& rho);
    Eigen::MatrixXd calc_cov_eff_T(float& sigma_eps, Eigen::MatrixXd& covar_w_inv);

// O_0
    Eigen::VectorXd calc_mean_eff_0( const Eigen::MatrixXd&  x_1, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_1, Eigen::VectorXd& beta, const Eigen::MatrixXd& S_0_inv,
                                    Eigen::VectorXd& mu_0, double& rho, double& sigma_0);
    Eigen::MatrixXd calc_cov_eff_0(Eigen::MatrixXd& covar_w_inv, const Eigen::MatrixXd& S_0_inv, double& rho, double& sigma_0);

// O_t
    Eigen::VectorXd calc_mean_eff_t(const Eigen::VectorXd& y_curr, const Eigen::MatrixXd& x_curr,const Eigen::MatrixXd& x_next, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_prev, Eigen::VectorXd& o_next,
                                           Eigen::VectorXd& beta, double& rho, double& sigma_eps);
// note: It is the exact same as for T, hence I will probably just use one declaration of the function
    Eigen::MatrixXd calc_cov_eff_t(double& sigma_eps, Eigen::MatrixXd& covar_w_inv);

    std::pair<double, double> calc_a_b_sigma_eps(double& a_prior, double& b_prior,
                                               const unsigned int& n,const  unsigned int& T , std::vector<Eigen::VectorXd>& y_store_vec, std::vector<Eigen::VectorXd>& o_store_vec);
    std::pair<double, double> calc_a_b_sigma_w(double& a_prior, double& b_prior,
                                             const unsigned int& n,const unsigned  int& T , std::vector<Eigen::MatrixXd>& x_store_vec,Eigen::MatrixXd& covar_w_inv ,
                                             std::vector<Eigen::VectorXd>& o_store_vec, Eigen::VectorXd& beta, double& rho);
    std::pair<double, double> calc_a_b_sigma_0(double& a_prior, double& b_prior,
                                             const unsigned int& n, const unsigned int& T, Eigen::MatrixXd& S_0_inv ,Eigen::VectorXd& o_0, Eigen::VectorXd& mu_0_prior);

    Eigen::VectorXd calc_mean_mu_0(Eigen::MatrixXd& S_0_inv, Eigen::VectorXd& o_0, double& sigma_0);
    Eigen::MatrixXd calc_cov_mu_0(Eigen::MatrixXd& S_0_inv, Eigen::VectorXd& o_0, double& sigma_0, double& sigma_mu_prior);

}
