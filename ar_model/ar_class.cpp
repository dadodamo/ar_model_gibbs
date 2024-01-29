#include "ar_class.h"
#include <boost/math/distributions/inverse_gamma.hpp>

ar_model::ar_model(std::vector<Eigen::VectorXd> &y_store,
                   std::vector<Eigen::MatrixXd> &x_store, std::vector<coord> &coord_vec,
                   std::vector<Eigen::VectorXd>& ot_store, Eigen::VectorXd& beta, Eigen::VectorXd& mu_0, double& rho,
                   double& sigma_eps_true, double& sigma_w_true, double& sigma_0_true, double& phi_true, double& nu) :
        y(y_store), X(x_store),
        coordinates(coord_vec),
        o_store_true(ot_store),
        beta_true(beta),
        mu_0_true(mu_0),
        rho_true(rho),
        sigma_eps_true(sigma_eps_true),
        sigma_w_true(sigma_w_true),
        sigma_0_true(sigma_0_true),
        phi_true(phi_true),
        nu(nu)
{
    T = (X).size();
    N = (X)[0].rows();
    p = (X)[0].cols();
};

void ar_model::init() {


    beta_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(p), Eigen::MatrixXd::Identity(p,p), use_cholesky, seed);
    o_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(N), Eigen::MatrixXd::Identity(N,N), use_cholesky, seed);
    mu0_sampler = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(N), mu0_sig_prior*Eigen::MatrixXd::Identity(N,N), use_cholesky, seed);
    rho_sampler = std::normal_distribution<double>(0,sqrt(rho_sig_prior));
    sig_eps_sampler = std::gamma_distribution<double>(ab_eps_prior.first, ab_eps_prior.second);
    sig_w_sampler = std::gamma_distribution<double>(ab_w_prior.first, ab_w_prior.second);
    sig_0_sampler = std::gamma_distribution<double> (ab_0_prior.first, ab_0_prior.second);
    phi_sampler = std::normal_distribution<double>(0,sqrt(phi_cand_var));

    phi = 1.;
    beta = beta_sampler.samples(1);
    rho = rho_sampler(generator);
    mu_0 = mu0_sampler.samples(1);
    sigma_eps = sig_eps_sampler(generator);
    sigma_w = sig_w_sampler(generator);
    sigma_0 = sig_0_sampler(generator);

//    o vector
//    initialize all spatial effects to zero at beginning, including index 0
    o_store = std::vector<Eigen::VectorXd>(T+1);
    for (int t = 0; t < T+1 ; ++t) {
        o_store[t] = Eigen::VectorXd::Zero(N);
    }


    coord_mat = eucl_dist_matrix(coordinates);
    matern_cov = calc_matern_mat(coord_mat, phi, nu);
    matern_inv = matern_cov.inverse();
    w_full_cov_inv =  matern_inv/sigma_w;
}


void ar_model::serialize() {
    proto::serialize(sample_stream);
}

void ar_model::standardize() {
    // std_mean_y
};

void ar_model::write_curr_state(){
    sampler_data::matrix* o_matrix = sample_stream.add_o();
    for (int i = 0; i < o_store.size(); ++i) {
        sampler_data::vector* o_vector = o_matrix->add_vec();
        std::vector<double> o_vec(o_store[0].begin(), o_store[0].end());
        o_vector->mutable_vec_value()->Add(o_vec.begin(), o_vec.end());
    }

    std::vector<double> mu0_vec(mu_0.begin(), mu_0.end());
    std::vector<double> beta_vec(beta.begin(), beta.end());

    sample_stream.add_rho(rho);
    sample_stream.add_phi(phi);
    sample_stream.add_sigma_0(sigma_0);
    sample_stream.add_sigma_eps(sigma_eps);
    sample_stream.add_sigma_w(sigma_w);
    sample_stream.add_mu_0()->mutable_vec_value()->Add(mu0_vec.begin(), mu0_vec.end());
    sample_stream.add_beta()->mutable_vec_value()->Add(beta_vec.begin(), beta_vec.end());

}

void ar_model::sample() {

//        o_store update

    //update 0 zero first
    Eigen::MatrixXd o_zero_update_cov = post::calc_cov_eff_0(w_full_cov_inv,matern_inv, rho, sigma_0);
    Eigen::VectorXd o_zero_update_mean = o_zero_update_cov* post::calc_mean_eff_0(X[0], w_full_cov_inv, o_store[1], beta, matern_inv, mu_0, rho, sigma_0);
    o_sampler.setMean(o_zero_update_mean);
    o_sampler.setCovar(o_zero_update_cov);
    o_store[0] = o_sampler.samples(1);

    //  now update the rest
    for (int t = 1; t < T; ++t) {
        Eigen::MatrixXd ot_update_cov = post::calc_cov_eff_t(sigma_eps,w_full_cov_inv, rho );
        Eigen::VectorXd ot_update_mean = ot_update_cov* post::calc_mean_eff_t(y[t-1],X[t-1], X[t],w_full_cov_inv, o_store[t-1], o_store[t+1], beta, rho, sigma_eps);
        o_sampler.setMean(ot_update_mean);
        o_sampler.setCovar(ot_update_cov);
        o_store[t] = o_sampler.samples(1);
    }
//        //update last in o_store
    Eigen::MatrixXd oT_update_cov = post::calc_cov_eff_T(sigma_eps,w_full_cov_inv);;
    Eigen::VectorXd oT_update_mean =    oT_update_cov*post::calc_mean_eff_T(y[T-1], X[T-1], w_full_cov_inv, o_store[T-1],beta, rho, sigma_eps);
    o_sampler.setMean(oT_update_mean);
    o_sampler.setCovar(oT_update_cov);
    o_store[T] = o_sampler.samples(1);

//         mu_0 update
    Eigen::MatrixXd mu0_update_cov = post::calc_cov_mu_0(matern_inv, sigma_0, mu0_sig_prior);
    Eigen::VectorXd mu0_update_mean = mu0_update_cov * post::calc_mean_mu_0(matern_inv, o_store[0], sigma_0);

    mu0_sampler.setMean(mu0_update_mean);
    mu0_sampler.setCovar(mu0_update_cov);
    mu_0 = mu0_sampler.samples(1);

    //beta update

    Eigen::MatrixXd beta_update_cov = post::calc_cov_beta(X,w_full_cov_inv, beta_sig_prior);

    Eigen::VectorXd beta_update_mean = beta_update_cov * post::calc_mean_beta(X, w_full_cov_inv, o_store, rho);
    beta_sampler.setCovar(beta_update_cov);
    beta_sampler.setMean(beta_update_mean);
    beta = beta_sampler.samples(1);

    //rho update

    double rho_update_cov = post::calc_var_rho(o_store, w_full_cov_inv, rho_sig_prior);
    double rho_update_mean = rho_update_cov * post::calc_mean_rho(X, w_full_cov_inv, o_store, beta);
    rho_sampler.param(std::normal_distribution<double>::param_type(rho_update_mean, sqrt(rho_update_cov) ));
    rho = rho_sampler(generator);

//         sigma calculation

    std::pair sig_eps_update_ab = post::calc_a_b_sigma_eps(ab_eps_prior.first, ab_eps_prior.second, N, T, y, o_store);
    sig_eps_sampler.param(std::gamma_distribution<double>::param_type(sig_eps_update_ab.first, 1./sig_eps_update_ab.second));
    sigma_eps = 1./sig_eps_sampler(generator);
//
    std::pair sig_w_update_ab = post::calc_a_b_sigma_w(ab_w_prior.first, ab_w_prior.second,N, T, X, matern_inv, o_store, beta, rho);
    sig_w_sampler.param(std::gamma_distribution<double>::param_type(sig_w_update_ab.first, 1/sig_w_update_ab.second));
    sigma_w = 1.0/sig_w_sampler(generator);

    std::pair sig_0_update_ab = post::calc_a_b_sigma_0(ab_0_prior.first, ab_0_prior.second,N, matern_inv,o_store[0], mu_0);
    sig_0_sampler.param(std::gamma_distribution<double>::param_type(sig_0_update_ab.first, 1./sig_0_update_ab.second));
    sigma_0 = 1./sig_0_sampler(generator);
//

    //phi MH step
    {
        phi_sampler.param(std::normal_distribution<double>::param_type(phi, sqrt(phi_cand_var)));
        double phi_cand = phi_sampler(generator);
        double u = unif(generator);
        if (phi_cand > 0) {
            double ratio = post::target_ratio_phi(phi, phi_cand, coord_mat, o_store, X, beta, mu_0, rho, sigma_w,
                                                  sigma_0, ab_phi_prior, nu);

            if (ratio > u || ratio > 1) {
                phi = phi_cand;
                matern_cov = calc_matern_mat(coord_mat, phi, nu);
                matern_inv = matern_cov.inverse();
                phi_accept_rate++;
            }
        }
    }
    //covariance update

    w_full_cov_inv =  matern_inv/sigma_w;
    iter_count++;
}