#include "ar_class.h"
#include <boost/math/distributions/inverse_gamma.hpp>

ar_model::ar_model(unsigned int n, unsigned int T, std::vector<Eigen::VectorXd> &y_store,
                   std::vector<Eigen::MatrixXd> &x_store, std::vector<coord> &coord_vec,
                   std::vector<Eigen::VectorXd>& ot_store, Eigen::VectorXd& beta, Eigen::VectorXd& mu_0, double& rho,
                   double& sigma_eps_true, double& sigma_w_true, double& sigma_0_true, double& phi_true, double& nu,
                   bool use_cholesky, u_int64_t seed) :
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
        nu(nu),
        n_iter(n), T(T),
        use_cholesky(use_cholesky), seed(seed)
{
    N = (X)[0].rows();
    p = (X)[0].cols();
};

void ar_model::sample() {
    std::mt19937 generator(seed);

    //proto buffers
    o_data::full_o_it o_stream;
    scalar::full_scalar_it rho_stream;
    vector::full_iter_vec beta_stream;
    vector::full_iter_vec mu_0_stream;
    scalar::full_scalar_it sig_0_stream;
    scalar::full_scalar_it sig_w_stream;
    scalar::full_scalar_it sig_eps_stream;
    scalar::full_scalar_it phi_stream;



    //Identity matrices


    //initialization DEBUG
//    Eigen::VectorXd beta = beta_true;
    std::vector<Eigen::VectorXd> o_store = o_store_true;
    Eigen::VectorXd mu_0 = mu_0_true;
    double rho = rho_true;
    phi = phi_true;
    double sigma_eps = sigma_eps_true;
    double sigma_w = sigma_w_true;
    double sigma_0 = sigma_0_true;

    //initialization of samplers

    Eigen::EigenMultivariateNormal<double> beta_sampler(Eigen::VectorXd::Zero(p), beta_sig_prior*Eigen::MatrixXd::Identity(p,p), use_cholesky, seed);
    Eigen::EigenMultivariateNormal<double> o_sampler(Eigen::VectorXd::Zero(N), Eigen::MatrixXd::Identity(N,N), use_cholesky, seed);
    Eigen::EigenMultivariateNormal<double> mu0_sampler(Eigen::VectorXd::Zero(N), mu0_sig_prior*Eigen::MatrixXd::Identity(N,N), use_cholesky, seed);
    std::normal_distribution<double> rho_sampler(0,sqrt(rho_sig_prior));
    std::gamma_distribution<double> sig_eps_sampler(ab_eps_prior.first, ab_eps_prior.second);
    std::gamma_distribution<double> sig_w_sampler(ab_w_prior.first, ab_w_prior.second);
    std::gamma_distribution<double> sig_0_sampler(ab_0_prior.first, ab_0_prior.second);
    std::normal_distribution<double> phi_sampler(0,sqrt(phi_cand_var));
    std::uniform_real_distribution<double> unif(0., 1.);


    //initialization

    //beta
    beta = beta_sampler.samples(1);
//
//    o vector
//    initialize all spatial effects to zero at beginning, including index 0
    o_store = std::vector<Eigen::VectorXd>(T+1);
    for (int t = 0; t < T+1 ; ++t) {
        o_store[t] = Eigen::VectorXd::Zero(N);
    }

//   / // rho
    rho = rho_sampler(generator);
    phi = 1;
//
//     mu_0
    mu_0 = mu0_sampler.samples(1);
//
    sigma_eps = sig_eps_sampler(generator);
    sigma_w = sig_w_sampler(generator);
    sigma_0 = sig_0_sampler(generator);
//
    Eigen::MatrixXd coord_mat = eucl_dist_matrix(coordinates);
    Eigen::MatrixXd matern_cov = calc_matern_mat(coord_mat, phi, nu);
    Eigen::MatrixXd matern_inv = matern_cov.inverse();
    Eigen::MatrixXd w_full_cov_inv =  matern_inv/sigma_w;
    // for n_iter loop: update o_store, then update beta, then update rho. write in proto files results after iteration

    for(int i = 0; i < n_iter; ++i) {


//        o_store update


        o_data::matrix* o_matrix= o_stream.add_m();
        o_matrix->set_iter(i);
//
        //update 0 zero first
        Eigen::MatrixXd o_zero_update_cov = post::calc_cov_eff_0(w_full_cov_inv,matern_inv, rho, sigma_0);
        Eigen::VectorXd o_zero_update_mean = o_zero_update_cov* post::calc_mean_eff_0(X[0], w_full_cov_inv, o_store[1], beta, matern_inv, mu_0, rho, sigma_0);
        o_sampler.setMean(o_zero_update_mean);
        o_sampler.setCovar(o_zero_update_cov);
        o_store[0] = o_sampler.samples(1);
//
//
//        // write in proto file
//
        o_data::vector* o_vector = o_matrix->add_vec();
        std::vector<double> o_vec(o_store[0].begin(), o_store[0].end());
        o_vector->mutable_vec_value()->Add(o_vec.begin(), o_vec.end());
        o_vector->set_t(0);
//
//        //  now update the rest
        for (int t = 1; t < T; ++t) {
            Eigen::MatrixXd ot_update_cov = post::calc_cov_eff_t(sigma_eps,w_full_cov_inv, rho );
            Eigen::VectorXd ot_update_mean = ot_update_cov* post::calc_mean_eff_t(y[t-1],X[t-1], X[t],w_full_cov_inv, o_store[t-1], o_store[t+1], beta, rho, sigma_eps);
            o_sampler.setMean(ot_update_mean);
            o_sampler.setCovar(ot_update_cov);
            o_store[t] = o_sampler.samples(1);
//            write in proto file
            o_data::vector* o_t_vector = o_matrix->add_vec(); //maybe not neccessary
            o_t_vector->mutable_vec_value()->Add(o_store[t].begin(), o_store[t].end());
            o_t_vector->set_t(t);
        }
//        //update last in o_store
        Eigen::MatrixXd oT_update_cov = post::calc_cov_eff_T(sigma_eps,w_full_cov_inv);;
        Eigen::VectorXd oT_update_mean =    oT_update_cov*post::calc_mean_eff_T(y[T-1], X[T-1], w_full_cov_inv, o_store[T-1],beta, rho, sigma_eps);
        o_sampler.setMean(oT_update_mean);
        o_sampler.setCovar(oT_update_cov);
        o_store[T] = o_sampler.samples(1);
//
        o_data::vector* o_T_vector = o_matrix->add_vec(); //maybe not neccessary
        o_T_vector->mutable_vec_value()->Add(o_store[T].begin(), o_store[T].end());
        o_T_vector->set_t(T);
//
//         mu_0 update
        Eigen::MatrixXd mu0_update_cov = post::calc_cov_mu_0(matern_inv, sigma_0, mu0_sig_prior);
        Eigen::VectorXd mu0_update_mean = mu0_update_cov * post::calc_mean_mu_0(matern_inv, o_store[0], sigma_0);

        mu0_sampler.setMean(mu0_update_mean);
        mu0_sampler.setCovar(mu0_update_cov);
        mu_0 = mu0_sampler.samples(1);
        //proto mu_0
        vector::vector* mu0_vector= mu_0_stream.add_vec_t();
        mu0_vector->set_iter(i);
        std::vector<double> mu0_vec(mu_0.begin(), mu_0.end());
        mu0_vector->mutable_vec_value()->Add(mu0_vec.begin(), mu0_vec.end());

        //beta update

        Eigen::MatrixXd beta_update_cov = post::calc_cov_beta(X,w_full_cov_inv, beta_sig_prior);

        Eigen::VectorXd beta_update_mean = beta_update_cov * post::calc_mean_beta(X, w_full_cov_inv, o_store, rho);
        beta_sampler.setCovar(beta_update_cov);
        beta_sampler.setMean(beta_update_mean);
        beta = beta_sampler.samples(1);

        //        proto beta
        vector::vector* beta_vector= beta_stream.add_vec_t();
        beta_vector->set_iter(i);
        std::vector<double> beta_vec(beta.begin(), beta.end());
        beta_vector->mutable_vec_value()->Add(beta_vec.begin(), beta_vec.end());


        //rho update

        double rho_update_cov = post::calc_var_rho(o_store, w_full_cov_inv, rho_sig_prior);
        double rho_update_mean = rho_update_cov * post::calc_mean_rho(X, w_full_cov_inv, o_store, beta);
        rho_sampler.param(std::normal_distribution<double>::param_type(rho_update_mean, sqrt(rho_update_cov) ));
        rho = rho_sampler(generator);

//      proto rho
        scalar::scalar* rho_scalar= rho_stream.add_scalar();
        rho_scalar->set_iter(i);
        rho_scalar->set_value(rho);


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
        scalar::scalar* sig_eps_scalar= sig_eps_stream.add_scalar();
        sig_eps_scalar->set_iter(i);
        sig_eps_scalar->set_value(sigma_eps);

        scalar::scalar* sig_w_scalar= sig_w_stream.add_scalar();
        sig_w_scalar->set_iter(i);
        sig_w_scalar->set_value(sigma_w);
//
        scalar::scalar* sig_0_scalar= sig_0_stream.add_scalar();
        sig_0_scalar->set_iter(i);
        sig_0_scalar->set_value(sigma_0);
//

        //phi MH step
        {
            phi_sampler.param(std::normal_distribution<double>::param_type(phi, sqrt(phi_cand_var)));
            double phi_cand = phi_sampler(generator);
            double u = unif(generator);
            std::cout<< "candidate phi: " << phi_cand << std::endl;
            if (phi_cand > 0) {
                double ratio = post::target_ratio_phi(phi, phi_cand, coord_mat, o_store, X, beta, mu_0, rho, sigma_w,
                                                      sigma_0, ab_phi_prior, nu);
                std::cout << "ratio: " << ratio << std::endl;
                std::cout << "u:  " << u << std::endl;
                if (ratio > u || ratio > 1) {
                    phi = phi_cand;
                    matern_cov = calc_matern_mat(coord_mat, phi, nu);
                    matern_inv = matern_cov.inverse();
                }
            }
        }
        std::cout << "phi: "<<phi <<std::endl;
        scalar::scalar* phi_scalar= phi_stream.add_scalar();
        phi_scalar->set_iter(i);
        phi_scalar->set_value(phi);

        //covariance update

        w_full_cov_inv =  matern_inv/sigma_w;

        std::cout<< "Iteration " << i << " finished" << std::endl;
    }

//    serializing streams
    proto::serialize_beta(beta_stream);
    proto::serialize_o(o_stream);
    proto::serialize_mu0(mu_0_stream);
    proto::serialize_rho(rho_stream);
    proto::serialize_phi(phi_stream);
    proto::serialize_sig_0(sig_0_stream);
    proto::serialize_sig_w(sig_w_stream);
    proto::serialize_sig_eps(sig_eps_stream);
}