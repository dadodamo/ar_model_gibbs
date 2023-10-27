#include "ar_class.h"

void ar_model::sample() const {
    std::mt19937 generator(seed);

    //proto buffers
    o_data::full_o_it o_stream;
    scalar::full_scalar_it rho_stream;
    vector::full_iter_vec beta_stream;
    vector::full_iter_vec mu_0_stream;
    scalar::full_scalar_it sig_0_stream;
    scalar::full_scalar_it sig_w_stream;
    scalar::full_scalar_it sig_eps_stream;

    // beta update vector
    // float rho
    // sigma_eps
    // sigma_w
    // sigma_0
    // mu_0
    //current iteration O store vector

    //UPDATES
    //1. Update big O
    //2. Update beta
    //3.  update rho
    //4. update sigma (for now fixed)
    //5. update mu_0 (for now fixed)

    // updating big O:
        // iterate over iteration store vector and update accordingly, last and first excluded
    //usage of sampler: Eigen::EigenMultivariateNormal<float> normal_sampler(mean_X, covar_X)
    //                  normal_sampler.samples(n) n no. of samples
    //Identity matrices
    Eigen::MatrixXd id_p = Eigen::VectorXd::Ones(p).asDiagonal();
    Eigen::MatrixXd id_N = Eigen::VectorXd::Ones(N).asDiagonal();


    //initialization
    Eigen::VectorXd beta(p);
//    Eigen::VectorXd beta = beta_true; // debug


    std::vector<Eigen::VectorXd> o_store(T+1);
//    std::vector<Eigen::VectorXd> o_store = ot; //debug


    Eigen::VectorXd mu_0(N);
    //    Eigen::VectorXd mu_0 = mu_0_true; // debug

    double rho;
//    double rho = rho_true; //debug

    double sigma_eps;
    double sigma_w;
    double sigma_0;

    //initialization of samplers
    Eigen::VectorXd zero_prior_p = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd zero_prior_N = Eigen::VectorXd::Zero(N);

    Eigen::EigenMultivariateNormal<double> beta_sampler(zero_prior_p, beta_sig_prior*id_p, use_cholesky, seed);
    Eigen::EigenMultivariateNormal<double> o_sampler(zero_prior_N, id_N, use_cholesky, seed);
    Eigen::EigenMultivariateNormal<double> mu0_sampler(zero_prior_N, id_N, use_cholesky, seed);
    std::normal_distribution<double> rho_sampler(rho_mean_prior,sqrt(rho_sig_prior));
    std::gamma_distribution<double> sig_eps_sampler(ab_eps_prior.first, ab_eps_prior.second);
    std::gamma_distribution<double> sig_w_sampler(ab_w_prior.first, ab_w_prior.second);
    std::gamma_distribution<double> sig_0_sampler(ab_0_prior.first, ab_0_prior.second);


    //initialization
    {
        //beta
        beta = beta_sampler.samples(1);

        //o vector
        //initialize all spatial effects to zero at beginning, including index 0
        for (int t = 0; t <= T; ++t) {
            o_store[t] = Eigen::VectorXd::Zero(N);
        }

        // rho
        rho = rho_sampler(generator);

        // mu_0
        mu_0 = mu0_sampler.samples(1);

        sigma_eps = sig_eps_sampler(generator);
        sigma_w = sig_w_sampler(generator);
        sigma_0 = sig_0_sampler(generator);

    }

    // for n_iter loop: update o_store, then update beta, then update rho. write in proto files results after iteration

    for(int i = 0; i < n_iter; ++i) {

        //update covar inverse
        Eigen::MatrixXd w_full_cov_inv = (1/sigma_w) * matern_inv;

        //o_store update

        o_data::matrix* o_matrix= o_stream.add_m();
        o_matrix->set_iter(i);
        //update 0 zero first
        Eigen::MatrixXd o_zero_update_cov = post::calc_cov_eff_0(w_full_cov_inv,matern_inv, rho, sigma_0);
        Eigen::VectorXd o_zero_update_mean = o_zero_update_cov
                                            * post::calc_mean_eff_0(X[0], w_full_cov_inv, o_store[1],
                                                                    beta, matern_inv, mu_0, rho, sigma_0);

        // insert calculations from posterior
        o_sampler.setMean(o_zero_update_mean);
        o_sampler.setCovar(o_zero_update_cov);
        o_store[0] = o_sampler.samples(1);

        // write in proto file

        o_data::vector* o_vector = o_matrix->add_vec();
        o_vector->mutable_vec_value()->Add(o_store[0].begin(), o_store[0].end());
        o_vector->set_t(0);

        //  now update the rest

        for (int t = 1; t < T-1; ++t) {
            Eigen::MatrixXd ot_update_cov = post::calc_cov_eff_t(sigma_eps,w_full_cov_inv );
            Eigen::VectorXd ot_update_mean = ot_update_cov
                                                * post::calc_mean_eff_t(y[t-1],X[t-1], X[t],
                                                                        w_full_cov_inv, o_store[t-1], o_store[t+1], beta, rho, sigma_eps);
            o_sampler.setMean(ot_update_mean);
            o_sampler.setCovar(ot_update_cov);

            o_store[t] = o_sampler.samples(1);
            //write in proto file
            o_data::vector* o_t_vector = o_matrix->add_vec(); //maybe not neccessary
            o_t_vector->mutable_vec_value()->Add(o_store[t].begin(), o_store[t].end());
            o_t_vector->set_t(t);
        }

        //update last in o_store
        Eigen::MatrixXd oT_update_cov = post::calc_cov_eff_t(sigma_eps,w_full_cov_inv );;
        Eigen::VectorXd oT_update_mean =    oT_update_cov
                                            *post::calc_mean_eff_T(y[T-1], X[T-1], w_full_cov_inv, o_store[T-1],
                                                               beta, sigma_eps, rho);
        o_sampler.setMean(oT_update_mean);
        o_sampler.setCovar(oT_update_cov);
        o_store[T] = o_sampler.samples(1);


        o_data::vector* o_T_vector = o_matrix->add_vec(); //maybe not neccessary
        o_T_vector->mutable_vec_value()->Add(o_store[T].begin(), o_store[T].end());
        o_T_vector->set_t(T);

        // mu_0 update

        Eigen::MatrixXd mu0_update_cov = post::calc_cov_mu_0(matern_inv, o_store[0], sigma_0, mu0_sig_prior);
        Eigen::VectorXd mu0_update_mean = mu0_update_cov * post::calc_mean_mu_0(matern_inv, o_store[0], sigma_0);

        mu0_sampler.setCovar(mu0_update_cov);
        mu0_sampler.setMean(mu0_update_mean);
        mu_0 = mu0_sampler.samples(1);

        //proto mu_0
        vector::vector* mu0_vector= mu_0_stream.add_vec_t();
        mu0_vector->set_iter(i);
        std::vector<double> mu0_vec(mu_0.data(), mu_0.data() + mu_0.rows() * mu_0.cols());
        mu0_vector->mutable_vec_value()->Add(mu0_vec.begin(), mu0_vec.end());

        //beta update

        Eigen::MatrixXd beta_update_cov = post::calc_cov_beta(X,w_full_cov_inv, beta_sig_prior);
        Eigen::VectorXd beta_update_mean = beta_update_cov * post::calc_mean_beta(X, w_full_cov_inv, o_store, rho);

        beta_sampler.setCovar(beta_update_cov);
        beta_sampler.setMean(beta_update_mean);
        beta = beta_sampler.samples(1);
        //proto beta
        vector::vector* beta_vector= beta_stream.add_vec_t();
        beta_vector->set_iter(i);
        std::vector<double> beta_vec(beta.data(), beta.data() + beta.rows() * beta.cols());
        beta_vector->mutable_vec_value()->Add(beta_vec.begin(), beta_vec.end());

        //rho update
        scalar::scalar* rho_scalar= rho_stream.add_scalar();
        rho_scalar->set_iter(i);

        double rho_update_cov = post::calc_var_rho(o_store, w_full_cov_inv, rho_sig_prior);
        double rho_update_mean = rho_update_cov * post::calc_mean_rho(X, w_full_cov_inv, o_store, beta);
        rho_sampler.param(std::normal_distribution<double>::param_type(rho_update_mean, sqrt(rho_update_cov) ));
        rho = rho_sampler(generator);

        rho_scalar->set_value(rho);


        // sigma calculation
        std::pair sig_eps_update_ab = post::calc_a_b_sigma_eps(ab_eps_prior.first, ab_eps_prior.second, N, T, y, o_store);
        sig_eps_sampler.param(std::gamma_distribution<double>::param_type(sig_eps_update_ab.first, sig_eps_update_ab.second));
        sigma_eps = sig_eps_sampler(generator);

        std::pair sig_w_update_ab = post::calc_a_b_sigma_w(ab_w_prior.first, ab_w_prior.second,N, T, X, matern_inv, o_store, beta, rho);
        sig_w_sampler.param(std::gamma_distribution<double>::param_type(sig_w_update_ab.first, sig_w_update_ab.second));
        sigma_w = sig_w_sampler(generator);

        std::pair sig_0_update_ab = post::calc_a_b_sigma_0(ab_0_prior.first, ab_0_prior.second,N, T, matern_inv,o_store[0], mu_0);
        sig_0_sampler.param(std::gamma_distribution<double>::param_type(sig_0_update_ab.first, sig_0_update_ab.second));
        sigma_0 = sig_0_sampler(generator);

        scalar::scalar* sig_eps_scalar= sig_eps_stream.add_scalar();
        sig_eps_scalar->set_iter(i);
        sig_eps_scalar->set_value(sigma_eps);

        scalar::scalar* sig_w_scalar= sig_w_stream.add_scalar();
        sig_w_scalar->set_iter(i);
        sig_w_scalar->set_value(sigma_w);

        scalar::scalar* sig_0_scalar= sig_0_stream.add_scalar();
        sig_0_scalar->set_iter(i);
        sig_0_scalar->set_value(sigma_0);



        std::cout<< "Iteration " << i << " finished" << std::endl;
    }

//    serializing streams
    std::string serialized_o;
    if (!o_stream.SerializeToString(&serialized_o)) {
        std::cerr << "Failed to write random effects data." << std::endl;
    }
    std::ofstream outputFile_o("o_serialized.bin", std::ios::binary);
    if (outputFile_o.is_open()) {
        outputFile_o.write(serialized_o.c_str(), serialized_o.size());
        outputFile_o.close();
        std::cout << "Serialized data written to o_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

    std::string serialized_beta;
    if (!beta_stream.SerializeToString(&serialized_beta)) {
        std::cerr << "Failed to write beta data." << std::endl;
    }
    std::ofstream outputFile_beta("beta_serialized.bin", std::ios::binary);
    if (outputFile_beta.is_open()) {
        outputFile_beta.write(serialized_beta.c_str(), serialized_beta.size());
        outputFile_beta.close();
        std::cout << "Serialized data written to beta_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
//
    std::string serialized_rho;
    if (!rho_stream.SerializeToString(&serialized_rho)) {
        std::cerr << "Failed to write rho data." << std::endl;
    }
    std::ofstream outputFile_rho("rho_serialized.bin", std::ios::binary);
    if (outputFile_rho.is_open()) {
        outputFile_rho.write(serialized_rho.c_str(), serialized_rho.size());
        outputFile_rho.close();
        std::cout << "Serialized data written to rho_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

    std::string serialized_mu0;
    if (!mu_0_stream.SerializeToString(&serialized_mu0)) {
        std::cerr << "Failed to write mu_0 data." << std::endl;
    }
    std::ofstream outputFile_mu0("mu0_serialized.bin", std::ios::binary);
    if (outputFile_mu0.is_open()) {
        outputFile_mu0.write(serialized_mu0.c_str(), serialized_mu0.size());
        outputFile_mu0.close();
        std::cout << "Serialized data written to mu0_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

    std::string serialized_sig_eps;
    if (!sig_eps_stream.SerializeToString(&serialized_sig_eps)) {
        std::cerr << "Failed to write sigma_eps data." << std::endl;
    }
    std::ofstream outputFile_sig_eps("sig_eps_serialized.bin", std::ios::binary);
    if (outputFile_sig_eps.is_open()) {
        outputFile_sig_eps.write(serialized_sig_eps.c_str(), serialized_sig_eps.size());
        outputFile_sig_eps.close();
        std::cout << "Serialized data written to sig_eps_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

    std::string serialized_sig_w;
    if (!sig_w_stream.SerializeToString(&serialized_sig_w)) {
        std::cerr << "Failed to write sigma_w data." << std::endl;
    }
    std::ofstream outputFile_sig_w("sig_w_serialized.bin", std::ios::binary);
    if (outputFile_sig_w.is_open()) {
        outputFile_sig_w.write(serialized_sig_w.c_str(), serialized_sig_w.size());
        outputFile_sig_w.close();
        std::cout << "Serialized data written to sig_w_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }

    std::string serialized_sig_0;
    if (!sig_0_stream.SerializeToString(&serialized_sig_0)) {
        std::cerr << "Failed to write sigma_0 data." << std::endl;
    }
    std::ofstream outputFile_sig_0("sig_0_serialized.bin", std::ios::binary);
    if (outputFile_sig_0.is_open()) {
        outputFile_sig_0.write(serialized_sig_0.c_str(), serialized_sig_0.size());
        outputFile_sig_0.close();
        std::cout << "Serialized data written to sig_0_serialized.bin" << std::endl;
    } else {
        std::cerr << "Error: Unable to open the output file." << std::endl;
    }
//
};