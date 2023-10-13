#include "ar_class.h"

void ar_model::sample() const {

    //proto buffers
    o_data::full_o_it o_stream;
    scalar::full_scalar_it rho_stream;
    vector::full_iter_vec beta_stream;
    vector::full_iter_vec mu_0_stream;


    std::mt19937 generator;
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
    Eigen::MatrixXd id_p = Eigen::VectorXd::Ones(this->p).asDiagonal();
    Eigen::MatrixXd id_N = Eigen::VectorXd::Ones(this->N).asDiagonal();


    //initialization
    Eigen::VectorXd beta(this->p);
    std::vector<Eigen::VectorXd> o_store(this->T+1);
    Eigen::VectorXd mu_0(this->N);
    double rho;
    double sigma_eps = 1.; //fixed for now
    double sigma_w = 1.;     // fixed for now
    double sigma_0 = 1.;      //fixed for now




    //initialization
    {
        //beta
        Eigen::VectorXd beta_prior = Eigen::VectorXd::Zero(this->p);
        Eigen::EigenMultivariateNormal<double> beta_sampler(beta_prior, this->beta_sig_prior*id_p);
        beta = beta_sampler.samples(1);

        //o vector
        Eigen::VectorXd o_zero_prior = Eigen::VectorXd::Zero(this->N);
        Eigen::EigenMultivariateNormal<double> o_zero_sampler(o_zero_prior, sigma_0*id_N);
        o_store[0] = o_zero_sampler.samples(1);
        //initialize all spatial effects to zero at beginning except for index 0
        for (int t = 1; t <= this->T; ++t) {
            o_store[t] = Eigen::VectorXd::Zero(N);
        }

        // rho
        std::normal_distribution<double> rho_sampler(this->rho_mean_prior,sqrt(this->rho_sig_prior));
        rho = rho_sampler(generator);

        //mu_0
        //beta
        Eigen::VectorXd mu_0_prior = Eigen::VectorXd::Zero(this->N);
        Eigen::EigenMultivariateNormal<double> mu_0_sampler(mu_0_prior, this->mu0_sig_prior*id_N);
        Eigen::VectorXd mu_0 = mu_0_sampler.samples(1);

    }


    // for n_iter loop: update o_store, then update beta, then update rho. write in proto files results after iteration
    for (int i = 0; i < this->n_iter; ++i) {

        //update covar inverse
        Eigen::MatrixXd w_full_cov_inv = (1/sigma_w) * this->matern_inv;




        //o_store update

        //update 0 zero first
        Eigen::MatrixXd o_zero_update_cov = post::calc_cov_eff_0(w_full_cov_inv,this->matern_inv, rho, sigma_0); // insert calculations from posterior
        Eigen::VectorXd o_zero_update_mean = o_zero_update_cov
                                            * post::calc_mean_eff_0((this->X)[0], w_full_cov_inv, o_store[1],
                                                                    beta, matern_inv, mu_0, rho, sigma_0); // insert calculations from posterior
        Eigen::EigenMultivariateNormal<double> o_zero_sampler(o_zero_update_mean, o_zero_update_cov);
        o_store[0] = o_zero_sampler.samples(1);

        // write in proto file
        o_data::matrix* o_matrix= o_stream.add_m();
        o_matrix->set_iter(i);
        o_data::vector* o_vector = o_matrix->add_vec();
        o_vector->mutable_vec_value()->Add(o_store[0].begin(), o_store[0].end());
        o_vector->set_t(0);



        // now update the rest

        for (int t = 1; t < this->T-1; ++t) {
            Eigen::MatrixXd ot_update_cov = post::calc_cov_eff_t(sigma_eps,w_full_cov_inv );
            Eigen::VectorXd ot_update_mean = ot_update_cov
                                                * post::calc_mean_eff_t((this->y)[t-1],(this->X)[t-1], (this->X)[t],
                                                                        w_full_cov_inv, o_store[t-1], o_store[t+1], beta, rho, sigma_eps);
            Eigen::EigenMultivariateNormal<double> ot_sampler(ot_update_mean, ot_update_cov);
            o_store[t] = ot_sampler.samples(1);
            //write in proto file
            o_data::vector* o_t_vector = o_matrix->add_vec(); //maybe not neccessary
            o_t_vector->mutable_vec_value()->Add(o_store[t].begin(), o_store[t].end());
            o_t_vector->set_t(t);
        }

        //update last in o_store
        Eigen::MatrixXd oT_update_cov = post::calc_cov_eff_t(sigma_eps,w_full_cov_inv );;
        Eigen::VectorXd oT_update_mean =    oT_update_cov
                                            *post::calc_mean_eff_T((this->y)[this->T -1], (this->X)[this->T -1], w_full_cov_inv, o_store[T-1],
                                                               beta, sigma_eps, rho);
        Eigen::EigenMultivariateNormal<double> ot_sampler(oT_update_mean, oT_update_cov);
        o_store[this->T-1] = ot_sampler.samples(1);


        o_data::vector* o_T_vector = o_matrix->add_vec(); //maybe not neccessary
        o_T_vector->mutable_vec_value()->Add(o_store[T].begin(), o_store[T].end());
        o_T_vector->set_t(T);

        //beta update
        vector::vector* beta_vector= beta_stream.add_vec_t();
        beta_vector->set_iter(i);
        Eigen::MatrixXd beta_update_cov = post::calc_cov_beta((this->X),w_full_cov_inv, this->beta_sig_prior);

        Eigen::VectorXd beta_update_mean = beta_update_cov * post::calc_mean_beta((this->X), w_full_cov_inv, o_store, rho);

        Eigen::EigenMultivariateNormal<double> beta_sampler(beta_update_mean, beta_update_cov);
        beta = beta_sampler.samples(1);

        //proto beta
        std::vector<double> beta_vec(beta.data(), beta.data() + beta.rows() * beta.cols());
        beta_vector->mutable_vec_value()->Add(beta_vec.begin(), beta_vec.end());

        //rho update
        scalar::scalar* rho_scalar= rho_stream.add_scalar();
        rho_scalar->set_iter(i);

        double rho_update_cov = post::calc_var_rho(o_store, w_full_cov_inv, this->rho_sig_prior);
        double rho_update_mean = rho_update_cov * post::calc_mean_rho(this->X, w_full_cov_inv, o_store, beta);
        std::normal_distribution<double> rho_sampler(rho_update_mean,rho_update_cov);
        rho = rho_sampler(generator);

        rho_scalar->set_value(rho);


        //rest later
        std::cout<< "Iteration " << i << " finished" << std::endl;
    }

    //serializing streams
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
};