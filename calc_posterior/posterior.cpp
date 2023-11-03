#include "posterior.h"

Eigen::VectorXd post::calc_mean_beta( const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, std::vector<Eigen::VectorXd>& ot_store_vec, double& rho){
        Eigen::VectorXd mean(x_store_vec[0].cols());
        for (int t = 1; t < x_store_vec.size(); t++) {
            mean += x_store_vec[t-1].transpose() * covar_w_inv * (ot_store_vec[t] - rho*ot_store_vec[t-1]);
        }
        return mean;
    };
    Eigen::MatrixXd post::calc_cov_beta(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, const double& sigma_beta_prior){
        auto id_mat = Eigen::MatrixXd::Identity(x_store_vec[0].cols(),x_store_vec[0].cols() );
        Eigen::MatrixXd temp(x_store_vec[0].cols(), x_store_vec[0].cols());
        //iteration starts at first element of x_store_vec which is t = 1
        for (int t = 0; t < x_store_vec.size()-1; ++t) {
            temp += x_store_vec[t].transpose() * covar_w_inv * x_store_vec[t];
        }
        temp +=  id_mat * 1/sigma_beta_prior;
        return temp.inverse();
    };

// rho
    double post::calc_mean_rho(const std::vector<Eigen::MatrixXd>& x_store_vec, Eigen::MatrixXd& covar_w_inv, std::vector<Eigen::VectorXd>& ot_store_vec, Eigen::VectorXd& beta){
        double mean = 0;
        for (int t = 1; t < ot_store_vec.size(); ++t) {
            mean += ot_store_vec[t-1].transpose() * covar_w_inv * (ot_store_vec[t] - x_store_vec[t-1]*beta);
        }
        return mean;
    };
    double post::calc_var_rho(std::vector<Eigen::VectorXd>& ot_store_vec, Eigen::MatrixXd& covar_w_inv,const double& sigma_rho_prior){
        double var = 0;
        for (int t = 1; t < ot_store_vec.size(); ++t) {
            var += ot_store_vec[t-1].transpose() * covar_w_inv * ot_store_vec[t-1];
        }
        return 1/(var + 1/sigma_rho_prior);
    };

// O_T
    Eigen::VectorXd post::calc_mean_eff_T( const Eigen::VectorXd& y_T, const Eigen::MatrixXd& x_T, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_prev,
                                           Eigen::VectorXd& beta, double& sigma_eps, double& rho){
        Eigen::VectorXd mean(y_T.rows());
        mean = (y_T * 1/sigma_eps); // + covar_w_inv*(rho*o_prev + x_T*beta);
        return mean;
    };
    Eigen::MatrixXd post::calc_cov_eff_T(double& sigma_eps, Eigen::MatrixXd& covar_w_inv){
        auto id_mat = Eigen::MatrixXd::Identity(covar_w_inv.cols(),covar_w_inv.cols() );
        Eigen::MatrixXd cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (id_mat * 1/sigma_eps + covar_w_inv);
        return cov.inverse();
    };

// O_0
    Eigen::VectorXd post::calc_mean_eff_0(const Eigen::MatrixXd&  x_1, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_1, Eigen::VectorXd& beta,const Eigen::MatrixXd& S_0_inv,
                                    Eigen::VectorXd& mu_0, double& rho, double& sigma_0){
        Eigen::VectorXd mean(covar_w_inv.cols());
        mean =  (1/sigma_0 * S_0_inv * mu_0) + (rho * covar_w_inv * (o_1 - x_1 * beta));
        return mean;
    };
    Eigen::MatrixXd post::calc_cov_eff_0(Eigen::MatrixXd& covar_w_inv,const Eigen::MatrixXd& S_0_inv, double& rho, double& sigma_0){
        Eigen::MatrixXd cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (pow(rho, 2) *  covar_w_inv + 1/sigma_0 * S_0_inv);
        return cov.inverse();
    };

// O_t
    Eigen::VectorXd post::calc_mean_eff_t(const Eigen::VectorXd& y_curr,const Eigen::MatrixXd& x_curr,const Eigen::MatrixXd& x_next, Eigen::MatrixXd& covar_w_inv, Eigen::VectorXd& o_prev, Eigen::VectorXd& o_next,
                                           Eigen::VectorXd& beta, double& rho, double& sigma_eps) {
        // important: o_store_vec stores only previous and next effect, i.e. O_t-1 and O_t+1
        Eigen::VectorXd mean(x_curr.rows());
        mean = y_curr * 1/sigma_eps + covar_w_inv * (x_curr*beta + rho  * (o_prev + o_next - x_next*beta));
        return mean;
    };
// note: It is the exact same as for T, hence I will probably just use one declaration of the function
    Eigen::MatrixXd post::calc_cov_eff_t(double& sigma_eps, Eigen::MatrixXd& covar_w_inv, double& rho){
        auto id_mat = Eigen::MatrixXd::Identity(covar_w_inv.cols(),covar_w_inv.cols());
        Eigen::MatrixXd cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (1/sigma_eps * id_mat + (1+pow(rho,2)) * covar_w_inv);
        return cov.inverse();
    };

    std::pair<double, double> post::calc_a_b_sigma_eps(const double& a_prior,const double& b_prior,
                                               const unsigned int& n,const  unsigned int& T ,const std::vector<Eigen::VectorXd>& y_store_vec, std::vector<Eigen::VectorXd>& o_store_vec){
        std::pair<double, double> param;
        param.first = a_prior + n*T /2;
        double temp = 0;
        for (int t = 1; t <= y_store_vec.size(); t++) {
            temp += (y_store_vec[t-1] - o_store_vec[t]).transpose() * (y_store_vec[t-1] - o_store_vec[t]);
        }
        param.second = b_prior +  temp/2;
        return param;
    };
    std::pair<double, double> post::calc_a_b_sigma_w(const double& a_prior,const double& b_prior,
                                             const unsigned int& n,const unsigned  int& T ,const std::vector<Eigen::MatrixXd>& x_store_vec, const Eigen::MatrixXd& matern_inv ,
                                             std::vector<Eigen::VectorXd>& o_store_vec, Eigen::VectorXd& beta, double& rho){
        std::pair<double, double> param;
        param.first = a_prior + 0.5 * n*T;
        double temp = 0;
        for (int t = 1; t <o_store_vec.size(); t++) {
            Eigen::VectorXd v = (o_store_vec[t] - rho *o_store_vec[t-1]  - x_store_vec[t-1] * beta);
            temp += v.transpose() * matern_inv * v;
        }
        param.second = b_prior + 0.5* temp;
        return param;
    };
    std::pair<double, double> post::calc_a_b_sigma_0(const double& a_prior, const double& b_prior,
                                             const unsigned int& n, const unsigned int& T, const Eigen::MatrixXd& S_0_inv ,Eigen::VectorXd& o_0,Eigen::VectorXd& mu_0){
        std::pair<double, double> param;
        param.first = a_prior + n/2;
        double temp = (o_0 - mu_0).transpose() * S_0_inv * (o_0 - mu_0);
        param.second = b_prior + temp/2;
        return param;
    };
    Eigen::VectorXd post::calc_mean_mu_0(const Eigen::MatrixXd& S_0_inv, Eigen::VectorXd& o_0, double& sigma_0){
        return S_0_inv * o_0 * 1/sigma_0;
    };
    Eigen::MatrixXd post::calc_cov_mu_0(const Eigen::MatrixXd& S_0_inv, double& sigma_0, const double& sigma_mu_prior){
        auto id_mat = Eigen::MatrixXd::Identity(S_0_inv.cols(),S_0_inv.cols());
        Eigen::MatrixXd cov = (1/sigma_0) * S_0_inv + (1/sigma_mu_prior) * id_mat;
        return cov.inverse();
    };


