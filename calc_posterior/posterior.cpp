#include "posterior.h"

Eigen::VectorXf post::calc_mean_beta( const std::vector<Eigen::MatrixXf>& x_store_vec, Eigen::MatrixXf& covar_w_inv, std::vector<Eigen::VectorXf>& ot_store_vec, float& rho){
        Eigen::VectorXf mean(x_store_vec[0].cols());
        for (int t = 1; t < x_store_vec.size(); t++) {
            mean += x_store_vec[t-1].transpose() * covar_w_inv * (ot_store_vec[t] - rho*ot_store_vec[t-1]);
        }
        return mean;
    };
    Eigen::MatrixXf post::calc_cov_beta(const std::vector<Eigen::MatrixXf>& x_store_vec, Eigen::MatrixXf& covar_w_inv, const float& sigma_beta_prior){
        auto id_mat = Eigen::MatrixXf::Identity(x_store_vec[0].cols(),x_store_vec[0].cols() );
        Eigen::MatrixXf temp(x_store_vec[0].cols(), x_store_vec[0].cols());
        //iteration starts at first element of x_store_vec which is t = 1
        for (int t = 0; t < x_store_vec.size()-1; ++t) {
            temp += x_store_vec[t].transpose() * covar_w_inv * x_store_vec[t];
        }
        temp +=  id_mat * 1/sigma_beta_prior;
        return temp.inverse();
    };

// rho
    float post::calc_mean_rho(const std::vector<Eigen::MatrixXf>& x_store_vec, Eigen::MatrixXf& covar_w_inv, std::vector<Eigen::VectorXf>& ot_store_vec, Eigen::VectorXf& beta){
        float mean = 0;
        for (int t = 1; t < x_store_vec.size(); ++t) {
            mean += ot_store_vec[t-1].transpose() * covar_w_inv * (ot_store_vec[t] - x_store_vec[t-1]*beta);
        }
        return mean;
    };
    float post::calc_var_rho(std::vector<Eigen::VectorXf>& ot_store_vec, Eigen::MatrixXf& covar_w_inv,const float& sigma_rho_prior){
        float var = 0;
        for (int t = 1; t < ot_store_vec.size(); ++t) {
            var += ot_store_vec[t-1].transpose() * covar_w_inv * ot_store_vec[t-1] + 1/sigma_rho_prior;
        }
        return 1/var;
    };

// O_T
    Eigen::VectorXf post::calc_mean_eff_T( const Eigen::VectorXf& y_T, const Eigen::MatrixXf& x_T, Eigen::MatrixXf& covar_w_inv, Eigen::VectorXf& o_prev,
                                           Eigen::VectorXf& beta, float& sigma_eps, float& rho){
        Eigen::VectorXf mean(y_T.rows());
        mean = (y_T * 1/sigma_eps); // + covar_w_inv*(rho*o_prev + x_T*beta);
        return mean;
    };
    Eigen::MatrixXf post::calc_cov_eff_T(float& sigma_eps, Eigen::MatrixXf& covar_w_inv){
        auto id_mat = Eigen::MatrixXf::Identity(covar_w_inv.cols(),covar_w_inv.cols() );
        Eigen::MatrixXf cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (id_mat * 1/sigma_eps + covar_w_inv);
        return cov.inverse();
    };

// O_0
    Eigen::VectorXf post::calc_mean_eff_0(const Eigen::MatrixXf&  x_1, Eigen::MatrixXf& covar_w_inv, Eigen::VectorXf& o_1, Eigen::VectorXf& beta,const Eigen::MatrixXf& S_0_inv,
                                    Eigen::VectorXf& mu_0, float& rho, float& sigma_0){
        Eigen::VectorXf mean(covar_w_inv.cols());
        mean =  (1/sigma_0 * S_0_inv * mu_0) + (rho * covar_w_inv * (o_1 - x_1 * beta));
        return mean;
    };
    Eigen::MatrixXf post::calc_cov_eff_0(Eigen::MatrixXf& covar_w_inv,const Eigen::MatrixXf& S_0_inv, float& rho, float& sigma_0){
        Eigen::MatrixXf cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (pow(rho, 2) *  covar_w_inv + S_0_inv);
        return cov.inverse();
    };

// O_t
    Eigen::VectorXf post::calc_mean_eff_t(const Eigen::VectorXf& y_curr,const Eigen::MatrixXf& x_curr,const Eigen::MatrixXf x_next, Eigen::MatrixXf& covar_w_inv, Eigen::VectorXf& o_prev, Eigen::VectorXf& o_next,
                                           Eigen::VectorXf& beta, float& rho, float& sigma_eps) {
        // important: o_store_vec stores only previous and next effect, i.e. O_t-1 and O_t+1
        Eigen::VectorXf mean(x_curr.rows());
        mean = y_curr * 1/sigma_eps + covar_w_inv * (x_curr*beta + rho  * (o_prev + o_next - x_next*beta));
        return mean;
    };
// note: It is the exact same as for T, hence I will probably just use one declaration of the function
    Eigen::MatrixXf post::calc_cov_eff_t(float& sigma_eps, Eigen::MatrixXf& covar_w_inv){
        auto id_mat = Eigen::MatrixXf::Identity(covar_w_inv.cols(),covar_w_inv.cols());
        Eigen::MatrixXf cov(covar_w_inv.cols(), covar_w_inv.cols());
        cov = (1/sigma_eps * id_mat + covar_w_inv);
        return cov.inverse();
    };

    std::pair<float, float> post::calc_a_b_sigma_eps(float& a_prior, float& b_prior,
                                               const unsigned int& n,const  unsigned int& T , std::vector<Eigen::VectorXf>& y_store_vec, std::vector<Eigen::VectorXf>& o_store_vec){
        std::pair<float, float> param;
        param.first = a_prior + n*T /2;
        float temp = 0;
        for (int t = 1; t <= y_store_vec.size(); t++) {
            temp += (y_store_vec[t-1] - o_store_vec[t]).transpose() * (y_store_vec[t-1] - o_store_vec[t]);
        }
        param.second = b_prior + 1/2 * temp;
        return param;
    };
    std::pair<float, float> post::calc_a_b_sigma_w(float& a_prior, float& b_prior,
                                             const unsigned int& n,const unsigned  int& T , std::vector<Eigen::MatrixXf>& x_store_vec,Eigen::MatrixXf& covar_w_inv ,
                                             std::vector<Eigen::VectorXf>& o_store_vec, Eigen::VectorXf& beta, float& rho){
        std::pair<float, float> param;
        param.first = a_prior + n*T /2;
        float temp = 0;
        for (int t = 1; t <= x_store_vec.size(); t++) {
            temp += (o_store_vec[t] - rho *o_store_vec[t-1]  - x_store_vec[t-1] * beta).transpose() * covar_w_inv * (o_store_vec[t] - rho *o_store_vec[t-1]  - x_store_vec[t-1] * beta);
        }
        param.second = b_prior + 1/2 * temp;
        return param;
    };
    std::pair<float, float> post::calc_a_b_sigma_0(float& a_prior, float& b_prior,
                                             const unsigned int& n, const unsigned int& T, Eigen::MatrixXf& S_0_inv ,Eigen::VectorXf& o_0, Eigen::VectorXf& mu_0_prior){
        std::pair<float,float> param;
        param.first = a_prior + n/2;
        param.second = b_prior + 1/2 * ( o_0 - mu_0_prior).transpose() * S_0_inv * (o_0 - mu_0_prior);
        return param;
    };
    Eigen::VectorXf post::calc_mean_mu_0(Eigen::MatrixXf& S_0_inv, Eigen::VectorXf& o_0, float& sigma_0){
        return S_0_inv * o_0 /sigma_0;
    };
    Eigen::MatrixXf post::calc_cov_mu_0(Eigen::MatrixXf& S_0_inv, Eigen::VectorXf& o_0, float& sigma_0, float& sigma_mu_prior){
        auto id_mat = Eigen::MatrixXf::Identity(S_0_inv.cols(),S_0_inv.cols());
        Eigen::MatrixXf cov = (1/sigma_0) * S_0_inv + 1/sigma_mu_prior * id_mat;
        return cov.inverse();
    };


