#include <cmath>
#include <vector>
#include <numeric> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "gsl/gsl_integration.h"


using namespace std;
using std::endl;

std::ostream &operator<<(std::ostream &out, const std::vector<double> &v);
std::ostream &operator<<(std::ostream &out, const std::vector<size_t> &v);
std::ostream &operator<<(std::ostream &out, const std::vector<std::vector<double>> &v);
std::ostream &operator<<(std::ostream &out, const std::vector<std::vector<size_t>> &v);

struct Params{
public:
    double n; // lenthg of x_vec
    double tau;
    double mu; // lognormal mean
    double sigma; // lognormal sd
    double x; // x value to be valualted
    std::vector<double>  x_vec; // vector of x_i's
    double iter;

    Params(double x, std::vector<double> x_vec, double tau, double mu, double sigma)
    {
        this->iter = 0;
        this->tau = tau;
        this->mu = mu;
        this->sigma = sigma;
        this->x = x;
        this->x_vec = x_vec;
        this->n = x_vec.size();
        return;
    }
};

double kernal(double x, void * params)
	{   
        Params p = * (Params *) params;
        double tau = p.tau;
        double mu = p.mu;
        double sigma = p.sigma;
        double y = p.x;
        std::vector<double> y_vec = p.x_vec;
        double n = y_vec.size();
        double iter = p.iter;
        
		return exp( - pow(y - y_vec[iter], 2) / 2 / tau / (x+1) - pow(log(x) - mu, 2) / 2 / pow(sigma, 2) ) / x / sqrt(x + 1);
        // return 0.0;
    }

double test_f(double x, void * params)
{
    return exp( - pow(x, 2) / 2) / sqrt( 2 * 3.1415926 ) ;
}


    double density_single(double x, std::vector<double> x_vec, double tau, double mu, double sigma)
    {
        double output = 0.0;
        Params params(x, x_vec, tau, mu, sigma);
        std::vector<double> int_vec(params.n);
        std::vector<double> error(params.n);
        // std::vector<double> integral_vec(this->n);

        gsl_function F;
        F.function = &kernal;
        F.params = &params;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(3000);
       
        for (size_t i = 0; i < params.n; i++){
            // std::cout << "density_single loop i " << i << endl;
            
            gsl_integration_qagiu (&F, 0, 0, 1e-6, 1000, w, &int_vec[i], &error[i]);

            params.iter += 1;
            F.params = &params;
        }
        output = std::accumulate(int_vec.begin(),int_vec.end(), 0.0) / params.n / 2 / 3.1415926 / sqrt(tau) / sigma;
        
        gsl_integration_workspace_free(w);
        return output;
    }
    
    double test_d(){
        gsl_function F;
        F.function = &test_f;
        // F.params = &;

        double output, error;

        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        gsl_integration_qagi (&F, 0, 1e-6, 1000, w, &output, &error);
        gsl_integration_workspace_free (w);
        return output;
    }

// double density_vec(std::vector<double> x_vec, std::vector<double> x_prior, double tau, bool take_log)
// {

//     size_t h;
//     size_t n = x_prior.size();
//     std::cout << "n " << n << endl;
//     size_t N = x_vec.size() + n + 1;
//     std::cout << "N " << N << endl;
//     std::vector<double> x_density(x_vec.size());
//     double eta_n, rho_h, nu_n1, nu_n2, nu_n, mu_n, sigma_n, temp, output, x;
//     eta_n = nu_n1 = nu_n2 = 1;

//     for (size_t h = 1; h <= n; h++)
//     {
//         std::cout << "x_prior, h " << h << endl;
//         rho_h = double (N-h) / double(N-h+1);
//         eta_n = eta_n / rho_h;
//         nu_n1 = nu_n1 * (2-rho_h) / pow(rho_h, 2);
//         nu_n2 = nu_n2 * pow(rho_h, -2);
//     }
//     nu_n = nu_n1 - nu_n2;
//     mu_n = 2*log(eta_n - 1) - 0.5*log(nu_n + pow(eta_n - 1, 2));
//     sigma_n = sqrt(log( 1 + nu_n / pow(eta_n - 1, 2)));
//     h = n;

//     std::cout << "rho_h " << rho_h << endl;
//     std::cout << "eta_n " << eta_n << endl;
//     std::cout << "nu_n1 " << nu_n1 << endl;
//     std::cout << "nu_n2 " << nu_n1 << endl;
//     std::cout << "nu_n " << nu_n << endl;
//     std::cout << "mu_n " << mu_n << endl;
//     std::cout << "sigma_n " << sigma_n << endl;

//     for (size_t i = 0; i < x_vec.size(); i++)
//     {
//         std::cout << "x_vec, i " << i << endl;
//         x = x_vec[i];
//         std::cout << "x " << x << endl;
//         temp = density_single(x, x_prior, tau, mu_n, sigma_n);
//         std::cout << "temp " << temp << endl;
//         x_density[i] = temp;

//         h++;
//         rho_h = double (N-h) / double(N-h+1);
//         eta_n *= 1 / rho_h;
//         nu_n1 *= (2-rho_h) / pow(rho_h, 2);
//         nu_n2 *= pow(rho_h, -2);
//         nu_n = nu_n1 - nu_n2;
//         mu_n = 2*log(eta_n - 1) - 0.5*log(nu_n + pow(eta_n - 1, 2));
//         sigma_n = sqrt(log( 1 + nu_n / pow(eta_n - 1, 2)));
//         x_prior.push_back(x);
//         std::cout << "h " << h << endl;
//         std::cout << "rho_h " << rho_h << endl;
//         std::cout << "eta_n " << eta_n << endl;
//         std::cout << "nu_n1 " << nu_n1 << endl;
//         std::cout << "nu_n2 " << nu_n1 << endl;
//         std::cout << "nu_n " << nu_n << endl;
//         std::cout << "mu_n " << mu_n << endl;
//         std::cout << "sigma_n " << sigma_n << endl;
            
//     }
//     if (take_log)
//     {
//         output = 0;
//         for (size_t i = 0; i < x_density.size(); i++) { output += log(x_density[i]);}
//         return output;
//     }
//     else
//     {
//         output = 0;
//         for (size_t i = 0; i < x_density.size(); i++) { output *= x_density[i];}
//         return output;
//     }
    
// }
