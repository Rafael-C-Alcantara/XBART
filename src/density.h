#include <cmath>
#include <vector>
#include <numeric> 
#include "gsl/gsl_integration.h"

using namespace std;

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

double kernal(double x, void * params);

double test_f(double x, void * params);

double density_single(double x, std::vector<double> x_vec, double tau, double mu, double sigma);

// double density_vec(std::vector<double> x_vec, std::vector<double> x_prior, double tau, bool take_log);

double test_d();

