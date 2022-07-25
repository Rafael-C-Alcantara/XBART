//////////////////////////////////////////////////////////////////////////////////////
// class to carry all intermediate data vectors, parameters across all functions
//////////////////////////////////////////////////////////////////////////////////////

#ifndef GUARD_fit_info_h
#define GUARD_fit_info_h

#include <ctime>
#include "common.h"
#include "utility.h"
#include <chrono>

class State
{
public:
    size_t dim_residual;         // residual size
    matrix<double> *residual_std; // a matrix to save all residuals

    // random number generators
    std::vector<double> prob;
    std::random_device rd;
    std::mt19937 gen;
    std::discrete_distribution<> d;

    // Splits
    matrix<double> split_count_all_tree;
    std::vector<double> split_count_all;
    std::vector<double> split_count_current_tree;
    std::vector<double> mtry_weight_current_tree;

    // for XBCF
    matrix<double> split_count_all_tree_ps;
    std::vector<double> split_count_all_ps;
    std::vector<double> mtry_weight_current_tree_ps;
    
    matrix<double> split_count_all_tree_trt;
    std::vector<double> split_count_all_trt;
    std::vector<double> mtry_weight_current_tree_trt;

    // mtry
    bool use_all = true;
    bool parallel = true;

    // fitinfo
    size_t n_min;
    size_t n_cutpoints;
    size_t p_categorical;
    size_t p_continuous;
    size_t p; // total number of variables = p_categorical + p_continuous
    size_t mtry;
    size_t n_y;                 // number of total data points in root node
    const double *X_std;        // pointer to original data
    std::vector<double> *y_std; // pointer to y data
    size_t max_depth;
    size_t num_trees;
    size_t num_sweeps;
    size_t burnin;
    bool sample_weights;
    double ini_var_yhat;

    matrix<size_t> *Xorder_std;

    // residual standard deviation
    double sigma;
    double sigma2; // sigma squared

    // paralization
    size_t nthread;

    // Logit Model
    // lambdas
    std::vector<std::vector<std::vector<double>>> *lambdas;
    std::vector<std::vector<std::vector<double>>> *lambdas_separate;

    // for continuous treatment XBCF
    matrix<double> *Z_std;
    std::vector<double> *tau_fit;
    std::vector<double> *mu_fit;
    bool treatment_flag;
    matrix<size_t> *Xorder_std_ps;
    matrix<size_t> *Xorder_std_trt;
    size_t p_ps;
    size_t p_trt;
    size_t p_categorical_ps;
    size_t p_categorical_trt;
    size_t p_continuous_ps;
    size_t p_continuous_trt;
    size_t mtry_ps;
    size_t mtry_trt;
    size_t num_trees_ps;
    size_t num_trees_trt;

    void update_sigma(double sigma)
    {
        this->sigma = sigma;
        this->sigma2 = pow(sigma, 2);
        return;
    }

    State(const double *Xpointer, matrix<size_t> &Xorder_std, size_t N, size_t p, size_t num_trees, size_t p_categorical, size_t p_continuous, bool set_random_seed, size_t random_seed, size_t n_min, size_t n_cutpoints, size_t mtry, const double *X_std, size_t num_sweeps, bool sample_weights, std::vector<double> *y_std, double sigma, size_t max_depth, double ini_var_yhat, size_t burnin, size_t dim_residual, size_t nthread)
    {

        // Init containers
        // initialize predictions_std at given value / number of trees
        this->residual_std = new matrix<double>();
        ini_matrix((*this->residual_std), N, dim_residual);

        // Random
        this->prob = std::vector<double>(2, 0.5);
        this->gen = std::mt19937(rd());
        if (set_random_seed)
        {
            gen.seed(random_seed);
        }
        this->d = std::discrete_distribution<>(prob.begin(), prob.end());

        // Splits
        ini_xinfo(this->split_count_all_tree, p, num_trees);

        this->split_count_all_tree_ps.resize(0);
        this->split_count_all_tree_trt.resize(0);
        this->split_count_current_tree = std::vector<double>(p, 0);
        this->mtry_weight_current_tree = std::vector<double>(p, 0);
        this->split_count_all = std::vector<double>(p, 0);
        this->sigma = sigma;

        this->n_min = n_min;
        this->n_cutpoints = n_cutpoints;
        this->p_categorical = p_categorical;
        this->p_continuous = p_continuous;
        this->mtry = mtry;
        this->X_std = X_std;
        this->p = p_categorical + p_continuous;
        this->n_y = N;
        this->num_trees = num_trees;
        this->num_sweeps = num_sweeps;
        this->sample_weights = sample_weights;
        this->y_std = y_std;
        this->max_depth = max_depth;
        this->burnin = burnin;
        this->ini_var_yhat = ini_var_yhat;
        this->Xorder_std = &Xorder_std;

        this->nthread = nthread;
        return;
    }

    void update_split_counts(size_t tree_ind)
    {
        mtry_weight_current_tree = mtry_weight_current_tree + split_count_current_tree;
        split_count_all_tree[tree_ind] = split_count_current_tree;
        return;
    }
};

class NormalState : public State
{
public:
    NormalState(const double *Xpointer, matrix<size_t> &Xorder_std, size_t N, size_t p, size_t num_trees, size_t p_categorical, size_t p_continuous, bool set_random_seed, size_t random_seed, size_t n_min, size_t n_cutpoints, size_t mtry, const double *X_std, size_t num_sweeps, bool sample_weights, std::vector<double> *y_std, double sigma, size_t max_depth, double ini_var_yhat, size_t burnin, size_t dim_residual, size_t nthread, bool parallel) : State(Xpointer, Xorder_std, N, p, num_trees, p_categorical, p_continuous, set_random_seed, random_seed, n_min, n_cutpoints, mtry, X_std, num_sweeps, sample_weights, y_std, sigma, max_depth, ini_var_yhat, burnin, dim_residual, nthread)
    {
        this->sigma = sigma;
        this->sigma2 = pow(sigma, 2);
        this->parallel = parallel;
    }
};

class LogitState : public State
{

    void ini_lambda(std::vector<std::vector<std::vector<double>>> &lambdas, size_t num_trees, size_t dim_residual)
    {
        // each tree has different number of theta vectors, each is of the size dim_residual (num classes)
        lambdas.resize(num_trees);
    }

    void ini_lambda_separate(std::vector<std::vector<std::vector<double>>> &lambdas, size_t num_trees, size_t dim_residual)
    {
        // each tree have (num classes) of lambda vectors
        lambdas.resize(num_trees);
        for (size_t i = 0; i < num_trees; i++)
        {
            lambdas[i].resize(dim_residual);
        }
    }

public:
    LogitState(const double *Xpointer, matrix<size_t> &Xorder_std, size_t N, size_t p, size_t num_trees, size_t p_categorical, size_t p_continuous, bool set_random_seed, size_t random_seed, size_t n_min, size_t n_cutpoints, size_t mtry, const double *X_std, size_t num_sweeps, bool sample_weights, std::vector<double> *y_std, double sigma, size_t max_depth, double ini_var_yhat, size_t burnin, size_t dim_residual, size_t nthread) : State(Xpointer, Xorder_std, N, p, num_trees, p_categorical, p_continuous, set_random_seed, random_seed, n_min, n_cutpoints, mtry, X_std, num_sweeps, sample_weights, y_std, sigma, max_depth, ini_var_yhat, burnin, dim_residual, nthread)
    {
        this->lambdas = new std::vector<std::vector<std::vector<double>>>();
        this->lambdas_separate = new std::vector<std::vector<std::vector<double>>>();
        ini_lambda((*this->lambdas), num_trees, dim_residual);
        ini_lambda_separate((*this->lambdas_separate), num_trees, dim_residual);
    }
};

class NormalLinearState : public State
{
public:
    NormalLinearState(matrix<double> *Z_std, const double *Xpointer_ps, const double *Xpointer_trt, matrix<size_t> &Xorder_std_ps, matrix<size_t> &Xorder_std_trt, size_t N, size_t p_ps, size_t p_trt, size_t num_trees_ps, size_t num_trees_trt, size_t p_categorical_ps, size_t p_categorical_trt, size_t p_continuous_ps, size_t p_continuous_trt, bool set_random_seed, size_t random_seed, size_t n_min, size_t n_cutpoints, size_t mtry_ps, size_t mtry_trt, const double *X_std, size_t num_sweeps, bool sample_weights, std::vector<double> *y_std, double sigma, size_t max_depth, double ini_var_yhat, size_t burnin, size_t dim_residual, size_t nthread, bool parallel) : State(Xpointer_ps, Xorder_std_ps, N, p_ps, num_trees_ps, p_categorical_ps, p_continuous_ps, set_random_seed, random_seed, n_min, n_cutpoints, mtry_ps, Xpointer_ps, num_sweeps, sample_weights, y_std, sigma, max_depth, ini_var_yhat, burnin, dim_residual, nthread)
    {
        ini_xinfo(this->split_count_all_tree_ps, p_ps, num_trees_ps);
        ini_xinfo(this->split_count_all_tree_trt, p_trt, num_trees_trt);
        this->split_count_all_ps = std::vector<double>(p_ps, 0);
        this->mtry_weight_current_tree_ps = std::vector<double>(p_ps, 0);
        this->split_count_all_trt = std::vector<double>(p_trt, 0);
        this->mtry_weight_current_tree_trt = std::vector<double>(p_trt, 0);
        this->Z_std = Z_std;
        this->sigma = sigma;
        this->sigma2 = pow(sigma, 2);
        this->parallel = parallel;
        this->tau_fit = (new std::vector<double>(N, 0));
        this->mu_fit = (new std::vector<double>(N, 0));
        this->Xorder_std_ps = &Xorder_std_ps;
        this->Xorder_std_trt = &Xorder_std_trt;
        this->p_ps = p_ps;
        this->p_trt = p_trt;
        this->p_categorical_ps = p_categorical_ps;
        this->p_categorical_trt = p_categorical_trt;
        this->p_continuous_ps = p_continuous_ps;
        this->p_continuous_trt = p_continuous_trt;
        this->mtry_ps = mtry_ps;
        this->mtry_trt = mtry_trt;
        this->num_trees_ps = num_trees_ps;
        this->num_trees_trt = num_trees_trt;
    }
};
#endif
