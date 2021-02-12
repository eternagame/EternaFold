//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi.cpp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

#include <InnerOptimizationWrapperViterbi.hpp>

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi::InnerOptimizationWrapperViterbi()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

InnerOptimizationWrapperViterbi::InnerOptimizationWrapperViterbi(OptimizationWrapper *optimizer,
                                                                 const std::vector<int> &units,
                                                                 const std::vector<double> &C) :
    SubgradientDescent<double>(1000,
                               optimizer->computation.ComputeParameterNormBound(units, optimizer->params.ExpandHyperparameters(C)),
                               optimizer->computation.ComputeGradientNormBound(units, optimizer->params.ExpandHyperparameters(C)),
                               Min(optimizer->params.ExpandHyperparameters(C))),
    optimizer(optimizer), units(units), C(C), bias(optimizer->params.GetNumParameters()), best_f(std::numeric_limits<double>::infinity())
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi::LoadBias()
//
// Load linear bias.
//////////////////////////////////////////////////////////////////////

void InnerOptimizationWrapperViterbi::LoadBias(const std::vector<double> &bias)
{
    this->bias = bias;
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi::ComputeFunction()
//
// Compute the regularized logloss using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

double InnerOptimizationWrapperViterbi::ComputeFunction(const std::vector<double> &x)
{
    std::vector<double> Ce = optimizer->params.ExpandHyperparameters(C);
    return optimizer->computation.ComputeFunction(units, x + optimizer->base_values, false) + 0.5 * DotProduct(Ce, x*x) + DotProduct(x, bias);  
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi::ComputeSubgradient()
//
// Compute the regularized logloss gradient using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

void InnerOptimizationWrapperViterbi::ComputeSubgradient(std::vector<double> &g, const std::vector<double> &x)
{
    std::vector<double> Ce = optimizer->params.ExpandHyperparameters(C);
    g = optimizer->computation.ComputeGradient(units, x + optimizer->base_values, false) + Ce * x + bias;
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi::Report()
//
// Routines for printing results and messages from the optimizer.
//////////////////////////////////////////////////////////////////////

void InnerOptimizationWrapperViterbi::Report(int iteration, double f, const std::vector<double> &x, const std::vector<double> &g,
                                             double norm_bound, double step_size)
{
    // write results to disk
    
    best_f = f;
    best_x = x;
    optimizer->params.WriteToFile(SPrintF("optimize.params.iter%d", iteration), best_x + optimizer->base_values);
    
    // write results to console
    
    std::vector<double> Ce = optimizer->params.ExpandHyperparameters(C);
    const double unregularized = f - 0.5 * DotProduct(Ce, x*x);
    optimizer->PrintMessage(SPrintF("Inner iteration %d: f = %lf (%lf), |x| = %lf, step = %lf, efficiency = %lf%%", 
                                    iteration, f, unregularized, Norm(x), step_size, optimizer->computation.GetEfficiency()));
}

void InnerOptimizationWrapperViterbi::Report(const std::string &s) 
{
    optimizer->PrintMessage(SPrintF("Inner message: %s", s.c_str()));
}

