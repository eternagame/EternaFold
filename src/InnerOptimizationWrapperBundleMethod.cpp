//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod.ipp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

#include "InnerOptimizationWrapperBundleMethod.hpp"

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod::InnerOptimizationWrapperBundleMethod()
//
// Constructor.
//////////////////////////////////////////////////////////////////////


InnerOptimizationWrapperBundleMethod<RealT>::InnerOptimizationWrapperBundleMethod(OptimizationWrapper *optimization_wrapper,
                                                                                            const std::vector<int> &units,
                                                                                            const std::vector<RealT> &C) :
    BundleMethod<double>(1000,C[1]),
    InnerOptimizationWrapper(optimization_wrapper, units, C)
{
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod::ComputeFunction()
//
// Compute the regularized logloss using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////


RealT InnerOptimizationWrapperBundleMethod<RealT>::ComputeFunction(const std::vector<RealT> &w)
{
  return this->optimization_wrapper->GetComputationWrapper().ComputeFunction(this->units, w, true, true, this->optimization_wrapper->GetOptions().GetRealValue("log_base"));
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod::ComputeSubgradient()
//
// Compute the regularized logloss gradient using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////


void InnerOptimizationWrapperBundleMethod<RealT>::ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &w)
{
  g = this->optimization_wrapper->GetComputationWrapper().ComputeGradient(this->units, w, true, true, this->optimization_wrapper->GetOptions().GetRealValue("log_base"));
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod::Report()
//
// Routines for printing results and messages from the optimizer.
//////////////////////////////////////////////////////////////////////


void InnerOptimizationWrapperBundleMethod<RealT>::Report(int iteration, const std::vector<RealT> &w, RealT f, const std::vector<RealT> &g, RealT norm_bound, RealT step_size)
{
    // write results to disk
    this->optimization_wrapper->GetParameterManager().WriteToFile(SPrintF("optimize.params.iter%d", iteration), w);
    
    // write results to console
    this->optimization_wrapper->PrintMessage(SPrintF("Inner iteration %d: f = %lf (%lf), |w| = %lf, |g| = %lf, norm bound = %lf, step = %lf, efficiency = %lf%%", 
                                                     iteration, double(f), double(f - RealT(0.5) * DotProduct(this->C, w*w)),
                                                     double(Norm(w)), double(step_size), double(Norm(g)), double(norm_bound),
                                                     double(this->optimization_wrapper->GetComputationEngine().GetEfficiency())));
}


void InnerOptimizationWrapperBundleMethod<RealT>::Report(const std::string &s) 
{
    this->optimization_wrapper->PrintMessage(SPrintF("Inner message: %s", s.c_str()));
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod::Minimize()
//
// Perform subgradient optimization.
//////////////////////////////////////////////////////////////////////


RealT InnerOptimizationWrapperBundleMethod<RealT>::Minimize(std::vector<RealT> &x0)
{
    return BundleMethod<RealT>::Minimize(x0);
}
