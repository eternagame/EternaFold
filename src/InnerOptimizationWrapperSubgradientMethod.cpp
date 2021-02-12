//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod.ipp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

#include "InnerOptimizationWrapperSubgradientMethod.hpp"

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod::InnerOptimizationWrapperSubgradientMethod()
//
// Constructor.
//////////////////////////////////////////////////////////////////////


InnerOptimizationWrapperSubgradientMethod::InnerOptimizationWrapperSubgradientMethod(OptimizationWrapper *optimization_wrapper,
                                                                                            const std::vector<int> &units,
                                                                                            const std::vector<RealT> &C) :
    SubgradientMethod(1000,
                              optimization_wrapper->GetComputationWrapper().ComputeSolutionNormBound(units, C, optimization_wrapper->GetOptions().GetRealValue("log_base")),
                              optimization_wrapper->GetComputationWrapper().ComputeGradientNormBound(units, C, optimization_wrapper->GetOptions().GetRealValue("log_base")),
                              Min(C)),
    InnerOptimizationWrapper(optimization_wrapper, units, C)
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod::ComputeFunction()
//
// Compute the regularized logloss using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////


RealT InnerOptimizationWrapperSubgradientMethod::ComputeFunction(const std::vector<RealT> &w)
{
    return this->optimization_wrapper->GetComputationWrapper().ComputeFunction(this->units, w, true, true, this->optimization_wrapper->GetOptions().GetRealValue("log_base"),
     this->optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
     this->optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
     this->optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"))
      + RealT(0.5) * DotProduct(this->C, w*w) + DotProduct(w, this->bias);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod::ComputeSubgradient()
//
// Compute the regularized logloss gradient using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////


void InnerOptimizationWrapperSubgradientMethod::ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &w)
{
    g = this->optimization_wrapper->GetComputationWrapper().ComputeGradient(this->units, w, true, true, this->optimization_wrapper->GetOptions().GetRealValue("log_base"),
     this->optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
     this->optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
     this->optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"))
      + this->C * w + this->bias;
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod::Report()
//
// Routines for printing results and messages from the optimizer.
//////////////////////////////////////////////////////////////////////


void InnerOptimizationWrapperSubgradientMethod::Report(int iteration, const std::vector<RealT> &w, RealT f, const std::vector<RealT> &g, RealT norm_bound, RealT step_size)
{
    // write results to disk
    this->optimization_wrapper->GetParameterManager().WriteToFile(SPrintF("optimize.params.iter%d", iteration), w);
    
    // write results to console
    this->optimization_wrapper->PrintMessage(SPrintF("Inner iteration %d: f = %lf (%lf), |w| = %lf, |g| = %lf, norm bound = %lf, step = %lf, efficiency = %lf%%", 
                                                     iteration, double(f), double(f - RealT(0.5) * DotProduct(this->C, w*w)),
                                                     double(Norm(w)), double(step_size), double(Norm(g)), double(norm_bound),
                                                     double(this->optimization_wrapper->GetComputationEngine().GetEfficiency())));
}


void InnerOptimizationWrapperSubgradientMethod::Report(const std::string &s) 
{
    this->optimization_wrapper->PrintMessage(SPrintF("Inner message: %s", s.c_str()));
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod::Minimize()
//
// Perform subgradient optimization.
//////////////////////////////////////////////////////////////////////


RealT InnerOptimizationWrapperSubgradientMethod::Minimize(std::vector<RealT> &x0)
{
    return SubgradientMethod::Minimize(x0);
}
