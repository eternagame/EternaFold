//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperLBFGS.hpp
//
// Inner optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPERLBFGS_HPP
#define INNEROPTIMIZATIONWRAPPERLBFGS_HPP

#include "InnerOptimizationWrapper.hpp"
#include "OptimizationWrapper.hpp"
#include "LBFGS.hpp"

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapperLBFGS
//////////////////////////////////////////////////////////////////////

class InnerOptimizationWrapperLBFGS : public LBFGS, public InnerOptimizationWrapper
{
    RealT log_base;
    RealT hyperparam_data;
    RealT kd_hyperparam_data;
    RealT lig_hyperparam_data;

public:
    InnerOptimizationWrapperLBFGS(OptimizationWrapper *optimization_wrapper,
                                  const std::vector<int> &units,
                                  const std::vector<RealT> &weights_initial,
                                  const std::vector<RealT> &C);
    
    RealT ComputeFunction(const std::vector<RealT> &x);
    void ComputeGradient(std::vector<RealT> &g, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, RealT step_size);
    void Report(const std::string &s);
    RealT Minimize(std::vector<RealT> &x0);
};

#endif
