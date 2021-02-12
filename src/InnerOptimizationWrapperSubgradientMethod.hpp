//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod.hpp
//
// Inner optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPERSUBGRADIENTMETHOD_HPP
#define INNEROPTIMIZATIONWRAPPERSUBGRADIENTMETHOD_HPP

#include "InnerOptimizationWrapper.hpp"
#include "OptimizationWrapper.hpp"
#include "SubgradientMethod.hpp"

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapperSubgradientMethod
//////////////////////////////////////////////////////////////////////

class InnerOptimizationWrapperSubgradientMethod : public SubgradientMethod, public InnerOptimizationWrapper
{
public:
    InnerOptimizationWrapperSubgradientMethod(OptimizationWrapper *optimization_wrapper,
                                               const std::vector<int> &units,
                                               const std::vector<RealT> &C);
    
    RealT ComputeFunction(const std::vector<RealT> &x);
    void ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, const std::vector<RealT> &g,
                RealT norm_bound, RealT step_size);
    void Report(const std::string &s);
    RealT Minimize(std::vector<RealT> &x0);
};

#endif
