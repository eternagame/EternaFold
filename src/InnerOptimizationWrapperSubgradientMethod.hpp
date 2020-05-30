//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperSubgradientMethod.hpp
//
// Inner optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPERSUBGRADIENTMETHOD_HPP
#define INNEROPTIMIZATIONWRAPPERSUBGRADIENTMETHOD_HPP

#include "OptimizationWrapper.hpp"
#include "SubgradientMethod.hpp"

template<class RealT>
class SubgradientMethod;

template<class RealT>
class OptimizationWrapper;

template<class RealT>
class InnerOptimizationWrapper;

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapperSubgradientMethod
//////////////////////////////////////////////////////////////////////

template<class RealT>
class InnerOptimizationWrapperSubgradientMethod : public SubgradientMethod<RealT>, public InnerOptimizationWrapper<RealT>
{
public:
    InnerOptimizationWrapperSubgradientMethod(OptimizationWrapper<RealT> *optimization_wrapper,
                                               const std::vector<int> &units,
                                               const std::vector<RealT> &C);
    
    RealT ComputeFunction(const std::vector<RealT> &x);
    void ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, const std::vector<RealT> &g,
                RealT norm_bound, RealT step_size);
    void Report(const std::string &s);
    RealT Minimize(std::vector<RealT> &x0);
};

#include "InnerOptimizationWrapperSubgradientMethod.ipp"

#endif
