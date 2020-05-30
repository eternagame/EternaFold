//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperBundleMethod.hpp
//
// Inner optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPERBUNDLEMETHOD_HPP
#define INNEROPTIMIZATIONWRAPPERBUNDLEMETHOD_HPP

#include "OptimizationWrapper.hpp"
#include "BundleMethod.hpp"

template<class RealT>
class BundleMethod;

template<class RealT>
class OptimizationWrapper;

template<class RealT>
class InnerOptimizationWrapper;

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapperBundleMethod
//////////////////////////////////////////////////////////////////////

template<class RealT>
class InnerOptimizationWrapperBundleMethod : public BundleMethod<RealT>, public InnerOptimizationWrapper<RealT>
{
public:
    InnerOptimizationWrapperBundleMethod(OptimizationWrapper<RealT> *optimization_wrapper,
                                               const std::vector<int> &units,
                                               const std::vector<RealT> &C);
    
    RealT ComputeFunction(const std::vector<RealT> &x);
    void ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, const std::vector<RealT> &g,
                RealT norm_bound, RealT step_size);
    void Report(const std::string &s);
    RealT Minimize(std::vector<RealT> &x0);
};

#include "InnerOptimizationWrapperBundleMethod.ipp"

#endif
