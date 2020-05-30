//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper.hpp
//
// Generic class for performing optimization of model with respect to
// a particular set of parameters.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPER_HPP
#define INNEROPTIMIZATIONWRAPPER_HPP

#include "OptimizationWrapper.hpp"
#include "LBFGS.hpp"

template<class RealT>
class OptimizationWrapper;

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapper
//////////////////////////////////////////////////////////////////////

template<class RealT>
class InnerOptimizationWrapper 
{
    
protected:
    
    OptimizationWrapper<RealT> *optimization_wrapper;
    const std::vector<int> units;
    const std::vector<RealT> weights_initial;
    const std::vector<RealT> C;
    std::vector<RealT> bias;

public:
    
    InnerOptimizationWrapper(OptimizationWrapper<RealT> *optimization_wrapper,
                             const std::vector<int> &units,
                             const std::vector<RealT> &weights_initial,
                             const std::vector<RealT> &C);

    InnerOptimizationWrapper(OptimizationWrapper<RealT> *optimization_wrapper,
                             const std::vector<int> &units,
                             const std::vector<RealT> &C);

    virtual ~InnerOptimizationWrapper();

    void LoadBias(const std::vector<RealT> &bias);
    
    virtual RealT Minimize(std::vector<RealT> &x0) = 0;
};

#include "InnerOptimizationWrapper.ipp"

#endif
