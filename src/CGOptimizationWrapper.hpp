//////////////////////////////////////////////////////////////////////
// CGOptimizationWrapper.hpp
//
// Conjugate gradient optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef CGOPTIMIZATIONWRAPPER_HPP
#define CGOPTIMIZATIONWRAPPER_HPP

#include "OptimizationWrapper.hpp"
#include "CGLinear.hpp"

template<class RealT>
class OptimizationWrapper;

//////////////////////////////////////////////////////////////////////
// class CGOptimizationWrapper
//////////////////////////////////////////////////////////////////////

template<class RealT>
class CGOptimizationWrapper : public CGLinear<RealT>
{
    OptimizationWrapper<RealT> *optimization_wrapper;
    const std::vector<int> units;
    const std::vector<RealT> w;
    const std::vector<RealT> C;
    
public:
    CGOptimizationWrapper(OptimizationWrapper<RealT> *optimizer,
                          const std::vector<int> &units,
                          const std::vector<RealT> &w,
                          const std::vector<RealT> &C);
    
    void ComputeAx(std::vector<RealT> &Ax, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, RealT step_size);
    void Report(const std::string &s);
};

#include "CGOptimizationWrapper.ipp"

#endif
