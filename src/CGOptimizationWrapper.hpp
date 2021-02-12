//////////////////////////////////////////////////////////////////////
// CGOptimizationWrapper.hpp
//
// Conjugate gradient optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef CGOPTIMIZATIONWRAPPER_HPP
#define CGOPTIMIZATIONWRAPPER_HPP

#include "OptimizationWrapper.hpp"
#include "CGLinear.hpp"

class OptimizationWrapper;

//////////////////////////////////////////////////////////////////////
// class CGOptimizationWrapper
//////////////////////////////////////////////////////////////////////

class CGOptimizationWrapper : public CGLinear
{
    OptimizationWrapper *optimization_wrapper;
    const std::vector<int> units;
    const std::vector<RealT> w;
    const std::vector<RealT> C;
    
public:
    CGOptimizationWrapper(OptimizationWrapper *optimizer,
                          const std::vector<int> &units,
                          const std::vector<RealT> &w,
                          const std::vector<RealT> &C);
    
    void ComputeAx(std::vector<RealT> &Ax, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, RealT step_size);
    void Report(const std::string &s);
};

#endif
