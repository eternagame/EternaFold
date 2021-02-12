//////////////////////////////////////////////////////////////////////
// OuterOptimizationWrapper.hpp
//
// Outer optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef OUTEROPTIMIZATIONWRAPPER_HPP
#define OUTEROPTIMIZATIONWRAPPER_HPP

#include "OptimizationWrapper.hpp"
#include "CGOptimizationWrapper.hpp"
#include "InnerOptimizationWrapper.hpp"
#include "OuterOptimizationWrapper.hpp"
#include "LBFGS.hpp"

//////////////////////////////////////////////////////////////////////
// class OuterOptimizationWrapper
//////////////////////////////////////////////////////////////////////

class OuterOptimizationWrapper : public LBFGS
{
    OptimizationWrapper *optimization_wrapper;
    const std::vector<RealT> initial_w;
    const std::vector<int> training;
    const std::vector<int> holdout;
    
public:
    OuterOptimizationWrapper(OptimizationWrapper *optimization_wrapper,
                             const std::vector<RealT> &initial_w,
                             const std::vector<int> &training,
                             const std::vector<int> &holdout);

    RealT ComputeFunction(const std::vector<RealT> &log_C);
    void ComputeGradient(std::vector<RealT> &g, const std::vector<RealT> &log_C);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, RealT step_size);
    void Report(const std::string &s);
};

#endif
