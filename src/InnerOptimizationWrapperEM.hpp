//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM.hpp
//
// Inner optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPEREM_HPP
#define INNEROPTIMIZATIONWRAPPEREM_HPP

#include "OptimizationWrapper.hpp"
#include "EM.hpp"
#include "GammaMLE.hpp"
#include "InnerOptimizationWrapper.hpp"

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapperEM
//////////////////////////////////////////////////////////////////////

class InnerOptimizationWrapperEM : public EM, public GammaMLE, public InnerOptimizationWrapper
{
    RealT log_base;
    
public:
    InnerOptimizationWrapperEM(OptimizationWrapper *optimization_wrapper,
                                  const std::vector<int> &units,
                                  const std::vector<RealT> &C);
    
    RealT ComputeFunction(const std::vector<RealT> &x);
    void ComputeGradient(std::vector<RealT> &g, const std::vector<RealT> &x);
    void Report(int iteration, const std::vector<RealT> &x, RealT f, const std::vector<RealT> &g, RealT step_size);
    void Report(const std::string &s);
    RealT OneStep(std::vector<RealT> &x0, int iter, std::vector<std::vector<bool> > config_params);
    RealT Minimize(std::vector<RealT> &x0) { Error("Not implemented"); return RealT(0); }

    RealT ComputeGammaMLEFunction(const std::vector<RealT> &x, int i, int j, int k, RealT scale, int which_data);
    void ComputeGammaMLEGradient(std::vector<RealT> &g, const std::vector<RealT> &x, int i, int j, int k, RealT scale, int which_data);
    int GetLogicalIndex(int i, int j, int k, int which_data);
    bool FindZerosInData(int i, int j, int which_data);
    void ComputeGammaMLEScalingFactor(std::vector<RealT> &g, const std::vector<RealT> &w, int i, int j, int which_data);

};
#endif
