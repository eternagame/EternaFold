//////////////////////////////////////////////////////////////////////
// SubgradientMethod.hpp
//
// This file contains an implementation of the subgradient
// optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef SUBGRADIENTMETHOD_HPP
#define SUBGRADIENTMETHOD_HPP

#include <vector>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// SubgradientMethod()
//
// Implementation of subgradient optimization routine.
//////////////////////////////////////////////////////////////////////

using RealT = double;
class SubgradientMethod
{
    const int   MAX_ITERATIONS;
    const RealT PARAMETER_NORM_BOUND;
    const RealT GRADIENT_NORM_BOUND;
    const RealT CURVATURE;

public:
    SubgradientMethod
    (
        const int     MAX_ITERATIONS                 = 1000,          // maximum number of iterations to run subgradient method
        const RealT   PARAMETER_NORM_BOUND           = RealT(1e-5),   // maximum parameter vector norm
        const RealT   GRADIENT_NORM_BOUND            = RealT(1e-5),   // maximum gradient vector norm
        const RealT   CURVATURE                      = RealT(0)       // strong convexity constant: f_t(v) >= f_t(w) + g'*(v - w) + 0.5 * CURVATURE * ||v - w||^2
    );
    
    virtual ~SubgradientMethod() {}
    
    RealT Minimize(std::vector<RealT> &x0);
    
    virtual RealT ComputeFunction(const std::vector<RealT> &x) = 0;
    virtual void ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &x) = 0;
    virtual void Report(int iteration, const std::vector<RealT> &x, RealT f, const std::vector<RealT> &g,
                        RealT norm_bound, RealT step_size) = 0;
    virtual void Report(const std::string &s) = 0;
};

#endif
