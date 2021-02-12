//////////////////////////////////////////////////////////////////////
// EM.hpp
//
// This file contains an implementation of a variant of EM.
// Specifically, the OneStep function will compute the ESS and then
// take a SINGLE gradient step in that direction; this is a G-EM algo.
//////////////////////////////////////////////////////////////////////

#ifndef EM_HPP
#define EM_HPP

#include <vector>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// EM()
//
// Implementation of EM optimization routine.
//////////////////////////////////////////////////////////////////////

using Real = double;
class EM
{
    const Real s0;
    const Real s1;
    
public:
    EM
    (
        const Real s0           = Real(1),     // step size is s0 / (1 + iter)^s1
        const Real s1           = Real(0.6)
    );
    
    virtual ~EM() {}
    
    Real OneStep(std::vector<Real> &x0, int iter);

    virtual double ComputeFunction(const std::vector<double> &x) = 0;
    virtual void ComputeGradient(std::vector<double> &g, const std::vector<double> &x) = 0;
    virtual void Report(int iteration, const std::vector<Real> &x, double f, const std::vector<Real> &g, double step_size) = 0;
    virtual void Report(const std::string &s) = 0;
};

#endif
