//////////////////////////////////////////////////////////////////////
// LBFGS.hpp
//
// This file contains an implementation of the standard
// limited-memory BFGS optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef LBFGS_HPP
#define LBFGS_HPP

#include <vector>
#include "Utilities.hpp"
#include "LineSearch.hpp"

//////////////////////////////////////////////////////////////////////
// LBFGS()
//
// Implementation of L-BFGS optimization routine.
//////////////////////////////////////////////////////////////////////

using Real = double;

class LBFGS : public LineSearch
{
    const int M;
    const Real TERMINATION_RATIO;
    const int MAX_ITERATIONS;
    const Real SMALL_STEP_RATIO;
    const int MAX_SMALL_STEPS;
    const Real MAX_STEP_NORM;
    
public:
    LBFGS
    (
        const int   M                           = 20,             // number of previous gradients to remember
        const Real  TERMINATION_RATIO           = Real(1e-5),     // required ratio of gradient norm to parameter norm for termination
        const int   MAX_ITERATIONS              = 1000,           // maximum number of iterations to run L-BFGS
        const Real  SMALL_STEP_RATIO            = Real(1e-5),     // ratio beneath which steps are considered "small"
        const int   MAX_SMALL_STEPS             = 3,              // maximum number of small steps before we quit
        const Real  MAX_STEP_NORM               = Real(1e10)      // maximum norm for a single step
    );
    
    virtual ~LBFGS() {}
    
    Real Minimize(std::vector<Real> &x0);

    virtual double ComputeFunction(const std::vector<double> &x) = 0;
    virtual void ComputeGradient(std::vector<double> &g, const std::vector<double> &x) = 0;
    virtual void Report(int iteration, const std::vector<double> &x, double f, double step_size) = 0;
    virtual void Report(const std::string &s) = 0;
};

#endif
