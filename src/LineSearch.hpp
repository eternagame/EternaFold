//////////////////////////////////////////////////////////////////////
// LineSearch.hpp
//
// This file contains an implementation of a safeguarded
// backtracking line search algorithm for cubic polynomial
// interpolation.
//////////////////////////////////////////////////////////////////////

#ifndef LINESEARCH_HPP
#define LINESEARCH_HPP

#include <vector>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// LineSearchFunctor()
//
// Implementation of backtracking line-search.  Returns a
// value of t, corresponding to the step size needed for
// updating
//
//    x <-- x + t*d
//
// As a side effect, this function also updates the value of
// the function f(x) to its new value f(x + t*d).
//////////////////////////////////////////////////////////////////////

using Real = double;

class LineSearch
{
    const Real T_INIT;
    const Real MU;
    const Real MIN_IMPROVEMENT_RATIO;
    const int MAX_EVALUATIONS;
    const Real GAMMA1;
    const Real GAMMA2;    
    
public:
    LineSearch
    (
        const Real T_INIT                = Real(1),         // initial step size
        const Real MU                    = Real(0.001),     // phi(t) <= phi(0) + MU * phi'(0) t      (sufficient decrease)
        const Real MIN_IMPROVEMENT_RATIO = Real(0.1),       // minimum proportion of overall improvement required to keep going
        const int MAX_EVALUATIONS        = 10,              // maximum number of function evaluations
        const Real GAMMA1                = Real(0.01),      // maximum step length shrinkage 
        const Real GAMMA2                = Real(0.8)        // minimum step length shrinkage 
    );

    virtual ~LineSearch() {}
    
    Real DoLineSearch
    (
        const std::vector<Real> &x,                        // initial parameters
        const Real f,                                      // initial function value
        const std::vector<Real> &g,                        // initial gradient vector
        const std::vector<Real> &d,                        // initial direction vector
        
        std::vector<Real> &new_x,                          // new parameters
        Real &new_f,                                       // new function value
        std::vector<Real> &new_g,                          // new gradient vector
        
        const Real T_MIN,                                  // minimum step size
        const Real T_MAX                                   // maximum step size
    );
    
    virtual double ComputeFunction(const std::vector<double> &x) = 0;
    virtual void ComputeGradient(std::vector<double> &g, const std::vector<double> &x) = 0;
};

#endif
