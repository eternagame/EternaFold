//////////////////////////////////////////////////////////////////////
// CGLinear.hpp
//
// This file contains an implementation of the conjugate gradient
// algorithm for solving linear systems.
//////////////////////////////////////////////////////////////////////

#ifndef CGLINEAR_HPP
#define CGLINEAR_HPP

#include <vector>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// CGLinear()
//
// Implementation of conjugate gradient for solving linear
// systems Ax = b, where A is symmetric positive definite.
//////////////////////////////////////////////////////////////////////

using Real = double;
class CGLinear {
public:
    Real Minimize
    (
        const std::vector<Real> &b,                              // right hand side
        std::vector<Real> &x,                                    // initial parameter vector
        
        const int    MAX_ITERATIONS              = 1000,         // maximum number of iterations to run CG
        const Real   SMALL_STEP_RATIO            = 0.001,        // ratio beneath which steps are considered "small"
        const int    MAX_SMALL_STEPS             = 5             // maximum number of small steps before we quit
    );

    virtual ~CGLinear() {}

    virtual void ComputeAx(std::vector<double> &Ax, const std::vector<double> &x) = 0;
    virtual void Report(int iteration, const std::vector<double> &x, double f, double step_size) = 0;
    virtual void Report(const std::string &s) = 0;
};

#endif
