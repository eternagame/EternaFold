//////////////////////////////////////////////////////////////////////
// EM.ipp
//
// This file contains an implementation of a variant of EM.
//////////////////////////////////////////////////////////////////////

#include <EM.hpp>

//////////////////////////////////////////////////////////////////////
// EM::EM()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

EM::EM
(
    const Real  s0,                            // step size is
    const Real  s1                             //   s0/(1+iter)^s1
) :
    s0(s0),
    s1(s1)
{}

//////////////////////////////////////////////////////////////////////
// EM::OneStep()
//
// Implementation of one E, one gradient step
//////////////////////////////////////////////////////////////////////

Real EM::OneStep(std::vector<Real> &x0, int iter)
{
    // initialize
    const int n = int(x0.size());
    std::vector<Real> g(n);
    double f = 0;

    Real stepsize = s0 * pow(1.0 + iter, -s1);   

    ComputeGradient(g, x0);
    f = ComputeFunction(x0); // this should be cached from the gradient
    
    // Report before taking gradient step
    Report(iter, x0, f, g, stepsize);

    x0 = x0 - stepsize*g; 

    // Ideally we would run ComputeFunction again, but that means more inference
    return f;
}   

