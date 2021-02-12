//////////////////////////////////////////////////////////////////////
// GammaMLE.hpp
//
// This file contains an implementation of MLE for a Gamma function
//////////////////////////////////////////////////////////////////////

#ifndef GammaMLE_HPP
#define GammaMLE_HPP

#include <vector>
#include <math.h>
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// GammaMLE()
//
// Implementation of MLE for Gamma function
//////////////////////////////////////////////////////////////////////

using Real = double;
class GammaMLE
{
    const int    MAX_ITER_K;
    const Real   THRESH_K;
    
public:
    GammaMLE
    (

        const int    MAX_ITER_K = 1000,
        const Real   THRESH_K = Real(1e-3)
    );
    
    virtual ~GammaMLE() {}
    
    Real Minimize(std::vector<Real> &x0, std::vector<std::vector<bool> > config_params, int which_data);
    Real Psi(Real k);
    Real PsiPrime(Real k);

    virtual double ComputeFunction(const std::vector<double> &x) = 0;
    virtual void ComputeGradient(std::vector<double> &g, const std::vector<double> &x) = 0;
    virtual void Report(int iteration, const std::vector<Real> &x, double f, const std::vector<Real> &g, double step_size) = 0;
    virtual void Report(const std::string &s) = 0;

    virtual double ComputeGammaMLEFunction(const std::vector<Real> &x, int i, int j, int k, Real scale, int which_data) = 0;
    virtual void ComputeGammaMLEGradient(std::vector<Real> &g, const std::vector<Real> &x, int i, int j, int k, Real scale, int which_data) = 0;
    virtual int GetLogicalIndex(int i, int j, int k, int which_data) = 0;
    virtual bool FindZerosInData(int i, int j, int which_data) = 0;
    virtual void ComputeGammaMLEScalingFactor(std::vector<Real> &g, const std::vector<Real> &w, int i, int j, int which_data) = 0;

};

#endif
