//////////////////////////////////////////////////////////////////////
// SubgradientMethod.hpp
//
// This file contains an implementation of the bundle
// optimization algorithm.
//////////////////////////////////////////////////////////////////////
#ifndef BUNDLEMETHOD_HPP
#define BUNDLEMETHOD_HPP

#include <vector>
#include "Utilities.hpp"

#include "utilities/sml.hpp"
#include "utilities/common.hpp"
#include "utilities/timer.hpp"
#include "utilities/configuration.hpp"
#include "utilities/bmrmexception.hpp"
#include "solver/bmrminnersolver/bmrminnersolver.hpp"
#include "solver/bmrminnersolver/l2n2_daifletcherpgm.hpp"
#include "solver/bmrminnersolver/l2n2_prloqo.hpp"
#include <fstream>
#include <sstream>
using namespace std;

//////////////////////////////////////////////////////////////////////
// BundleMethod()
//
// Implementation of bundle optimization routine.
//////////////////////////////////////////////////////////////////////

class BundleMethod
{
  const int   MAX_ITERATIONS;
  double lambda;


public:
    virtual ~BundleMethod() {}
  BundleMethod
    (
        const int     MAX_ITERATIONS                 = 1000,          // maximum number of iterations to run subgradient method
        const double  lambda                         = 1
    );

    RealT Minimize(std::vector<RealT> &x0);
    
    virtual RealT ComputeFunction(const std::vector<RealT> &x) = 0;
    virtual void ComputeSubgradient(std::vector<RealT> &g, const std::vector<RealT> &x) = 0;
    virtual void Report(int iteration, const std::vector<RealT> &x, RealT f, const std::vector<RealT> &g,
                        RealT norm_bound, RealT step_size) = 0;
    virtual void Report(const std::string &s) = 0;

};

#include "BundleMethod.ipp"

#endif

