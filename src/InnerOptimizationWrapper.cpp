//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper.ipp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

#include "InnerOptimizationWrapper.hpp"

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper::InnerOptimizationWrapper()
//
// Constructor.
//////////////////////////////////////////////////////////////////////


InnerOptimizationWrapper::InnerOptimizationWrapper(OptimizationWrapper *optimization_wrapper,
                                                          const std::vector<int> &units,
                                                          const std::vector<RealT> &weights_initial,
                                                          const std::vector<RealT> &C) :
    optimization_wrapper(optimization_wrapper),
    units(units),
    weights_initial(weights_initial),
    C(C),
    bias(optimization_wrapper->GetParameterManager().GetNumLogicalParameters())
{}


InnerOptimizationWrapper::InnerOptimizationWrapper(OptimizationWrapper *optimization_wrapper,
                                                          const std::vector<int> &units,
                                                          const std::vector<RealT> &C) :
    optimization_wrapper(optimization_wrapper),
    units(units),
    C(C),
    bias(optimization_wrapper->GetParameterManager().GetNumLogicalParameters())
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper::~InnerOptimizationWrapper()
//
// Destructor.
//////////////////////////////////////////////////////////////////////


InnerOptimizationWrapper::~InnerOptimizationWrapper()
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper::LoadBias()
//
// Load linear bias.
//////////////////////////////////////////////////////////////////////


void InnerOptimizationWrapper::LoadBias(const std::vector<RealT> &bias)
{
    this->bias = bias;
}
