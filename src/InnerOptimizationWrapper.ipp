//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper.ipp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper<RealT>::InnerOptimizationWrapper()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
InnerOptimizationWrapper<RealT>::InnerOptimizationWrapper(OptimizationWrapper<RealT> *optimization_wrapper,
                                                          const std::vector<int> &units,
                                                          const std::vector<RealT> &weights_initial,
                                                          const std::vector<RealT> &C) :
    optimization_wrapper(optimization_wrapper),
    units(units),
    weights_initial(weights_initial),
    C(C),
    bias(optimization_wrapper->GetParameterManager().GetNumLogicalParameters())
{}

template<class RealT>
InnerOptimizationWrapper<RealT>::InnerOptimizationWrapper(OptimizationWrapper<RealT> *optimization_wrapper,
                                                          const std::vector<int> &units,
                                                          const std::vector<RealT> &C) :
    optimization_wrapper(optimization_wrapper),
    units(units),
    C(C),
    bias(optimization_wrapper->GetParameterManager().GetNumLogicalParameters())
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper<RealT>::~InnerOptimizationWrapper()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
InnerOptimizationWrapper<RealT>::~InnerOptimizationWrapper()
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapper<RealT>::LoadBias()
//
// Load linear bias.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapper<RealT>::LoadBias(const std::vector<RealT> &bias)
{
    this->bias = bias;
}
