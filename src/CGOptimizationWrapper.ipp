//////////////////////////////////////////////////////////////////////
// CGOptimizationWrapper.ipp
//
// CG optimization code.
//////////////////////////////////////////////////////////////////////

#include "CGOptimizationWrapper.hpp"

//////////////////////////////////////////////////////////////////////
// CGOptimizationWrapper<RealT>::CGOptimizationWrapper()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
CGOptimizationWrapper<RealT>::CGOptimizationWrapper(OptimizationWrapper<RealT> *optimization_wrapper,
                                                    const std::vector<int> &units,
                                             const std::vector<RealT> &w,
                                                    const std::vector<RealT> &C) :
    CGLinear<RealT>(), optimization_wrapper(optimization_wrapper), units(units), w(w), C(C)
{}

//////////////////////////////////////////////////////////////////////
// CGOptimizationWrapper<RealT>::ComputeAx()
//
// Compute Hessian-vector product.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void CGOptimizationWrapper<RealT>::ComputeAx(std::vector<RealT> &Ax, const std::vector<RealT> &x)
{
    std::vector<RealT> Ce = optimization_wrapper->GetParameterManager().ExpandParameterGroupValues(C);
    Ax = optimization_wrapper->GetComputationWrapper().ComputeHessianVectorProduct(units, w, x, true, optimization_wrapper->GetOptions().GetRealValue("log_base")) + Ce * x;
}

//////////////////////////////////////////////////////////////////////
// CGOptimizationWrapper<RealT>::Report()
//
// Provide progress report on CG algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void CGOptimizationWrapper<RealT>::Report(int iteration, const std::vector<RealT> &x, RealT f, RealT step_size)
{
    optimization_wrapper->PrintMessage(SPrintF("CG iteration %d: f = %lf, |x| = %lf, step = %lf, efficiency = %lf%%", 
                                               iteration, double(f), double(Norm(x)), double(step_size),
                                               double(optimization_wrapper->GetComputationEngine().GetEfficiency())));
}

template<class RealT>
void CGOptimizationWrapper<RealT>::Report(const std::string &s)
{
    optimization_wrapper->PrintMessage(SPrintF("CG message: %s", s.c_str()));
}
