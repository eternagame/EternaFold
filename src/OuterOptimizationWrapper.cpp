//////////////////////////////////////////////////////////////////////
// OuterOptimizationWrapper.cpp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

#include "OuterOptimizationWrapper.hpp"

//////////////////////////////////////////////////////////////////////
// OuterOptimizationWrapper::OuterOptimizationWrapper()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

OuterOptimizationWrapper::OuterOptimizationWrapper(OptimizationWrapper *optimization_wrapper,
                                                          const std::vector<RealT> &initial_w,
                                                          const std::vector<int> &training,
                                                          const std::vector<int> &holdout):
    LBFGS(20,1e-5,100,1e-5,3,1),
    optimization_wrapper(optimization_wrapper),
    initial_w(initial_w),
    training(training),
    holdout(holdout)
{}

//////////////////////////////////////////////////////////////////////
// OuterOptimizationWrapper::ComputeFunction()
//
// Compute function for outer iteration.
//////////////////////////////////////////////////////////////////////

RealT OuterOptimizationWrapper::ComputeFunction(const std::vector<RealT> &log_C)
{
    std::ostringstream oss;
    oss << "Computing outer function using C = " << Exp(log_C);
    optimization_wrapper->PrintMessage(oss.str());
    optimization_wrapper->Indent();
    
    // w = solution of OPT1 for current C
    std::vector<RealT> w = initial_w;
    optimization_wrapper->PrintMessage("Solving OPT1...");
    optimization_wrapper->Indent();
    optimization_wrapper->Train(training, w, w, Exp(log_C)); //HKWS added another w
    optimization_wrapper->Unindent();
    
    // compute holdout logloss
    RealT ret = optimization_wrapper->GetComputationWrapper().ComputeFunction(holdout, w, false, true, optimization_wrapper->GetOptions().GetRealValue("log_base"),
     optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
     optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
     optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"));
    
    optimization_wrapper->Unindent();
    optimization_wrapper->PrintMessage(SPrintF("Finished outer function: %lf", ret));
    return ret;
}

//////////////////////////////////////////////////////////////////////
// OuterOptimizationWrapper::ComputeGradient()
//
// Compute the regularized logloss gradient using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

void OuterOptimizationWrapper::ComputeGradient(std::vector<RealT> &g, const std::vector<RealT> &log_C)
{
    const std::vector<RealT> C = Exp(log_C);

    std::ostringstream oss;
    oss << "Computing outer gradient using C = " << C;
    optimization_wrapper->PrintMessage(oss.str());
    optimization_wrapper->Indent();
    
    // show current set of hyperparameters
    optimization_wrapper->PrintMessage("Current hyperparameters:");
    optimization_wrapper->Indent();
    const std::vector<ParameterGroup> &groups = optimization_wrapper->GetParameterManager().GetParameterGroups();
    for (size_t i = 0; i < groups.size(); i++)
        optimization_wrapper->PrintMessage(SPrintF("Hyperparameter group %d (%s): %lf", i+1, groups[i].name.c_str(), C[i]));
    optimization_wrapper->Unindent();
    
    // w = solution of OPT1 for current C
    std::vector<RealT> w = initial_w;
    optimization_wrapper->PrintMessage("Solving OPT1...");
    optimization_wrapper->Indent();
    optimization_wrapper->Train(training, w, w, C);
    optimization_wrapper->Unindent();
    
    // compute holdout logloss
    std::vector<RealT> holdout_gradient = optimization_wrapper->GetComputationWrapper().ComputeGradient(holdout, w, false, true, optimization_wrapper->GetOptions().GetRealValue("log_base"),
        optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
        optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
        optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"));
    
    // solve linear system
    CGOptimizationWrapper cg_linear(optimization_wrapper, training, w, C);
    std::vector<RealT> x(holdout_gradient.size());
    
    optimization_wrapper->PrintMessage("Solving linear system...");
    optimization_wrapper->Indent();
    cg_linear.Minimize(holdout_gradient,x);
    optimization_wrapper->Unindent();
    
    // form "B" matrix    
    const std::vector<RealT> log_C_grad = Exp(log_C);
    std::vector<std::vector<RealT> > B(x.size(), std::vector<RealT>(optimization_wrapper->GetParameterManager().GetNumParameterGroups()));
    for (size_t i = 0; i < groups.size(); i++)
        for (int j = groups[i].begin; j < groups[i].end; j++)
            B[j][i] = w[j] * log_C_grad[i];
    
    // compute gradient
    g.clear();
    g.resize(log_C.size());
    for (size_t i = 0; i < B.size(); i++)
        g -= x[i] * B[i];
    
    optimization_wrapper->Unindent();
    optimization_wrapper->PrintMessage(SPrintF("Finished outer gradient: norm = %lf", Norm(g)));
}

//////////////////////////////////////////////////////////////////////
// OuterOptimizationWrapper::Report()
//
// Routines for printing results and messages from the optimization_wrapper.
//////////////////////////////////////////////////////////////////////

void OuterOptimizationWrapper::Report(int iteration, const std::vector<RealT> &log_C, RealT f, RealT step_size)
{
    std::ostringstream oss;
    oss << "Outer iteration " << iteration << ": holdout f = " << f << ", C = " << Exp(log_C)
        << ", step length = " << step_size << ", efficiency = " << optimization_wrapper->GetComputationEngine().GetEfficiency() << "%";
    optimization_wrapper->PrintMessage(oss.str());
}

void OuterOptimizationWrapper::Report(const std::string &s) 
{
    optimization_wrapper->PrintMessage(SPrintF("Outer message: %s", s.c_str()));
}

