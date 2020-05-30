//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM.ipp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::InnerOptimizationWrapperEM()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
InnerOptimizationWrapperEM<RealT>::InnerOptimizationWrapperEM(OptimizationWrapper<RealT> *optimization_wrapper,
                                                                    const std::vector<int> &units,
                                                                    const std::vector<RealT> &C) :
    EM<RealT>(0.00001, 0),
    GammaMLE<RealT>(1000,1e-3),
    InnerOptimizationWrapper<RealT>(optimization_wrapper, units, C),
    log_base(optimization_wrapper->GetOptions().GetRealValue("log_base"))
{}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::ComputeFunction()
//
// Compute the regularized logloss using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT InnerOptimizationWrapperEM<RealT>::ComputeFunction(const std::vector<RealT> &w)
{
    return this->optimization_wrapper->GetComputationWrapper().ComputeFunction(this->units, w, false, true, log_base, 
        this->optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
        this->optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
        this->optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"))
         + RealT(0.5) * DotProduct(this->C, w*w) + DotProduct(w, this->bias);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::ComputeGradient()
//
// Compute the regularized logloss gradient using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapperEM<RealT>::ComputeGradient(std::vector<RealT> &g, const std::vector<RealT> &w)
{
    g = this->optimization_wrapper->GetComputationWrapper().ComputeGradient(this->units, w, false, true, log_base,
        this->optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
        this->optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
        this->optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"))
      + this->C * w + this->bias;
}


// Functions for learning the evidence CPD

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::ComputeGammaMLEFunction()
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT InnerOptimizationWrapperEM<RealT>::ComputeGammaMLEFunction(const std::vector<RealT> &w, int i, int j, int k, RealT scale, int which_data)
{
    return this->optimization_wrapper->GetComputationWrapper().ComputeGammaMLEFunction(this->units, w, false, true, log_base, i, j, k, scale, which_data);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::ComputeGammaMLEGradient()
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapperEM<RealT>::ComputeGammaMLEGradient(std::vector<RealT> &g, const std::vector<RealT> &w, int i, int j, int k, RealT scale, int which_data)
{
    g = this->optimization_wrapper->GetComputationWrapper().ComputeGammaMLEGradient(this->units, w, false, true, log_base, i, j, k, scale, which_data);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::ComputeGammaMLEScalingFactor()
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapperEM<RealT>::ComputeGammaMLEScalingFactor(std::vector<RealT> &g, const std::vector<RealT> &w, int i, int j, int which_data)
{
    g = this->optimization_wrapper->GetComputationWrapper().ComputeGammaMLEScalingFactor(this->units, w, i, j, which_data);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::GetLogicalIndex()
//
// return the logical index in w to which this parameter corresponds
//////////////////////////////////////////////////////////////////////

template<class RealT>
int InnerOptimizationWrapperEM<RealT>::GetLogicalIndex(int i, int j, int k, int which_data)
{
    // dim0: k or theta and dim1: A,C,G,T and dim2: paired (0),unpaired (1)
    std::pair<RealT,RealT>* param = this->optimization_wrapper->GetComputationWrapper().GetInferenceEngine().GetLogScoreEvidence(i,j,k,which_data);
    return this->optimization_wrapper->GetComputationWrapper().GetParameterManager().GetLogicalIndex(param);

}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::FindZerosInData()
//
// returns true if there are sequences that have 0-data positions corresponding to that CPD
//////////////////////////////////////////////////////////////////////
template<class RealT>
bool InnerOptimizationWrapperEM<RealT>::FindZerosInData(int i, int j, int which_data)
{
    // dim1: A,C,G,T and dim2: paired (0),unpaired (1)
    return this->optimization_wrapper->GetComputationWrapper().FindZerosInData(this->units, i, j, which_data);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::Report()
//
// Routines for printing results and messages from the optimizer.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapperEM<RealT>::Report(int iteration, const std::vector<RealT> &w, RealT f, const std::vector<RealT> &g, RealT step_size)
{
    // write results to disk
    this->optimization_wrapper->GetParameterManager().WriteToFile(SPrintF("optimize.params.iter%d", iteration), w);
    
    // write results to console
    this->optimization_wrapper->PrintMessage(SPrintF("Inner iteration %d: f = %lf (%lf) [%lf], |w| = %lf, |g| = %lf, step = %lf, efficiency = %lf%%", 
                                                     iteration, double(f), double(f - RealT(0.5) * DotProduct(this->C, w*w)),
                                                     double(this->optimization_wrapper->GetComputationWrapper().ComputeFunction(this->units, w, true, true, log_base,
                                                      this->optimization_wrapper->GetOptions().GetRealValue("hyperparam_data"),
                                                      this->optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data"),
                                                      this->optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"))
                                                       + RealT(0.5) * DotProduct(this->C, w*w) + DotProduct(w, this->bias)),
                                                     double(Norm(w)), double(Norm(g)), double(step_size),
                                                     double(this->optimization_wrapper->GetComputationEngine().GetEfficiency())));
}

template<class RealT>
void InnerOptimizationWrapperEM<RealT>::Report(const std::string &s) 
{
    this->optimization_wrapper->PrintMessage(SPrintF("Inner message: %s", s.c_str()));
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperEM::OneStep()
//
// Perform EM optimization.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT InnerOptimizationWrapperEM<RealT>::OneStep(std::vector<RealT> &x0, int iter, std::vector<std::vector<bool> > config_params)
{

    RealT result = 0;
    int num_data_sources = this->optimization_wrapper->GetOptions().GetIntValue("num_data_sources");
    for (int i = 0; i < num_data_sources; i++)
    {
        result = result + GammaMLE<RealT>::Minimize(x0,config_params,i);
    }

    result = result + EM<RealT>::OneStep(x0, iter);

    return result;
        
}
