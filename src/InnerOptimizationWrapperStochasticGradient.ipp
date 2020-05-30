//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient.ipp
//
// Implementation of functors needed for optimization.
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::InnerOptimizationWrapperStochasticGradient()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
InnerOptimizationWrapperStochasticGradient<RealT>::InnerOptimizationWrapperStochasticGradient(OptimizationWrapper<RealT> *optimization_wrapper,
                                                                                              const std::vector<int> &units,
                                                                                              const std::vector<RealT> &C
                                                                                              ) :
    InnerOptimizationWrapper<RealT>(optimization_wrapper, units, C),
    log_base(optimization_wrapper->GetOptions().GetRealValue("log_base")),
    batch_size(optimization_wrapper->GetOptions().GetIntValue("batch_size")),
    MAX_ITERATIONS(optimization_wrapper->GetOptions().GetIntValue("train_max_iter")),
    s0(optimization_wrapper->GetOptions().GetRealValue("s0")),
    s1(optimization_wrapper->GetOptions().GetRealValue("s1")),
    hyperparam_data(optimization_wrapper->GetOptions().GetRealValue("hyperparam_data")),
    kd_hyperparam_data(optimization_wrapper->GetOptions().GetRealValue("kd_hyperparam_data")),
    lig_hyperparam_data(optimization_wrapper->GetOptions().GetRealValue("lig_hyperparam_data"))

{
    std::srand(GetSystemTime());
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::ComputeFunction()
//
// Compute the regularized logloss using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT InnerOptimizationWrapperStochasticGradient<RealT>::ComputeFunction(const std::vector<RealT> &w)
{
    return this->optimization_wrapper->GetComputationWrapper().ComputeFunctionSE(this->units, w, false, true, log_base, hyperparam_data, kd_hyperparam_data,lig_hyperparam_data) + RealT(0.5) * DotProduct(this->C, w*w) + DotProduct(w, this->bias);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::ComputeGradient()
//
// Compute the regularized logloss gradient using a particular
// parameter set and fixed regularization hyperparameters.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapperStochasticGradient<RealT>::ComputeGradient(std::vector<RealT> &g, const std::vector<RealT> &w, const int batch_size)
{
    const int num_examples = this->units.size();
    std::vector<int> units;

    if (batch_size == 0 || batch_size >= num_examples) {
        units = this->units;
    } else {
        for (int i = 0; i < batch_size; i++) {
           units.push_back(RandInt(num_examples));
        }
    }
    // TODO - do we need to rescale the gradient by 1/num_examples?
    g = this->optimization_wrapper->GetComputationWrapper().ComputeGradientSE(units, w, false, true, log_base, hyperparam_data, kd_hyperparam_data,lig_hyperparam_data) + this->C * w + this->bias;
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::GetLogicalIndex()
//
// return the logical index in w to which this parameter corresponds
//////////////////////////////////////////////////////////////////////

template<class RealT>
int InnerOptimizationWrapperStochasticGradient<RealT>::GetLogicalIndex(int i, int j, int k, int which_data)
{
    // dim0: k or theta and dim1: A,C,G,T and dim2: paired (0),unpaired (1)
//    const std::pair<RealT,RealT> param = this->optimization_wrapper->GetComputationWrapper().GetInferenceEngine().GetScoreEvidence()[i][j][k];
//    this->optimization_wrapper->GetComputationWrapper().GetInferenceEngine().ComputePosterior();  // return the logical index in w to which this 
    std::pair<RealT,RealT>* param = this->optimization_wrapper->GetComputationWrapper().GetInferenceEngine().GetLogScoreEvidence(i,j,k,which_data);
    return this->optimization_wrapper->GetComputationWrapper().GetParameterManager().GetLogicalIndex(param);

}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::FindZerosInData()
//
// returns true if there are sequences that have 0-data positions corresponding to that CPD
//////////////////////////////////////////////////////////////////////
template<class RealT>
bool InnerOptimizationWrapperStochasticGradient<RealT>::FindZerosInData(int i, int j, int which_data)
{
    // dim1: A,C,G,T and dim2: paired (0),unpaired (1)
    return this->optimization_wrapper->GetComputationWrapper().FindZerosInData(this->units, i, j, which_data);
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::Report()
//
// Routines for printing results and messages from the optimizer.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InnerOptimizationWrapperStochasticGradient<RealT>::Report(int iteration, const std::vector<RealT> &w, RealT step_size)
{
    // write results to disk
    this->optimization_wrapper->GetParameterManager().WriteToFile(SPrintF("optimize.params.iter%d", iteration), w);

    std::vector<RealT> g;
    ComputeGradient(g, w, 0);
    RealT f = ComputeFunction(w);

    // write results to console
    this->optimization_wrapper->PrintMessage(SPrintF("Inner iteration %d: f = %lf (%lf), |w| = %lf, |g| = %lf, step = %lf, efficiency = %lf%%", 
                                                     iteration, double(f), double(f - RealT(0.5) * DotProduct(this->C, w*w)),
                                                     double(Norm(w)), double(Norm(g)), double(step_size),
                                                     double(this->optimization_wrapper->GetComputationEngine().GetEfficiency())));
}

template<class RealT>
void InnerOptimizationWrapperStochasticGradient<RealT>::Report(const std::string &s)
{
    this->optimization_wrapper->PrintMessage(SPrintF("Inner message: %s", s.c_str()));
}

//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperStochasticGradient::Minimize()
//
// Perform StochasticGradient optimization.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT InnerOptimizationWrapperStochasticGradient<RealT>::Minimize(std::vector<RealT> &x0)
{
    RealT result = 0;
    std::vector<RealT> g;
    int next_report_iter = 1;

    for (int iter = 1; iter <= MAX_ITERATIONS; iter++) {
        ComputeGradient(g, x0, batch_size);

        RealT stepsize = s0 / pow(1.0 + iter, s1);
        x0 -= stepsize*g;

        // TODO: project gamma parameters to positive reals
        if (iter == next_report_iter || iter == MAX_ITERATIONS) {
            Report(iter, x0, stepsize);
            next_report_iter *= 2;
        }
    }

    return ComputeFunction(x0);
    /*
    RealT f = RealT(1e20);
    for (log_base = 1; log_base < RealT(1e6); log_base *= 2)
    {
        this->optimization_wrapper->PrintMessage(SPrintF("Inner optimization using LOG_BASE = %lf", double(log_base)));
        this->optimization_wrapper->Indent();
        f = EM<RealT>::Minimize(x0);        
        this->optimization_wrapper->Unindent();
    }
    return f;
    */
}
