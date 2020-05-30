//////////////////////////////////////////////////////////////////////
// OptimizationWrapper.ipp
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::OptimizationWrapper()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
OptimizationWrapper<RealT>::OptimizationWrapper(ComputationWrapper<RealT> &computation_wrapper) :
    computation_wrapper(computation_wrapper),
    indent(0)
{
    logfile.open("optimize.log");
    if (logfile.fail()) Error("Could not open log file for writing.");
}

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::~OptimizationWrapper()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
OptimizationWrapper<RealT>::~OptimizationWrapper()
{
    logfile.close();
}

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::Indent()
// OptimizationWrapper<RealT>::Unindent()
// OptimizationWrapper<RealT>::PrintMessage()
//
// Print indented message.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void OptimizationWrapper<RealT>::Indent() { indent++; }

template<class RealT>
void OptimizationWrapper<RealT>::Unindent() { indent--; Assert(indent >= 0, "Cannot unindent!"); }

template<class RealT>
void OptimizationWrapper<RealT>::PrintMessage(const std::string &s)
{
    for (int i = 0; i < indent; i++) std::cerr << "    ";
    for (int i = 0; i < indent; i++) logfile << "    ";
    std::cerr << s << std::endl;
    logfile << s << std::endl;
}

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::Train()
//
// Run optimization algorithm with fixed regularization
// constants.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT OptimizationWrapper<RealT>::Train(const std::vector<int> &units,
                                        std::vector<RealT> &w,
                                        std::vector<RealT> &weights_initial,
                                        const std::vector<RealT> &C)
{
    static std::vector<int> cached_units;
    static std::vector<RealT> cached_initial_w;
    static std::vector<RealT> cached_C;
    static std::vector<RealT> cached_learned_w;
    static RealT cached_f;

    if (cached_units != units ||
        cached_initial_w != w ||
        cached_C != C)
    {
        cached_units = units;
        cached_initial_w = w;
        cached_C = C;
        cached_learned_w = w;
        
        WriteProgressMessage("Starting training...");

#if STOCHASTIC_GRADIENT
        Error("Not yet implemented.");
#else

        const std::vector<RealT> Ce = GetParameterManager().ExpandParameterGroupValues(C);
        const RealT log_base = RealT(GetOptions().GetRealValue("log_base"));
        const RealT hyperparam_data = RealT(GetOptions().GetRealValue("hyperparam_data"));
        const RealT kd_hyperparam_data = RealT(GetOptions().GetRealValue("kd_hyperparam_data"));
        const RealT lig_hyperparam_data = RealT(GetOptions().GetRealValue("lig_hyperparam_data"));



        if (GetOptions().GetBoolValue("viterbi_parsing"))
        {
            std::vector<RealT> bias(w.size());

            
            
#if SMOOTH_MAX_MARGIN
            InnerOptimizationWrapperLBFGS<RealT> inner_optimization_wrapper(this, units, Ce);
#else
#if BMRM_AVAILABLE
            InnerOptimizationWrapperBundleMethod<RealT> inner_optimization_wrapper(this, units, Ce);
#else
            InnerOptimizationWrapperSubgradientMethod<RealT> inner_optimization_wrapper(this, units, Ce);
            PrintMessage("BMRM not available, so defaulting to subgradient algorithm.");
#endif
#endif
                
            for (int i = 0; i < NUM_CCCP_STEPS; i++)
            {
                PrintMessage(SPrintF("Starting inner loop (pass %d)...", i));
                if (i > 0) bias = -RealT(NONCONVEX_MULTIPLIER) * computation_wrapper.ComputeGradient(units, cached_learned_w, true, false, log_base, 
                    hyperparam_data,kd_hyperparam_data,lig_hyperparam_data);

                std::cerr << bias << std::endl;
                inner_optimization_wrapper.LoadBias(bias);
                cached_f = inner_optimization_wrapper.Minimize(cached_learned_w);
                GetParameterManager().WriteToFile(SPrintF("optimize.params.stage%d", i+1), cached_learned_w);
                
                RealT loss = computation_wrapper.ComputeLoss(units, cached_learned_w, log_base);
                PrintMessage(SPrintF("Current loss: %lf", loss));
                if (RealT(NONCONVEX_MULTIPLIER) == RealT(0)) break;
            }

        }
        else
        {
            InnerOptimizationWrapperLBFGS<RealT> inner_optimization_wrapper(this, units, weights_initial, Ce);
            
            cached_f = inner_optimization_wrapper.Minimize(cached_learned_w);
        }

#endif
    }
    else
    {
        PrintMessage ("Using cached result from Train()...");
    }
    
    w = cached_learned_w;
    return cached_f;
}

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::TrainSGD()
//
// Run optimization algorithm with fixed regularization
// constants.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT OptimizationWrapper<RealT>::TrainSGD(const std::vector<int> &units,
                                           std::vector<RealT> &w,
                                           const std::vector<RealT> &C)
{
    static std::vector<int> cached_units;
    static std::vector<RealT> cached_initial_w;
    static std::vector<RealT> cached_C;
    static std::vector<RealT> cached_learned_w;
    static RealT cached_f;

    if (cached_units != units ||
        cached_initial_w != w ||
        cached_C != C)
    {
        cached_units = units;
        cached_initial_w = w;
        cached_C = C;
        cached_learned_w = w;

        WriteProgressMessage("Starting SGD training...");

#if STOCHASTIC_GRADIENT
        Error("Not yet implemented.");
#else

        const std::vector<RealT> Ce = GetParameterManager().ExpandParameterGroupValues(C);
        // const RealT log_base = RealT(GetOptions().GetRealValue("log_base"));

        if (GetOptions().GetBoolValue("viterbi_parsing"))
        {
            Error("Hard-EM (viterbi-parsing) training not implemented.");
            std::vector<RealT> bias(w.size());
        }
        else
        {
            InnerOptimizationWrapperStochasticGradient<RealT> inner_optimization_wrapper(this, units, Ce);

            // Not sure if this is necessary
            /*
            std::vector<std::vector<bool> > config_params;
            for (int i = 0; i < num_data_sources; i++)
            {
                config_params.push_back(std::vector<bool>(M*2, false));
                for (int j = 0; j < M; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        config_params[i][j*2 + k] = inner_optimization_wrapper.FindZerosInData(j,k,i);
                        if (config_params[i][j*2 + k])
                            std::cout << "GammaMLE::Using Method of Moments since there are 0-counts for data source " << i << " CPD(" << j << "," << k << ")" << std::endl;
                   }
               }
            }
            */

            cached_f = inner_optimization_wrapper.Minimize(cached_learned_w);
        }

#endif
    }
    else
    {
        PrintMessage ("Using cached result from TrainSGD()...");
    }

    w = cached_learned_w;
    return cached_f;
}

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::TrainEM()
//
// Run optimization algorithm with fixed regularization
// constants.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT OptimizationWrapper<RealT>::TrainEM(const std::vector<int> &units,
                                        std::vector<RealT> &w,
                                        const std::vector<RealT> &C, const int train_max_iter)
{
    static std::vector<int> cached_units;
    static std::vector<RealT> cached_initial_w;
    static std::vector<RealT> cached_C;
    static std::vector<RealT> cached_learned_w;
    static RealT cached_f;

    if (cached_units != units ||
        cached_initial_w != w ||
        cached_C != C)
    {
        cached_units = units;
        cached_initial_w = w;
        cached_C = C;
        cached_learned_w = w;
        
        WriteProgressMessage("Starting training...");

#if STOCHASTIC_GRADIENT
        Error("Not yet implemented.");
#else

        const std::vector<RealT> Ce = GetParameterManager().ExpandParameterGroupValues(C);
        // const RealT log_base = RealT(GetOptions().GetRealValue("log_base"));
        
        if (GetOptions().GetBoolValue("viterbi_parsing"))
        {
            Error("Hard-EM (viterbi-parsing) training not implemented.");
            std::vector<RealT> bias(w.size());            
        }
        else
        {
            InnerOptimizationWrapperEM<RealT> inner_optimization_wrapper(this, units, Ce);            
            RealT old_f = 1e100;
            const double opt_tol = 1e-3;
            const int MAX_ITER = train_max_iter;
//            const int MAX_ITER = 1000;

            std::cout << "Using MAX_ITER = " << train_max_iter << std::endl;

            int num_data_sources = GetOptions().GetIntValue("num_data_sources");
            std::vector<std::vector<bool> > config_params;
            for (int i = 0; i < num_data_sources; i++)
            {
                config_params.push_back(std::vector<bool>(M*2, false));
                for (int j = 0; j < M; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        config_params[i][j*2 + k] = inner_optimization_wrapper.FindZerosInData(j,k,i);
                        if (config_params[i][j*2 + k])
                            std::cout << "GammaMLE::Using Method of Moments since there are 0-counts for data source " << i << " CPD(" << j << "," << k << ")" << std::endl;
                   }
               }
            }
            // Main EM loop: do full E, one gradient (M) step each iteration.           
            for (int i = 1; i <= MAX_ITER; i++) {
                cached_f = inner_optimization_wrapper.OneStep(cached_learned_w, i, config_params);

                // Terminate when the function (marginal likelihood) doesn't change much
                if (fabs(old_f - cached_f) < opt_tol) break;
                old_f = cached_f;
            }
        }

#endif
    }
    else
    {
        PrintMessage ("Using cached result from TrainEM()...");
    }
    
    w = cached_learned_w;
    return cached_f;
}


//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::LearnHyperparameters()
//
// Use holdout cross validation in order to estimate
// regularization constants.
//////////////////////////////////////////////////////////////////////

#if HYPERPARAMETER_GRID_SEARCH
template<class RealT>
void OptimizationWrapper<RealT>::LearnHyperparameters(std::vector<int> units,
                                                      std::vector<RealT> &w)
{
    // split data into training and holdout sets

    if (GetOptions().GetBoolValue("shuffle")){
    std::random_shuffle(units.begin(), units.end());
    }
    
    const RealT holdout_ratio = GetOptions().GetRealValue("holdout_ratio");
    const std::vector<int> holdout(units.begin(), units.begin() + int(units.size() * holdout_ratio));
    const std::vector<int> training(units.begin() + int(units.size() * holdout_ratio), units.end());

    const RealT log_base = RealT(GetOptions().GetRealValue("log_base"));
    const RealT hyperparam_data = RealT(GetOptions().GetRealValue("hyperparam_data"));
    const RealT kd_hyperparam_data = RealT(GetOptions().GetRealValue("kd_hyperparam_data"));
    const RealT lig_hyperparam_data = RealT(GetOptions().GetRealValue("lig_hyperparam_data"));


    if (training.size() == 0 || holdout.size() == 0) 
        Error("Not enough training samples for cross-validation.");

    // do hyperparameter optimization
    PrintMessage("Starting hyperparameter optimization...");
    Indent();
    
    PrintMessage("List of hyperparameters:");
    Indent();
    const std::vector<ParameterGroup> &groups = GetParameterManager().GetParameterGroups();
    for (size_t i = 0; i < groups.size(); i++)
        PrintMessage(SPrintF("Parameter group %d: %s", i+1, groups[i].name.c_str()));
    Unindent();

    RealT best_C = 0, best_holdout_loss = 1e20;
    std::vector<RealT> C = std::vector<RealT>(GetParameterManager().GetNumParameterGroups());

    std::vector<RealT> w0(w);

    // perform cross-validation
    for (int k = -5; k <= 10; k++)
    {
        // perform training
        std::fill(C.begin(), C.end(), Pow(2.0, RealT(k)));
        PrintMessage(SPrintF("Performing optimization using C = %lf", C[0]));
        Indent();
        std::vector<RealT> x(w);
        const RealT f = Train(training, x, w0, C);
        Unindent();

        // compute holdout loss
#if CROSS_VALIDATE_USING_LOGLOSS
        if (GetOptions().GetBoolValue("viterbi_parsing")) Error("Cannot use logloss for cross validation if Viterbi parsing.");
        RealT loss = computation_wrapper.ComputeFunction(holdout, x, false, false, log_base, hyperparam_data, kd_hyperparam_data, lig_hyperparam_data);
#else
        RealT loss = computation_wrapper.ComputeLoss(holdout, x, true, log_base);
#endif
        
        PrintMessage(SPrintF("Using C = %lf, regularized training loss = %lf, holdout loss = %lf", double(C[0]), double(f), double(loss)));
        
        if (loss < best_holdout_loss)
        {
            best_holdout_loss = loss;
            best_C = C[0];
        }
    }

    Unindent();
    PrintMessage(SPrintF("Chose C = %lf, holdout loss = %lf", best_C, best_holdout_loss));
    std::fill(C.begin(), C.end(), best_C / (1.0 - holdout_ratio));
    
    // now, retrain on all data
    PrintMessage("Retraining on entire training set...");
    Indent();
    Train(units, w, w0, C);
    Unindent();
}
#endif


//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::LearnHyperparametersEM()
//
// Use holdout cross validation in order to estimate
// regularization constants.
//////////////////////////////////////////////////////////////////////

#if HYPERPARAMETER_GRID_SEARCH
template<class RealT>
void OptimizationWrapper<RealT>::LearnHyperparametersEM(std::vector<int> units,
                                                      std::vector<RealT> &w, const int train_max_iter)
{
    // split data into training and holdout sets
    if (GetOptions().GetBoolValue("shuffle")){
    std::random_shuffle(units.begin(), units.end());
    }

    const RealT holdout_ratio = GetOptions().GetRealValue("holdout_ratio");
    const std::vector<int> holdout(units.begin(), units.begin() + int(units.size() * holdout_ratio));
    const std::vector<int> training(units.begin() + int(units.size() * holdout_ratio), units.end());

    const RealT log_base = RealT(GetOptions().GetRealValue("log_base"));
    const RealT hyperparam_data = RealT(GetOptions().GetRealValue("hyperparam_data"));
    const RealT kd_hyperparam_data = RealT(GetOptions().GetRealValue("kd_hyperparam_data"));
    const RealT lig_hyperparam_data = RealT(GetOptions().GetRealValue("lig_hyperparam_data"));


    if (training.size() == 0 || holdout.size() == 0) 
        Error("Not enough training samples for cross-validation.");

    // do hyperparameter optimization
    PrintMessage("Starting hyperparameter optimization...");
    Indent();
    
    PrintMessage("List of hyperparameters:");
    Indent();
    const std::vector<ParameterGroup> &groups = GetParameterManager().GetParameterGroups();
    for (size_t i = 0; i < groups.size(); i++)
        PrintMessage(SPrintF("Parameter group %d: %s", i+1, groups[i].name.c_str()));
    Unindent();

    RealT best_C = 0, best_holdout_loss = 1e20;
    std::vector<RealT> C = std::vector<RealT>(GetParameterManager().GetNumParameterGroups());

    // perform cross-validation
    for (int k = -5; k <= 10; k++)
    {
        // perform training
        std::fill(C.begin(), C.end(), Pow(2.0, RealT(k)));
        PrintMessage(SPrintF("Performing optimization using C = %lf", C[0]));
        Indent();
        std::vector<RealT> x(w);
        const RealT f = TrainEM(training, x, C, train_max_iter);
        Unindent();

        // compute holdout loss
#if CROSS_VALIDATE_USING_LOGLOSS
        if (GetOptions().GetBoolValue("viterbi_parsing")) Error("Cannot use logloss for cross validation if Viterbi parsing.");
        RealT loss = computation_wrapper.ComputeFunction(holdout, x, false, false, log_base, hyperparam_data, kd_hyperparam_data, lig_hyperparam_data);
#else
        RealT loss = computation_wrapper.ComputeLoss(holdout, x, true, log_base);
#endif
        
        PrintMessage(SPrintF("Using C = %lf, regularized training loss = %lf, holdout loss = %lf", double(C[0]), double(f), double(loss)));
        
        if (loss < best_holdout_loss)
        {
            best_holdout_loss = loss;
            best_C = C[0];
        }
    }

    Unindent();
    PrintMessage(SPrintF("Chose C = %lf, holdout loss = %lf", best_C, best_holdout_loss));
    std::fill(C.begin(), C.end(), best_C / (1.0 - holdout_ratio));
    
    // now, retrain on all data
    PrintMessage("Retraining on entire training set...");
    Indent();
    TrainEM(units, w, C, train_max_iter);
    Unindent();
}
#endif


//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::LearnHyperparameters()
//
// Use gradient-based holdout cross-validation in order estimate
// regularization constants.
//////////////////////////////////////////////////////////////////////

#if HYPERPARAMETER_GRADIENT_OPTIMIZATION
template<class RealT>
void OptimizationWrapper<RealT>::LearnHyperparameters(std::vector<int> units,
                                                      std::vector<RealT> &w)
{
    // split data into training and holdout sets
    if (GetOptions().GetBoolValue("shuffle")){
    std::random_shuffle(units.begin(), units.end());
    }
    
    const RealT holdout_ratio = GetOptions().GetRealValue("holdout_ratio");
    const std::vector<int> holdout(units.begin(), units.begin() + int(units.size() * holdout_ratio));
    const std::vector<int> training(units.begin() + int(units.size() * holdout_ratio), units.end());

    if (training.size() == 0 || holdout.size() == 0) 
        Error("Not enough training samples for cross-validation.");

    // do hyperparameter optimization
    PrintMessage("Starting hyperparameter optimization...");
    Indent();
    
    PrintMessage("List of hyperparameters:");
    Indent();
    const std::vector<ParameterGroup> &groups = GetParameterManager().GetParameterGroups();
    for (size_t i = 0; i < groups.size(); i++)
        PrintMessage(SPrintF("Parameter group %d: %s", i+1, groups[i].name.c_str()));
    Unindent();

    std::vector<RealT> log_C;

    if (GetOptions().GetStringValue("opt2_regularization_weights") != ""){

    std::cout << GetOptions().GetStringValue("opt2_regularization_weights") << std::endl;
        // Read in initial regularization weights to use from file

    std::ifstream regfile(GetOptions().GetStringValue("opt2_regularization_weights").c_str());
    if (regfile.fail()) Error("Could not open regfile for reading.");

    std::string line;

    while (std::getline(regfile, line))
    {
        float value;
        std::stringstream ss(line);

        while (ss >> value)
        {
            log_C.push_back(value);
        }
    }

    //TODO: write check that length log_C = NumParameterGroups

    } else { //start OPT2 optimization with parameters from scratch
    std::cout << "Starting OPT2 from scratch" << std::endl;
    log_C = std::vector<RealT>(GetParameterManager().GetNumParameterGroups(), RealT(INITIAL_LOG_C)); //HKWS THIS
    }

    if (GetOptions().GetBoolValue("viterbi_parsing"))
    {
        Error("Not yet implemented.");
    }
    else
    {
        OuterOptimizationWrapper<RealT> outer_optimization_wrapper(this, w, training, holdout);
        outer_optimization_wrapper.Minimize(log_C);

        std::ofstream outfile("optimize.opt2wts");
            if (outfile.fail()) Error("Could not open file for writing opt2wts.");
            for (size_t i = 0; i < log_C.size(); i++)
                outfile << std::setprecision(10) << log_C[i] << std::endl;
            outfile.close();
        GetParameterManager().WriteToFile(SPrintF("optimize.opt2wts.iter0"), log_C);


    }
    
    Unindent();
    std::ostringstream oss;
    const std::vector<RealT> C = Exp(log_C);
    oss << "Chose hyperparameters, C = " << C;
    PrintMessage(oss.str());
    
    std::vector<RealT> w0(w);

    // Now, retrain on all data
    PrintMessage("Retraining on entire training set...");
    Indent();
    Train(units, w, w0, C);
    Unindent();
}
#endif

//////////////////////////////////////////////////////////////////////
// OptimizationWrapper<RealT>::LearnHyperparameters()
//
// Use Bayesian hyperparameter selection algorithm in order to
// estimate regularization constants.
//////////////////////////////////////////////////////////////////////

#if HYPERPARAMETER_MAJORIZATION_MINIMIZATION
template<class RealT>
void OptimizationWrapper<RealT>::LearnHyperparameters(std::vector<int> units,
                                               std::vector<RealT> &w,
                                               RealT holdout_ratio,
                                               bool toggle_viterbi)
{
    // do hyperparameter optimization
    
    PrintMessage("Starting hyperparameter optimization...");
    Indent();
    
    PrintMessage("List of hyperparameters:");
    Indent();
    const std::vector<HyperparameterGroup> &groups = params.GetHyperparameterGroups();
    for (size_t i = 0; i < groups.size(); i++)
        PrintMessage(SPrintF("Hyperparameter group %d: %s", i+1, groups[i].name.c_str()));
    Unindent();
    
    std::vector<RealT> C = std::vector<RealT>(params.GetNumHyperparameterGroups(), 1);
    
    // iterative relinearization

    for (int iters = 0; iters < NUM_ITERATIVE_RELINEARIZATION_STEPS; iters++)
    {
        // show current set of hyperparameters
        
        PrintMessage("Current hyperparameters:");
        Indent();
        const std::vector<HyperparameterGroup> &groups = params.GetHyperparameterGroups();
        for (size_t i = 0; i < groups.size(); i++)
            PrintMessage(SPrintF("Hyperparameter group %d (%s): %lf", i+1, groups[i].name.c_str(), C[i]));
        Unindent();

        // perform training

        std::ostringstream oss;
        const std::vector<RealT> Ce = params.ExpandHyperparameters(C);
        oss << "Performing optimization using C = " << C;
        PrintMessage(oss.str());
        Indent();
        std::vector<RealT> x(w);
        const RealT f = Train(units, x, C, toggle_viterbi);
        Unindent();

        // compute new hyperparameters

        for (size_t g = 0; g < groups.size(); g++)
        {
            RealT numerator = (groups[g].end - groups[g].begin + 1.0) / 2.0;
            RealT denominator = RealT(MM_SMOOTHING);
            for (int i = groups[g].begin; i < groups[g].end; i++)
                denominator += 0.5 * x[i] * x[i];
            C[g] = numerator / denominator;
        }

        // adjust for Viterbi mode

        if (toggle_viterbi)
        {
            const RealT loss = f - 0.5 * DotProduct(Ce, x*x);
            const RealT loss_multiplier = RealT(units.size()) / (RealT(MM_SMOOTHING) + loss);
            C /= loss_multiplier;
        }
    }
    
    Unindent();
    std::ostringstream oss;
    oss << "Chose hyperparameters, C = " << C;
    PrintMessage(oss.str());
    
    // now, retrain on all data
    
    PrintMessage("Retraining on entire training set...");
    Indent();
    Train(units, w, C, toggle_viterbi);
    Unindent();
}
#endif
