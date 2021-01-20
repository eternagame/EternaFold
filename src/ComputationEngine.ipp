//////////////////////////////////////////////////////////////////////
// ComputationEngine.cpp
//////////////////////////////////////////////////////////////////////

#include "ComputationEngine.hpp"

//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputationEngine()
// ComputationEngine::~ComputationEngine()
//
// Constructor and destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
ComputationEngine<RealT>::ComputationEngine(const Options &options,
                                            const std::vector<FileDescription> &descriptions,
                                            InferenceEngine<RealT> &inference_engine,
                                            ParameterManager<RealT> &parameter_manager) :
    DistributedComputation<RealT, SharedInfo<RealT>, NonSharedInfo>(options.GetBoolValue("verbose_output")),
    options(options),
    descriptions(descriptions),
    inference_engine(inference_engine),
    parameter_manager(parameter_manager)
{ }

template<class RealT>
ComputationEngine<RealT>::~ComputationEngine()
{}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::DoComputation()
//
// Decide what type of computation needs to be done and then
// pass the work on to the appropriate routine.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::DoComputation(std::vector<RealT> &result, 
                                             const SharedInfo<RealT> &shared,
                                             const NonSharedInfo &nonshared)
{
    switch (shared.command)
    {
        case CHECK_PARSABILITY:
            CheckParsability(result, nonshared);
            break;
        case COMPUTE_SOLUTION_NORM_BOUND:
            ComputeSolutionNormBound(result, shared, nonshared);
            break;
        case COMPUTE_GRADIENT_NORM_BOUND:
            ComputeGradientNormBound(result, nonshared);
            break;
        case COMPUTE_LOSS:
            ComputeLoss(result, shared, nonshared);
            break;
        case COMPUTE_FUNCTION:
            ComputeFunctionAndGradient(result, shared, nonshared, false);
            break;
        case COMPUTE_GRADIENT:
            ComputeFunctionAndGradient(result, shared, nonshared, true);
            break;
        case COMPUTE_MSTEP_FUNCTION:
            ComputeMStepFunctionAndGradient(result, shared, nonshared, false);
            break;
        case COMPUTE_MSTEP_GRADIENT:
            ComputeMStepFunctionAndGradient(result, shared, nonshared, true);
            break;
        case COMPUTE_GAMMAMLE_FUNCTION:
            ComputeGammaMLEFunctionAndGradient(result, shared, nonshared, false);
            break;
        case COMPUTE_GAMMAMLE_GRADIENT:
            ComputeGammaMLEFunctionAndGradient(result, shared, nonshared, true);
            break;
        case COMPUTE_GAMMAMLE_SCALING_FACTOR:
            ComputeGammaMLEScalingFactor(result, shared, nonshared);
            break;
        case CHECK_ZEROS_IN_DATA:
            CheckZerosInData(result, shared, nonshared);
            break;
        case COMPUTE_FUNCTION_SE:
            ComputeFunctionAndGradientSE(result, shared, nonshared, false);
            break;
        case COMPUTE_GRADIENT_SE:
            ComputeFunctionAndGradientSE(result, shared, nonshared, true);
            break;
        case COMPUTE_HV:
            ComputeHessianVectorProduct(result, shared, nonshared);
            break;
        case PREDICT:
            Predict(result, shared, nonshared);
            break;
        case PREDICT_FOLDCHANGE:
            PredictFoldChange(result, shared, nonshared);
            break;
        case SAMPLE:
            Sample(result, shared, nonshared);
            break;
        case REVI:
            RunREVI(result, shared, nonshared);
            break;
        case TEST_ENERGIES:
            TestEnergies(result,shared,nonshared);
            break;
        default: 
            Assert(false, "Unknown command type.");
            break;
    }
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::CheckParsability()
//
// Check to see if a sequence is parsable or not.  Return a
// vector with a "0" in the appropriate spot indicating that a
// file is not parsable.
//////////////////////////////////////////////////////////////////////

template <class RealT>
void ComputationEngine<RealT>::CheckParsability(std::vector<RealT> &result, 
                                                const NonSharedInfo &nonshared)
{

    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // conditional inference
    inference_engine.LoadValues(std::vector<RealT>(parameter_manager.GetNumLogicalParameters()));
    inference_engine.UseConstraints(sstruct.GetMapping());
    inference_engine.UpdateEvidenceStructures();

    inference_engine.ComputeViterbi();
    RealT conditional_score = inference_engine.GetViterbiScore();

    // check for bad parse
    result.clear();
    result.resize(descriptions.size());
    result[nonshared.index] = (conditional_score < RealT(NEG_INF/2) ? 0 : 1);
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeSolutionNormBound()
//
// Compute the max entropy and loss possible for an example.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeSolutionNormBound(std::vector<RealT> &result, 
                                                        const SharedInfo<RealT> &shared,
                                                        const NonSharedInfo &nonshared)
{

    RealT max_entropy = RealT(0);
    RealT max_loss = RealT(0);

    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // load parameters
    const std::vector<RealT> w(parameter_manager.GetNumLogicalParameters(), RealT(0));
    inference_engine.LoadValues(w);
    inference_engine.UpdateEvidenceStructures();

    // perform computation
#if !SMOOTH_MAX_MARGIN
    if (!options.GetBoolValue("viterbi_parsing"))
#endif
    {
        inference_engine.ComputeInside();
        max_entropy += inference_engine.ComputeLogPartitionCoefficient();
    }
        
#if defined(HAMMING_LOSS)
    inference_engine.UseLoss(sstruct.GetMapping(), RealT(HAMMING_LOSS));
    inference_engine.ComputeViterbi();
    max_loss += inference_engine.GetViterbiScore();
#endif

    result.clear();
    result.resize(descriptions.size());
    result[nonshared.index] = max_entropy / shared.log_base + max_loss;

    result *= RealT(descriptions[nonshared.index].weight);
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeGradientNormBound()
//
// Compute the max L1 norm for the features of an example.
// Return a vector with this value in the appropriate spot for
// this example.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeGradientNormBound(std::vector<RealT> &result,
                                                        const NonSharedInfo &nonshared)
{
    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // load parameters
    const std::vector<RealT> w(parameter_manager.GetNumLogicalParameters(), RealT(1));
    inference_engine.LoadValues(w);
    inference_engine.UpdateEvidenceStructures();

    // perform inference
    inference_engine.ComputeViterbi();
    const RealT max_L1_norm = inference_engine.GetViterbiScore();

    result.clear();
    result.resize(descriptions.size());
    result[nonshared.index] = max_L1_norm;
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeLoss()
//
// Return a vector containing a single entry with the loss value.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeLoss(std::vector<RealT> &result, 
                                           const SharedInfo<RealT> &shared,
                                           const NonSharedInfo &nonshared)
{
    std::cout << "ComputeLoss" << std::endl;
    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);
    inference_engine.UpdateEvidenceStructures();

    // perform inference
    SStruct *solution;
    if (options.GetBoolValue("viterbi_parsing"))
    {
        inference_engine.ComputeViterbi();
        solution = new SStruct(sstruct);
        solution->SetMapping(inference_engine.PredictPairingsViterbi());
    }
    else
    {
        inference_engine.ComputeInside();
        inference_engine.ComputeOutside();
        inference_engine.ComputePosterior();
        solution = new SStruct(sstruct);
        solution->SetMapping(inference_engine.PredictPairingsPosterior(shared.gamma));
    }

    // compute loss
    if (!shared.use_loss) Error("Must be using loss function in order to compute loss.");
#if defined(HAMMING_LOSS)
    inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
#endif
    inference_engine.LoadValues(std::vector<RealT>(w.size()));
    inference_engine.UseConstraints(solution->GetMapping());
    inference_engine.UpdateEvidenceStructures();
    inference_engine.ComputeViterbi();

    delete solution;

    result.clear();
    result.push_back(inference_engine.GetViterbiScore());

    result *= RealT(descriptions[nonshared.index].weight);
    result.back() /= shared.log_base;
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeMStepFunctionAndGradient();
//
// Return a vector containing the gradient and function value for the
// M-step in EM training. Specifically, uses the expected sufficient 
// statistics for the feature counts under the evidence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeMStepFunctionAndGradient(std::vector<RealT> &result, 
                                                          const SharedInfo<RealT> &shared,
                                                          const NonSharedInfo &nonshared,
                                                          bool need_gradient)
{
    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);
    inference_engine.UpdateEvidenceStructures();

#if defined(HAMMING_LOSS)
    Error("HAMMING_LOSS not implemented within EM training");
    if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
#endif
    
    // unconditional inference - unchanged in M-step optimization
    RealT unconditional_score;
    std::vector<RealT> unconditional_counts;

    if (shared.use_nonsmooth)
    {
        Error("Viterbi training not supported within EM training");
        // To support this, need to add a new function ComputeViterbiEM which does MAP with evidence potential
        inference_engine.ComputeViterbi();
        unconditional_score = inference_engine.GetViterbiScore();
        if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
        inference_engine.ComputeInside();
        unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
        if (need_gradient)
        {
            inference_engine.ComputeOutside();
            unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
        }
    }

    // conditional inference - need to use ESS (if you don't know the labels)
    RealT conditional_score;
    std::vector<RealT> conditional_counts;
    
    if (shared.use_nonsmooth)
    {
        Error("Viterbi training not supported within EM training");
        inference_engine.ComputeViterbi();
        conditional_score = inference_engine.GetViterbiScore();
        if (need_gradient) conditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
	// If we don't know the true structure

	if (!sstruct.HasStruct()) {
            // Compute the ESS - the function value should be dot(w, ESS)
            
	    // No constraints used here if we don't know what the actual structure is.
            inference_engine.ComputeInsideESS();
            inference_engine.ComputeOutsideESS();
            conditional_counts = inference_engine.ComputeFeatureCountExpectationsESS();
            conditional_score = DotProduct(w, conditional_counts);

	} else {

            // Otherwise use the true feature counts

	    // Clamp to true
            inference_engine.UseConstraints(sstruct.GetMapping());
            inference_engine.UpdateEvidenceStructures();

            inference_engine.ComputeInside();
            conditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                conditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
	} 
    }

    result.clear();

    // compute subgradient
    if (need_gradient) result = unconditional_counts - conditional_counts;
    
    // compute function value
    Assert(conditional_score <= unconditional_score, "Conditional score cannot exceed unconditional score.");
    result.push_back(unconditional_score - conditional_score);
    
    // check for bad parse
    if (conditional_score < RealT(NEG_INF/2))
    {
        std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
        fill(result.begin(), result.end(), RealT(0));
        return;
    }

    if (NONCONVEX_MULTIPLIER != 0)
    {
        Error("Nonconvex training not supported within EM training");
      
#if STOCHASTIC_GRADIENT  // TODO: need to check whether this is correct
        if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), RealT(0));
        
        // unconditional counts
        inference_engine.UseMapping(std::vector<int>(sstruct.GetLength() + 1, UNKNOWN));  // TODO: call after this: inference_engine.UpdateEvidenceStructures(); 
        if (shared.use_nonsmooth)
        {
            inference_engine.ComputeViterbi();
            unconditional_score = inference_engine.GetViterbiScore();
            if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
        }
        else
        {
            inference_engine.ComputeInside();
            unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
        
        // conditional counts
        inference_engine.UseMapping(sstruct.GetMapping());
        if (shared.use_nonsmooth)
        {
            inference_engine.ComputeViterbi();
            unconditional_score = inference_engine.GetViterbiScore();
            if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
        }
        else
        {
            inference_engine.ComputeInside();
            unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
        
        std::vector<RealT> result2;
        
        // compute subgradient
        if (need_gradient) result2 = unconditional_counts - conditional_counts;
        
        // compute function value
        Assert(conditional_score <= unconditional_score, "Conditional score cannot exceed unconditional score.");
        result2.push_back(unconditional_score - conditional_score);
        
        // check for bad parse
        if (conditional_score < RealT(NEG_INF/2))
        {
            std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
            fill(result.begin(), result.end(), 0);
            return;
        }
        
        result -= NONCONVEX_MULTIPLIER * result2;
#endif
    }

    // avoid precision problems
    if (result.back() < 0)
    {
        if (result.back() < -1e-6)
        {
            std::cerr << "Encountered negative function value for " << descriptions[nonshared.index].input_filename << ": " << result.back() << std::endl;
            parameter_manager.WriteToFile(SPrintF("neg_params.%s", GetBaseName(descriptions[nonshared.index].input_filename).c_str()), w);
            exit(0);
        }
        std::fill(result.begin(), result.end(), RealT(0));
        return;
    }

    result *= RealT(descriptions[nonshared.index].weight);
    result.back() /= shared.log_base;
}


//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeGammaMLEFunctionAndGradient();
//
// Return a vector containing the gradient and function value for the
// M-step in EM training. Specifically, uses the expected sufficient 
// statistics for the feature counts under the evidence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeGammaMLEFunctionAndGradient(std::vector<RealT> &result, 
                                                          const SharedInfo<RealT> &shared,
                                                          const NonSharedInfo &nonshared,
                                                          bool need_gradient)
{
    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;

    int which_data = shared.which_data;

    // ignore structures that have no evidence for this dataset
    if (!sstruct.HasEvidence(which_data))
    {
        result.clear();
        if (need_gradient) {
            result.push_back(RealT(0));  // sum d
            result.push_back(RealT(0));  // sumlog d (for MLE) or sum d^2 (for MM)
            result.push_back(RealT(0));  // num examples
        }
        result.push_back(RealT(0));  // ll
        return;
    }

    inference_engine.LoadSequence(sstruct);
    
    // set the true structure if it exists
    inference_engine.UseConstraints(sstruct.GetMapping());

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);
#if defined(HAMMING_LOSS)
    Error("HAMMING_LOSS not implemented within EM training");
//    if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
#endif

    inference_engine.UpdateEvidenceStructures();


    RealT update_gammamle_sssum;
    RealT update_gammamle_sssumlog;
    RealT update_gammamle_num_examples;
    RealT update_gammamle_ssq;
    RealT ll = 0;
    std::vector<RealT> stats(3);

    int j = shared.id_base; 
    int k = shared.id_pairing; 
    bool areZeros = (bool)shared.areZeros; 
    RealT scale = (RealT)shared.evidence_data_scale;
    std::vector<int> evidence_cpd_id;
    evidence_cpd_id.push_back(j);
    evidence_cpd_id.push_back(k);
    evidence_cpd_id.push_back(areZeros);
    
    // Note: we adjust the suff stats using the scale parameter as follows:
    // sssum => ssum / scale
    // sssumlog => ssumlog - N * log scale
    // esssumlog => essumlog - (log scale) * ssq

    RealT current_k = exp(w[parameter_manager.GetLogicalIndex(inference_engine.GetLogScoreEvidence(0,j,k,which_data))]);
    RealT current_theta = exp(w[parameter_manager.GetLogicalIndex(inference_engine.GetLogScoreEvidence(1,j,k,which_data))]);

    if (shared.use_nonsmooth)
    {
        Error("Viterbi training not supported within EM training");
        //        inference_engine.ComputeViterbi();
//        conditional_score = inference_engine.GetViterbiScore();
//        if (need_gradient) conditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
	if (!sstruct.HasStruct()) {

            // Compute expected suff stats and use these in the update equations
            inference_engine.ComputeInsideESS();
            inference_engine.ComputeOutsideESS();
            inference_engine.ComputePosteriorESS();

  	    // If we don't know the true structure
	    stats = inference_engine.ComputeGammaMLEESS(evidence_cpd_id, !areZeros, !areZeros, which_data);  
            update_gammamle_num_examples = inference_engine.GetNumExamplesSeq(evidence_cpd_id,false, which_data);  

            // adjust suff stats as noted above
            update_gammamle_sssum = stats[0];
            update_gammamle_ssq = stats[2];
            update_gammamle_sssumlog = stats[1];

            if (!areZeros)
                ll = (current_k-1)*update_gammamle_sssumlog - update_gammamle_sssum / current_theta - update_gammamle_num_examples * current_k * log(current_theta) - update_gammamle_num_examples * lgamma(current_k);
            else {
                // in the LL calculation, ignore the 0-counts and use MLE SS
                stats = inference_engine.ComputeGammaMLEESS(evidence_cpd_id, true, true, which_data);  // either: true, true (so don't ignore zeros and use MLE SS) OR false, false (so ignore zeros and use method of
		RealT update_gammamle_sssumlog_nozeros = stats[1];
                RealT update_gammamle_num_examples_nozeros = (int)inference_engine.GetNumExamplesSeq(evidence_cpd_id,true, which_data);  // ignore the zeros
                ll = (current_k-1)*update_gammamle_sssumlog_nozeros - update_gammamle_sssum / current_theta - update_gammamle_num_examples_nozeros * current_k * log(current_theta) - update_gammamle_num_examples_nozeros * lgamma(current_k);
            }

            // adjust suff stats as noted above
            update_gammamle_sssum = update_gammamle_sssum / scale;
            update_gammamle_sssumlog = update_gammamle_sssumlog - update_gammamle_ssq * log(scale);

	} else {

  	    // If we know the true structure
	    stats = inference_engine.ComputeGammaMLESS(evidence_cpd_id, !areZeros, !areZeros, which_data);  
	    update_gammamle_num_examples = inference_engine.GetNumExamplesSeqPairing(evidence_cpd_id,false, which_data);  

            // adjust suff stats as noted above
            update_gammamle_sssum = stats[0];
            update_gammamle_sssumlog = stats[1];

            if (!areZeros)
                ll = (current_k-1)*update_gammamle_sssumlog - update_gammamle_sssum / current_theta - update_gammamle_num_examples * current_k * log(current_theta) - update_gammamle_num_examples * lgamma(current_k);
            else {
                // in the LL calculation, ignore the 0-counts and use MLE SS
                stats = inference_engine.ComputeGammaMLESS(evidence_cpd_id, true, true, which_data);  
		RealT update_gammamle_sssumlog_nozeros = stats[1];
                RealT update_gammamle_num_examples_nozeros = (int)inference_engine.GetNumExamplesSeqPairing(evidence_cpd_id,true, which_data);  // ignore the zeros
                ll = (current_k-1)*update_gammamle_sssumlog_nozeros - update_gammamle_sssum / current_theta - update_gammamle_num_examples_nozeros * current_k * log(current_theta) - update_gammamle_num_examples_nozeros * lgamma(current_k);
            }

            // adjust suff stats as noted above
            update_gammamle_sssum = update_gammamle_sssum / scale;
            update_gammamle_sssumlog = update_gammamle_sssumlog - update_gammamle_num_examples * log(scale);

	}
    }

    result.clear();

    if (need_gradient) {
         result.push_back(update_gammamle_sssum);
         result.push_back(update_gammamle_sssumlog);
         result.push_back(update_gammamle_num_examples);
    }

    // compute function value
    result.push_back(ll);
    
/*
    // check for bad parse
    if (conditional_score < RealT(NEG_INF/2))
    {
        std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
        fill(result.begin(), result.end(), RealT(0));
        return;
    }

    if (NONCONVEX_MULTIPLIER != 0)
    {
        Error("Nonconvex training not supported within EM training");
#if STOCHASTIC_GRADIENT
#endif
    }

    // avoid precision problems
    if (result.back() < 0)
    {
        if (result.back() < -1e-6)
        {
            std::cerr << "Encountered negative function value for " << descriptions[nonshared.index].input_filename << ": " << result.back() << std::endl;
            parameter_manager.WriteToFile(SPrintF("neg_params.%s", GetBaseName(descriptions[nonshared.index].input_filename).c_str()), w);
            exit(0);
        }
        std::fill(result.begin(), result.end(), RealT(0));
        return;
    }
*/

    result *= RealT(descriptions[nonshared.index].weight);

}


//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeGammaMLEScalingFactor();
//
// Return a vector containing the gradient and function value for the
// M-step in EM training. Specifically, uses the expected sufficient 
// statistics for the feature counts under the evidence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeGammaMLEScalingFactor(std::vector<RealT> &result, 
                                                          const SharedInfo<RealT> &shared,
                                                          const NonSharedInfo &nonshared)
{

    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    int which_data = shared.which_data;

    // ignore structures that have no evidence for this dataset
    if (!sstruct.HasEvidence(which_data))
    {
        result.clear();
        result.push_back(0);  // sum
        result.push_back(0);  // num_examples
        return;
    }

    inference_engine.LoadSequence(sstruct);
    
    // set the true structure if it exists
    inference_engine.UseConstraints(sstruct.GetMapping());

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    inference_engine.UpdateEvidenceStructures();

    RealT update_gammamle_sssum;
    RealT update_gammamle_num_examples;

    int j = shared.id_base;
    int k = shared.id_pairing;
    std::vector<int> evidence_cpd_id;
    evidence_cpd_id.push_back(j);
    evidence_cpd_id.push_back(k);

    if (!sstruct.HasStruct()) {

        // If we don't know the true structure
        inference_engine.ComputeInsideESS();
        inference_engine.ComputeOutsideESS();
        inference_engine.ComputePosteriorESS();

	update_gammamle_sssum = inference_engine.ComputeGammaMLESum(evidence_cpd_id, true, true, which_data);  // ignore pairing since don't know it, use posterior
        update_gammamle_num_examples = inference_engine.GetNumExamplesSeq(evidence_cpd_id,false, which_data);  // don't ignore zeros
    } else {

        // If we don't know the true structure
	update_gammamle_sssum = inference_engine.ComputeGammaMLESum(evidence_cpd_id, false, false, which_data);  // don't ignore pairing since we know it, don't use posterior
        update_gammamle_num_examples = inference_engine.GetNumExamplesSeqPairing(evidence_cpd_id,false, which_data);  // don't ignore zeros
    }

    result.clear();
    result.push_back(update_gammamle_sssum);
    result.push_back(update_gammamle_num_examples);
}


template<class RealT>
void ComputationEngine<RealT>::CheckZerosInData(std::vector<RealT> &result, 
                                                          const SharedInfo<RealT> &shared,
                                                          const NonSharedInfo &nonshared)
{

    int which_data = shared.which_data;

    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // conditional inference
    inference_engine.LoadValues(std::vector<RealT>(parameter_manager.GetNumLogicalParameters()));
    inference_engine.UseConstraints(sstruct.GetMapping());

    result.clear();
    result.resize(descriptions.size());
    // ignore structures that only have ground truth

    // ignore structures that have no evidence for this dataset
    if (!sstruct.HasEvidence(which_data))
    {
        result[nonshared.index] = 0;
        return;
    }

    inference_engine.UpdateEvidenceStructures(which_data);

    int areZeros;
    if (!sstruct.HasStruct()) {
        areZeros = inference_engine.AreZerosInSeq(shared.id_base, shared.which_data);
    } else {
        areZeros = inference_engine.AreZerosInSeqPairing(shared.id_base, shared.id_pairing, shared.which_data);
    }

    result[nonshared.index] = (RealT)areZeros;
}


//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeFunctionAndGradient();
//
// Return a vector containing the gradient and function value.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeFunctionAndGradient(std::vector<RealT> &result, 
                                                          const SharedInfo<RealT> &shared,
                                                          const NonSharedInfo &nonshared,
                                                          bool need_gradient)
{

    // load training example
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);
#if defined(HAMMING_LOSS)
    if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
#endif
    
    // unconditional inference
    RealT unconditional_score;
    std::vector<RealT> unconditional_counts;

    if (shared.use_nonsmooth)
    {
        inference_engine.ComputeViterbi();
        unconditional_score = inference_engine.GetViterbiScore();
        if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
        inference_engine.ComputeInside();
        unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
        if (need_gradient)
        {
            inference_engine.ComputeOutside();
            unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
        }
    }

    // conditional inference
    RealT conditional_score;
    std::vector<RealT> conditional_counts;

    inference_engine.UseConstraints(sstruct.GetMapping());
    if (shared.use_nonsmooth)
    {
        inference_engine.ComputeViterbi();
        conditional_score = inference_engine.GetViterbiScore();
        if (need_gradient) conditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
        inference_engine.ComputeInside();
        conditional_score = inference_engine.ComputeLogPartitionCoefficient();
        if (need_gradient)
        {
            inference_engine.ComputeOutside();
            conditional_counts = inference_engine.ComputeFeatureCountExpectations();
        }
    }

    result.clear();

    // compute subgradient
    if (need_gradient) result = unconditional_counts - conditional_counts;
    
    // compute function value
    Assert(conditional_score <= unconditional_score, "Conditional score cannot exceed unconditional score.");
    result.push_back(unconditional_score - conditional_score);
    
    // check for bad parse
    if (conditional_score < RealT(NEG_INF/2))
    {
        std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
        fill(result.begin(), result.end(), RealT(0));
        return;
    }

    if (NONCONVEX_MULTIPLIER != 0)
    {
        
#if STOCHASTIC_GRADIENT
        if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), RealT(0));
        
        // unconditional counts
        inference_engine.UseMapping(std::vector<int>(sstruct.GetLength() + 1, UNKNOWN));
        if (shared.use_nonsmooth)
        {
            inference_engine.ComputeViterbi();
            unconditional_score = inference_engine.GetViterbiScore();
            if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
        }
        else
        {
            inference_engine.ComputeInside();
            unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
        
        // conditional counts
        inference_engine.UseMapping(sstruct.GetMapping());
        if (shared.use_nonsmooth)
        {
            inference_engine.ComputeViterbi();
            unconditional_score = inference_engine.GetViterbiScore();
            if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
        }
        else
        {
            inference_engine.ComputeInside();
            unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
        
        std::vector<RealT> result2;
        
        // compute subgradient
        if (need_gradient) result2 = unconditional_counts - conditional_counts;
        
        // compute function value
        Assert(conditional_score <= unconditional_score, "Conditional score cannot exceed unconditional score.");
        result2.push_back(unconditional_score - conditional_score);
        
        // check for bad parse
        if (conditional_score < RealT(NEG_INF/2))
        {
            std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
            fill(result.begin(), result.end(), 0);
            return;
        }
        
        result -= NONCONVEX_MULTIPLIER * result2;
#endif
    }

    // avoid precision problems
    if (result.back() < 0)
    {
        if (result.back() < -1e-6)
        {
            std::cerr << "Encountered negative function value for " << descriptions[nonshared.index].input_filename << ": " << result.back() << std::endl;
            parameter_manager.WriteToFile(SPrintF("neg_params.%s", GetBaseName(descriptions[nonshared.index].input_filename).c_str()), w);
            exit(0);
        }
        std::fill(result.begin(), result.end(), RealT(0));
        return;
    }

    result *= RealT(descriptions[nonshared.index].weight);
    result.back() /= shared.log_base;
}


//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeFunctionAndGradientSE();
//
// Return a vector containing the gradient and function value for the
// M-step in EM training. Specifically, uses the expected sufficient 
// statistics for the feature counts under the evidence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeFunctionAndGradientSE(std::vector<RealT> &result, 
                                                            const SharedInfo<RealT> &shared,
                                                            const NonSharedInfo &nonshared,
                                                            bool need_gradient)
{

    RealT unconditional_score_ms2_ref;
    std::vector<RealT> unconditional_counts_ms2_ref;
    RealT conditional_score_ms2_ref;
    std::vector<RealT> conditional_counts_ms2_ref;

    // load training example to check
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);
    inference_engine.UpdateEvidenceStructures();

    if (sstruct.HasKD()){
        // Calculate values for MS2 hairpin reference

        std::string MS2_seq = "ACAUGAGGAUCACCCAUGU";
        std::string MS2_cst = "(((((.((....)))))))"; //contrafold seq parsing uses dots for unpaired ('?' is unspecified)
        SStruct sstruct_0;
        sstruct_0.LoadAPI(MS2_seq, MS2_cst);

        inference_engine.LoadSequence(sstruct_0); //changed from sstruct_0

        inference_engine.ComputeInside();
        unconditional_score_ms2_ref = inference_engine.ComputeLogPartitionCoefficient();

        if (need_gradient)
        {
        inference_engine.ComputeOutside();
        unconditional_counts_ms2_ref = inference_engine.ComputeFeatureCountExpectations();
        }

        inference_engine.UseConstraints(sstruct_0.GetMapping());
        inference_engine.UpdateEvidenceStructures();

        inference_engine.ComputeInside();
        conditional_score_ms2_ref = inference_engine.ComputeLogPartitionCoefficient();

        if (need_gradient)
        {
        inference_engine.ComputeOutside();
        conditional_counts_ms2_ref = inference_engine.ComputeFeatureCountExpectations();
        }

    }

    // load training example
    inference_engine.LoadSequence(sstruct);

#if defined(HAMMING_LOSS)
    Error("HAMMING_LOSS not implemented within EM training");
    if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), shared.log_base * RealT(HAMMING_LOSS));
#endif
    
    // unconditional inference - only required if we want the gradient for evidence-only data
    RealT unconditional_score = 0;
    std::vector<RealT> unconditional_counts;

    if (shared.use_nonsmooth)
    {
        Error("Viterbi training not supported within EM training");
        // To support this, need to add a new function ComputeViterbiEM which does MAP with evidence potential
        inference_engine.ComputeViterbi();
        unconditional_score = inference_engine.GetViterbiScore();
        if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
        inference_engine.ComputeInside();
        // We still need to compute Z for the Contrafold distribution even when
        // we have evidence because we need it to normalize the Contrafold
        // potentials to get log P(d | x) as the log-partition function of Q(y).
        unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
        //std::cout << "line 966 uncond score " << unconditional_score << std::endl;
        if (need_gradient)
        {
            inference_engine.ComputeOutside();
            unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
        }
    }

    // conditional inference - need to use ESS (if you don't know the labels)
    RealT conditional_score;
    std::vector<RealT> conditional_counts;
    RealT conditional_score2;
    std::vector<RealT> conditional_counts2;
    RealT conditional_score3;
    std::vector<RealT> conditional_counts3;

    if (shared.use_nonsmooth)
    {
        Error("Viterbi training not supported within EM training");
        inference_engine.ComputeViterbi();
        conditional_score = inference_engine.GetViterbiScore();
        if (need_gradient) conditional_counts = inference_engine.ComputeViterbiFeatureCounts();
    }
    else
    {
        if (sstruct.HasKD()){

            // MS2 aptamer
            inference_engine.UseConstraints(sstruct.GetMapping());
            inference_engine.UpdateEvidenceStructures();

            inference_engine.ComputeInside();
            conditional_score = inference_engine.ComputeLogPartitionCoefficient();

            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                conditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
                if (options.GetBoolValue("train_with_ligand_data"))
                {
                // lig aptamer
                // 
                inference_engine.UseConstraints(sstruct.GetMapping2());

                inference_engine.ComputeInside();
                conditional_score2 = inference_engine.ComputeLogPartitionCoefficient();

                if (need_gradient)
                {
                    inference_engine.ComputeOutside();
                    conditional_counts2 = inference_engine.ComputeFeatureCountExpectations();
                }
                // MS2 & lig aptamer
                //
                inference_engine.UseConstraints(sstruct.GetMapping3());

                inference_engine.ComputeInside();
                conditional_score3 = inference_engine.ComputeLogPartitionCoefficient();
                if (need_gradient)
                {
                    inference_engine.ComputeOutside();
                    conditional_counts3 = inference_engine.ComputeFeatureCountExpectations();
                }
            }

        // If we don't know the true structure
        } else {
        if (!sstruct.HasStruct()) {
            // Compute the ESS - the function value is log Q = log sum_y P(y,d|x)
            // No constraints used here if we don't know what the actual structure is.
            inference_engine.ComputeInsideESS();
            conditional_score = inference_engine.ComputeLogPartitionCoefficientESS();

            if (need_gradient) {
                inference_engine.ComputeOutsideESS();
                inference_engine.ComputePosteriorESS(); // needed for Gamma ESS
                conditional_counts = inference_engine.ComputeFeatureCountExpectationsESS();
            }
        } else {
            // Otherwise use the true feature counts
            // Clamp to true
            inference_engine.UseConstraints(sstruct.GetMapping());
            inference_engine.UpdateEvidenceStructures();

            inference_engine.ComputeInside();
            conditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                conditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
    }
    }

    result.clear();

    // compute gradient
    if (need_gradient){
        if (sstruct.HasKD()){ // HKWS ed.
        const std::vector<double> &kd_data = sstruct.GetKdData(); // list of log_kd_no_lig, log_kd_w_lig, ligand_bonus
        RealT log_kd_no_lig = kd_data[0] - 0.72; // 0.72 taken as min kd value in ribologic dataset
        RealT log_kd_w_lig = kd_data[1] - 0.72; // 0.72 taken as min kd value in ribologic dataset
        RealT ligand_bonus = kd_data[2];
        //RealT ligand_bonus = options.GetRealValue("ligand_bonus");
        RealT HKWS_gradient_mode = options.GetRealValue("HKWS_gradient_mode");
        RealT pred_kd_no_lig = unconditional_score - conditional_score - unconditional_score_ms2_ref + conditional_score_ms2_ref; //jan24 update

        result = unconditional_counts - conditional_counts - unconditional_counts_ms2_ref + conditional_counts_ms2_ref; //jan24

        for(size_t i = 0; i < result.size(); i++){
            if (HKWS_gradient_mode == 0){
                result[i] *= pred_kd_no_lig - log_kd_no_lig;
            }
            else if (HKWS_gradient_mode == 1){
                result[i] *= log_kd_no_lig - pred_kd_no_lig; //no good from scratch
            }
            else if (HKWS_gradient_mode == 2){
        result[i] *= pow(pow((pred_kd_no_lig - log_kd_no_lig),2),0.5); //equals mode=0 from scratch on test, but no good
            }
        }

        //junk for debugging
        // std::cout << "unconditional_score" << unconditional_score << std::endl;
        // std::cout << "conditional_score  " << conditional_score   << std::endl;
        // std::cout << "pred_kd_no_lig" << pred_kd_no_lig << std::endl;
        // std::cout << "log_kd_no_lig" << log_kd_no_lig << std::endl;
        
        if (options.GetBoolValue("train_with_ligand_data"))
        {
        RealT pred_kd_w_lig = log(exp(unconditional_score)+ligand_bonus*exp(conditional_score2))
                              -log(exp(conditional_score)+ligand_bonus*exp(conditional_score3))
                               - unconditional_score_ms2_ref + conditional_score_ms2_ref;

        RealT prefactor_1 = exp(unconditional_score)/(exp(unconditional_score)+ligand_bonus*exp(conditional_score2));

        RealT prefactor_2 = (ligand_bonus*exp(conditional_score2))/(exp(unconditional_score)+ligand_bonus*exp(conditional_score2));

        RealT prefactor_3 = exp(conditional_score)/(exp(conditional_score)+ligand_bonus*exp(conditional_score3));

        RealT prefactor_4 = (ligand_bonus*exp(conditional_score3))/(exp(conditional_score)+ligand_bonus*exp(conditional_score3));

        // std::cout << "unconditional_score" << unconditional_score << std::endl;
        // std::cout << "conditional_score  " << conditional_score   << std::endl;
        // std::cout << "conditional_score2 " << conditional_score2  << std::endl;
        // std::cout << "conditional_score3 " << conditional_score3  << std::endl;
        // std::cout << "prefactor_1        " << prefactor_1         << std::endl;
        // std::cout << "prefactor_2        " << prefactor_2         << std::endl;
        // std::cout << "prefactor_3        " << prefactor_3         << std::endl;
        // std::cout << "prefactor_4        " << prefactor_4         << std::endl;

        // std::cout << "pred_kd_no_lig" << pred_kd_no_lig << std::endl;
        // std::cout << "log_kd_no_lig" << log_kd_no_lig << std::endl;
        // std::cout << "pred_kd_w_lig" << pred_kd_w_lig << std::endl;
        // std::cout << "log_kd_w_lig" << log_kd_w_lig << std::endl;

              for(size_t i = 0; i < result.size(); i++){

                // implemented with lig/kd multiplier because entire `result` vector will get multiplied at end by shared.kd_hyperparam_data

                result[i] += RealT(shared.lig_hyperparam_data)/RealT(shared.kd_hyperparam_data)*(pred_kd_w_lig - log_kd_w_lig)*
                  (prefactor_1 * unconditional_counts[i]
                 + prefactor_2 * conditional_counts2[i]
                 - prefactor_3 * conditional_counts[i] 
                 - prefactor_4 * conditional_counts3[i]
                 - unconditional_counts_ms2_ref[i] + conditional_counts_ms2_ref[i]);
      }
        }

    }
        else //structural data
{
    result = unconditional_counts - conditional_counts;  
     //std::cout << "struct" << unconditional_score << " " << conditional_score << std::endl;
}      
    } 

    RealT function_value = 0;

    // Compute Gamma CPD (Expected) Sufficient Stats
    int num_data_sources = options.GetIntValue("num_data_sources");

    for (int dataset_id = 0; dataset_id < num_data_sources; dataset_id++) {
        for (int i = 0; i < M; i++) { // nucleotide
            for (int j = 0; j < 2; j++) { // base pairing
                RealT sum_d = 0;
                RealT sum_log_d = 0;
                RealT N = 0;

                std::vector<int> evidence_cpd_id;
                evidence_cpd_id.push_back(i);
                evidence_cpd_id.push_back(j);

                if (sstruct.HasEvidence(dataset_id)) {
                    std::vector<RealT> stats(3, 0);

                    if (sstruct.HasStruct()) {
                        stats = inference_engine.ComputeGammaMLESS(evidence_cpd_id, true, true, dataset_id);
                    } else if (need_gradient) {
                        stats = inference_engine.ComputeGammaMLEESS(evidence_cpd_id, true, true, dataset_id);
                    }
                    sum_d = stats[0];
                    sum_log_d = stats[1];
                    N = stats[2];
                }

                int index_k = parameter_manager.GetLogicalIndex(
                    inference_engine.GetLogScoreEvidence(0, i, j, dataset_id)); // i: seq, j: pr or unpaired, ds=0
                int index_theta = parameter_manager.GetLogicalIndex(
                    inference_engine.GetLogScoreEvidence(1, i, j, dataset_id));

                RealT k = exp(w[index_k]);
                RealT log_theta = w[index_theta];
                RealT theta_inv = exp(-log_theta);

                // Add to the function_value (only if we have struct + evidence)
                if (sstruct.HasStruct()) {
                    function_value -= (k-1)*sum_log_d - sum_d*theta_inv - N*k*log_theta - N*lgamma(k);
                }

                // Compute gradient (use implicit differentiation for actual parameters log(k) and log(theta))
                if (need_gradient) {
                    RealT grad_k = -(sum_log_d - N*log_theta - N*Psi(k)) * k;
                    RealT grad_theta = -(sum_d*theta_inv - N*k);

                    result[index_k] = grad_k;
                    result[index_theta] = grad_theta;
                }

            }
        }

    }

   // compute function value
    RealT function_value_nonevidence = 0;
    if (sstruct.HasKD()){
        const std::vector<double> &kd_data = sstruct.GetKdData();
        RealT log_kd_no_lig = kd_data[0] - 0.72; // 0.72 taken as min kd value in ribologic dataset
        RealT log_kd_w_lig = kd_data[1] - 0.72; // 0.72 taken as min kd value in ribologic dataset 
        RealT ligand_bonus = options.GetRealValue("ligand_bonus");
        RealT pred_kd_no_lig = unconditional_score - conditional_score;

        function_value_nonevidence = pow((log_kd_no_lig - pred_kd_no_lig),2);

        if (options.GetBoolValue("train_with_ligand_data")) {
        RealT pred_kd_w_lig = log(exp(unconditional_score)+ligand_bonus*exp(conditional_score2)) - log(exp(conditional_score)+ligand_bonus*exp(conditional_score3));

        // implemented with lig/kd multiplier because entire `result` vector will get multiplied at end by shared.kd_hyperparam_data  -hkws
        function_value_nonevidence += RealT(shared.lig_hyperparam_data)/RealT(shared.kd_hyperparam_data)*pow((log_kd_w_lig - pred_kd_w_lig),2);
        }

    function_value += function_value_nonevidence;

    }
    else{
        function_value_nonevidence = unconditional_score - conditional_score;
    function_value += function_value_nonevidence;
    }

    if (sstruct.HasStruct()) {
        Assert(conditional_score <= unconditional_score, "Conditional score cannot exceed unconditional score.");
    }

    result.push_back(function_value);

    // check for bad parse
    if (conditional_score < RealT(NEG_INF/2))
    {
        std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
        fill(result.begin(), result.end(), RealT(0));
        return;
    }

    if (NONCONVEX_MULTIPLIER != 0)
    {
        Error("Nonconvex training not supported within EM training");
      
#if STOCHASTIC_GRADIENT  // TODO: need to check whether this is correct
        if (shared.use_loss) inference_engine.UseLoss(sstruct.GetMapping(), RealT(0));
        
        // unconditional counts
        inference_engine.UseMapping(std::vector<int>(sstruct.GetLength() + 1, UNKNOWN));  
        if (shared.use_nonsmooth)
        {
            inference_engine.ComputeViterbi();
            unconditional_score = inference_engine.GetViterbiScore();
            if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
        }
        else
        {
            inference_engine.ComputeInside();
            unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
        
        // conditional counts
        inference_engine.UseMapping(sstruct.GetMapping());
        if (shared.use_nonsmooth)
        {
            inference_engine.ComputeViterbi();
            unconditional_score = inference_engine.GetViterbiScore();
            if (need_gradient) unconditional_counts = inference_engine.ComputeViterbiFeatureCounts();
        }
        else
        {
            inference_engine.ComputeInside();
            unconditional_score = inference_engine.ComputeLogPartitionCoefficient();
            if (need_gradient)
            {
                inference_engine.ComputeOutside();
                unconditional_counts = inference_engine.ComputeFeatureCountExpectations();
            }
        }
        
        std::vector<RealT> result2;
        
        // compute subgradient
        if (need_gradient) result2 = unconditional_counts - conditional_counts;
        
        // compute function value
        Assert(conditional_score <= unconditional_score, "Conditional score cannot exceed unconditional score.");
        result2.push_back(unconditional_score - conditional_score);
        
        // check for bad parse
        if (conditional_score < RealT(NEG_INF/2))
        {
            std::cerr << "Unexpected bad parse for file: " << descriptions[nonshared.index].input_filename << std::endl;
            fill(result.begin(), result.end(), 0);
            return;
        }
        
        result -= NONCONVEX_MULTIPLIER * result2;
#endif
    }

    // avoid precision problems
    if (function_value_nonevidence < 0)
    {
        if (function_value_nonevidence < -1e-6 && sstruct.HasStruct())
        {
            std::cerr << "Encountered negative function value for " << descriptions[nonshared.index].input_filename << ": " << result.back() << std::endl;
            parameter_manager.WriteToFile(SPrintF("neg_params.%s", GetBaseName(descriptions[nonshared.index].input_filename).c_str()), w);
            exit(0);
        }
        std::fill(result.begin(), result.end(), RealT(0));
        return;
    }

    result *= RealT(descriptions[nonshared.index].weight);
    result.back() /= shared.log_base;

    // if this is a data-only example, multiply by the hyperparameter for data
    if (!sstruct.HasStruct()) {
    	result *= RealT(shared.hyperparam_data);
    }
 if (sstruct.HasKD()){ //HKWS ed.
     result *= RealT(shared.kd_hyperparam_data);
 }
}


//////////////////////////////////////////////////////////////////////
// ComputationEngine::ComputeHessianVectorProduct()
//
// Return a vector containing Hv.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::ComputeHessianVectorProduct(std::vector<RealT> &result, 
                                                           const SharedInfo<RealT> &shared,
                                                           const NonSharedInfo &nonshared)
{
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    const std::vector<RealT> v(shared.v, shared.v + parameter_manager.GetNumLogicalParameters());

    if (options.GetBoolValue("viterbi_parsing"))
    {
        Error("Should not use Hessian-vector products with Viterbi parsing.");
    }
    
    const RealT EPSILON = RealT(1e-8);
    SharedInfo<RealT> shared_temp(shared);
    std::vector<RealT> result2;

    for (size_t i = 0; i < parameter_manager.GetNumLogicalParameters(); i++)
        shared_temp.w[i] = shared.w[i] + EPSILON * v[i];
    ComputeFunctionAndGradient(result, shared_temp, nonshared, true);
    
    for (size_t i = 0; i < parameter_manager.GetNumLogicalParameters(); i++)
        shared_temp.w[i] = shared.w[i] - EPSILON * v[i];
    ComputeFunctionAndGradient(result2, shared_temp, nonshared, true);
    
    result = (result - result2) / (RealT(2) * EPSILON);
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::TestEnergies()
//
// Print energies for a single sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::TestEnergies(std::vector<RealT> &result, 
                                       const SharedInfo<RealT> &shared,
                                       const NonSharedInfo &nonshared)
{
    result.clear();

    // load sequence, with constraints if necessary
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);
    if (options.GetBoolValue("use_constraints")) inference_engine.UseConstraints(sstruct.GetMapping());

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    inference_engine.UpdateEvidenceStructures();

    // perform inference
    SStruct *solution;

    inference_engine.ComputeViterbi();

    if (options.GetBoolValue("use_evidence")){
        inference_engine.GetViterbiFeaturesESS();
    }
    else{
            inference_engine.GetViterbiFeatures();
    }

    std::cout << "Viterbi score for \"" << descriptions[nonshared.index].input_filename << "\": " 
            << inference_engine.GetViterbiScore() << std::endl;

    solution = new SStruct(sstruct);
    solution->SetMapping(inference_engine.PredictPairingsViterbi());

    WriteProgressMessage("");
    solution->WriteParens(std::cout);

    delete solution;
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::REVI()
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::RunREVI(std::vector<RealT> &result, 
                                       const SharedInfo<RealT> &shared,
                                       const NonSharedInfo &nonshared)
{
    result.clear();


    // perform inference (repurposing chem mapping potentials as pairwise potentials)

    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);
    inference_engine.InitializeREVIVec();

    if (options.GetBoolValue("use_constraints")) inference_engine.UseConstraints(sstruct.GetMapping());
    // load parameters

    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    int L = sstruct.GetLength();
    int n_samples = 1000;
    int n_iters = 100;
    double last_loss_0 = 100000000;
    double last_loss_1 = 100000000;
    double last_loss_2 = 100000000;
    double step_size = 0.1;

    RealT sigma = shared.sigma;

    for (int k=1; k<=n_iters; k++){

    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    inference_engine.ComputeInsideESS();

    std::vector<RealT> unp_i(L,0.0);
    std::vector<RealT> pr_i(L,0.0);

    std::vector<std::vector<RealT> > p_j_unp_given_i_unp;
    std::vector<std::vector<RealT> > p_j_paired_given_i_paired;

    for (int i=0; i< L; i++)
    {    
        std::vector<RealT> vec(L,0.0);
        p_j_unp_given_i_unp.push_back(vec);
        p_j_paired_given_i_paired.push_back(vec);
    }

    for (int sample_ind = 0; sample_ind < n_samples; sample_ind++)
    { 

    int seed = rand();
    inference_engine.InitRand(seed);
    std::vector<int> map = inference_engine.PredictPairingsStochasticTracebackESS();

    // count occurences from samples
    for (int i=0; i< L; i++)
        {
            if (map[i+1] == 0){
                unp_i[i] += 1;

                for (int j=0; j<L; j++){
                    if (map[j+1] == 0)
                    {
                        p_j_unp_given_i_unp[i][j] += 1;
                    }
                }
        }
        else{ // base is paired
            pr_i[i] += 1;
             for (int j=0; j<L; j++){
                        if (map[j+1] != 0)
                        {
                            p_j_paired_given_i_paired[i][j] += 1;
                        }
                    }   
        }
    }

}
    // normalize probability counts for p_unpaired vector

    for (int i=0; i< L; i++)
        {
            if (unp_i[i] > 0){

                for (int j=0; j<L; j++){
                    p_j_unp_given_i_unp[i][j] /= unp_i[i]; 
                }
            }
            unp_i[i] /= n_samples;
    }

    // normalize probability counts for p_paired vector
    for (int i=0; i< L; i++)
        {
            if (pr_i[i] > 0){

                for (int j=0; j<L; j++){
                    p_j_paired_given_i_paired[i][j] /= pr_i[i]; 
                }
            }
            pr_i[i] /= n_samples;
    }

    std::vector<RealT> error = inference_engine.GetREVIError(unp_i);
    std::vector<std::vector< double> > currREVIvec_up = inference_engine.GetREVIvec_up();
    std::vector<std::vector< double> > currREVIvec_pr = inference_engine.GetREVIvec_pr();

    std::vector<RealT> update_unp(L,0.0);
    std::vector<RealT> update_pr(L,0.0);

    // get p_i, p_i_given_j
    // write gradient

    double curr_loss = 0.0;

    // get gradient for unpaired potentials
    for (int j=0; j< L; j++){

     for (int i=0; i< L; i++)
        {    
        update_unp[j] += sigma*step_size*error[i]*unp_i[i]*(unp_i[j] - p_j_unp_given_i_unp[i][j]);
        update_pr[j] += sigma*step_size*error[i]*pr_i[i]*(pr_i[j] - p_j_paired_given_i_paired[i][j]);

            }
        update_unp[j] += -2*step_size*currREVIvec_up[0][j];
        update_pr[j] += -2*step_size*currREVIvec_pr[0][j]; //UPDATE to be curr REVI vec for paired

        }

    for (int i=0; i< L; i++)
        {   
        curr_loss += pow(currREVIvec_up[0][i],2) + sigma*pow(error[i],2) + pow(currREVIvec_pr[0][i],2); // UPDATE to be curr REVI vec for paired
        }     


    if ((abs(last_loss_2 - curr_loss) < 0.1) || k==n_iters){
        // We're done, sample and print structures one more time

    SStruct *solution;

    if (options.GetStringValue("output_bpseq_destination") != "")
    {
        solution = new SStruct(sstruct);
        //solution->SetMapping(inference_engine.PredictPairingsViterbi());

        const std::string filename = MakeOutputFilename(descriptions[nonshared.index].input_filename,
                                                        options.GetStringValue("output_bpseq_destination"),
                                                        options.GetRealValue("gamma") < 0,
                                                        shared.gamma);
        std::ofstream outfile(filename.c_str());
        if (outfile.fail()) Error("Unable to open output bpseq file '%s' for writing.", filename.c_str());
        solution->WriteBPPSEQ(outfile, inference_engine.GetREVIvec_up());
        outfile.close();
        delete solution;
    }

    if (options.GetStringValue("output_posteriors_destination") != "")
    {
        const std::string filename = MakeOutputFilename(descriptions[nonshared.index].input_filename,
                                                        options.GetStringValue("output_posteriors_destination"),
                                                        options.GetRealValue("gamma") < 0,
                                                        shared.gamma);

        inference_engine.ComputeOutsideESS();
        inference_engine.ComputePosteriorESS();

        RealT *posterior = inference_engine.GetPosterior(options.GetRealValue("output_posteriors_cutoff"));
        SparseMatrix<RealT> sparse(posterior, sstruct.GetLength()+1, RealT(0));
        delete [] posterior;
        std::ofstream outfile(filename.c_str());
        if (outfile.fail()) Error("Unable to open output posteriors file '%s' for writing.", filename.c_str());
        sparse.PrintSparseBPSEQ(outfile, sstruct.GetSequences()[0]);
        outfile.close();
    }

        for (int sample_ind = 0; sample_ind < 100; sample_ind++)
            { 
            solution = new SStruct(sstruct);
            int seed = rand();
            inference_engine.InitRand(seed);

            solution->SetMapping(inference_engine.PredictPairingsStochasticTracebackESS());
            solution->WriteParensOnly(std::cout);
            delete solution;
            }

        break;
    } else {
        // continue
        std::cerr << k << " " << curr_loss << std::endl;
        std::cerr << "p_i " << unp_i << std::endl; 
        std::cerr << "err " << error << std::endl;
        last_loss_2 = last_loss_1;
        last_loss_1 = last_loss_0;
        last_loss_0 = curr_loss;

        step_size /= 2;

        inference_engine.UpdateREVIVec(update_unp, update_pr);
    
    }


}
}
//////////////////////////////////////////////////////////////////////
// ComputationEngine::Predict()
//
// Predict structure of a single sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::Predict(std::vector<RealT> &result, 
                                       const SharedInfo<RealT> &shared,
                                       const NonSharedInfo &nonshared)
{
    result.clear();
    
    // load sequence, with constraints if necessary
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);
    if (options.GetBoolValue("use_constraints")) inference_engine.UseConstraints(sstruct.GetMapping());

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    inference_engine.UpdateEvidenceStructures();

    // perform inference
    SStruct *solution;
    if (options.GetBoolValue("viterbi_parsing"))
    {
        if (options.GetBoolValue("use_evidence"))
            Error("Viterbi parsing is not supported with evidence yet");
        // Basically, add a ComputeViterbiESS and then call it to support this.

        inference_engine.ComputeViterbi();
        inference_engine.ComputeViterbiFeatureCounts(); // HKWS why is this here?

        if (options.GetBoolValue("partition_function_only"))
        {
            std::cout << "Viterbi score for \"" << descriptions[nonshared.index].input_filename << "\": " 
                      << inference_engine.GetViterbiScore() << std::endl;
            return;
        }
        solution = new SStruct(sstruct);
        solution->SetMapping(inference_engine.PredictPairingsViterbi());
    }
    else
    {
        if (options.GetBoolValue("use_evidence")) 
        {

            inference_engine.ComputeInsideESS();
            if (options.GetBoolValue("partition_function_only"))
            {
                std::cout << "Log partition coefficient for \"" << descriptions[nonshared.index].input_filename << "\": " 
                          << inference_engine.ComputeLogPartitionCoefficientESS() << std::endl;
                return;
            }
            inference_engine.ComputeOutsideESS();
            inference_engine.ComputePosteriorESS();

        } 
        else
        {

            inference_engine.ComputeInside();
            if (options.GetBoolValue("partition_function_only"))
            {
                std::cout << "Log partition coefficient for \"" << descriptions[nonshared.index].input_filename << "\": " 
                          << inference_engine.ComputeLogPartitionCoefficient() << std::endl;
                return;
            }
            inference_engine.ComputeOutside();
            inference_engine.ComputePosterior();

        }

        solution = new SStruct(sstruct);
        if (options.GetBoolValue("centroid_estimator")) {
            std::cout << "Predicting using centroid estimator." << std::endl;
            solution->SetMapping(inference_engine.PredictPairingsPosteriorCentroid(shared.gamma));
        } else {
            std::cout << "Predicting using MEA estimator." << std::endl;
            solution->SetMapping(inference_engine.PredictPairingsPosterior(shared.gamma));
        }
    }

    // write output
    if (options.GetStringValue("output_parens_destination") != "")
    {
        const std::string filename = MakeOutputFilename(descriptions[nonshared.index].input_filename,
                                                        options.GetStringValue("output_parens_destination"),
                                                        options.GetRealValue("gamma") < 0,
                                                        shared.gamma);
        std::ofstream outfile(filename.c_str());
        if (outfile.fail()) Error("Unable to open output parens file '%s' for writing.", filename.c_str());
        solution->WriteParens(outfile);
        outfile.close();
    }
  
    if (options.GetStringValue("output_bpseq_destination") != "")
    {
        const std::string filename = MakeOutputFilename(descriptions[nonshared.index].input_filename,
                                                        options.GetStringValue("output_bpseq_destination"),
                                                        options.GetRealValue("gamma") < 0,
                                                        shared.gamma);
        std::ofstream outfile(filename.c_str());
        if (outfile.fail()) Error("Unable to open output bpseq file '%s' for writing.", filename.c_str());
        solution->WriteBPSEQ(outfile);
        outfile.close();
    }
    
    if (options.GetStringValue("output_posteriors_destination") != "")
    {
        const std::string filename = MakeOutputFilename(descriptions[nonshared.index].input_filename,
                                                        options.GetStringValue("output_posteriors_destination"),
                                                        options.GetRealValue("gamma") < 0,
                                                        shared.gamma);
        RealT *posterior = inference_engine.GetPosterior(options.GetRealValue("output_posteriors_cutoff"));
        SparseMatrix<RealT> sparse(posterior, sstruct.GetLength()+1, RealT(0));
        delete [] posterior;
        std::ofstream outfile(filename.c_str());
        if (outfile.fail()) Error("Unable to open output posteriors file '%s' for writing.", filename.c_str());
        sparse.PrintSparseBPSEQ(outfile, sstruct.GetSequences()[0]);
        outfile.close();
    }
    
    if (options.GetStringValue("output_parens_destination") == "" &&
        options.GetStringValue("output_bpseq_destination") == "" &&
        options.GetStringValue("output_posteriors_destination") == "")
    {
        WriteProgressMessage("");
        solution->WriteParens(std::cout);
    }
    
    delete solution;
}

//////////////////////////////////////////////////////////////////////
// ComputationEngine::Sample()
//
// Sample structures from predicted ensemble.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::Sample(std::vector<RealT> &result, 
                                       const SharedInfo<RealT> &shared,
                                       const NonSharedInfo &nonshared)
{
    result.clear();
    
    // load sequence, with constraints if necessary
    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);
    if (options.GetBoolValue("use_constraints")) inference_engine.UseConstraints(sstruct.GetMapping());

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    inference_engine.UpdateEvidenceStructures();

    // perform inference
    SStruct *solution;

    int N = options.GetIntValue("nsamples");

    //inference_engine.ComputePosterior(); // do we need this

    if (options.GetBoolValue("use_evidence"))
    {
        inference_engine.ComputeInsideESS();

    for (int sample_ind = 0; sample_ind < N; sample_ind++) { 

    solution = new SStruct(sstruct);
    inference_engine.InitRand(sample_ind);

    solution->SetMapping(inference_engine.PredictPairingsStochasticTracebackESS());
    solution->WriteParensOnly(std::cout);
    delete solution;
    }

    } else {
        inference_engine.ComputeInside();
    
    for (int sample_ind = 0; sample_ind < N; sample_ind++) { 

    solution = new SStruct(sstruct);
    inference_engine.InitRand(sample_ind);

    solution->SetMapping(inference_engine.PredictPairingsStochasticTraceback());
    solution->WriteParensOnly(std::cout);
    delete solution;
    }
    }

}


//////////////////////////////////////////////////////////////////////
// ComputationEngine::PredictFoldChange() HKWS
//
// Predict output affinities and fold change of a riboswitch.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void ComputationEngine<RealT>::PredictFoldChange(std::vector<RealT> &result, 
                                       const SharedInfo<RealT> &shared,
                                       const NonSharedInfo &nonshared)
{
    result.clear();

    // load parameters
    const std::vector<RealT> w(shared.w, shared.w + parameter_manager.GetNumLogicalParameters());
    inference_engine.LoadValues(w * shared.log_base);

    RealT unconditional_score_ms2_ref;
    RealT conditional_score_ms2_ref;

    std::string MS2_seq = "ACAUGAGGAUCACCCAUGU";
    std::string MS2_cst = "(((((.((....)))))))"; //contrafold seq parsing uses dots for unpaired ('?' is unspecified)
    SStruct sstruct_0;
    sstruct_0.LoadAPI(MS2_seq, MS2_cst);

    inference_engine.LoadSequence(sstruct_0); //changed from sstruct_0

    inference_engine.ComputeInside();
    unconditional_score_ms2_ref = inference_engine.ComputeLogPartitionCoefficient();

    inference_engine.UseConstraints(sstruct_0.GetMapping());
    inference_engine.UpdateEvidenceStructures();

    inference_engine.ComputeInside();
    conditional_score_ms2_ref = inference_engine.ComputeLogPartitionCoefficient();

    //RealT ligand_bonus = options.GetRealValue("ligand_bonus");

    const SStruct &sstruct = descriptions[nonshared.index].sstruct;
    inference_engine.LoadSequence(sstruct);

    const std::vector<double> &kd_data = sstruct.GetKdData(); // list of log_kd_no_lig, log_kd_w_lig, ligand_bonus
    RealT log_kd_no_lig = kd_data[0] - 0.72; // 0.72 taken as min kd value in ribologic dataset
    RealT log_kd_w_lig = kd_data[1] - 0.72; // 0.72 taken as min kd value in ribologic dataset
    RealT ligand_bonus = kd_data[2]; 
    RealT unconditional_score;
    RealT conditional_score;
    RealT conditional_score2;
    RealT conditional_score3;

    // unconditional

    inference_engine.ComputeInside();
    unconditional_score = inference_engine.ComputeLogPartitionCoefficient();

    // MS2 aptamer
    
    inference_engine.UseConstraints(sstruct.GetMapping());
    inference_engine.UpdateEvidenceStructures();
    inference_engine.ComputeInside();
    conditional_score = inference_engine.ComputeLogPartitionCoefficient();

    // ligand aptamer
    
    inference_engine.UseConstraints(sstruct.GetMapping2());
    inference_engine.ComputeInside();
    conditional_score2 = inference_engine.ComputeLogPartitionCoefficient();

    // MS2 and ligand aptamer
    
    inference_engine.UseConstraints(sstruct.GetMapping3());
    inference_engine.ComputeInside();
    conditional_score3 = inference_engine.ComputeLogPartitionCoefficient();

    RealT pred_kd_w_lig = log(exp(unconditional_score)+ligand_bonus*exp(conditional_score2)) 
                            - log(exp(conditional_score)+ligand_bonus*exp(conditional_score3))
                            - unconditional_score_ms2_ref + conditional_score_ms2_ref;

    RealT pred_kd_no_lig = unconditional_score - conditional_score - unconditional_score_ms2_ref + conditional_score_ms2_ref;

    const std::string filename = MakeOutputFilename(descriptions[nonshared.index].input_filename,
                                                    options.GetStringValue("output_bpseq_destination"),
                                                    options.GetRealValue("gamma") < 0,
                                                    shared.gamma);
    std::ofstream outfile(filename.c_str());
    if (outfile.fail()) Error("Unable to open output parens file '%s' for writing.", filename.c_str());
    outfile << descriptions[nonshared.index].input_filename << "\t" << SPrintF("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
    unconditional_score,conditional_score,conditional_score2,conditional_score3,pred_kd_no_lig,log_kd_no_lig,pred_kd_w_lig,log_kd_w_lig);
    outfile.close();

}
//////////////////////////////////////////////////////////////////////
// ComputationEngine::MakeOutputFilename()
//
// Decide on output filename, if any.  The arguments to this function
// consist of (1) a boolean variable indicating whether the output
// destination should be treated as the name of an output directory
// (and the output filename is chosen to match the input file) or
// whether the output destination should be interpreted as the output
// filename; (2) the name of the input file to be processed; and (3)
// the supplied output destination.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::string ComputationEngine<RealT>::MakeOutputFilename(const std::string &input_filename,
                                                         const std::string &output_destination,
                                                         const bool cross_validation,
                                                         const RealT gamma) const 
{
    if (output_destination == "") return "";

    const std::string dir_name = GetDirName(output_destination);
    const std::string base_name = GetBaseName(output_destination);

    const std::string prefix = (dir_name != "" ? (dir_name + DIR_SEPARATOR_CHAR) : std::string(""));
    
    // check if output directory required
    if (descriptions.size() > 1)
    {
        if (cross_validation)
        {
            return SPrintF("%s%s%c%s.gamma=%lf%c%s",
                           prefix.c_str(),
                           base_name.c_str(),
                           DIR_SEPARATOR_CHAR,
                           base_name.c_str(),
                           double(gamma),
                           DIR_SEPARATOR_CHAR,
                           GetBaseName(input_filename).c_str());
        }
        return SPrintF("%s%s%c%s",
                       prefix.c_str(),
                       base_name.c_str(),
                       DIR_SEPARATOR_CHAR,
                       GetBaseName(input_filename).c_str());
    }
    
    if (cross_validation)
    {
        return SPrintF("%s%s%c%s.gamma=%lf",
                       prefix.c_str(),
                       base_name.c_str(),
                       DIR_SEPARATOR_CHAR,
                       base_name.c_str(),
                       double(gamma));
    }
    return SPrintF("%s%s",
                   prefix.c_str(),
                   base_name.c_str());
}
