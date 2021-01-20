#include "Options.hpp"
#include "FileDescription.hpp"
#include "ComputationEngine.hpp"
#include "Config.hpp"
#include "Utilities.hpp"
#include "ComputationWrapper.hpp"
#include "InferenceEngine.hpp"
#include "ParameterManager.hpp"
#include "OptimizationWrapper.hpp"
#include "SStruct.hpp"
#include "Defaults.ipp"
#include "contrafold_api.hpp"

void InitOptions(Options &options)
{
    // register default options
    options.SetStringValue("training_mode", "");
    options.SetBoolValue("verbose_output", false);
    options.SetRealValue("log_base", 1.0);
    options.SetBoolValue("viterbi_parsing", false);
    options.SetBoolValue("allow_noncomplementary", false);

    options.SetStringValue("parameter_filename", "");
    options.SetBoolValue("use_constraints", false);
    options.SetBoolValue("centroid_estimator", false);
    options.SetBoolValue("use_evidence", false);
    //options.SetRealValue("gamma", GAMMA_DEFAULT);
    options.SetStringValue("output_parens_destination", "");
    options.SetStringValue("output_bpseq_destination", "");
    options.SetRealValue("output_posteriors_cutoff", 0);
    options.SetStringValue("output_posteriors_destination", "");
    options.SetBoolValue("partition_function_only", false);

    options.SetBoolValue("gradient_sanity_check", false);
    options.SetRealValue("holdout_ratio", 0);
    //options.SetRealValue("regularization_coefficient", REGULARIZATION_DEFAULT);

    options.SetStringValue("train_examplefile", "");
    //options.SetIntValue("train_max_iter", TRAIN_MAX_ITER_DEFAULT);
    options.SetStringValue("train_initweights_filename", "");
    options.SetStringValue("train_priorweights_filename", "");
    options.SetIntValue("num_data_sources",0);
    options.SetIntValue("batch_size", 1);
    options.SetRealValue("s0", 0.0001);
    options.SetRealValue("s1", 0);
    //options.SetRealValue("hyperparam_data",HYPERPARAM_DATA_DEFAULT);
}

float pfunc(char* seq, char* c)
{
    std::string sequence(seq);
    std::string constraints(c);
    Options options;
    InitOptions(options);
    options.SetBoolValue("use_constraints", true);
    options.SetBoolValue("partition_function_only", true);
    options.SetBoolValue("allow_noncomplementary", false);
    std::string constr;

    if (constraints.compare("?") == 0){
    	constr.assign(sequence.length(),'?');
    }
    else{
    	constr.assign(constraints);
    }

	std::cout << "constraints: " << constr << std::endl;

    ParameterManager<float> parameter_manager;
    InferenceEngine<float> inference_engine(options.GetBoolValue("allow_noncomplementary"),0, options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);

    SStruct sstruct;
    sstruct.LoadAPI(sequence, constr);
    inference_engine.LoadSequence(sstruct);
    inference_engine.UseConstraints(sstruct.GetMapping());

	std::vector<FileDescription> descriptions;

    ComputationEngine<float> computation_engine(options, descriptions, inference_engine, parameter_manager);

    inference_engine.ComputeInside();
    float Z = inference_engine.ComputeLogPartitionCoefficient();
    inference_engine.ComputeOutside();
    inference_engine.ComputePosterior();

    return Z;
}


char* predict_struct(char* seq, char* c)
{
    std::string sequence(seq);
    std::string constraints(c);
    Options options;
    InitOptions(options);
    options.SetBoolValue("use_constraints", true);
    options.SetBoolValue("partition_function_only", true);
    options.SetBoolValue("allow_noncomplementary", false);
    std::string constr;

    if (constraints.compare("?") == 0){
    	constr.assign(sequence.length(),'?');
    }
    else{
    	constr.assign(constraints);
    }

	std::cout << "constraints: " << constr << std::endl;

    ParameterManager<float> parameter_manager;
    InferenceEngine<float> inference_engine(options.GetBoolValue("allow_noncomplementary"),0, options.GetRealValue("kappa"));
    inference_engine.RegisterParameters(parameter_manager);

    SStruct sstruct;
    sstruct.LoadAPI(sequence, constr);
    inference_engine.LoadSequence(sstruct);
    inference_engine.UseConstraints(sstruct.GetMapping());

	std::vector<FileDescription> descriptions;

    ComputationEngine<float> computation_engine(options, descriptions, inference_engine, parameter_manager);

	SStruct *solution;
    solution = new SStruct(sstruct);

    //inference_engine.ComputeViterbi();
    inference_engine.ComputeInside();
    inference_engine.ComputeOutside();
    inference_engine.ComputePosterior();

    //solution->SetMapping(inference_engine.PredictPairingsViterbi());
    solution->SetMapping(inference_engine.PredictPairingsPosterior(6));

    std::string structure = solution->ReturnMappingAsParens();

    char* outpt = new char[structure.length() + 1];
    strcpy(outpt, structure.c_str());

    return outpt;
}

float c_fold(char* sequence, char* structure){
	structure = predict_struct(sequence);
	float logZ = pfunc(sequence);
	float logZ_MEA = pfunc(sequence, structure);
	return logZ_MEA - logZ;
}

float c_energy_of_structure(char* sequence, char* structure){
	float logZ = pfunc(sequence);
	float logZ_MEA = pfunc(sequence, structure);
	return logZ_MEA - logZ;
}
