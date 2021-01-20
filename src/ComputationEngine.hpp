///////////////////////////////////////////////////////////////////////
// ComputationEngine.hpp
//
// This class provides an implementation of the DoComputation()
// routine needed by the DistributedComputation class.
//////////////////////////////////////////////////////////////////////

#ifndef COMPUTATIONENGINE_HPP
#define COMPUTATIONENGINE_HPP

#include "Config.hpp"
#include "Options.hpp"
#include "LogSpace.hpp"
#include "Utilities.hpp"
#include "SparseMatrix.hpp"
#include "SStruct.hpp"
#include "InferenceEngine.hpp"
#include "DistributedComputation.hpp"
#include "FileDescription.hpp"
#include <vector>
#include <cstdlib>

//////////////////////////////////////////////////////////////////////
// struct SharedInfo
//
// Class for storing information shared between processing units.  In
// particular, this includes:
//
//    values = all parameter values
//////////////////////////////////////////////////////////////////////

template<class RealT>
struct SharedInfo
{
    int command;
    RealT w[SHARED_PARAMETER_SIZE];
    RealT v[SHARED_PARAMETER_SIZE];
    bool use_nonsmooth;
    bool use_loss;
    RealT gamma;
    RealT sigma;
    RealT log_base;
//    std::vector<int> evidence_cpd_id;
    RealT evidence_data_scale;

    int id_base;
    int id_pairing;
    int which_data;
    int areZeros;

    RealT hyperparam_data;
    RealT kd_hyperparam_data;
    RealT lig_hyperparam_data;
};

//////////////////////////////////////////////////////////////////////
// struct NonSharedInfo
//
// Class for storing information unique to each processing unit.  In
// particular, this includes:
//
//    command = type of command to be executed
//    id = index of the work unit to be processed
//////////////////////////////////////////////////////////////////////

enum ProcessingType
{ 
    CHECK_PARSABILITY,
    COMPUTE_SOLUTION_NORM_BOUND,
    COMPUTE_GRADIENT_NORM_BOUND,
    COMPUTE_LOSS,
    COMPUTE_FUNCTION,
    COMPUTE_GRADIENT,
    COMPUTE_MSTEP_FUNCTION,
    COMPUTE_MSTEP_GRADIENT,
    COMPUTE_GAMMAMLE_FUNCTION,
    COMPUTE_GAMMAMLE_GRADIENT,
    COMPUTE_GAMMAMLE_SCALING_FACTOR,
    COMPUTE_FUNCTION_SE,
    COMPUTE_GRADIENT_SE,
    CHECK_ZEROS_IN_DATA,
    COMPUTE_HV,
    PREDICT,
    PREDICT_FOLDCHANGE,
    SAMPLE,
    TEST_ENERGIES,
    REVI,
};

struct NonSharedInfo
{
    int index;
};

//////////////////////////////////////////////////////////////////////
// class ComputationEngine
//
// Wrapper class for DistributedComputation.
//////////////////////////////////////////////////////////////////////

template<class RealT>
class ComputationEngine : public DistributedComputation<RealT, SharedInfo<RealT>, NonSharedInfo>
{
    const Options &options;
    const std::vector<FileDescription> &descriptions;
    InferenceEngine<RealT> &inference_engine;
    ParameterManager<RealT> &parameter_manager;

    std::string MakeOutputFilename(const std::string &input_filename,
                                   const std::string &output_destination,
                                   const bool cross_validation,
                                   const RealT gamma) const;

public:
    
    // constructor, destructor
    ComputationEngine(const Options &options,
                      const std::vector<FileDescription> &descriptions,
                      InferenceEngine<RealT> &inference_engine,
                      ParameterManager<RealT> &parameter_manager);
    ~ComputationEngine();

    // routine for performing an individual work unit
    void DoComputation(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);

    // methods to act on individual work units
    void CheckParsability(std::vector<RealT> &result, const NonSharedInfo &nonshared);
    void ComputeSolutionNormBound(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void ComputeGradientNormBound(std::vector<RealT> &result, const NonSharedInfo &nonshared);
    void ComputeLoss(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void ComputeFunctionAndGradient(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared, bool need_gradient);
    void ComputeMStepFunctionAndGradient(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared, bool need_gradient);
    void ComputeGammaMLEFunctionAndGradient(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared, bool need_gradient);
    void ComputeHessianVectorProduct(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void Predict(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void PredictFoldChange(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void Sample(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void RunREVI(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void TestEnergies(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void CheckZerosInData(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void ComputeGammaMLEScalingFactor(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared);
    void ComputeFunctionAndGradientSE(std::vector<RealT> &result, const SharedInfo<RealT> &shared, const NonSharedInfo &nonshared, bool need_gradient);

    // getters
    const Options &GetOptions() const { return options; }
    const std::vector<FileDescription> &GetDescriptions() const { return descriptions; }
    InferenceEngine<RealT> &GetInferenceEngine() { return inference_engine; }
    ParameterManager<RealT> &GetParameterManager() { return parameter_manager; }
};

#include "ComputationEngine.ipp"

#endif
