//////////////////////////////////////////////////////////////////////
// ComputationWrapper.hpp
//
// This class provides a wrapper around the Computation class that
// provides a framework for translating basic queries into the format
// needed by the Computation class.  This class also provides caching
// facilities for preventing redundant computations.
//////////////////////////////////////////////////////////////////////

#ifndef COMPUTATIONWRAPPER_HPP
#define COMPUTATIONWRAPPER_HPP

#include "ComputationEngine.hpp"

//////////////////////////////////////////////////////////////////////
// class ComputationWrapper
//
// Wrapper class for Computation.
//////////////////////////////////////////////////////////////////////

class ComputationWrapper
{
    ComputationEngine &computation_engine;

    SharedInfo shared_info;
    std::vector<NonSharedInfo> nonshared_info;
    
    // the following member variables are used to "cache" work to
    // ensure that it is not repeated unnecessarily
    bool cached_toggle_use_nonsmooth;
    bool cached_toggle_use_loss;
    std::vector<int> cached_units;
    std::vector<double> cached_w;
    std::vector<double> cached_function;
    std::vector<double> cached_gradient;

    std::vector<double> cached_function_gammamle;
    std::vector<double> cached_gradient_gammamle;
    
public:
    
    // constructor, destructor
    ComputationWrapper(ComputationEngine &computation_engine);
    ~ComputationWrapper();

    // retrieve list of work units
    std::vector<int> GetAllUnits() const;
    
    // methods to act on vectors of work units
    std::vector<int> FilterNonparsable(const std::vector<int> &units);
    RealT ComputeSolutionNormBound(const std::vector<int> &units, const std::vector<RealT> &C, RealT log_base);
    RealT ComputeGradientNormBound(const std::vector<int> &units, const std::vector<RealT> &C, RealT log_base);
    void Predict(const std::vector<int> &units, const std::vector<RealT> &w, RealT gamma, RealT log_base);
    void PredictFoldChange(const std::vector<int> &units, const std::vector<RealT> &w, RealT gamma, RealT log_base);
    void Sample(const std::vector<int> &units, const std::vector<RealT> &w, RealT gamma, RealT log_base);
    void RunREVI(const std::vector<int> &units, const std::vector<RealT> &w, RealT gamma, RealT log_base, RealT sigma);
    RealT ComputeLoss(const std::vector<int> &units, const std::vector<RealT> &w, RealT log_base);
    RealT ComputeFunction(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base, RealT hyperparam_data, RealT kd_hyperparam_data, RealT lig_hyperparam_data);
    std::vector<RealT> ComputeGradient(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base, RealT hyperparam_data, RealT kd_hyperparam_data, RealT lig_hyperparam_data);
    RealT ComputeEMFunction(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base);
    std::vector<RealT> ComputeEMGradient(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base);
    RealT ComputeFunctionSE(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base, RealT hyperparam_data, RealT kd_hyperparam_data, RealT lig_hyperparam_data);
    std::vector<RealT> ComputeGradientSE(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base, RealT hyperparam_data, RealT kd_hyperparam_data, RealT lig_hyperparam_data);
    RealT ComputeGammaMLEFunction(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base, int evidence_cpd_id1, int evidence_cpd_id2, int evidence_cpd_id3, RealT evidence_data_scale, int which_data);
    std::vector<RealT> ComputeGammaMLEGradient(const std::vector<int> &units, const std::vector<RealT> &w, bool toggle_use_nonsmooth, bool toggle_use_loss, RealT log_base, int evidence_cpd_id1, int evidence_cpd_id2, int evidence_cpd_id3, RealT evidence_data_scale, int which_data);

    bool FindZerosInData(const std::vector<int> &units, int evidence_cpd_id1, int evidence_cpd_id2, int which_data);
    std::vector<RealT> ComputeGammaMLEScalingFactor(const std::vector<int> &units, const std::vector<RealT> &w, int evidence_cpd_id1, int evidence_cpd_id2, int which_data);

    std::vector<RealT> ComputeHessianVectorProduct(const std::vector<int> &units, const std::vector<RealT> &w, const std::vector<RealT> &v, bool toggle_use_loss, RealT log_base);
    
    // for debugging
    void SanityCheckGradient(const std::vector<int> &units, const std::vector<RealT> &w);
    void TestEnergies(const std::vector<int> &units, const std::vector<RealT> &w, RealT gamma, RealT log_base);

    // getters
    const Options &GetOptions() const { return computation_engine.GetOptions(); }
    const std::vector<FileDescription> &GetDescriptions() const { return computation_engine.GetDescriptions(); }
    InferenceEngine &GetInferenceEngine() { return computation_engine.GetInferenceEngine(); }
    ParameterManager &GetParameterManager() { return computation_engine.GetParameterManager(); }
    ComputationEngine &GetComputationEngine() { return computation_engine; }
};

#endif
