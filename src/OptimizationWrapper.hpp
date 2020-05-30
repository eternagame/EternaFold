//////////////////////////////////////////////////////////////////////
// OptimizationWrapper.hpp
//////////////////////////////////////////////////////////////////////

#ifndef OPTIMIZATIONWRAPPER_HPP
#define OPTIMIZATIONWRAPPER_HPP

#include "Config.hpp"
#include "Utilities.hpp"
#include "ComputationWrapper.hpp"
#include "CGOptimizationWrapper.hpp"
#include "InnerOptimizationWrapper.hpp"
#include "InnerOptimizationWrapperLBFGS.hpp"
#include "InnerOptimizationWrapperSubgradientMethod.hpp"
#include "InnerOptimizationWrapperEM.hpp"
#include "InnerOptimizationWrapperStochasticGradient.hpp"
#if BMRM_AVAILABLE
#include "InnerOptimizationWrapperBundleMethod.hpp"
#endif
#include "OuterOptimizationWrapper.hpp"

//////////////////////////////////////////////////////////////////////
// class OptimizationWrapper
//
// Wrapper class for performing optimization.
//////////////////////////////////////////////////////////////////////

template<class RealT>
class OptimizationWrapper
{
    ComputationWrapper<RealT> &computation_wrapper;
    std::ofstream logfile;
    int indent;
    
public:
    
    OptimizationWrapper(ComputationWrapper<RealT> &computation_wrapper);
    ~OptimizationWrapper();
    
    RealT Train(const std::vector<int> &units, std::vector<RealT> &w, std::vector<RealT> &w0, const std::vector<RealT> &C);
    RealT TrainEM(const std::vector<int> &units, std::vector<RealT> &w, const std::vector<RealT> &C, const int train_max_iter);
    RealT TrainSGD(const std::vector<int> &units, std::vector<RealT> &w, const std::vector<RealT> &C);
   
    void LearnHyperparameters(std::vector<int> units, std::vector<RealT> &values);
    void LearnHyperparametersEM(std::vector<int> units, std::vector<RealT> &values, const int train_max_iter);

    void Indent();
    void Unindent();
    void PrintMessage(const std::string &s);
    
    // getters
    const Options &GetOptions() const { return computation_wrapper.GetOptions(); }
    const std::vector<FileDescription> &GetDescriptions() const { return computation_wrapper.GetDescriptions(); }
    InferenceEngine<RealT> &GetInferenceEngine() { return computation_wrapper.GetInferenceEngine(); }
    ParameterManager<RealT> &GetParameterManager() { return computation_wrapper.GetParameterManager(); }
    ComputationEngine<RealT> &GetComputationEngine() { return computation_wrapper.GetComputationEngine(); }
    ComputationWrapper<RealT> &GetComputationWrapper() { return computation_wrapper; }
};

#include "OptimizationWrapper.ipp"

#endif
