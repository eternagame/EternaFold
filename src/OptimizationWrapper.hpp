//////////////////////////////////////////////////////////////////////
// OptimizationWrapper.hpp
//////////////////////////////////////////////////////////////////////

#ifndef OPTIMIZATIONWRAPPER_HPP
#define OPTIMIZATIONWRAPPER_HPP

#include <Config.hpp>
#include <Utilities.hpp>
#include <ComputationWrapper.hpp>

//////////////////////////////////////////////////////////////////////
// class OptimizationWrapper
//
// Wrapper class for performing optimization.
//////////////////////////////////////////////////////////////////////

class OptimizationWrapper
{
    ComputationWrapper &computation_wrapper;
    std::ofstream logfile;
    int indent;
    
public:
    
    OptimizationWrapper(ComputationWrapper &computation_wrapper);
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
    InferenceEngine &GetInferenceEngine() { return computation_wrapper.GetInferenceEngine(); }
    ParameterManager &GetParameterManager() { return computation_wrapper.GetParameterManager(); }
    ComputationEngine &GetComputationEngine() { return computation_wrapper.GetComputationEngine(); }
    ComputationWrapper &GetComputationWrapper() { return computation_wrapper; }
};

#endif
