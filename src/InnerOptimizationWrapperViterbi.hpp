//////////////////////////////////////////////////////////////////////
// InnerOptimizationWrapperViterbi.hpp
//
// Inner optimization algorithm.
//////////////////////////////////////////////////////////////////////

#ifndef INNEROPTIMIZATIONWRAPPERVITERBI_HPP
#define INNEROPTIMIZATIONWRAPPERVITERBI_HPP

#include <Computation.hpp>
#include <OptimizationWrapper.hpp>
#include <InnerOptimizationWrapperViterbi.hpp>
#include <SubgradientDescent.hpp>

//////////////////////////////////////////////////////////////////////
// class InnerOptimizationWrapperViterbi
//////////////////////////////////////////////////////////////////////


class InnerOptimizationWrapperViterbi : public SubgradientDescent<double>
{
    OptimizationWrapper *optimizer;
    const std::vector<int> units;
    const std::vector<double> C;
    std::vector<double> best_x;
    std::vector<double> bias;
    double best_f;  

public:
    InnerOptimizationWrapperViterbi(OptimizationWrapper *optimizer,
                                    const std::vector<int> &units,
                                    const std::vector<double> &C);

    void LoadBias(const std::vector<double> &bias);
    double ComputeFunction(const std::vector<double> &x);
    void ComputeSubgradient(std::vector<double> &g, const std::vector<double> &x);
    void Report(int iteration, double f, const std::vector<double> &x, const std::vector<double> &g,
                double norm_bound, double step_size);
    void Report(const std::string &s);
};

#endif
