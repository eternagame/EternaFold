//////////////////////////////////////////////////////////////////////
// Subgradientmethod.ipp
//
// This file contains an implementation of an improved
// Pegasos style subgradient optimization algorithm.
//////////////////////////////////////////////////////////////////////

#include "SubgradientMethod.hpp"
#include <Utilities.ipp>

#define DOUBLY_ADAPTIVE           0
#define PROXIMAL_ADAPTIVE         1
#define BARTLETT_ADAPTIVE         0
#define DYNAMIC_STEPSIZE          0

const double LOWER_BOUND = 0.0;
const double TOLERANCE = 1e-4;
const double ACCEPTANCE_RATIO = 0.15;
const double IMPROVEMENT_FRACTION = 0.75;
const double MULTIPLIER = 0.95;
const int MAX_INNER_STEPS = 75;
const int MAX_NON_IMPROVEMENT_STEPS = 20;

//////////////////////////////////////////////////////////////////////
// SubgradientMethod::SubgradientMethod()
//
// Constructor.
//////////////////////////////////////////////////////////////////////


SubgradientMethod::SubgradientMethod
(
    const int     MAX_ITERATIONS,                     // maximum number of iterations to run subgradient method
    const RealT   PARAMETER_NORM_BOUND,               // maximum parameter vector norm
    const RealT   GRADIENT_NORM_BOUND,                // maximum gradient vector norm
    const RealT   CURVATURE                           // strong convexity constant: f_t(v) >= f_t(w) + g'*(v - w) + 0.5 * CURVATURE * ||v - w||^2
) :
    MAX_ITERATIONS(MAX_ITERATIONS),
    PARAMETER_NORM_BOUND(PARAMETER_NORM_BOUND),
    GRADIENT_NORM_BOUND(GRADIENT_NORM_BOUND),
    CURVATURE(CURVATURE)
{}

//////////////////////////////////////////////////////////////////////
// SubgradientMethod::Minimize()
//
// Implementation of an improved Pegasos style subgradient
// optimization algorithm.
//////////////////////////////////////////////////////////////////////


RealT SubgradientMethod::Minimize(std::vector<RealT> &x)
{
    std::vector<RealT> g;
    ComputeSubgradient(g, x);
    RealT f = ComputeFunction(x);

#if DOUBLY_ADAPTIVE || PROXIMAL_ADAPTIVE || BARTLETT_ADAPTIVE
    RealT bound = 0;
    RealT sigma_sum = 0;
    RealT tau_sum = 0;
#endif
#if DOUBLY_ADAPTIVE
    RealT gamma_sum = 0;
#endif

    // check early termination criteria
    
    if (f >= RealT(1e20))
    {
        Report(SPrintF("Termination before optimization: function value too big (%lf > %lf)", f, 1e20));
        return f;
    }

    // keep track of best parameter vector

    RealT best_f = f;
    std::vector<RealT> best_x = x;
    std::vector<RealT> best_g = g;

#if DYNAMIC_STEPSIZE
    RealT target_value = std::max(LOWER_BOUND, f - DotProduct(g, g) / 2.0);
    RealT outer_acceptance_interval = std::max((best_f - target_value) * ACCEPTANCE_RATIO, TOLERANCE);
    RealT path_length = 0;
    int inner_counter = 0;
    int non_improvement_counter = 0;
    
    RealT delta = TOLERANCE;
    int failure_count = 0;
#endif

    // run optimization algorithm

    for (int epoch = 1; epoch <= MAX_ITERATIONS; epoch++)
    {
        // compute learning rate

#if DOUBLY_ADAPTIVE
        RealT At = 0.5 * Pow(PARAMETER_NORM_BOUND + Norm(x), RealT(2));
        RealT Bt = 0.5 * PARAMETER_NORM_BOUND * PARAMETER_NORM_BOUND + 0.5 * DotProduct(x, x);

        sigma_sum += CURVATURE;
        RealT sum = sigma_sum + tau_sum + gamma_sum;
        RealT tau = (-sum + Sqrt(sum * sum + 4 * (1 + At/Bt) * DotProduct(g, g) / At)) / (2 * (1 + At/Bt));
        RealT gamma = tau * At/Bt;
        tau_sum += tau;
        gamma_sum += gamma;

        RealT eta = 1.0 / (sigma_sum + tau_sum + gamma_sum);
        bound += tau * At + gamma * Bt + DotProduct(g, g) * eta;
        g += gamma * x;
#endif

#if PROXIMAL_ADAPTIVE
        sigma_sum += CURVATURE;
        RealT tau = 0.5 * (-(sigma_sum + tau_sum) + Sqrt(Pow(sigma_sum + tau_sum, RealT(2)) + RealT(4) * DotProduct(g, g) / Pow(PARAMETER_NORM_BOUND + Norm(x), RealT(2))));
        tau_sum += tau;
        RealT eta = 1.0 / (sigma_sum + tau_sum);
        bound += 0.5 * tau * Pow(PARAMETER_NORM_BOUND + Norm(x), RealT(2)) + 0.5 * DotProduct(g, g) / (sigma_sum + tau_sum);
#endif

#if BARTLETT_ADAPTIVE
        sigma_sum += CURVATURE;
        RealT tau = 0.5 * (-(sigma_sum + tau_sum) + Sqrt(Pow(sigma_sum + tau_sum, RealT(2)) + RealT(8) * DotProduct(g, g) / (3*PARAMETER_NORM_BOUND*PARAMETER_NORM_BOUND)));
        tau_sum += tau;

        RealT eta = 1.0 / (sigma_sum + tau_sum);
        bound += 0.5 * tau * PARAMETER_NORM_BOUND * PARAMETER_NORM_BOUND + 0.5 * Pow(Norm(g) + tau*PARAMETER_NORM_BOUND, 2.0) / (sigma_sum + tau_sum);
#endif
        
#if DYNAMIC_STEPSIZE
        RealT eta = MULTIPLIER * (f - target_value) / std::max(DotProduct(g, g), 1e-10);
#endif
        // take a step

        x -= eta * g;

        // project back to ball

        RealT norm = Norm(x);
        if (norm > PARAMETER_NORM_BOUND)
        {
            x *= PARAMETER_NORM_BOUND / norm;
        }

        // compute new subgradient and function
        
        ComputeSubgradient(g, x);
        f = ComputeFunction(x);

#if DYNAMIC_STEPSIZE
        ++inner_counter;
        
        if (f < best_f)
        {
            non_improvement_counter = 0;
            path_length += best_f - f;
            
            // outer loop success

            if (f <= target_value + outer_acceptance_interval)
            {
                target_value = std::max(LOWER_BOUND, f - outer_acceptance_interval - IMPROVEMENT_FRACTION * path_length);
                outer_acceptance_interval = std::max((f - target_value) * ACCEPTANCE_RATIO, TOLERANCE);
                path_length = 0;
                Report(SPrintF("Outer loop success after %d inner iterations: new target value = %lf, best f = %lf", inner_counter, double(target_value), double(best_f)));
                inner_counter = 0;
            }

            // outer loop failure
            
            else if (inner_counter >= MAX_INNER_STEPS)
            {
                target_value = std::max(LOWER_BOUND, (f - outer_acceptance_interval + target_value) / 2);
                outer_acceptance_interval = std::max((f - target_value) * ACCEPTANCE_RATIO, TOLERANCE);
                path_length = 0;
                Report(SPrintF("Outer loop failure after %d inner iterations: new target value = %lf, best f = %lf", inner_counter, double(target_value), double(best_f)));
                inner_counter = 0;
            }
        }
        else
        {
            ++non_improvement_counter;

            // outer loop failure

            if (inner_counter >= MAX_INNER_STEPS || non_improvement_counter >= MAX_NON_IMPROVEMENT_STEPS)
            {
                target_value = std::max(LOWER_BOUND, (best_f - outer_acceptance_interval + target_value) / 2);
                outer_acceptance_interval = std::max((best_f - target_value) * ACCEPTANCE_RATIO, TOLERANCE);
                path_length = 0;
                Report(SPrintF("Outer loop failure after %d inner iterations: new target value = %lf, best f = %lf", inner_counter, double(target_value), double(best_f)));
                inner_counter = 0;
                non_improvement_counter = 0;
            }
        }
#endif

        // update best parameter values

        if (Norm(x) != 0 && (f < best_f || Norm(best_x) == 0))
        {
            best_f = f;
            best_g = g;
            best_x = x;
        }

        // print updates

        const int update_frequency = std::max(1, MAX_ITERATIONS / 100);
        if (epoch % update_frequency == 0)
        {
            Report(epoch, best_x, best_f, best_g, PARAMETER_NORM_BOUND, eta);
        }

        // check convergence criteria

        if (epoch >= MAX_ITERATIONS)
        {
            Report("Termination condition: maximum number of iterations reached");
            break; 
        }
    }

    Report(SPrintF("Cumulative regret bound: %lf", double(bound)));

    x = best_x;
    return best_f;
}
