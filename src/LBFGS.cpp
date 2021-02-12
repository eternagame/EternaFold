//////////////////////////////////////////////////////////////////////
// LBFGS.ipp
//
// This file contains an implementation of the standard
// limited-memory BFGS optimization algorithm.
//////////////////////////////////////////////////////////////////////

#include "LBFGS.hpp"
#include "LineSearch.hpp"
#include <Utilities.ipp>

//////////////////////////////////////////////////////////////////////
// LBFGS::LBFGS()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

LBFGS::LBFGS
(
    const int    M,                                          // number of previous gradients to remember
    const Real   TERMINATION_RATIO,                          // required ratio of gradient norm to parameter norm for termination
    const int    MAX_ITERATIONS,                             // maximum number of iterations to run L-BFGS
    const Real   SMALL_STEP_RATIO,                           // ratio beneath which steps are considered "small"
    const int    MAX_SMALL_STEPS,                            // maximum number of small steps before we quit
    const Real   MAX_STEP_NORM                               // maximum norm for a single step
) :
    LineSearch(),
    M(M),
    TERMINATION_RATIO(TERMINATION_RATIO),
    MAX_ITERATIONS(MAX_ITERATIONS),
    SMALL_STEP_RATIO(SMALL_STEP_RATIO),
    MAX_SMALL_STEPS(MAX_SMALL_STEPS),
    MAX_STEP_NORM(MAX_STEP_NORM)
{}

//////////////////////////////////////////////////////////////////////
// LBFGS::Minimize()
//
// Implementation of L-BFGS optimization routine.
//////////////////////////////////////////////////////////////////////

Real LBFGS::Minimize(std::vector<Real> &x0)
{
    // initialize

    const int n = int(x0.size());
    std::vector<Real> f(2);
    std::vector<Real> gamma(2);
    std::vector<std::vector<Real> > x(2, std::vector<Real>(n));
    std::vector<std::vector<Real> > g(2, std::vector<Real>(n));
    std::vector<std::vector<Real> > s(M, std::vector<Real>(n));            
    std::vector<std::vector<Real> > y(M, std::vector<Real>(n));
    std::vector<Real> rho(M);
    Real gradient_ratio;
    Real f0;

    // check for termination criteria at beginning
    
    x[0] = x0;
    f[0] = f0 = ComputeFunction(x[0]);
    
    if (f[0] > Real(1e20))
    {
        Report(SPrintF("Termination before optimization: function value too big (%lf > %lf)", f[0], 1e20));
        return f[0];
    }

    ComputeGradient(g[0], x[0]);
    gradient_ratio = Norm(g[0]) / std::max(Real(1), Norm(x[0]));
    if (gradient_ratio < TERMINATION_RATIO)
    {
        Report(SPrintF("Termination before optimization: gradient vector small (%lf < %lf)", gradient_ratio, TERMINATION_RATIO));
        return f[0];
    }

    // initial scaling

    gamma[0] = Real(1) / Norm(g[0]);

    // report initial iteration
    
    Report(0, x[0], f[0], 0);

    // main loop

    bool progress_made = false;
    int num_consecutive_small_steps = 0;
    int k = 0;
       
    while (true)
    {
        // compute search direction, d = -H[k] g[k]

        std::vector<Real> d(-g[k%2]);
        std::vector<Real> a(M);
        for (int i = k-1; i >= k-M; i--)
        {
            a[(i+M)%M] = rho[(i+M)%M] * DotProduct(s[(i+M)%M], d);
            d -= a[(i+M)%M] * y[(i+M)%M];
        }
        
        d *= gamma[k%2];

        for (int i = k-M; i <= k-1; i++)
        {
            Real b = rho[(i+M)%M] * DotProduct(y[(i+M)%M], d);
            d += (a[(i+M)%M] - b) * s[(i+M)%M];
        }

        // perform line search, update f, and take step

        Real step = this->DoLineSearch(x[k%2], f[k%2], g[k%2], d,
                                       x[(k+1)%2], f[(k+1)%2], g[(k+1)%2],
                                       Real(0), std::min(Real(10), MAX_STEP_NORM / std::max(Real(1), Norm(d))));
        
        Report(k+1, x[(k+1)%2], f[(k+1)%2], step);
        
        // check termination conditions 
        
        if (k+1 >= MAX_ITERATIONS)
        {
            Report("Termination condition: maximum number of iterations reached");
            break; 
        }
        
        // check gradient termination condition
        
        gradient_ratio = Norm(g[(k+1)%2]) / std::max(Real(1), Norm(x[(k+1)%2]));
        if (gradient_ratio < TERMINATION_RATIO)
        {
            Report(SPrintF("Termination condition: gradient vector small (%lf < %lf)", gradient_ratio, TERMINATION_RATIO));
            break;
        }

        // heuristics for detecting slow progress (needed for large-scale
        // problems due to floating-point precision problems in gradient
        // computation)

        // check for slow progress
        
        if (step == Real(0))
            num_consecutive_small_steps = MAX_SMALL_STEPS;
        else if ((f[k%2] - f[(k+1)%2]) / std::max(Real(1), f0 - f[(k+1)%2]) < SMALL_STEP_RATIO)
            num_consecutive_small_steps++;
        else
        {
            num_consecutive_small_steps = 0;
            progress_made = true;
        }
        
        // if we're making slow progress
        
        if (num_consecutive_small_steps == MAX_SMALL_STEPS)
        {
            // give us a second chance if we made some
            // progress since the last restart
            
            if (M > 0 && progress_made)
            {
                progress_made = false;
                num_consecutive_small_steps = 0;
                Report("Restart: Too many consecutive small steps");

                for (int i = 0; i < M; i++)
                {
                    std::fill(s[i].begin(), s[i].end(), Real(0));
                    std::fill(y[i].begin(), y[i].end(), Real(0));
                    rho[i] = Real(0);
                }
            }
            else
            {
                Report("Termination: Too many consecutive small steps");
                break;
            }
        }

        // update iterates

        s[k%M] = x[(k+1)%2] - x[k%2];
        y[k%M] = g[(k+1)%2] - g[k%2];
        rho[k%M] = Real(1) / DotProduct(y[k%M], s[k%M]);

        // skip update if non-positive-definite Hessian update
        // (setting all of these quantities to zero is equivalent
        // to skipping the update, based on the BFGS recursions)

        if (!std::isfinite(rho[k%M]) || rho[k%M] <= Real(0))
        {
            std::fill(s[k%M].begin(), s[k%M].end(), Real(0));
            std::fill(y[k%M].begin(), y[k%M].end(), Real(0));
            rho[k%M] = Real(0);
        }

        // update scaling factor

        gamma[(k+1)%2] = DotProduct(s[(k-1+M)%M], y[(k-1+M)%M]) / DotProduct(y[(k-1+M)%M], y[(k-1+M)%M]);
        if (!std::isfinite(gamma[(k+1)%2]))
        {
            gamma[(k+1)%2] = gamma[k%2];
        }
        
        ++k;
    }

    x0 = x[(k+1)%2];
    return f[(k+1)%2];
}

