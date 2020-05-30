//////////////////////////////////////////////////////////////////////
// LineSearch.ipp
//
// This file contains an implementation of a safeguarded
// backtracking line search algorithm for cubic polynomial
// interpolation.
//////////////////////////////////////////////////////////////////////

template<class Real>
LineSearch<Real>::LineSearch
(
    const Real T_INIT,                                  // initial step size
    const Real MU,                                      // phi(t) <= phi(0) + MU * phi'(0) t      (sufficient decrease)
    const Real MIN_IMPROVEMENT_RATIO,                   // minimum proportion of overall improvement required to keep going
    const int MAX_EVALUATIONS,                          // maximum number of function evaluations
    const Real GAMMA1,                                  // maximum step length shrinkage 
    const Real GAMMA2                                   // minimum step length shrinkage 
) :
    T_INIT(T_INIT),
    MU(MU),
    MIN_IMPROVEMENT_RATIO(MIN_IMPROVEMENT_RATIO),
    MAX_EVALUATIONS(MAX_EVALUATIONS),
    GAMMA1(GAMMA1),
    GAMMA2(GAMMA2)
{}

//////////////////////////////////////////////////////////////////////
// UpdateQuoc()
//
// UpdateQuoc best line search result seen so far.
//////////////////////////////////////////////////////////////////////

#define UpdateQuoc(t_new, f_new) \
{ \
    if (f_new < f_best) \
    { \
        f_best = f_new; \
        t_best = t_new; \
    } \
    if (f_new <= f + MU * t_new * dot_prod) \
        sufficient_decrease = true; \
}

//////////////////////////////////////////////////////////////////////
// DoLineSearch()
//
// Implementation of backtracking line-search.  Returns a
// value of t, corresponding to the step size needed for
// updating
//
//    x <-- x + t*d
//
// As a side effect, this function also updates the value of
// the function f(x) to its new value f(x + t*d).
//////////////////////////////////////////////////////////////////////

template<class Real>
Real LineSearch<Real>::DoLineSearch
(
    const std::vector<Real> &x,                        // initial parameters
    const Real f,                                      // initial function value
    const std::vector<Real> &g,                        // initial gradient vector
    const std::vector<Real> &d,                        // initial direction vector
    
    std::vector<Real> &new_x,                          // new parameters
    Real &new_f,                                       // new function value
    std::vector<Real> &new_g,                          // new gradient vector
    
    const Real T_MIN,                                  // minimum step size
    const Real T_MAX                                   // maximum step size
)
{
    Assert(T_MIN <= T_MAX, "Line search called with T_MIN > T_MAX.");
    const Real dot_prod = DotProduct(d, g);
    bool sufficient_decrease = false;

    // try initial point

    Real t_best = Real(0), t_last = T_INIT, t_prev = Real(0);
    t_last = std::max(T_MIN, t_last);
    t_last = std::min(T_MAX, t_last);
    Real f_best = f, f_last = ComputeFunction(x + t_last * d), f_prev = Real(0);
    UpdateQuoc(t_last, f_last);
    
    for (int iteration = 2; iteration <= MAX_EVALUATIONS; ++iteration)
    {
        // termination criteria
        
        if (sufficient_decrease && iteration > 2 && (f_prev - f_last) / std::max(Real(1), f - f_best) < MIN_IMPROVEMENT_RATIO) break;

        Real t_new;
        
        if (iteration == 2)
        {
            
            // fit using a quadratic: at^2 + bt + c
            //
            // This function must be equal to f at t=0, f_last at t=t_last,
            // and must match the directional derivative of f at t=0.
            
            Real a = ((f_last - f) / t_last - dot_prod) / t_last;
            Real b = dot_prod;
            
            t_new = -b / (Real(2)*a);
            
        }
        else
        {
            
            // fit using a cubic: at^3 + bt^2 + ct + d
            //
            // This function must be equal to f at t=0, f_last at t=t_last,
            // f_prev at t=t_prev, and must match the directional derivative 
            // of f at t=0.
            //
            // c and d are easily obtained.  The remaining two terms amount
            // to solving a linear system of the form
            //
            //   [ t_prev^3  t_prev^2 ] [ a ] = [ f_prev - d - c t_prev ]
            //   [ t_last^3  t_last^2 ] [ b ]   [ f_last - d - c t_last ]
            
            Real c = dot_prod;
            Real d = f;
            
            Real f_l = f_last - d - c * t_last;
            Real f_p = f_prev - d - c * t_prev;
            
            Real a = (f_l / (t_last * t_last) - f_p / (t_prev * t_prev)) / (t_last - t_prev);
            Real b = (f_l / (t_last * t_last * t_last) - f_p / (t_prev * t_prev * t_prev)) / (Real(1) / t_last - Real(1) / t_prev);
            
            Real A = Real(3)*a;
            Real B = Real(2)*b;
            Real C = c;
            
            t_new = (-B + sqrt(B*B - Real(4)*A*C)) / (Real(2)*A);  // pick the left root
        }
        
        // use safe-guarding: clip t to a safe range and evaluate function
        
        if (!std::isfinite(t_new)) t_new = t_last;

        Real lower_bound = T_MIN;
        Real upper_bound = T_MAX;

        if (iteration > 3)
        {
            lower_bound = std::max(lower_bound, GAMMA1 * t_last);
            upper_bound = std::min(upper_bound, GAMMA2 * t_last);
        }
        
        if (lower_bound > upper_bound) break;
        
        t_new = std::max(lower_bound, t_new);
        t_new = std::min(upper_bound, t_new);
        
        // now, move to this point and update iterates
        
        Real f_new = ComputeFunction(x + t_new * d);
        UpdateQuoc(t_new, f_new);
        t_prev = t_last; f_prev = f_last;
        t_last = t_new; f_last = f_new;
    }

    new_f = f_best;
    new_x = x + t_best * d;
    ComputeGradient(new_g, new_x);
    return t_best;    
}
