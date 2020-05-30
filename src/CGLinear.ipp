//////////////////////////////////////////////////////////////////////
// CGLinear.ipp
//
// This file contains an implementation of the conjugate gradient
// algorithm for solving linear systems.
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// CGLinear()
//
// Implementation of conjugate gradient for solving linear
// systems Ax = b, where A is symmetric positive definite.
//////////////////////////////////////////////////////////////////////

template<class Real>
Real CGLinear<Real>::Minimize
(
     const std::vector<Real> &b,                              // right hand side
     std::vector<Real> &x,                                    // initial parameter vector
     
     const int    MAX_ITERATIONS,                             // maximum number of iterations to run CG
     const Real   SMALL_STEP_RATIO,                           // ratio beneath which steps are considered "small"
     const int    MAX_SMALL_STEPS                             // maximum number of small steps before we quit
)
{
    std::vector<Real> Ax;
    ComputeAx(Ax, x);
    std::vector<Real> r = b - Ax;
    std::vector<Real> d = r;
    Real rTr = DotProduct(r,r);
    Real f = DotProduct(x, 0.5 * Ax - b);
    
    Real best_f = f;
    std::vector<Real> best_x = x;
    Real best_rTr = rTr;
    
    int num_consecutive_small_steps = 0;
    bool progress_made = false;
    
    // report initial iteration
    
    Report(0, x, f, 0);
    
    for (int iteration = 1; iteration <= MAX_ITERATIONS; ++iteration){
        
        // compute step size
        
        std::vector<Real> Ad;
        ComputeAx(Ad, d);
        Real alpha = rTr / DotProduct(d,Ad);
        
        // update x and r
        
        x += alpha * d;

        // to prevent loss of precision
        
        if (iteration % 10 == 0)
        {
            ComputeAx(Ax, x);
            r = b - Ax;
        }
        else
        {
            r -= alpha * Ad;
            Ax = b - r;
        }
        
        // update direction
        
        Real rpTrp = rTr;
        rTr = DotProduct(r,r);
        d = r + (rTr / rpTrp) * d;
        
        // update function value
        
        f = DotProduct(x, 0.5 * Ax - b);
        Report(iteration, x, f, alpha);
        
        // note if we're making progress slowly
        
        if ((best_f - f) / Abs(best_f) < SMALL_STEP_RATIO)
        {
            num_consecutive_small_steps++;
        }
        else
        {
            num_consecutive_small_steps = 0;
            progress_made = true;
        }    
        
        if (f < best_f)
        {
            best_f = f;
            best_x = x;
            best_rTr = rTr;
        }
        
        // prevent increasing steps

        if (DotProduct(d, r) < 0)
        {
            d = r;
        }
        
        // if we're making slow progress
        
        if (num_consecutive_small_steps == MAX_SMALL_STEPS)
        {
            // give us a second chance if we made some
            // progress since the last restart
            
            if (progress_made)
            {
                progress_made = false;
                num_consecutive_small_steps = 0;
                Report("Restart: Too many consecutive small steps");
                d = r;
            }
            else
            {
                Report("Termination: Too many consecutive small steps");
                break;
            }
        }	
    }

    x = best_x;
    return Sqrt(best_rTr / DotProduct(b,b));
}

