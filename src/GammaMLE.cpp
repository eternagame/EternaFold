//////////////////////////////////////////////////////////////////////
// GammaMLE.ipp
//
// This file contains an implementation of MLE for a Gamma function
//////////////////////////////////////////////////////////////////////

#include <Config.hpp>
#include "GammaMLE.hpp"

//////////////////////////////////////////////////////////////////////
// GammaMLE::GammaMLE()
//
// Constructor.
//////////////////////////////////////////////////////////////////////


GammaMLE::GammaMLE
(
    const int    MAX_ITER_K,
    const Real   THRESH_K
) :
    MAX_ITER_K(MAX_ITER_K), THRESH_K(THRESH_K)
{}


Real GammaMLE::Psi(Real z)
{
    // Based on: http://en.wikipedia.org/wiki/Digamma_function#Computation_.26_approximation

    Real x = 0;
    if (fabs(z) < 5)
    {
        x = z;
        return Psi(x + 5) - 1/x - 1/(x + 1) - 1/(x + 2) - 1/(x + 3) - 1/(x + 4);
    }

    x = 1/(z*z);
    return log(z) - 1/(2 * z) + ((x * (-1/12 + ((x * (1/120 + ((x * (-1/252 + ((x * (1/240 + ((x * (-1/132 + ((x * (691/32760 + ((x * (-1/12 + (3617 * x)/8160)))))))))))))))))))));

}


Real GammaMLE::PsiPrime(Real z)
{
    Real x = 0;
    if (fabs(z) < 5)
    {
        x = z;
        return PsiPrime(x + 5) + 1/(x*x) + 1/((x + 1)*(x + 1)) + 1/((x + 2)*(x + 2)) + 1/((x + 3)*(x + 3)) + 1/((x + 4)*(x + 4));
    }

    x = 1/(z*z);
    return 1/(2 * z*z) + (1 + (x * (1/6 + (x * (-1/30 + (x * (1/42 + (x * (-1/30 + (x * (5/66 + (x * (-691/2370 + (x * (7/6 -(3617 * x)/510)))))))))))))))/z;

}


//////////////////////////////////////////////////////////////////////
// GammaMLE::Minimize()
//
// Implementation of L-BFGS optimization routine.
//////////////////////////////////////////////////////////////////////

Real GammaMLE::Minimize(std::vector<Real> &x0, std::vector<std::vector<bool> > config_params, int which_data)
{
    // theta = (1/k*N) * sum{x_i}
    // optimize k using Newton-Raphson update

    Real s = 0;

    Real current_theta = 0;
    Real current_k = 0;
    Real current_k_inv = 0;
    Real gammamle_ss_sum = 0;
    Real gammamle_ss_sumlog = 0;
    Real num_examples = 0;
    bool use_MM;

    int iter2 = 0;
    Real old_k2 = 0;
    Real diff_k2 = 0;

    Real ll = 0;

    std::vector<Real> g(3);
    std::vector<Real> g2(2);
    Real scale = 0;

    int index_k = 0;
    int index_theta = 0;

    // Based on: http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf

    for (int i = 0; i < M; i++)  // nucleotide
    {
        for (int j = 0; j < 2; j++)  // base pairing
        {
            index_k = GetLogicalIndex(0,i,j,which_data);
            index_theta = GetLogicalIndex(1,i,j,which_data);

            current_k = exp(x0[index_k]);
            current_theta = exp(x0[index_theta]);

            use_MM = config_params[which_data][i*2 +j];

            scale = 1;

            ComputeGammaMLEGradient(g,x0,i,j,(int)use_MM,scale,which_data);
            ll += ComputeGammaMLEFunction(x0,i,j,(int)use_MM,scale,which_data);  // this should be cached from gradient

            gammamle_ss_sum = g[0];
            gammamle_ss_sumlog = g[1];
            num_examples = g[2];  // find N

            if (use_MM)
            {
                s = gammamle_ss_sumlog / (num_examples-1) - (gammamle_ss_sum/num_examples)*(gammamle_ss_sum/num_examples)*num_examples/(num_examples-1);  // variance = std^2 = sum {d^2} / n-1 - mu^2 * n/(n-1)
                current_k = (gammamle_ss_sum/num_examples)*(gammamle_ss_sum/num_examples)/s;
                current_theta = s/(gammamle_ss_sum/num_examples);
            }
            else
            {
                // find k via the Newton-Raphson update
                s = log(gammamle_ss_sum / num_examples) - gammamle_ss_sumlog / num_examples;
                iter2 = 0;
                while (iter2++ < MAX_ITER_K)
                {
                    old_k2 = current_k;

                    current_k_inv = 1/current_k + (log(current_k) - Psi(current_k) - s)/(current_k*current_k*(1/current_k - PsiPrime(current_k)));
                    current_k = 1/current_k_inv;

                    diff_k2 = fabs(current_k - old_k2);
                    if (diff_k2 < THRESH_K)
                        break;
                }
                // find theta
                current_theta = gammamle_ss_sum / (num_examples * current_k);
            }


            x0[index_k] = log(current_k);
            x0[index_theta] = log(current_theta);

        }
    }

    return ll;
}

