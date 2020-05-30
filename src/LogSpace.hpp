//////////////////////////////////////////////////////////////////////
// LogSpace.hpp
//
// Routines for dealing with numbers in log space.
//////////////////////////////////////////////////////////////////////

#ifndef LOGSPACE_HPP
#define LOGSPACE_HPP

#include "Utilities.hpp"

#define NEG_INF -2e20

//////////////////////////////////////////////////////////////////////
// Fast_Exp
//
// Fast exponentiation using Chebyshev approximating polynomials,
// optimized for negative inputs.
//////////////////////////////////////////////////////////////////////

inline double Fast_Exp(double x)
{
    if (x <= double(NEG_INF/2)) return 0;
    return exp(x);
}

inline float Fast_Exp(float x)
{
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.
    
    if (x < float(-2.4915033807))
    {
        if (x < float(-5.8622823336))
        {
            if (x < float(-9.91152))
                return float(0);
            return ((float(0.0000803850)*x+float(0.0021627428))*x+float(0.0194708555))*x+float(0.0588080014);
        }
        if (x < float(-3.8396630909))
            return ((float(0.0013889414)*x+float(0.0244676474))*x+float(0.1471290604))*x+float(0.3042757740);
        return ((float(0.0072335607)*x+float(0.0906002677))*x+float(0.3983111356))*x+float(0.6245959221);
    }
    if (x < float(-0.6725053211))
    {
        if (x < float(-1.4805375919))
            return ((float(0.0232410351)*x+float(0.2085645908))*x+float(0.6906367911))*x+float(0.8682322329);
        return ((float(0.0573782771)*x+float(0.3580258429))*x+float(0.9121133217))*x+float(0.9793091728);
    }
    if (x < float(0))
        return ((float(0.1199175927)*x+float(0.4815668234))*x+float(0.9975991939))*x+float(0.9999505077);
    return (x > float(46.052) ? float(1e20) : expf(x));
}

//////////////////////////////////////////////////////////////////////
// Fast_LogExpPlusOne()
//
// Fast log(exp(x)+1) using Chebyshev approximating polynomials.
//////////////////////////////////////////////////////////////////////

inline double Fast_LogExpPlusOne(double x)
{
    Assert(double(0) <= x && x <= double(30), "Argument out-of-range.");
    return log (exp(x) + double(1));
}

inline float Fast_LogExpPlusOne(float x){
  
    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.
    
    Assert(float(0.0000000000) <= x && x <= float(11.8624794162), "Argument out-of-range.");
    if (x < float(3.3792499610))
    {
        if (x < float(1.6320158198))
        {
            if (x < float(0.6615367791))
                return ((float(-0.0065591595)*x+float(0.1276442762))*x+float(0.4996554598))*x+float(0.6931542306);
            return ((float(-0.0155157557)*x+float(0.1446775699))*x+float(0.4882939746))*x+float(0.6958092989);
        }
        if (x < float(2.4912588184))
            return ((float(-0.0128909247)*x+float(0.1301028251))*x+float(0.5150398748))*x+float(0.6795585882);
        return ((float(-0.0072142647)*x+float(0.0877540853))*x+float(0.6208708362))*x+float(0.5909675829);
    }
    if (x < float(5.7890710412))
    {
        if (x < float(4.4261691294))
            return ((float(-0.0031455354)*x+float(0.0467229449))*x+float(0.7592532310))*x+float(0.4348794399);
        return ((float(-0.0010110698)*x+float(0.0185943421))*x+float(0.8831730747))*x+float(0.2523695427);
    }
    if (x < float(7.8162726752))
        return ((float(-0.0001962780)*x+float(0.0046084408))*x+float(0.9634431978))*x+float(0.0983148903);
    return ((float(-0.0000113994)*x+float(0.0003734731))*x+float(0.9959107193))*x+float(0.0149855051);

    /*
    // Bounds for tolerance of 9.99e-05: (0, 9.21129)
    // Approximating interval: (0, 1.40131) --> ((T(-0.0118287252)*x+T(0.1342168806))*x+T(0.4976005362))*x+T(0.6932470806);
    // Approximating interval: (1.40131, 3.06792) --> ((T(-0.0117040733)*x+T(0.1232945547))*x+T(0.5276092444))*x+T(0.6721240615);
    // Approximating interval: (3.06792, 5.15409) --> ((T(-0.0027005983)*x+T(0.0419040665))*x+T(0.7762991688))*x+T(0.4152395732);
    // Approximating interval: (5.15409, 9.21129) --> ((T(-0.0001617326)*x+T(0.0040111354))*x+T(0.9666890441))*x+T(0.0929363811);
    // 4 polynomials needed.
    
    Assert(float(0.0000000000) <= x && x <= float(9.2112909219), "Argument out-of-range.");
    if (x < float(3.0679202382))
    {
        if (x < float(1.4013117629))
            return ((float(-0.0118287252)*x+float(0.1342168806))*x+float(0.4976005362))*x+float(0.6932470806);
        return ((float(-0.0117040733)*x+float(0.1232945547))*x+float(0.5276092444))*x+float(0.6721240615);
    }
    if (x < float(5.1540922927))
        return ((float(-0.0027005983)*x+float(0.0419040665))*x+float(0.7762991688))*x+float(0.4152395732);
    return ((float(-0.0001617326)*x+float(0.0040111354))*x+float(0.9666890441))*x+float(0.0929363811);
    */
}

//////////////////////////////////////////////////////////////////////
// Fast_LogExpMinusOne()
//
// Fast log(exp(x)-1) using Chebyshev approximating polynomials.
//////////////////////////////////////////////////////////////////////

inline double Fast_LogExpMinusOne(double x)
{
    Assert(double(0) <= x && x <= double(30), "Argument out-of-range.");
    return log(exp(x) - double(1));
}

inline float Fast_LogExpMinusOne(float x)
{
    // Bounds for tolerance of 9.07e-06: (0.01, 11.6105)
    // Approximating interval: (0.01, 0.0159005) --> (((T(-9371727.3239750639)*x+T(645158.5209300558))*x+T(-18614.2673037550))*x+T(316.6449790062))*x+T(-6.4566212567);
    // Approximating interval: (0.0159005, 0.0252825) --> (((T(-1466149.1313003162)*x+T(160485.3209227881))*x+T(-7362.4729488413))*x+T(199.3272540294))*x+T(-5.9928567315);
    // Approximating interval: (0.0252825, 0.0402005) --> (((T(-229370.0200164427)*x+T(39921.2557091097))*x+T(-2912.0525520632))*x+T(125.5447724820))*x+T(-5.5290922070);
    // Approximating interval: (0.0402005, 0.0639207) --> (((T(-35883.5301576035)*x+T(9930.5447448950))*x+T(-1151.7784915849))*x+T(79.1421117756))*x+T(-5.0653276832);
    // Approximating interval: (0.0639207, 0.101637) --> (((T(-5613.7580302035)*x+T(2470.2559349975))*x+T(-455.5375751385))*x+T(49.9589387475))*x+T(-4.6015631591);
    // Approximating interval: (0.101637, 0.161608) --> (((T(-878.2382307251)*x+T(614.4843546811))*x+T(-180.1535174345))*x+T(31.6053018071))*x+T(-4.1377986348);
    // Approximating interval: (0.161608, 0.256964) --> (((T(-137.3952630777)*x+T(152.8550193085))*x+T(-71.2309051609))*x+T(20.0624860440))*x+T(-3.6740341093);
    // Approximating interval: (0.256964, 0.408586) --> (((T(-21.4949296304)*x+T(38.0231865170))*x+T(-28.1487599617))*x+T(12.8030746284))*x+T(-3.2102695618);
    // Approximating interval: (0.408586, 0.649671) --> (((T(-3.3630233327)*x+T(9.4583683684))*x+T(-11.1084741802))*x+T(8.2375430535))*x+T(-2.7465046655);
    // Approximating interval: (0.649671, 1.033) --> (((T(-0.5263704166)*x+T(2.3527028558))*x+T(-4.3684611166))*x+T(5.3661953909))*x+T(-2.2827344486);
    // Approximating interval: (1.033, 1.64244) --> (((T(-0.0825320025)*x+T(0.5848902657))*x+T(-1.7021855216))*x+T(3.5600909460))*x+T(-1.8188898105);
    // Approximating interval: (1.64244, 2.61014) --> (((T(-0.0129828609)*x+T(0.1444031861))*x+T(-0.6457761629))*x+T(2.4222036000))*x+T(-1.3542078422);
    // Approximating interval: (2.61014, 4.13932) --> (((T(-0.0019714155)*x+T(0.0335196695))*x+T(-0.2229718131))*x+T(1.6981586065))*x+T(-0.8841398010);
    // Approximating interval: (4.13932, 6.61007) --> (((T(-0.0002180722)*x+T(0.0055578600))*x+T(-0.0541456822))*x+T(1.2404314137))*x+T(-0.4137049114);
    // Approximating interval: (6.61007, 11.6105) --> (((T(-0.0000062593)*x+T(0.0002549731))*x+T(-0.0039028514))*x+T(1.0266538999))*x+T(-0.0686856567);
    // 15 polynomials needed.
    
    Assert(float(0) <= x && x <= float(11.6105428289), "Argument out-of-range.");
    if (x < float(0.2569641966))
    {
        if (x < float(0.0402004692))
        {
            if (x < float(0.0159004851))
            {
                if (x < float(0.01))
                    return Log (Exp(x) - float(1));
                return (((float(-9371727.3239750639)*x+float(645158.5209300558))*x+float(-18614.2673037550))*x+float(316.6449790062))*x+float(-6.4566212567);
            }
            if (x < float(0.0252825426))
                return (((float(-1466149.1313003162)*x+float(160485.3209227881))*x+float(-7362.4729488413))*x+float(199.3272540294))*x+float(-5.9928567315);
            return (((float(-229370.0200164427)*x+float(39921.2557091097))*x+float(-2912.0525520632))*x+float(125.5447724820))*x+float(-5.5290922070);
        }
        if (x < float(0.1016370074))
        {
            if (x < float(0.0639206961))
                return (((float(-35883.5301576035)*x+float(9930.5447448950))*x+float(-1151.7784915849))*x+float(79.1421117756))*x+float(-5.0653276832);
            return (((float(-5613.7580302035)*x+float(2470.2559349975))*x+float(-455.5375751385))*x+float(49.9589387475))*x+float(-4.6015631591);
        }
        if (x < float(0.1616077721))
            return (((float(-878.2382307251)*x+float(614.4843546811))*x+float(-180.1535174345))*x+float(31.6053018071))*x+float(-4.1377986348);
        return (((float(-137.3952630777)*x+float(152.8550193085))*x+float(-71.2309051609))*x+float(20.0624860440))*x+float(-3.6740341093);
    }
    if (x < float(1.6424387600))
    {
        if (x < float(0.6496706424))
        {
            if (x < float(0.4085855305))
                return (((float(-21.4949296304)*x+float(38.0231865170))*x+float(-28.1487599617))*x+float(12.8030746284))*x+float(-3.2102695618);
            return (((float(-3.3630233327)*x+float(9.4583683684))*x+float(-11.1084741802))*x+float(8.2375430535))*x+float(-2.7465046655);
        }
        if (x < float(1.0330037540))
            return (((float(-0.5263704166)*x+float(2.3527028558))*x+float(-4.3684611166))*x+float(5.3661953909))*x+float(-2.2827344486);
        return (((float(-0.0825320025)*x+float(0.5848902657))*x+float(-1.7021855216))*x+float(3.5600909460))*x+float(-1.8188898105);
    }
    if (x < float(4.1393216929))
    {
        if (x < float(2.6101444897))
            return (((float(-0.0129828609)*x+float(0.1444031861))*x+float(-0.6457761629))*x+float(2.4222036000))*x+float(-1.3542078422);
        return (((float(-0.0019714155)*x+float(0.0335196695))*x+float(-0.2229718131))*x+float(1.6981586065))*x+float(-0.8841398010);
    }
    if (x < float(6.6100708779))
        return (((float(-0.0002180722)*x+float(0.0055578600))*x+float(-0.0541456822))*x+float(1.2404314137))*x+float(-0.4137049114);
    return (((float(-0.0000062593)*x+float(0.0002549731))*x+float(-0.0039028514))*x+float(1.0266538999))*x+float(-0.0686856567);
}

//////////////////////////////////////////////////////////////////////
// Fast_LogAdd()
// Fast_LogPlusEquals()
//
// Compute log(exp(x)+exp(y)).
//////////////////////////////////////////////////////////////////////

inline double Fast_LogAdd(double x, double y)
{
    if (x < y) std::swap (x, y);
    if (y <= double(NEG_INF/2) || x-y >= double(30)) return x;
    return Fast_LogExpPlusOne(x-y) + y;
}

inline float Fast_LogAdd(float x, float y)
{
    if (x < y) std::swap (x, y);
    if (y <= float(NEG_INF/2) || x-y >= float(11.8624794162)) return x;
    return Fast_LogExpPlusOne(x-y) + y;
}

inline void Fast_LogPlusEquals (double &x, double y)
{
    if (x < y) std::swap (x, y);
    if (y > double(NEG_INF/2) && x-y < double(30))
        x = Fast_LogExpPlusOne(x-y) + y;
}

inline void Fast_LogPlusEquals (float &x, float y)
{
    if (x < y) std::swap (x, y);
    if (y > float(NEG_INF/2) && x-y < float(11.8624794162))
        x = Fast_LogExpPlusOne(x-y) + y;
}

//////////////////////////////////////////////////////////////////////
// Fast_LogSubtract()
// Fast_LogMinusEquals()
//
// Compute log(exp(x)-exp(y)).
//////////////////////////////////////////////////////////////////////

inline double Fast_LogSubtract (double x, double y)
{
    Assert(x > y, "Cannot represent negative numbers in log space.");
    if (y <= double(NEG_INF/2) || x-y >= double(30)) return x;
    return Fast_LogExpMinusOne(x-y) + y;
}

inline float Fast_LogSubtract (float x, float y)
{
    Assert(x > y, "Cannot represent negative numbers in log space.");
    if (y <= float(NEG_INF/2) || x-y >= float(11.6105428289)) return x;
    return Fast_LogExpMinusOne(x-y) + y;
}

inline void Fast_LogMinusEquals (double &x, double y)
{
    Assert(x > y, "Cannot represent negative numbers in log space.");
    if (y > double(NEG_INF/2) && x-y < double(30))
        x = Fast_LogExpMinusOne(x-y) + y;
}

inline void Fast_LogMinusEquals (float &x, float y)
{
    Assert(x > y, "Cannot represent negative numbers in log space.");
    if (y > float(NEG_INF/2) && x-y < float(11.6105428289))
        x = Fast_LogExpMinusOne(x-y) + y;
}

#endif
