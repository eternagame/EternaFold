//////////////////////////////////////////////////////////////////////
// Utilities.ipp
//////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename T3>
inline triple<T1,T2,T3>::triple() :
    first(), second(), third() 
{}

template<typename T1, typename T2, typename T3>
inline triple<T1,T2,T3>::triple(const T1 &first, const T2 &second, const T3 &third) :
    first(first), second(second), third(third) 
{}

template<typename T1, typename T2, typename T3>
inline triple<T1,T2,T3>::triple(const triple<T1,T2,T3> &rhs) :
    first(rhs.first), second(rhs.second), third(rhs.third) 
{}

template<typename T1, typename T2, typename T3>
inline bool operator==(const triple<T1,T2,T3> &x,
                       const triple<T1,T2,T3> &y)
{
    return
        x.first == y.first &&
        x.second == y.second &&
        x.third == y.third;
}

template<typename T1, typename T2, typename T3>
inline bool operator<(const triple<T1,T2,T3> &x,
                      const triple<T1,T2,T3> &y)
{
    return 
        x.first < y.first ||
        !(y.first < x.first) &&
        (x.second < y.second ||
         !(y.second < x.second) &&
         x.third < y.third);
}

template<typename T1, typename T2, typename T3>
inline bool operator!=(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y)
{
    return !(x == y);
}

template<typename T1, typename T2, typename T3>
inline bool operator>(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y)
{
    return y < x;
}

template<typename T1, typename T2, typename T3>
inline bool operator<=(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y)
{
    return !(y < x);
}

template<typename T1, typename T2, typename T3>
inline bool operator>=(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y)
{
    return !(x < y);
}

template<typename T1, typename T2, typename T3>
inline triple<T1,T2,T3> make_triple(T1 first, T2 second, T3 third)
{
    return triple<T1,T2,T3>(first, second, third);
}

template<typename T1, typename T2>
std::ostream &operator<<(std::ostream &out, const std::pair<T1,T2> &x)
{
    out << '(' << x.first << ',' << x.second << ')';
    return out;
}

template<typename T1, typename T2, typename T3>
std::ostream &operator<<(std::ostream &out, const triple<T1,T2,T3> &x)
{
    out << '(' << x.first << ',' << x.second << ',' << x.third << ')';
    return out;
}

#ifndef sqrtl
#define sqrtl(x) (static_cast<long double>(sqrt(double(x))))
#endif
#ifndef expl
#define expl(x) (static_cast<long double>(exp(double(x))))
#endif
#ifndef logl
#define logl(x) (static_cast<long double>(log(double(x))))
#endif
#ifndef powl
#define powl(x,y) (static_cast<long double>(pow(double(x),double(y))))
#endif
#ifndef tanhl
#define tanhl(x) (static_cast<long double>(tanh(double(x))))
#endif
#ifndef floorl
#define floorl(x) (static_cast<long double>(floor(double(x))))
#endif
#ifndef ceill
#define ceill(x) (static_cast<long double>(ceil(double(x))))
#endif

template<> inline float Sqrt(const float x) { return sqrtf(x); }
template<> inline double Sqrt(const double x) { return sqrt(x); }
template<> inline long double Sqrt(const long double x) { return sqrtl(x); }

template<> inline float Exp(const float x) { return expf(x); }
template<> inline double Exp(const double x) { return exp(x); }
template<> inline long double Exp(const long double x) { return expl(x); }

template<> inline float Log(const float x) { return logf(x); }
template<> inline double Log(const double x) { return log(x); }
template<> inline long double Log(const long double x) { return logl(x); }

template<> inline float Pow(const float x, const float p) { return powf(x,p); }
template<> inline double Pow(const double x, const double p) { return pow(x,p); }
template<> inline long double Pow(const long double x, const long double p) { return powl(x,p); }

template<> inline float Tanh(const float x) { return tanhf(x); }
template<> inline double Tanh(const double x) { return tanh(x); }
template<> inline long double Tanh(const long double x) { return tanhl(x); }

template<> inline float Floor(const float x) { return floorf(x); }
template<> inline double Floor(const double x) { return floor(x); }
template<> inline long double Floor(const long double x) { return floorl(x); }

template<> inline float Ceil(const float x) { return ceilf(x); }
template<> inline double Ceil(const double x) { return ceil(x); }
template<> inline long double Ceil(const long double x) { return ceill(x); }

inline int RandInt(const int n) {
    // generate a pseudo random int uniformly distributed in the inclusive range 0 - (n-1)
    // invariant: std::srand() has been called
    int r = n;
    do  r = ( std::rand() / double(RAND_MAX) ) * n ; while( r == n ) ;
    return r;
}


/* Digamma function adapted from: http://web.science.mq.edu.au/~mjohnson/code/digamma.c
 * We include the comments from the original source code below:
 *
 * digamma.c
 *
 * Mark Johnson, 2nd September 2007
 *
 * Computes the Ψ(x) or digamma function, i.e., the derivative of the
 * log gamma function, using a series expansion.
 *
 * Warning:  I'm not a numerical analyst, so I may have made errors here!
 *
 * The parameters of the series were computed using the Maple symbolic
 * algebra program as follows:
 *
 * series(Psi(x+1/2), x=infinity, 21);
 *
 * which produces:
 *
 *  ln(x)+1/(24*x^2)-7/960/x^4+31/8064/x^6-127/30720/x^8+511/67584/x^10-1414477/67092480/x^12+8191/98304/x^14-118518239/267386880/x^16+5749691557/1882718208/x^18-91546277357/3460300800/x^20+O(1/(x^21))
 *
 * It looks as if the terms in this expansion *diverge* as the powers
 * get larger.  However, for large x, the x^-n term will dominate.
 *
 * I used Maple to examine the difference between this series and
 * Digamma(x+1/2) over the range 7 < x < 20, and discovered that the
 * difference is less that 1e-8 if the terms up to x^-8 are included.
 * This determined the power used in the code here.  Of course,
 * Maple uses some kind of approximation to calculate Digamma,
 * so all I've really done here is find the smallest power that produces
 * the Maple approximation; still, that should be good enough for our
 * purposes.
 *
 * This expansion is accurate for x > 7; we use the recurrence
 *
 * digamma(x) = digamma(x+1) - 1/x
 *
 * to make x larger than 7.
 */
template<class RealT>
RealT Psi(RealT x)
{
    RealT result = 0, xx, xx2, xx4;
    Assert(x > 0, "Psi function only works for positive x.");
    for ( ; x < 7; ++x)
        result -= 1/x;
    x -= 1.0/2.0;
    xx = 1.0/x;
    xx2 = xx*xx;
    xx4 = xx2*xx2;
    result += log(x)+(1./24.)*xx2-(7.0/960.0)*xx4+(31.0/8064.0)*xx4*xx2-(127.0/30720.0)*xx4*xx4;
    return result;
}

template<typename T> inline T Clip(const T x, const T lower, const T upper)
{
    return std::min(std::max(x, lower), upper);
}

template<typename T> T DotProduct(const std::vector<T> &x, const std::vector<T> &y)
{
    T ret = 0;
    for (std::size_t i = 0; i < x.size(); i++) ret += x[i] * y[i];
    return ret;
}

template<typename T>
T Abs(const T x)
{
    return x < 0 ? -x : x;
}

template<typename T>
T Sign(const T x)
{
    return x < 0 ? -1 : 0 < x ? 1 : 0;
}

template<typename T>
T Norm(const std::vector<T> &x)
{ 
    return Sqrt(DotProduct(x,x)); 
}

template<typename T>
std::vector<T> Sqrt(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Sqrt(ret[i]);
    return ret;
}

template<typename T>
std::vector<T> Exp(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Exp(ret[i]);
    return ret;
}

template<typename T>
std::vector<T> Log(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Log(ret[i]);
    return ret;
}

template<typename T>
std::vector<T> Pow(const std::vector<T> &x, const T p)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Pow(ret[i],p);
    return ret;
}

template<typename T>
std::vector<T> Tanh(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Tanh(ret[i]);
    return ret;
}

template<typename T>
std::vector<T> Abs(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Abs(ret[i]);
    return ret;
}

template<typename T>
std::vector<T> Sign(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Sign(ret[i]);
    return ret;
}

template<typename T, typename P>
std::vector<T> Test(const std::vector<T> &x, P pred)
{
    std::vector<T> ret(x.size());
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = (pred(x[i]) ? T(1) : T(0));
    return ret;
}

template<typename T>
T Sum(const std::vector<T> &x)
{
    T ret = 0;
    for (std::size_t i = 0; i < x.size(); i++) ret += x[i];
    return ret;
}

template<typename T>
T Prod(const std::vector<T> &x)
{
    T ret = 1;
    for (std::size_t i = 0; i < x.size(); i++) ret *= x[i];
    return ret;
}

template<typename T> 
const std::vector<T> Min(const std::vector<T> &x, const T &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = std::min(ret[i], y);
    return ret;
}

template<typename T> 
const std::vector<T> Max(const std::vector<T> &x, const T &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = std::max(ret[i], y);
    return ret;
}

template<typename T> 
const std::vector<T> Clip(const std::vector<T> &x, const T &lower, const T &upper)
{
    Assert(lower <= upper, "Invalid clipping range.");
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = Clip(x[i], lower, upper);
    return ret;
}

template<typename T> 
const std::vector<T> Min(const T &x, const std::vector<T> &y)
{
    std::vector<T> ret(y);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = std::min(x, ret[i]);
    return ret;
}

template<typename T> 
const std::vector<T> Max(const T &x, const std::vector<T> &y)
{
    std::vector<T> ret(y);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = std::max(x, ret[i]);
    return ret;
}

template<typename T>
T Min(const std::vector<T> &x)
{
    T ret = x[0];
    for (std::size_t i = 1; i < x.size(); i++) if (x[i] < ret) ret = x[i];
    return ret;
}

template<typename T>
T Max(const std::vector<T> &x)
{
    T ret = x[0];
    for (std::size_t i = 1; i < x.size(); i++) if (ret < x[i]) ret = x[i];
    return ret;
}

template<typename T>
int ArgMin(const std::vector<T> &x)
{
    int ret = 0;
    for (std::size_t i = 1; i < x.size(); i++) if (x[i] < x[ret]) ret = int(i);
    return ret;
}

template<typename T>
int ArgMax(const std::vector<T> &x)
{
    int ret = 0;
    for (std::size_t i = 1; i < x.size(); i++) if (x[ret] < x[i]) ret = int(i);
    return ret;
}

template<typename T> 
const std::vector<T> operator-(const std::vector<T> &x)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = -ret[i];
    return ret;
}

template<typename T> 
const std::vector<T> operator*(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] *= y[i];
    return ret;
}

template<typename T> 
const std::vector<T> operator/(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] /= y[i];
    return ret;
}

template<typename T> 
const std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] += y[i];
    return ret;
}

template<typename T> 
const std::vector<T> operator-(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] -= y[i];
    return ret;
}

template<typename T> 
const std::vector<T> operator*(const std::vector<T> &x, const T &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] *= y;
    return ret;
}

template<typename T> 
const std::vector<T> operator/(const std::vector<T> &x, const T &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] /= y;
    return ret;
}

template<typename T> 
const std::vector<T> operator+(const std::vector<T> &x, const T &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] += y;
    return ret;
}

template<typename T> 
const std::vector<T> operator-(const std::vector<T> &x, const T &y)
{
    std::vector<T> ret(x);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] -= y;
    return ret;
}

template<typename T> 
const std::vector<T> operator*(const T &x, const std::vector<T> &y)
{
    std::vector<T> ret(y);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] *= x;
    return ret;
}

template<typename T> 
const std::vector<T> operator/(const T &x, const std::vector<T> &y)
{
    std::vector<T> ret(y);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = x / ret[i];
    return ret;
}

template<typename T> 
const std::vector<T> operator+(const T &x, const std::vector<T> &y)
{
    std::vector<T> ret(y);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] += x;
    return ret;
}

template<typename T> 
const std::vector<T> operator-(const T &x, const std::vector<T> &y)
{
    std::vector<T> ret(y);
    for (std::size_t i = 0; i < ret.size(); i++) ret[i] = x - ret[i];
    return ret;
}

template<typename T> 
std::vector<T> &operator*=(std::vector<T> &x, const std::vector<T> &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] *= y[i];
    return x;
}

template<typename T> 
std::vector<T> &operator/=(std::vector<T> &x, const std::vector<T> &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] /= y[i];
    return x;
}

template<typename T> 
std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] += y[i];
    return x;
}

template<typename T> 
std::vector<T> &operator-=(std::vector<T> &x, const std::vector<T> &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] -= y[i];
    return x;
}

template<typename T> 
std::vector<T> &operator*=(std::vector<T> &x, const T &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] *= y;
    return x;
}

template<typename T> 
std::vector<T> &operator/=(std::vector<T> &x, const T &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] /= y;
    return x;
}

template<typename T> 
std::vector<T> &operator+=(std::vector<T> &x, const T &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] += y;
    return x;
}

template<typename T> 
std::vector<T> &operator-=(std::vector<T> &x, const T &y)
{
    for (std::size_t i = 0; i < x.size(); i++) x[i] -= y;
    return x;
}

template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &x)
{
    out << '[';
    for (std::size_t i = 0; i < x.size(); i++)
        out << (i > 0 ? " " : "") << x[i];
    out << ']';
    return out;
}

template<typename T, typename U>
std::vector<T> ConvertVector(const std::vector<U> &x)
{
    std::vector<T> ret;
    ret.reserve(x.size());
    for (size_t i = 0; i < x.size(); i++)
        ret.push_back(x[i]);
    return ret;
}

template<typename T>
std::vector<T> Concatenate(const std::vector<T> &u, const std::vector<T> &v)
{
    std::vector<T> ret = u;
    ret.insert(ret.end(), v.begin(), v.end());
    return ret;
}

template<typename T>
std::vector<T> Transpose(const std::vector<T> &m, const int rows, const int cols)
{
    Assert(rows * cols == int(m.size()), "Dimension mismatch.");
    
    std::vector<T> ret(m.size());
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            ret[j*rows+i] = m[i*cols+j];
    
    return ret;
}

//////////////////////////////////////////////////////////////////////
// ExpandMatrix()
//
// Expand matrix by adding blank rows and columns.
//////////////////////////////////////////////////////////////////////

template<class T>
std::vector<T> ExpandMatrix(const std::vector<T> &mat,
                            const int new_rows,
                            const int new_cols,
                            const std::vector<int> &positions_rows,
                            const std::vector<int> &positions_cols)
{
    Assert(new_rows >= 0, "Invalid length.");
    Assert(new_cols >= 0, "Invalid length.");
    Assert(positions_rows.size() > 0, "positions_rows should contain at least one element.");
    Assert(positions_cols.size() > 0, "positions_cols should contain at least one element.");
    Assert(positions_rows.size() * positions_cols.size() == mat.size(), "Dimension mismatch.");
    
    const int cols = int(positions_cols.size());
    std::vector<T> res(new_rows * new_cols);
    
    for (size_t i = 0; i < positions_rows.size(); i++)
    {
        for (size_t j = 0; j < positions_cols.size(); j++)
        {
            Assert(0 <= positions_rows[i] && positions_rows[i] < new_rows, "Index out-of-range.");
            Assert(0 <= positions_cols[j] && positions_cols[j] < new_cols, "Index out-of-range.");
            res[positions_rows[i] * new_cols + positions_cols[j]] = mat[i * cols + j];
        }
    }
    
    return res;  
}

//////////////////////////////////////////////////////////////////////
// ExpandVector()
//
// Expand vector by adding blank entries.
//////////////////////////////////////////////////////////////////////

template<class T>
std::vector<T> ExpandVector(const std::vector<T> &v,
                            const int new_length,
                            const std::vector<int> &positions)
{
    Assert(new_length > 0, "Invalid length.");
    Assert(positions.size() > 0, "positions should contain at least one element.");
    Assert(positions.size() == v.size(), "Dimension mismatch.");
    
    std::vector<T> res(new_length);
    
    for (size_t i = 0; i < positions.size(); i++)
    {
        Assert(0 <= positions[i] && positions[i] < new_length, "Index out-of-range.");
        res[positions[i]] = v[i];
    }
    
    return res;  
}

//////////////////////////////////////////////////////////////////////
// IsComplementary()
//
// Check if a base-pairing is one of AU, CG, or GU.
//////////////////////////////////////////////////////////////////////

inline bool IsComplementary(char c, char d)
{
    if ('a' <= c && c <= 'z') c += 'A' - 'a';
    if ('a' <= d && d <= 'z') d += 'A' - 'a';
    
    return 
        (c == 'A' && d == 'U') ||
        (c == 'U' && d == 'A') ||
        (c == 'C' && d == 'G') ||
        (c == 'G' && d == 'C') ||
        (c == 'G' && d == 'U') ||
        (c == 'U' && d == 'G');
}
