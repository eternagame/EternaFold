//////////////////////////////////////////////////////////////////////
// Utilities.hpp
//////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <sys/time.h>
#include <vector>
#include <unistd.h>

typedef unsigned char BYTE;
const char DIR_SEPARATOR_CHAR = '/';

// necessary macros
#define _level_2_str(s) _level_1_str(s)
#define _level_1_str(s) #s
#define __LINESTR__ _level_2_str(__LINE__)
#define __FL__ "Assertion failed in file \"" __FILE__ "\", line " __LINESTR__
int _ASSERT_FAILED(const char *fmt, ...);

// print error message for a user error and terminate
void Error(const char *fmt, ...);

// print error message for a user error but do not terminate
void Warning(const char *fmt, ...);

// print error message for a failed assertion and terminate
#ifdef NDEBUG
#define Assert(test,fmt,...)
#else
#define Assert(test,fmt,...) (test ? 0 : _ASSERT_FAILED(__FL__ ": " fmt "\n", ## __VA_ARGS__))
#endif

// attempt to parse a number from the character string given; return
// true only if no parsing error occurs.
bool ConvertToNumber(const std::string &s, int &val);
bool ConvertToNumber(const std::string &s, unsigned int &val);
bool ConvertToNumber(const std::string &s, long int &val);
bool ConvertToNumber(const std::string &s, unsigned long int &val);
#if defined(LLONG_MIN) && defined(LLONG_MAX)
bool ConvertToNumber(const std::string &s, long long int &val);
#endif
#if defined(ULLONG_MAX)
bool ConvertToNumber(const std::string &s, unsigned long long int &val);
#endif
bool ConvertToNumber(const std::string &s, float &val);
bool ConvertToNumber(const std::string &s, double &val);
bool ConvertToNumber(const std::string &s, long double &val);

// convert lowercase/uppercase letters in a string to uppercase/lowercase
std::string ConvertToUpperCase(const std::string &s);
std::string ConvertToLowerCase(const std::string &s);

// remove whitespace from either end of a string
std::string Trim(const std::string &s);

// remove gap characters from a string
std::string RemoveGaps(const std::string &s);

// same as sprintf but returns a string
std::string SPrintF(const char *fmt, ...);

// write progress message to console (stderr) and return to 
// beginning of line; wipes out any previous message on current line
void WriteProgressMessage(const std::string &message);

// retrieve system time in seconds past the Epoch
double GetSystemTime();

// make a directory if one doesn't exist
void MakeDirectory(const std::string &directory);

// make temporary directory
std::string MakeTempDirectory();

// return an array whose ith element is the index of the ith
// letter in the input string.
std::vector<int> GetSequencePositions(const std::string &s);

// return an array from positions in a gapped sequence to
// positions in the ungapped sequence
std::vector<int> GetSequenceMapping(const std::string &s);

// indicator function
inline int Ind(bool condition){ return condition ? 1 : 0; }

// struct triple
template<typename T1, typename T2, typename T3>
struct triple {
    T1 first;
    T2 second;
    T3 third;

    // constructors
    triple();
    triple(const T1 &first, const T2 &second, const T3 &third);
    triple(const triple &rhs);
};

// comparators
template<typename T1, typename T2, typename T3> inline bool operator==(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y);
template<typename T1, typename T2, typename T3> inline bool operator<(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y);
template<typename T1, typename T2, typename T3> inline bool operator!=(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y);
template<typename T1, typename T2, typename T3> inline bool operator>(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y);
template<typename T1, typename T2, typename T3> inline bool operator<=(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y);
template<typename T1, typename T2, typename T3> inline bool operator>=(const triple<T1,T2,T3> &x, const triple<T1,T2,T3> &y);

// utility function for making triples
template<typename T1, typename T2, typename T3> inline triple<T1,T2,T3> make_triple(T1 first, T2 second, T3 third);

// printing pairs and triples
template<typename T1, typename T2> std::ostream &operator<<(std::ostream &out, const std::pair<T1,T2> &x);
template<typename T1, typename T2, typename T3> std::ostream &operator<<(std::ostream &out, const triple<T1,T2,T3> &x);

// math operators
template<typename T> T Sqrt(const T x);
template<typename T> T Exp(const T x);
template<typename T> T Log(const T x);
template<typename T> T Pow(const T x, const T p);
template<typename T> T Tanh(const T x);
template<typename T> T Floor(const T x);
template<typename T> T Ceil(const T x);
template<typename T> T Abs(const T x);
template<typename T> T Sign(const T x);
template<typename T> T Clip(const T x, const T lower, const T upper);

// standard linear algebra
template<typename T> T DotProduct(const std::vector<T> &x, const std::vector<T> &y);
template<typename T> T Norm(const std::vector<T> &x);
template<typename T> std::vector<T> Sqrt(const std::vector<T> &x);
template<typename T> std::vector<T> Exp(const std::vector<T> &x);
template<typename T> std::vector<T> Log(const std::vector<T> &x);
template<typename T> std::vector<T> Pow(const std::vector<T> &x, const T p);
template<typename T> std::vector<T> Tanh(const std::vector<T> &x);
template<typename T> std::vector<T> Abs(const std::vector<T> &x);
template<typename T> std::vector<T> Sign(const std::vector<T> &x);
template<typename T, typename P> std::vector<T> Test(const std::vector<T> &x, P pred);
template<typename T> T Sum(const std::vector<T> &x);
template<typename T> T Prod(const std::vector<T> &x);
template<typename T> const std::vector<T> Min(const std::vector<T> &x, const T &y);
template<typename T> const std::vector<T> Max(const std::vector<T> &x, const T &y);
template<typename T> const std::vector<T> Clip(const std::vector<T> &x, const T &lower, const T &upper);
template<typename T> const std::vector<T> Min(const T &x, const std::vector<T> &y);
template<typename T> const std::vector<T> Max(const T &x, const std::vector<T> &y);
template<typename T> T Min(const std::vector<T> &x);
template<typename T> T Max(const std::vector<T> &x);
template<typename T> int ArgMin(const std::vector<T> &x);
template<typename T> int ArgMax(const std::vector<T> &x);
template<typename T> const std::vector<T> operator-(const std::vector<T> &x);
template<typename T> const std::vector<T> operator*(const std::vector<T> &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator/(const std::vector<T> &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator-(const std::vector<T> &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator*(const std::vector<T> &x, const T &y);
template<typename T> const std::vector<T> operator/(const std::vector<T> &x, const T &y);
template<typename T> const std::vector<T> operator+(const std::vector<T> &x, const T &y);
template<typename T> const std::vector<T> operator-(const std::vector<T> &x, const T &y);
template<typename T> const std::vector<T> operator*(const T &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator/(const T &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator+(const T &x, const std::vector<T> &y);
template<typename T> const std::vector<T> operator-(const T &x, const std::vector<T> &y);
template<typename T> std::vector<T> &operator*=(std::vector<T> &x, const std::vector<T> &y);
template<typename T> std::vector<T> &operator/=(std::vector<T> &x, const std::vector<T> &y);
template<typename T> std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y);
template<typename T> std::vector<T> &operator-=(std::vector<T> &x, const std::vector<T> &y);
template<typename T> std::vector<T> &operator*=(std::vector<T> &x, const T &y);
template<typename T> std::vector<T> &operator/=(std::vector<T> &x, const T &y);
template<typename T> std::vector<T> &operator+=(std::vector<T> &x, const T &y);
template<typename T> std::vector<T> &operator-=(std::vector<T> &x, const T &y);
template<typename T> std::ostream &operator<<(std::ostream &out, const std::vector<T> &x);
template<typename T, typename U> std::vector<T> ConvertVector(const std::vector<U> &x);
template<typename T> std::vector<T> Concatenate(const std::vector<T> &u, const std::vector<T> &v);
template<typename T> std::vector<T> Transpose(const std::vector<T> &m, const int rows, const int cols);

// expand matrix by adding blank rows and columns
template<class T>
std::vector<T> ExpandMatrix(const std::vector<T> &mat,
                            const int new_rows,
                            const int new_cols,
                            const std::vector<int> &positions_rows,
                            const std::vector<int> &positions_cols);

// expand vector by adding blank entries
template<class T>
std::vector<T> ExpandVector(const std::vector<T> &v,
                            const int new_length,
                            const std::vector<int> &positions);

// check if two nucleotides are complementary (AU, CG, GU)
inline bool IsComplementary(char c, char d);

// retrieve directory and basename for a given full path
std::string GetDirName(const std::string &path);
std::string GetBaseName(const std::string &path);


#include "Utilities.ipp"

#endif
