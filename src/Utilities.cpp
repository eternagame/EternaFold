//////////////////////////////////////////////////////////////////////
// Utilities.cpp
//////////////////////////////////////////////////////////////////////

#include "Utilities.hpp"
#include "Utilities.ipp"
#include <unistd.h>

bool toggle_error = false;

//////////////////////////////////////////////////////////////////////
// _ASSERT_FAILED()
// 
// Print error message for a failed assertion and terminate.
//////////////////////////////////////////////////////////////////////

int _ASSERT_FAILED(const char *fmt, ...)
{
    if (toggle_error) return 0;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    abort();
    return 0;
}

//////////////////////////////////////////////////////////////////////
// Error()
// 
// Print error message for a user error and terminate.
//////////////////////////////////////////////////////////////////////

void Error(const char *fmt, ...)
{
    fprintf(stderr, "ERROR: ");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    toggle_error = true;
    exit(1);
}

//////////////////////////////////////////////////////////////////////
// Warning()
// 
// Print warning message for a user error without terminating.
//////////////////////////////////////////////////////////////////////

void Warning(const char *fmt, ...)
{
    fprintf(stderr, "WARNING: ");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
}

//////////////////////////////////////////////////////////////////////
// ConvertToNumber()
//
// Attempts to parse a number from the character string given.
// Returns true only if no parsing error occurs.  
//////////////////////////////////////////////////////////////////////

bool ConvertToNumber(const std::string &s, int &val)
{
    char *end_ptr;
    long int temp;
    errno = 0;
    temp = strtol(s.c_str(), &end_ptr, 10);
    if (errno == ERANGE && (temp == LONG_MIN || temp == LONG_MAX)) return false;
    if (temp == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    if (temp < INT_MIN || temp > INT_MAX) return false;
    val = int(temp);
    return true;
}

bool ConvertToNumber(const std::string &s, unsigned int &val)
{
    char *end_ptr;
    unsigned long int temp;
    errno = 0;
    temp = strtol(s.c_str(), &end_ptr, 10);
    if (errno == ERANGE && temp == ULONG_MAX) return false;
    if (temp == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    if (temp > UINT_MAX) return false;
    val = int(temp);
    return true;
}

bool ConvertToNumber(const std::string &s, long int &val)
{
    char *end_ptr;
    errno = 0;
    val = strtol(s.c_str(), &end_ptr, 10);
    if (errno == ERANGE && (val == LONG_MIN || val == LONG_MAX)) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}

bool ConvertToNumber(const std::string &s, unsigned long int &val)
{
    char *end_ptr;
    errno = 0;
    val = strtoul(s.c_str(), &end_ptr, 10);
    if (errno == ERANGE && val == ULONG_MAX) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}

#if defined(LLONG_MIN) && defined(LLONG_MAX)

bool ConvertToNumber(const std::string &s, long long int &val)
{
    char *end_ptr;
    errno = 0;
    val = strtoll(s.c_str(), &end_ptr, 10);
    if (errno == ERANGE && (val == LLONG_MIN || val == LLONG_MAX)) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}

#endif

#if defined(ULLONG_MAX)

bool ConvertToNumber(const std::string &s, unsigned long long int &val)
{
    char *end_ptr;
    errno = 0;
    val = strtoull(s.c_str(), &end_ptr, 10);
    if (errno == ERANGE && val == ULLONG_MAX) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}


#endif

bool ConvertToNumber(const std::string &s, float &val)
{
    char *end_ptr;
    errno = 0;
    val = strtof(s.c_str(), &end_ptr);
    if (errno == ERANGE && (val == HUGE_VALF || val == HUGE_VALF)) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}

bool ConvertToNumber(const std::string &s, double &val)
{
    char *end_ptr;
    errno = 0;
    val = strtod(s.c_str(), &end_ptr);
    if (errno == ERANGE && (val == HUGE_VAL || val == -HUGE_VAL)) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}

bool ConvertToNumber(const std::string &s, long double &val)
{
    char *end_ptr;
    errno = 0;
#ifdef strtold
    val = strtold(s.c_str(), &end_ptr);
#else
    val = double(strtod(s.c_str(), &end_ptr));
#endif
    if (errno == ERANGE && (val == HUGE_VALL || val == -HUGE_VALL)) return false;
    if (val == 0 && errno != 0) return false;
    if (s.c_str() == end_ptr) return false;
    return true;
}

//////////////////////////////////////////////////////////////////////
// ConvertToUpperCase()
//
// Converts lowercase letters in a string to uppercase.
//////////////////////////////////////////////////////////////////////

std::string ConvertToUpperCase(const std::string &s)
{
    std::string t(s);
    for (size_t i = 0; i < t.length(); i++) 
        t[i] = toupper(t[i]);
    return t;
}

//////////////////////////////////////////////////////////////////////
// ConvertToLowerCase()
//
// Converts uppercase letters in a string to lowercase.
//////////////////////////////////////////////////////////////////////

std::string ConvertToLowerCase(const std::string &s)
{
    std::string t(s);
    for (size_t i = 0; i < t.length(); i++) 
        t[i] = tolower(t[i]);
    return t;
}

////////////////////////////////////////////////////////////
// Trim()
//
// Remove whitespace from either end of a string.
////////////////////////////////////////////////////////////

std::string Trim(const std::string &s){
    int left = 0, right = int(s.length());
    while (left < right && std::isspace(s[left])) left++;
    while (left < right && std::isspace(s[right-1])) right--;
    return s.substr(left,right-left);
}

//////////////////////////////////////////////////////////////////////
// RemoveGaps()
// 
// Remove gap characters from a string.
//////////////////////////////////////////////////////////////////////

std::string RemoveGaps(const std::string &s)
{
    std::string ret;

    for (size_t i = 0; i < s.length(); i++)
    {
        if (s[i] != '.' && s[i] != '-')
            ret.push_back(s[i]);
    }

    return ret;
}

//////////////////////////////////////////////////////////////////////
// SPrintF()
//
// Same as sprintf but returns a string.
//////////////////////////////////////////////////////////////////////

std::string SPrintF(const char *fmt, ...){
    int buf_size = 1024;
    char *buffer = new char[buf_size];
    Assert(buffer, "Failed to allocate memory.");
    
    while (true)
    {
        // print text to buffer
        
        va_list ap;
        va_start(ap, fmt);
        int num_written = vsnprintf(buffer, buf_size, fmt, ap);
        va_end(ap);
        
        // double size of buffer if needed
        
        if (num_written >= buf_size)
        {
            char *temp = new char[buf_size*2];
            Assert(temp, "Failed to allocate memory.");
            memcpy(temp, buffer, sizeof(char) * buf_size);
            delete [] buffer;
            buffer = temp;
            buf_size *= 2;      
        }
        else
        {
            break;
        }
    }
    
    // return text
    
    const std::string s(buffer);
    delete [] buffer;
    return s;
}

//////////////////////////////////////////////////////////////////////
// WriteProgressMessage()
//
// Write progress message to console (stderr) and return to 
// beginning of line.  Wipes out any previous message.
//////////////////////////////////////////////////////////////////////

void WriteProgressMessage(const std::string &message)
{
    static int old_length = 0;
    std::cerr << '\r' << message;
    for (int i = message.length(); i < old_length; i++) std::cerr << ' ';
    std::cerr << '\r';
    old_length = int(message.length());
}

//////////////////////////////////////////////////////////////////////
// GetSystemTime()
//
// Retrieve system time in seconds past the Epoch.
//////////////////////////////////////////////////////////////////////

double GetSystemTime()
{
    timeval t;
    if (gettimeofday(&t, NULL) != 0) Error("Failed to obtain system time.");
    return t.tv_sec + 1e-6 * t.tv_usec;
}

//////////////////////////////////////////////////////////////////////
// MakeDirectory()
//
// Make a directory if one doesn't exist.
//////////////////////////////////////////////////////////////////////

void MakeDirectory(const std::string &directory)
{
    if (directory != "" && system(("mkdir -p " + directory).c_str()))
        Error(SPrintF("Unable to create directory \"%s\"", directory.c_str()).c_str());
}

//////////////////////////////////////////////////////////////////////
// MakeTempDirectory()
//
// Make a temporary directory that will automatically be deleted
// once the program is complete.
//////////////////////////////////////////////////////////////////////

std::string MakeTempDirectory()
{
    char *temp_dir_name = new char[10000];
    Assert(temp_dir_name, "Failed to allocate memory.");
    strcpy(temp_dir_name, "temp_XXXXXX");
    char *ret = mkdtemp(temp_dir_name);

    // check for error

    if (ret == NULL)
    {
        delete [] temp_dir_name;
        Error("Unable to create temp directory!");
        return "";
    }

    // return name of temporary directory

    std::string res(temp_dir_name);
    delete [] temp_dir_name;
    return res;
}

//////////////////////////////////////////////////////////////////////
// GetSequencePositions()
//
// Return an array whose ith element is the index of the ith
// letter in the input string.  The 0th element is always zero,
// since there is no 0-th character in the input string.
//////////////////////////////////////////////////////////////////////

std::vector<int> GetSequencePositions(const std::string &s)
{
    std::vector<int> ret(1);
    for (size_t i = 0; i < s.length(); i++)
        if (isalpha(s[i])) ret.push_back(i);
    return ret;
}

//////////////////////////////////////////////////////////////////////
// GetSequenceMapping()
//
// Return an array mapping from positions in a gapped sequence to
// indices of letteres in an ungapped sequence.  Positions in
// the array corresponding to gaps are given a mapping of 0.
//////////////////////////////////////////////////////////////////////

std::vector<int> GetSequenceMapping(const std::string &s)
{
    std::vector<int> ret(1);
    int ct = 0;
    for (size_t i = 1; i < s.length(); i++)
    {
        if (isalpha(s[i]))
        {
            ++ct;
            ret.push_back(ct);
        }
        else
            ret.push_back(0);
    }
    return ret;
}

//////////////////////////////////////////////////////////////////////
// GetDirName()
// GetBaseName()
//
// Retrieve directory and base file name for a given full path.
//////////////////////////////////////////////////////////////////////

std::string GetDirName(const std::string &filename)
{
    const std::string::size_type separator_pos = filename.find_last_of(DIR_SEPARATOR_CHAR);
    std::string dir_name = ((separator_pos == std::string::npos) ? std::string("") : filename.substr(0, separator_pos));
    while (dir_name.length() > 0 && dir_name[dir_name.length() - 1] == DIR_SEPARATOR_CHAR) dir_name = dir_name.substr(0, dir_name.length() - 1);
    return dir_name;
}

std::string GetBaseName(const std::string &filename)
{
    const std::string::size_type separator_pos = filename.find_last_of(DIR_SEPARATOR_CHAR);
    return ((separator_pos == std::string::npos) ? filename : filename.substr(separator_pos + 1));
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