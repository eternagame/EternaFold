//////////////////////////////////////////////////////////////////////
// Options.hpp
//////////////////////////////////////////////////////////////////////

#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// class Options
//////////////////////////////////////////////////////////////////////

class Options
{
    std::map<std::string, std::string> mapping;
    
public:

    // constructor, copy constructor, and assignment operator
    Options();
    Options(const Options &rhs);
    Options &operator=(const Options &rhs);

    // getters and setters
    bool GetBoolValue(const std::string &key) const;
    void SetBoolValue(const std::string &key, bool value);
    
    int GetIntValue(const std::string &key) const;
    void SetIntValue(const std::string &key, int value);
    
    double GetRealValue(const std::string &key) const;
    void SetRealValue(const std::string &key, double value);
    
    const std::string &GetStringValue(const std::string &key) const;
    void SetStringValue(const std::string &key, const std::string &value);
};

#endif
