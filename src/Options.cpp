//////////////////////////////////////////////////////////////////////
// Options.cpp
//////////////////////////////////////////////////////////////////////

#include "Options.hpp"

//////////////////////////////////////////////////////////////////////
// Options::Options()
// Options::Options()
// Options::operator=()
//
// Constructor, copy constructor, and assignment operator
//////////////////////////////////////////////////////////////////////

Options::Options()
{}

Options::Options(const Options &rhs) :
    mapping(rhs.mapping)
{}

Options &Options::operator=(const Options &rhs)
{
    if (this != &rhs)
    {
        mapping = rhs.mapping;
    }
    return *this;
}


//////////////////////////////////////////////////////////////////////
// Options::GetBoolValue()
//
// Get boolean value.
//////////////////////////////////////////////////////////////////////

bool Options::GetBoolValue(const std::string &key) const
{
    std::map<std::string, std::string>::const_iterator iter = mapping.find(key);
    if (iter == mapping.end()) Error("Requested key '%s' not found in configuration.", key.c_str());
    if (iter->second != "true" && iter->second != "false") Error("Failed to parse boolean value in configuration.");
    return (iter->second == "true");
}

//////////////////////////////////////////////////////////////////////
// Options::SetBoolValue()
//
// Set boolean value.
//////////////////////////////////////////////////////////////////////

void Options::SetBoolValue(const std::string &key, bool value)
{
    mapping[key] = (value ? "true" : "false");
}

//////////////////////////////////////////////////////////////////////
// Options::GetIntValue()
//
// Get integer value.
//////////////////////////////////////////////////////////////////////

int Options::GetIntValue(const std::string &key) const
{
    int value = 0;
    std::map<std::string, std::string>::const_iterator iter = mapping.find(key);
    if (iter == mapping.end()) Error("Requested key '%s' not found in configuration.", key.c_str());
    if (!ConvertToNumber(iter->second, value)) Error("Failed to parse integer value in configuration.");
    return value;
}

//////////////////////////////////////////////////////////////////////
// Options::SetIntValue()
//
// Set integer value.
//////////////////////////////////////////////////////////////////////

void Options::SetIntValue(const std::string &key, int value)
{
    mapping[key] = SPrintF("%d", value);
}

//////////////////////////////////////////////////////////////////////
// Options::GetRealValue()
//
// Get real value.
//////////////////////////////////////////////////////////////////////

double Options::GetRealValue(const std::string &key) const
{
    double value = 0;
    std::map<std::string, std::string>::const_iterator iter = mapping.find(key);
    if (iter == mapping.end()) Error("Requested key '%s' not found in configuration.", key.c_str());
    if (!ConvertToNumber(iter->second, value)) Error("Failed to parse real value in configuration.");
    return value;
}

//////////////////////////////////////////////////////////////////////
// Options::SetRealValue()
//
// Set real value.
//////////////////////////////////////////////////////////////////////

void Options::SetRealValue(const std::string &key, double value)
{
    mapping[key] = SPrintF("%lf", value);
}

//////////////////////////////////////////////////////////////////////
// Options::GetStringValue()
//
// Get string value.
//////////////////////////////////////////////////////////////////////

const std::string &Options::GetStringValue(const std::string &key) const
{
    std::map<std::string, std::string>::const_iterator iter = mapping.find(key);
    if (iter == mapping.end()) Error("Requested key '%s' not found in configuration.", key.c_str());
    return iter->second;
}

//////////////////////////////////////////////////////////////////////
// Options::SetStringValue()
//
// Set string value.
//////////////////////////////////////////////////////////////////////

void Options::SetStringValue(const std::string &key, const std::string &value)
{
    mapping[key] = value;
}
