//////////////////////////////////////////////////////////////////////
// ParameterManager.ipp
//////////////////////////////////////////////////////////////////////

#include "ParameterManager.hpp"

//////////////////////////////////////////////////////////////////////
// ParameterGroup::ParameterGroup()
// ParameterGroup::operator=()
// 
// Constructors and assignment operator.
//////////////////////////////////////////////////////////////////////

ParameterGroup::ParameterGroup() {}

ParameterGroup::ParameterGroup(const std::string &name, int begin, int end) :
    name(name),
    begin(begin),
    end(end)
{
    Assert(begin <= end, "Inconsistent begin and end indices.");
}

ParameterGroup::ParameterGroup(const ParameterGroup &rhs) :
    name(rhs.name),
    begin(rhs.begin),
    end(rhs.end)
{
    Assert(begin <= end, "Inconsistent begin and end indices.");
}

ParameterGroup &ParameterGroup::operator=(const ParameterGroup &rhs)
{
    if (this != &rhs)
    {
        name = rhs.name;
        begin = rhs.begin;
        end = rhs.end;
        Assert(begin <= end, "Inconsistent begin and end indices.");
    }
    return *this;
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::ParameterManager()
//
// Constructor.
//////////////////////////////////////////////////////////////////////


ParameterManager::ParameterManager() :
    names(),
    groups(),
    physical_to_logical(),
    logical_to_physical(),
    logical_name_to_index()
{}

//////////////////////////////////////////////////////////////////////
// ParameterManager::ParameterManager()
//
// Copy constructor.
//////////////////////////////////////////////////////////////////////


ParameterManager::ParameterManager(const ParameterManager &rhs) :
    names(rhs.names),
    groups(rhs.groups),
    physical_to_logical(rhs.physical_to_logical),
    logical_to_physical(rhs.logical_to_physical),
    logical_name_to_index(rhs.logical_name_to_index)
{}

//////////////////////////////////////////////////////////////////////
// ParameterManager::operator=()
//
// Assignment operator.
//////////////////////////////////////////////////////////////////////


ParameterManager &ParameterManager::operator=(const ParameterManager &rhs)
{
    if (this != &rhs)
    {
        names = rhs.names;
        groups = rhs.groups;
        physical_to_logical = rhs.physical_to_logical;
        logical_to_physical = rhs.logical_to_physical;
        logical_name_to_index = rhs.logical_name_to_index;
    }
    return *this;
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::~ParameterManager()
//
// Destructor.
//////////////////////////////////////////////////////////////////////


ParameterManager::~ParameterManager()
{}

//////////////////////////////////////////////////////////////////////
// ParameterManager::ClearParameters()
//
// Clear parameter manager.
//////////////////////////////////////////////////////////////////////


void ParameterManager::ClearParameters()
{
    names.clear();
    groups.clear();
    physical_to_logical.clear();
    logical_to_physical.clear();
    logical_name_to_index.clear();
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::AddParameterGroup()
//
// Mark the beginning of a new parameter group.
//////////////////////////////////////////////////////////////////////


void ParameterManager::AddParameterGroup(const std::string &name)
{
    groups.push_back(ParameterGroup(name, int(names.size()), int(names.size())));
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::AddParameterMapping()
//
// Create a mapping from a physical parameter to a logical parameter.
//////////////////////////////////////////////////////////////////////


void ParameterManager::AddParameterMapping(const std::string &logical_name, std::pair<RealT, RealT> *physical_ptr)
{
    // check if the logical parameter name has been seen before
    std::map<std::string,int>::iterator iter = logical_name_to_index.find(logical_name);
    if (iter == logical_name_to_index.end())
    {
        // if not, add it
        iter = logical_name_to_index.insert(std::make_pair(logical_name, int(names.size()))).first;
        names.push_back(logical_name);
        logical_to_physical.push_back(std::vector<std::pair<RealT, RealT> *>());
        ++(groups.back().end);
    }

    // save mapping from physical parameter pointer to logical index, and vice versa
    physical_to_logical[physical_ptr] = iter->second;
    logical_to_physical[iter->second].push_back(physical_ptr);

}

//////////////////////////////////////////////////////////////////////
// ParameterManager::ReadFromFile()
//
// Read parameters from file.
//////////////////////////////////////////////////////////////////////


void ParameterManager::ReadFromFile(const std::string &filename, std::vector<RealT> &values)
{
    std::map<std::string, RealT> params;
    RealT value;
    std::string name;
    std::string s;
    
    // read parameter file
    std::ifstream infile(filename.c_str());
    if (infile.fail()) Error(("Could not open file \"" + filename + "\" for reading.").c_str());
    
    while (getline(infile, s))
    {
        // skip blank lines and comments
        if (s.length() == 0 || s[0] == '#') continue;

        // read parameter names and values
        std::istringstream iss(s);
        if (iss >> name >> value)
        {
            if (params.find(name) != params.end())
                Error("Parameter file contains a duplicate parameter: %s", name.c_str());
            params[name] = value;
        }
    }
    infile.close();

    // convert read parameters to vector format
    values.clear();
    values.resize(names.size());
    for (size_t i = 0; i < names.size(); i++)
    {
        typename std::map<std::string, RealT>::iterator iter = params.find(names[i]);
        if (iter == params.end()) Error("Parameter file missing parameter: %s", names[i].c_str());
        values[i] = iter->second;
        params.erase(iter);
    }

    // print an error message for extra parameters
    for (typename std::map<std::string, RealT>::iterator iter = params.begin(); iter != params.end(); ++iter)
        Warning("Parameter file contains extra parameter: %s", iter->first.c_str());
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::WriteToFile()
//
// Write parameters to file.
//////////////////////////////////////////////////////////////////////


void ParameterManager::WriteToFile(const std::string &filename, const std::vector<RealT> &values)
{
    if (values.size() != names.size()) Error("Incorrect number of parameters.");
    std::ofstream outfile(filename.c_str());
    if (outfile.fail()) Error(("Could not open file \"" + filename + "\" for writing.").c_str());
    for (size_t i = 0; i < values.size(); i++)
        outfile << names[i] << " " << std::setprecision(10) << values[i] << std::endl;
    outfile.close();
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::ExpandParameterGroupValues()
//
// Expand a vector of group values.
//////////////////////////////////////////////////////////////////////


const std::vector<RealT> ParameterManager::ExpandParameterGroupValues(const std::vector<RealT> &values) const
{
    std::vector<RealT> expanded;
    std::cout << "here" << values.size() << " " << groups.size() << std::endl;
    if (values.size() != groups.size()) Error("Incorrect number of hyperparametrs.");
    for (size_t i = 0; i < groups.size(); i++)
        for (int j = groups[i].begin; j < groups[i].end; j++)
            expanded.push_back(values[i]);
    return expanded;
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::GetPhysicalParameters()
//
// Retrieve the vector of physical parameters associated with a
// particular logical parameter.
//////////////////////////////////////////////////////////////////////


std::vector<std::pair<RealT, RealT> *> ParameterManager::GetPhysicalParameters(int logical_index) const
{
    if (logical_index < 0 || logical_index >= int(names.size()))
        Error("Requested for invalid logical parameter index: %d", logical_index);
    return logical_to_physical[logical_index];
}

//////////////////////////////////////////////////////////////////////
// ParameterManager::GetLogicalIndex()
//
// Retrieve the logical index for a particular physical parameter.
//////////////////////////////////////////////////////////////////////


int ParameterManager::GetLogicalIndex(std::pair<RealT, RealT> *physical_ptr) const
{
    typename std::map<std::pair<RealT, RealT> *, int>::const_iterator iter = physical_to_logical.find(physical_ptr);
    if (iter == physical_to_logical.end()) Error("Request for unknown physical parameter.");
    return iter->second;
}
