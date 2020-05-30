//////////////////////////////////////////////////////////////////////
// ParameterManager.hpp
//
// This class is used for parameter management for models in which
// 
//    (1) the number of logical (free) parameters in the model is not
//        the same as the number of physical parameters stored in
//        memory
//
//    (2) certain sets of logical parameters are bound together in
//        parameter groups.
//
// To accomplish this, the ParameterManager class provides a set of
// routines for constructing a mapping between logical parameters and
// physical parameters:
//
//    AddParameterGroup(group_name)
//    AddParameterMapping(logical_name, physical_ptr)
//
// -------------------------------------------------------------------
//
// The ParameterManager class provides the following functionality:
//
//    GetPhysicalParameters(logical_index)
//
//       return a vector containing pointers to all physical
//       parameters which map to a particular logical_index (or all
//       physical parameters if logical_index == -1)
//
//    GetLogicalIndex(physical_ptr)
//
//       return the index of the logical parameter corresponding to
//       physical_ptr
//
//    GetNumPhysicalParameters()
//
//       return the number of physical parameters
//       
//    GetNumLogicalParameters()
//
//       return the number of logical parameters
//       
//////////////////////////////////////////////////////////////////////

#ifndef PARAMETERMANAGER_HPP
#define PARAMETERMANAGER_HPP

#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// struct ParameterGroup
//
// Contains the indices of all parameters corresponding to a
// particular parameter group.  In particular, "begin" is the
// index of the first parameter belonging to the parameter group,
// and "end" is the index of the element *after* the last parameter
// belonging to the group.
//////////////////////////////////////////////////////////////////////

struct ParameterGroup
{
    std::string name;
    int begin, end;
    
    ParameterGroup();
    ParameterGroup(const std::string &name, int begin, int end);
    ParameterGroup(const ParameterGroup &rhs);  
    ParameterGroup &operator=(const ParameterGroup &rhs);
};

//////////////////////////////////////////////////////////////////////
// class ParameterManager
//////////////////////////////////////////////////////////////////////

template<class RealT>
class ParameterManager
{
    std::vector<std::string> names;
    std::vector<ParameterGroup> groups;
    std::map<std::pair<RealT, RealT> *, int> physical_to_logical;
    std::vector<std::vector<std::pair<RealT, RealT> *> > logical_to_physical;
    std::map<std::string, int> logical_name_to_index;
    
public:

    // constructors, assignment operator, and destructor
    ParameterManager();
    ParameterManager(const ParameterManager &rhs);
    ParameterManager &operator= (const ParameterManager &rhs);
    virtual ~ParameterManager();

    // routines for adding new parameters and parameter groups
    void ClearParameters();
    void AddParameterGroup(const std::string &name);
    void AddParameterMapping(const std::string &logical_name, std::pair<RealT, RealT> *physical_ptr);

    // file input/output
    void ReadFromFile(const std::string &filename, std::vector<RealT> &values);
    void WriteToFile(const std::string &filename, const std::vector<RealT> &values);
    
    // expand a vector of values for each parameter group
    const std::vector<RealT> ExpandParameterGroupValues(const std::vector<RealT> &values) const;

    // retrieve physical parameters corresponding to a particular logical index
    std::vector<std::pair<RealT, RealT> *> GetPhysicalParameters(int logical_index) const;

    // return logical index corresponding to a physical parameter
    int GetLogicalIndex(std::pair<RealT, RealT> *physical_ptr) const;

    // simple getters
    const std::vector<std::string> GetNames() const { return names; }
    const std::vector<ParameterGroup> &GetParameterGroups() const { return groups; }
    size_t GetNumParameterGroups() const { return groups.size(); }
    size_t GetNumPhysicalParameters() const { return physical_to_logical.size(); }
    size_t GetNumLogicalParameters() const { return logical_to_physical.size(); }
};

#include "ParameterManager.ipp"

#endif
