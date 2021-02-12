//////////////////////////////////////////////////////////////////////
// FileDescription.hpp
//
// Contains a description of a file for processing.
//////////////////////////////////////////////////////////////////////

#ifndef FILEDESCRIPTION_HPP
#define FILEDESCRIPTION_HPP

#include <string>
#include "Config.hpp"
#include "SStruct.hpp"

//////////////////////////////////////////////////////////////////////
// struct FileDescription
//////////////////////////////////////////////////////////////////////

struct FileDescription
{
    SStruct sstruct;
    std::string input_filename;
    int size;
    double weight;
    
    // constructors, assignment operator, destructor
    FileDescription(const std::string &input_filename="",
                    const bool allow_noncomplementary=false,
                    const int num_data_sources=0);
    FileDescription(const FileDescription &rhs);
    FileDescription &operator=(const FileDescription &rhs);
    ~FileDescription();

    // comparator for sorting by decreasing size
    bool operator<(const FileDescription &rhs) const;
};

#endif
