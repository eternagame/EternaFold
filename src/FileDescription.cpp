//////////////////////////////////////////////////////////////////////
// FileDescription.cpp
//////////////////////////////////////////////////////////////////////

#include "FileDescription.hpp"
#include <Utilities.ipp>

//////////////////////////////////////////////////////////////////////
// FileDescription::~FileDescription()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

FileDescription::~FileDescription()
{}

//////////////////////////////////////////////////////////////////////
// FileDescription::FileDescription()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

FileDescription::FileDescription(const std::string &input_filename,
                                 const bool allow_noncomplementary,
                                 const int num_data_sources) :
    sstruct(input_filename,num_data_sources),
    input_filename(input_filename),
    size(int(Pow(double(sstruct.GetLength()), 3.0))),
    weight(1.0)
{
    if (!allow_noncomplementary)
        sstruct.RemoveNoncomplementaryPairings();
#if PROFILE
    if (sstruct.GetNumSequences() == 1)
        Warning("Using --profile mode with only one input sequence.");
#else
    if (sstruct.GetNumSequences() > 1)
        Warning("Folding multiple input sequences without --profile mode enabled.");
#endif
 //   std::cout << "kd data in FileDesc " << sstruct.GetKdData() << std::endl;
}

//////////////////////////////////////////////////////////////////////
// FileDescription::FileDescription()
//
// Copy constructor.
//////////////////////////////////////////////////////////////////////

FileDescription::FileDescription(const FileDescription &rhs) :
    sstruct(rhs.sstruct),
    input_filename(rhs.input_filename), 
    size(rhs.size),
    weight(rhs.weight)
{}

//////////////////////////////////////////////////////////////////////
// FileDescription::operator=()
//
// Assignment operator.
//////////////////////////////////////////////////////////////////////

FileDescription &FileDescription::operator=(const FileDescription &rhs)
{
    if (this != &rhs)
    {
        sstruct = rhs.sstruct;
        input_filename = rhs.input_filename;
        size = rhs.size;
        weight = rhs.weight;
    }
    return *this;
}

//////////////////////////////////////////////////////////////////////
// FileDescription::operator<()
//
// Comparator used to sort by decreasing size.
//////////////////////////////////////////////////////////////////////

bool FileDescription::operator<(const FileDescription &rhs) const 
{
    return size > rhs.size;
}

