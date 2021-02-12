//////////////////////////////////////////////////////////////////////
// SStruct.hpp
//
// This is a class for reading and writing of RNA secondary
// structures.  The file formats supported include
// 
//     (1) BPSEQ
//     (2) FASTA
//     (3) plain text (raw)
//////////////////////////////////////////////////////////////////////

#ifndef SSTRUCT_HPP
#define SSTRUCT_HPP

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Config.hpp" // why?
#include "Utilities.hpp"

//////////////////////////////////////////////////////////////////////
// class SStruct
//////////////////////////////////////////////////////////////////////

class SStruct
{
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    std::vector<int> mapping;
    std::vector<int> mapping2;
    std::vector<int> mapping3;
    std::vector<double> kd_data;
    std::vector<std::vector<double> > unpaired_potentials;
    bool has_struct;
    bool has_kd;
    bool has_evidence;
    int num_data_sources;
    std::vector<bool> which_evidence;

    // automatic file format detection
    int AnalyzeFormat(const std::string &filename) const;

    // load file of a particular file format
    void LoadFASTA(const std::string &filename);
    void LoadRAW(const std::string &filename);
    void LoadBPSEQ(const std::string &filename);
    void LoadBPP2SEQ(const std::string &filename);
    void LoadBPP2TSEQ(const std::string &filename);
    void LoadBPSEQR(const std::string &filename);

    // perform standard character conversions for RNA sequence and structures
    std::string FilterSequence(std::string sequence) const;
    std::string FilterParens(std::string sequence) const;

    // convert a pseudoknot-free parenthesized structure to a mapping and back
    std::vector<int> ConvertParensToMapping(const std::string &parens) const;
    std::string ConvertMappingToParens(const std::vector<int> &mapping) const;

    // check that a (possibly pseudoknotted) mapping is valid
    void ValidateMapping(const std::vector<int> &mapping) const;
    
public:

    // integer constants used to identify nucleotides which are either
    // unpaired or whose pairing is not known
    static const int UNPAIRED;
    static const int UNKNOWN;
    
    // constant used to identify nucleotides which have unknown probabilities
    // for being unpaired
    static const double UNKNOWN_POTENTIAL;
    
    // constructor and destructor
    SStruct();
    SStruct(const std::string &filename);
    SStruct(const std::string &filename, const int num_data_sources);
    SStruct(const SStruct &rhs);
    ~SStruct();

    // load sequence and struture from file
    void Load(const std::string &filename);

    // load sequence directly in API
    void LoadAPI(const std::string sequence, const std::string structure);

    // assignment operator
    const SStruct& operator=(const SStruct &rhs);

    // check for pseudoknots
    bool ContainsPseudoknots() const;

    // remove noncomplementary base-pairs
    void RemoveNoncomplementaryPairings(const int seq = 0);

    // output in various formats
    void WriteBPSEQ(std::ostream &outfile, const int seq = 0) const;
    void WriteParens(std::ostream &outfile) const;
    void WriteParensOnly(std::ostream &outfile) const;
    void WriteBPPSEQ(std::ostream &outfile, std::vector<std::vector<double> > potential, const int seq = 0);

    // compute alignment percent identity
    double ComputePercentIdentity() const;

    // compute position-based sequence weights
    std::vector<double> ComputePositionBasedSequenceWeights() const;

    // set mapping
    void SetMapping(const std::vector<int> &mapping);

    //////////////////////////////////////////////////////////////////////
    // Getters
    //////////////////////////////////////////////////////////////////////
    
    const std::vector<std::string> &GetNames() const { return names; }
    const std::vector<std::string> &GetSequences() const { return sequences; }
    const std::vector<int> &GetMapping() const { return mapping; }
    const std::vector<int> &GetMapping2() const { return mapping2; }
    const std::vector<int> &GetMapping3() const { return mapping3; }

    const std::vector<double> &GetKdData() const { return kd_data; }
    const std::vector<double> &GetUnpairedPotentials(int which) const { return unpaired_potentials[which]; }
    const std::vector<double> &GetPairedPotentials(int which) const { return unpaired_potentials[which]; }
    // Note: paired_potentials and unpaired_potentials arrays are the same here; the logic to handle them differently is in InferenceEngine
    int GetLength() const { return int(mapping.size())-1; }
    int GetNumSequences() const { return int(sequences.size()); }
    bool HasStruct() const { return has_struct; }
    bool HasKD() const { return has_kd; }
    bool HasEvidence(int which_data) const { return which_evidence[which_data]; }
    bool HasEvidence() const { return has_evidence; }
    const std::string ReturnMappingAsParens();
};

#endif
