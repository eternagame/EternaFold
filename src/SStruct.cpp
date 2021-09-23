//////////////////////////////////////////////////////////////////////
// SStruct.cpp
//////////////////////////////////////////////////////////////////////

#include "SStruct.hpp"

enum FileFormat
{ 
    FileFormat_FASTA,
    FileFormat_BPSEQ,
    FileFormat_BPP2SEQ,
    FileFormat_BPP2TSEQ,
    FileFormat_RAW,
    FileFormat_BPSEQR,
    FileFormat_BPSEQ_hp,
    FileFormat_UNKNOWN
};                  

const int SStruct::UNPAIRED = 0;
const int SStruct::UNKNOWN = -1;

const double SStruct::UNKNOWN_POTENTIAL = -1;

const double THRESH_NO_DATA = 1e-5;

//////////////////////////////////////////////////////////////////////
// SStruct::SStruct()
//
// Constructor.
//////////////////////////////////////////////////////////////////////

SStruct::SStruct()
{ 
    has_struct = false;
    has_kd = false;
    has_log_p_struct = false;
    has_evidence = false;
    num_data_sources = 1;
}

//////////////////////////////////////////////////////////////////////
// SStruct::SStruct()
//
// Create object from file.
//////////////////////////////////////////////////////////////////////

SStruct::SStruct(const std::string &filename, const int num_data_sources)
{
   has_struct = false;
   has_kd = false;
   has_evidence = false;
   has_log_p_struct = false;
   this->num_data_sources = num_data_sources;
   Load(filename);
}

//////////////////////////////////////////////////////////////////////
// SStruct::SStruct()
//
// Create object from file.
//////////////////////////////////////////////////////////////////////

SStruct::SStruct(const std::string &filename)
{
   has_struct = false;
   has_kd = false;
   has_log_p_struct = false;
   has_evidence = false;
   Load(filename);
   num_data_sources = 1;
}

//////////////////////////////////////////////////////////////////////
// SStruct::Load()
//
// Load from file.  Attempt to detect the format of the file
// automatically.
//////////////////////////////////////////////////////////////////////

void SStruct::Load(const std::string &filename)
{
    // auto-detect file format and load file
    switch (AnalyzeFormat(filename))
    {
        case FileFormat_FASTA: LoadFASTA(filename); break;
        case FileFormat_RAW: LoadRAW(filename); break;
        case FileFormat_BPSEQ: LoadBPSEQ(filename); break;
        case FileFormat_BPP2TSEQ: LoadBPP2TSEQ(filename); break;
        case FileFormat_BPP2SEQ: LoadBPP2SEQ(filename); break;
        case FileFormat_BPSEQR: LoadBPSEQR(filename);break;
        case FileFormat_BPSEQ_hp: LoadBPSEQ_hp(filename);break;

        default: Error("Unable to determine file type.");
    }

    // perform character conversions
    for (size_t i = 0; i < sequences.size(); i++)
        sequences[i] = FilterSequence(sequences[i]);

    // error-checking
    ValidateMapping(mapping);
    
    //ValidateMapping(mapping2);
    //ValidateMapping(mapping3);

    //std::cout << "kd_data in SStruct::Load " << kd_data << std::endl;
}

//////////////////////////////////////////////////////////////////////
// SStruct::AnalyzeFormat()
//
// Determine file format.
//////////////////////////////////////////////////////////////////////

int SStruct::AnalyzeFormat(const std::string &filename) const
{
    std::ifstream data(filename.c_str());
    if (data.fail()) Error(("Unable to open input file: " + filename).c_str());
    
    // look for first non-blank line
    std::string s;

    while (std::getline(data, s))
        if (s.length() > 0) break;
    
    // analyze to determine file format
    FileFormat format;
    if (s[0] == '>')
        format = FileFormat_FASTA;

    else if (s[0] == 'k')
        format = FileFormat_BPSEQR;
    else if (s[0] == 'n')
        format = FileFormat_BPSEQ_hp;
    else
    {
        std::istringstream iss(s);
        int number;
        std::string i, c, j;
        if ((iss >> i >> c >> j)) 
        { 
            if (ConvertToNumber(i, number) && c.length() == 1 && ConvertToNumber(j, number))
                format = FileFormat_BPSEQ;
            // there is a number following indicating the number of data sources
            else if (j.substr(0,1) == "e")
                format = FileFormat_BPP2SEQ;
            // there is a number following indicating the number of data sources
            else if (j.substr(0,1) == "t")
                format = FileFormat_BPP2TSEQ;
            else
                format = FileFormat_UNKNOWN;
        }
        else
            format = FileFormat_RAW;
    }
    
    data.close();

    return format;
}

//////////////////////////////////////////////////////////////////////
// SStruct::LoadFASTA()
//
// Create object from a FASTA file.  Assumes that the data file has a
// FASTA format.  Optionally, a parenthesized base-pairing structure
// may be provided as one of the sequences in the file.
//////////////////////////////////////////////////////////////////////

void SStruct::LoadFASTA(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // open file for reading
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());

    // process sequences
    std::string s;
    while (std::getline(data, s))
    {
        s = Trim(s);
        if (s.length() == 0) continue;

        // check for MFA header
        if (s[0] == '>')
        {
            names.push_back(s.substr(1));
            sequences.push_back("@");
        }

        // otherwise process sequence data
        else
        {
            if (sequences.size() == 0) Error("Expected header for FASTA file: %s", filename.c_str());
            for (size_t i = 0; i < s.length(); i++)
            {
                if (isspace(s[i])) continue;
                sequences.back() += s[i];
            }
        }
    }
    
    // sanity-checks
    if (sequences.size() == 0) Error("No sequences read.");
    if (sequences[0].length() == 1) Error("Zero-length sequence read.");
    for (size_t i = 1; i < sequences.size(); i++)
        if (sequences[i].length() != sequences[0].length())
            Error("Not all sequences have the same length.");

    // determine if any of the sequences could be a consensus sequence
    bool consensus_found = false;
    size_t i = 0;
    while (i < sequences.size())
    {
        // check for alphabetic characters
        bool is_consensus = true;
        for (size_t j = 1; is_consensus && j < sequences[i].length(); j++)
            if (isalpha(sequences[i][j])) is_consensus = false;

        // extract consensus mapping
        if (is_consensus)
        {
            if (consensus_found)
                Error("More than one consensus base-pairing structure found.");
            else
            {
                mapping = ConvertParensToMapping(FilterParens(sequences[i]));
                sequences.erase(sequences.begin() + i);
                names.erase(names.begin() + i);
                consensus_found = true;
                continue;
            }
        }
        i++;
    }

    // supply empty mapping if none found
    if (!consensus_found)
    {
        mapping = std::vector<int>(sequences[0].length(), UNKNOWN);
    } else {
        has_struct = true;
    }
    
    // initialize unknown unpairedness potentials
    for (int i = 0; i < num_data_sources; i++) {
        unpaired_potentials.push_back(std::vector<double>(sequences[0].length(), UNKNOWN_POTENTIAL));
    }
    has_evidence = false;
    which_evidence.resize(num_data_sources,false);
}

//////////////////////////////////////////////////////////////////////
// SStruct::LoadRAW()
//
// Create object from raw (unformatted) file.  Assumes that exactly
// one sequence is provided, with no secondary structure.
//////////////////////////////////////////////////////////////////////

void SStruct::LoadRAW(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back(filename);
    sequences.push_back("@");

    // open file for reading
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());
    
    // now retrieve sequence data    
    std::string s;
    while (std::getline(data, s))
    {
        for (size_t i = 0; i < s.length(); i++)
        {
            if (isspace(s[i])) continue;
            sequences.back() += s[i];

        }
    }

    // sanity-checks
    if (sequences[0].length() == 1) Error("Zero-length sequence read.");

    // initialize empty secondary structure
    mapping.resize(sequences[0].length(), UNKNOWN);
    
    // initialize unknown unpairedness potentials
    for (int i = 0; i < num_data_sources; i++) {
        unpaired_potentials.push_back(std::vector<double>(sequences[0].length(), UNKNOWN_POTENTIAL));
    }
    which_evidence.resize(num_data_sources,false);
    has_evidence = false;
}

//////////////////////////////////////////////////////////////////////
// SStruct::LoadBPSEQ()
//
// Create object from BPSEQ file.  Assumes that exactly one sequence
// is provided.  Base-pairings in the file may contain pseudoknots.
// Unpaired nucleotides should be annotated with base-pairing '0', and
// nucleotides with no known pairing should be annotated with
// base-pairing '-1'.
//////////////////////////////////////////////////////////////////////

void SStruct::LoadBPSEQ(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back(filename);
    sequences.push_back("@");
    mapping.push_back(UNKNOWN);

    // open file
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());

    // process file
    std::string token;
    int row = 0;
    while (data >> token)
    {
        // read row        
        int index = 0;
        if (!ConvertToNumber(token, index)) Error("Could not read row number: %s", filename.c_str());
        if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
        if (index != row+1) Error("Rows of BPSEQ file must occur in increasing order: %s", filename.c_str());
        row = index;

        // read sequence letter
        if (!(data >> token)) Error("Expected sequence letter after row number: %s", filename.c_str());
        if (token.length() != 1) Error("Expected sequence letter after row number: %s", filename.c_str());      
        char ch = token[0];

        // read mapping        
        int maps_to = 0;
        if (!(data >> token)) Error("Expected mapping after sequence letter: %s", filename.c_str());
        if (!ConvertToNumber(token, maps_to)) Error("Could not read matching row number: %s", filename.c_str());
        if (maps_to < -1) Error("Matching row numbers must be greater than or equal to -1: %s", filename.c_str());

        sequences.back().push_back(ch);
        mapping.push_back(maps_to);
    }
    
    has_struct = true;

    // initialize unknown unpairedness potentials
    for (int i = 0; i < num_data_sources; i++) {
        unpaired_potentials.push_back(std::vector<double>(sequences[0].length(), UNKNOWN_POTENTIAL));
    }
    which_evidence.resize(num_data_sources,false);
    has_evidence = false;
}

const std::string SStruct::ReturnMappingAsParens(){
    Assert(!ContainsPseudoknots(), "Should not attempt to convert a mapping with pseudoknots.");
    std::string parens = "";

    for (int i = 1; i < int(mapping.size()); i++)
    {
        if (mapping[i] == UNKNOWN)
            parens += "?";
        else if (mapping[i] == UNPAIRED)
            parens += ".";
        else if (mapping[i] > i)
            parens += "(";
        else if (mapping[i] < i)
            parens += ")";
        else
            Error("Invalid structure.");
    }

    return parens;
}

//////////////////////////////////////////////////////////////////////
// SStruct::LoadBPSEQR() HKWS
//
// Create object from BPSEQ file with one line at top of kd values.
//  Assumes that exactly one sequence
// is provided.  Base-pairings in the file may contain pseudoknots.
// Unpaired nucleotides should be annotated with base-pairing '0', and
// nucleotides with no known pairing should be annotated with
// base-pairing '-1'.
//////////////////////////////////////////////////////////////////////

void SStruct::LoadBPSEQR(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<double>().swap(kd_data);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back(filename);
    sequences.push_back("@");
    mapping.push_back(UNKNOWN);
    mapping2.push_back(UNKNOWN);
    mapping3.push_back(UNKNOWN);


    // open file
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());

    // process file
    std::string token;
    int row = 0;
    double kd = 0.0;
    data >> token;
    if (!ConvertToNumber(token.substr(1), kd)) Error("Could not read log_kd_no_lig from %s", filename.c_str());
    kd_data.push_back(kd);

    data >> token;
    if (!ConvertToNumber(token, kd)) Error("Could not read log_kd_w_lig from %s", filename.c_str());
    kd_data.push_back(kd);

    data >> token;
    if (!ConvertToNumber(token, kd)) Error("Could not read ligand_bonus from %s", filename.c_str());
    kd_data.push_back(kd);


    //std::cout << "loading kd data " << kd << std::endl;
    //std::cout << kd_data << std::endl;
    
    while (data >> token)
    {
        // read row 
        int index = 0;
        if (!ConvertToNumber(token, index)) Error("Could not read row number: %s", filename.c_str());
        if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
        if (index != row+1) Error("Rows of BPSEQ file must occur in increasing order: %s", filename.c_str());
        row = index;

        // read sequence letter
        if (!(data >> token)) Error("Expected sequence letter after row number: %s", filename.c_str());
        if (token.length() != 1) Error("Expected sequence letter after row number: %s", filename.c_str());      
        char ch = token[0];

        // read mapping
        int maps_to = 0;
        if (!(data >> token)) Error("Expected mapping after sequence letter: %s", filename.c_str());
        if (!ConvertToNumber(token, maps_to)) Error("Could not read matching row number: %s", filename.c_str());
        if (maps_to < -1) Error("Matching row numbers must be greater than or equal to -1: %s", filename.c_str());

        // read mapping 2 
        int maps_to2 = 0;
        if (!(data >> token)) Error("Expected mapping after sequence letter: %s", filename.c_str());
        if (!ConvertToNumber(token, maps_to2)) Error("Could not read matching row number: %s", filename.c_str());
        if (maps_to2 < -1) Error("Matching row numbers must be greater than or equal to -1: %s", filename.c_str());

        // read mapping 3
        int maps_to3 = 0;
        if (!(data >> token)) Error("Expected mapping after sequence letter: %s", filename.c_str());
        if (!ConvertToNumber(token, maps_to3)) Error("Could not read matching row number: %s", filename.c_str());
        if (maps_to3 < -1) Error("Matching row numbers must be greater than or equal to -1: %s", filename.c_str());

        sequences.back().push_back(ch);
        mapping.push_back(maps_to);
        mapping2.push_back(maps_to2);
        mapping3.push_back(maps_to3);

    }
    
    has_kd = true;
    has_struct = true;

    // std::cout << "mapping" <<std::endl;
    // std::cout << mapping <<std::endl;
    // std::cout << "mapping2" <<std::endl;
    // std::cout << mapping2 <<std::endl;
    // std::cout << "mapping3" <<std::endl;
    // std::cout << mapping3 <<std::endl;
    
    // initialize unknown unpairedness potentials
    for (int i = 0; i < num_data_sources; i++) {
        unpaired_potentials.push_back(std::vector<double>(sequences[0].length(), UNKNOWN_POTENTIAL));
    }
    which_evidence.resize(num_data_sources,false);
    has_evidence = false;
}

//////////////////////////////////////////////////////////////////////
// SStruct::LoadBPSEQ_hp() HKWS, July 2021
//
// Create object from BPSEQ file with one line at top, which contains one value reflecting prob(Struct) from melt data.

//////////////////////////////////////////////////////////////////////

void SStruct::LoadBPSEQ_hp(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<double>().swap(log_p_data);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back(filename);
    sequences.push_back("@");
    mapping.push_back(UNKNOWN);

    // open file
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());

    // process file
    std::string token;
    int row = 0;
    double num = 0.0;
    data >> token;
    if (!ConvertToNumber(token.substr(1), num)) Error("Could not read log_p_struct from %s", filename.c_str());
    log_p_data.push_back(num);

    //std::cout << "loading p struct data " << log_p_data << std::endl;
    
    while (data >> token)
    {
        // read row 
        int index = 0;
        if (!ConvertToNumber(token, index)) Error("Could not read row number: %s", filename.c_str());
        if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
        if (index != row+1) Error("Rows of BPSEQ file must occur in increasing order: %s", filename.c_str());
        row = index;

        // read sequence letter
        if (!(data >> token)) Error("Expected sequence letter after row number: %s", filename.c_str());
        if (token.length() != 1) Error("Expected sequence letter after row number: %s", filename.c_str());      
        char ch = token[0];

        // read mapping
        int maps_to = 0;
        if (!(data >> token)) Error("Expected mapping after sequence letter: %s", filename.c_str());
        if (!ConvertToNumber(token, maps_to)) Error("Could not read matching row number: %s", filename.c_str());
        if (maps_to < -1) Error("Matching row numbers must be greater than or equal to -1: %s", filename.c_str());

        sequences.back().push_back(ch);
        mapping.push_back(maps_to);
    }
    
    has_log_p_struct = true;
    has_struct = true;

    //std::cout << "mapping" <<std::endl;
    //std::cout << mapping <<std::endl;

    // initialize unknown unpairedness potentials // probably don't need these for non-structure-probing data? 
    for (int i = 0; i < num_data_sources; i++) {
        unpaired_potentials.push_back(std::vector<double>(sequences[0].length(), UNKNOWN_POTENTIAL));
    }
    which_evidence.resize(num_data_sources,false);
    has_evidence = false;
}

////////////////////////////////////////////////////////////////////////////
// SStruct::LoadAPI()
//
// Create object from std::string sequence and std::string structure, to be 
// used in contrafold API developed for solving inverse folding problem.
////////////////////////////////////////////////////////////////////////////

void SStruct::LoadAPI(const std::string sequence, const std::string structure)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back("api");
    sequences.push_back("@" + sequence);

    // sanity-checks
    if (sequences[0].length() == 1) Error("Zero-length sequence read.");

    // initialize empty secondary structure
    mapping.resize(sequences[0].length(), UNKNOWN);

    mapping = ConvertParensToMapping(FilterParens("@" + structure));

    // initialize unknown unpairedness potentials
    for (int i = 0; i < num_data_sources; i++) {
        unpaired_potentials.push_back(std::vector<double>(sequences[0].length(), UNKNOWN_POTENTIAL));
    }
    which_evidence.resize(num_data_sources,false);
    has_evidence = false;
}

////////////////////////////////////////////////////////////////////////////
// SStruct::LoadBPP2SEQ()
//
// Create object from BPP2SEQ file.  Assumes that exactly one sequence
// is provided. Potentials from different data sources without true pairing.
// Potentials for the unpairedness of a nucleotide should be positive.
////////////////////////////////////////////////////////////////////////////

void SStruct::LoadBPP2SEQ(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back(filename);
    sequences.push_back("@");
    mapping.push_back(UNKNOWN);

    // open file
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());

    // process file
    std::string token;
    int row = 0;

    which_evidence.resize(num_data_sources,false);
    int num_data_sources_local = 0;
    for (int i = 0; i < num_data_sources; i++) {
        std::vector<double> unpaired_potentials_current;
        unpaired_potentials.push_back(unpaired_potentials_current);
    }
    
    while (data >> token)
    {
        // read row        
        int index = 0;
        if (!ConvertToNumber(token, index)) Error("Could not read row number: %s", filename.c_str());
        if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
        if (index != row+1) Error("Rows of BPSEQ file must occur in increasing order: %s", filename.c_str());
        row = index;

        // read sequence letter
        if (!(data >> token)) Error("Expected sequence letter after row number: %s", filename.c_str());
        if (token.length() != 1) Error("Expected sequence letter after row number: %s", filename.c_str());      
        char ch = token[0];
        
        // read the "e" letter
        if (!(data >> token)) Error("Expected 'e' after sequence letter: %s", filename.c_str());
        if (token.substr(0,1) != "e") Error("Expected 'e' after sequence letter: %s", filename.c_str());

        bool success = ConvertToNumber(token.substr(1,token.length()),num_data_sources_local);
	if (!success)
            Error("Number of data sources must be an integer!");
        if (num_data_sources != num_data_sources_local)
            Error("Number of data sources at command line (%d) do not match number of data sources in file (%d)!",num_data_sources,num_data_sources_local);

        for (int i = 0; i < num_data_sources_local; i++)
        {
            // read probing data
            double potential;
            if (!(data >> token)) Error("Expected unpaired potential after sequence letter: %s", filename.c_str());
            if (!ConvertToNumber(token, potential)) Error("Could not read unpaired potential: %s", filename.c_str());
            
            unpaired_potentials[i].push_back(potential);

            // evidence for dataset only if saw a non-zero potential
            if (potential > THRESH_NO_DATA)
                which_evidence[i] = true;
        }

        sequences.back().push_back(ch);
        
    }

    // initialize empty secondary structure
    mapping.resize(sequences[0].length(), UNKNOWN);

    int i;
    for (i = 0; i < num_data_sources; i++)
    {
        if (which_evidence[i])
            break;
    }
    if (i >= num_data_sources)
       Error("The data potentials appear to all be null. If no evidence, use a different format for the input file.");

    has_evidence = true;
}

////////////////////////////////////////////////////////////////////////////////
// SStruct::LoadBPP2TSEQ()
//
// Create object from BPP2TSEQ file.  Assumes that exactly one sequence
// is provided. Potentials from different data sources followed by true pairing.
// Potentials for the unpairedness of a nucleotide should be positive.
////////////////////////////////////////////////////////////////////////////////

void SStruct::LoadBPP2TSEQ(const std::string &filename)
{
    // clear any previous data
    std::vector<std::string>().swap(names);
    std::vector<std::string>().swap(sequences);
    std::vector<int>().swap(mapping);
    std::vector<std::vector<double> >().swap(unpaired_potentials);
    std::vector<bool>().swap(which_evidence);

    // initialize
    names.push_back(filename);
    sequences.push_back("@");
    mapping.push_back(UNKNOWN);

    // open file
    std::ifstream data(filename.c_str());
    if (data.fail()) Error("Unable to open input file: %s", filename.c_str());

    // process file
    std::string token;
    int row = 0;

    which_evidence.resize(num_data_sources,false);
    int num_data_sources_local = 0;
    for (int i = 0; i < num_data_sources; i++) {
        std::vector<double> unpaired_potentials_current;
        unpaired_potentials.push_back(unpaired_potentials_current);
    }

    while (data >> token)
    {
        // read row        
        int index = 0;
        if (!ConvertToNumber(token, index)) Error("Could not read row number: %s", filename.c_str());
        if (index <= 0) Error("Row numbers must be positive: %s", filename.c_str());
        if (index != row+1) Error("Rows of BPSEQ file must occur in increasing order: %s", filename.c_str());
        row = index;

        // read sequence letter
        if (!(data >> token)) Error("Expected sequence letter after row number: %s", filename.c_str());
        if (token.length() != 1) Error("Expected sequence letter after row number: %s", filename.c_str());      
        char ch = token[0];

        // read the "t" letter
        if (!(data >> token)) Error("Expected 't' after sequence letter: %s", filename.c_str());
        if (token.substr(0,1) != "t") Error("Expected 't' after sequence letter: %s", filename.c_str());

        bool success = ConvertToNumber(token.substr(1,token.length()),num_data_sources_local);
	if (!success)
            Error("Number of data sources must be an integer!");
        if (num_data_sources != num_data_sources_local)
            Error("Number of data sources at command line (%d) do not match number of data sources in file (%d)!",num_data_sources,num_data_sources_local);

        for (int i = 0; i < num_data_sources_local; i++)
        {
            // read probing data
            double potential;
            if (!(data >> token)) Error("Expected unpaired potential after sequence letter: %s", filename.c_str());
            if (!ConvertToNumber(token, potential)) Error("Could not read unpaired potential: %s", filename.c_str());
            
            unpaired_potentials[i].push_back(potential);

            // evidence for dataset only if saw a non-zero potential
            if (potential > THRESH_NO_DATA)
                which_evidence[i] = true;
        }
        
        // read mapping        
        int maps_to = 0;
        if (!(data >> token)) Error("Expected mapping after sequence letter: %s", filename.c_str());
        if (!ConvertToNumber(token, maps_to)) Error("Could not read matching row number: %s", filename.c_str());
        if (maps_to < -1) Error("Matching row numbers must be greater than or equal to -1: %s", filename.c_str());

        sequences.back().push_back(ch);
        mapping.push_back(maps_to);
    }

    int i;
    for (i = 0; i < num_data_sources; i++)
    {
        if (which_evidence[i])
            break;
    }
    if (i >= num_data_sources)
       Error("The data potentials appear to all be null. If no evidence, use a different format for the input file.");

    has_struct = true;
    has_evidence = true;
}

//////////////////////////////////////////////////////////////////////
// SStruct::SStruct()
//
// Create object from existing SStruct object.
//////////////////////////////////////////////////////////////////////

SStruct::SStruct(const SStruct &rhs) :
    names(rhs.names),
    sequences(rhs.sequences),
    mapping(rhs.mapping),
    mapping2(rhs.mapping2),
    mapping3(rhs.mapping3),

    kd_data(rhs.kd_data),
    log_p_data(rhs.log_p_data),
    unpaired_potentials(rhs.unpaired_potentials),
    has_struct(rhs.has_struct),
    has_kd(rhs.has_kd),
    has_log_p_struct(rhs.has_log_p_struct),
    has_evidence(rhs.has_evidence),
    num_data_sources(rhs.num_data_sources),
    which_evidence(rhs.which_evidence)
{}

//////////////////////////////////////////////////////////////////////
// SStruct::~SStruct()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

SStruct::~SStruct()
{}

//////////////////////////////////////////////////////////////////////
// SStruct::operator=
//
// Assignment operator.
//////////////////////////////////////////////////////////////////////

const SStruct &SStruct::operator=(const SStruct &rhs)
{
    if (this != &rhs)
    {
        names = rhs.names;
        sequences = rhs.sequences;
        mapping = rhs.mapping;
        mapping2 = rhs.mapping2;
        mapping3 = rhs.mapping3;

        kd_data = rhs.kd_data;
        log_p_data = rhs.log_p_data;
        unpaired_potentials = rhs.unpaired_potentials;
        which_evidence = rhs.which_evidence;
        has_struct = rhs.has_struct;
        has_kd = rhs.has_kd;
        has_log_p_struct = rhs.has_log_p_struct;
        has_evidence = rhs.has_evidence;
        num_data_sources = rhs.num_data_sources;
    }

    return *this;
}

//////////////////////////////////////////////////////////////////////
// SStruct::FilterSequence()
//
// Perform character conversions to put the RNA sequence in a standard
// format.
//////////////////////////////////////////////////////////////////////

std::string SStruct::FilterSequence(std::string sequence) const
{
    if (sequence[0] != '@') Error("Improperly formatted sequence.");
    
    for (size_t i = 1; i < sequence.length(); i++)
    {
        bool uppercase = isupper(sequence[i]);
        char c = tolower(sequence[i]);

        switch (c)
        {
            case '.': c = '-'; break;
            case 't': c = 'u'; break;
            case '-': case 'a': case 'c': case 'g': case 'u': break;
            default:
                if (isalpha(c))
                    c = 'n';
                else 
                    Error("Unexpected character '%c' in sequence.", c);
                break;
        }
        
        if (uppercase) c = toupper(c);
        sequence[i] = c;
    }
    
    return sequence;
}

//////////////////////////////////////////////////////////////////////
// SStruct::FilterParens()
//
// Perform character conversions as needed.
//////////////////////////////////////////////////////////////////////

std::string SStruct::FilterParens(std::string sequence) const
{
    if (sequence[0] != '@') Error("Improperly formatted sequence.");
    
    for (size_t i = 1; i < sequence.length(); i++)
    {
        switch (sequence[i])
        {
            case '-': sequence[i] = '.'; break;
            case '?': case '.': case '(': case ')': break;
            default: Error("Unexpected character '%c' in parenthesized structure.", sequence[i]);
        }
    }

    return sequence;
}

//////////////////////////////////////////////////////////////////////
// SStruct::ConvertParensToMapping()
//
// Convert a parenthesized string to a mapping.  No pseudoknots
// allowed.
//////////////////////////////////////////////////////////////////////

std::vector<int> SStruct::ConvertParensToMapping(const std::string &parens) const
{
    std::vector<int> mapping(parens.length(), UNKNOWN);
    std::vector<int> stack;
   
    Assert(parens[0] == '@', "Invalid parenthesized string.");
    for (int i = 1; i < int(parens.length()); i++)
    {
        switch (parens[i])
        {
            case '?': break;
            case '.': mapping[i] = UNPAIRED; break;
            case '(': stack.push_back(i); break;
            case ')':
                if (stack.size() == 0) Error("Parentheses mismatch.");
                mapping[i] = stack.back();
                mapping[stack.back()] = i;
                stack.pop_back();
                break;
            default:
                Error("Unexpected character '%c' in parenthesized structure.", parens[i]);
        }
    }
    if (stack.size() != 0) Error("Parentheses mismatch.");

    return mapping;
}

//////////////////////////////////////////////////////////////////////
// SStruct::ConvertMappingToParens()
//
// Convert a mapping to a parenthesized string.  No pseudoknots
// allowed.
//////////////////////////////////////////////////////////////////////

std::string SStruct::ConvertMappingToParens(const std::vector<int> &mapping) const
{
    Assert(!ContainsPseudoknots(), "Should not attempt to convert a mapping with pseudoknots.");
    std::string parens = "@";

    for (int i = 1; i < int(mapping.size()); i++)
    {
        if (mapping[i] == UNKNOWN)
            parens += "?";
        else if (mapping[i] == UNPAIRED)
            parens += ".";
        else if (mapping[i] > i)
            parens += "(";
        else if (mapping[i] < i)
            parens += ")";
        else
            Error("Invalid structure.");
    }

    return parens;
}

//////////////////////////////////////////////////////////////////////
// SStruct::ValidateMapping()
//
// Check that a std::vector<int> represents a valid secondary
// structure mapping.  Pseudoknots are allowed.
//////////////////////////////////////////////////////////////////////

void SStruct::ValidateMapping(const std::vector<int> &mapping) const
{
    if (mapping.size() == 0 || mapping[0] != UNKNOWN) Error("Invalid mapping.");
    for (int i = 1; i < int(mapping.size()); i++)
    {
        if (mapping[i] == UNPAIRED || mapping[i] == UNKNOWN)
            continue;
        if (mapping[i] < 1 || mapping[i] >= int(mapping.size()))
            Error("Position %d of sequence maps to invalid position.", i);
        if (mapping[mapping[i]] != i)
            Error("Positions %d and %d of sequence do not map to each other.", i, mapping[i]);
        if (mapping[i] == i)
            Error("Position %d of sequence maps to itself.", i);
    }
}

//////////////////////////////////////////////////////////////////////
// SStruct::ContainsPseudoknots()
//
// Check if secondary structure contains pseudoknots.
//////////////////////////////////////////////////////////////////////

bool SStruct::ContainsPseudoknots() const
{
    std::vector<int> stack;

    for (int i = 1; i < int(mapping.size()); i++)
    {
        if (mapping[i] == UNPAIRED || mapping[i] == UNKNOWN)
            continue;
        if (mapping[i] > i)
            stack.push_back(i);
        else if (mapping[i] < i)
        {
            if (stack.back() == mapping[i])
                stack.pop_back();
            else
                return true;
        }
        else
            Error("Invalid structure: positions may not map to themselves.");
    }
    if (stack.size() != 0) Error("Invalid structure: bad pairings found.");

    return false;
}

//////////////////////////////////////////////////////////////////////
// SStruct::RemoveNoncomplementaryPairings()
//
// Remove all non-{AU,CG,GU} pairings from mapping.
//////////////////////////////////////////////////////////////////////

void SStruct::RemoveNoncomplementaryPairings(const int seq)
{
    if (seq < 0 || seq >= int(sequences.size())) Error("Reference to invalid sequence.");
    Assert(sequences[seq].length() == mapping.size(), "Inconsistent lengths.");
    
    for (int i = 1; i < int(mapping.size()); i++)
    {
        if (mapping[i] > i && !IsComplementary(sequences[seq][i], sequences[seq][mapping[i]]))
        {
            mapping[mapping[i]] = UNPAIRED;
            mapping[i] = UNPAIRED;
        }
    }
}

//////////////////////////////////////////////////////////////////////
// SStruct::WriteBPSEQ()
//
// Write sequence in BPSEQ format.  The BPSEQ format can only handle
// single sequences, so it will only print out the sequence "seq".
//////////////////////////////////////////////////////////////////////

void SStruct::WriteBPSEQ(std::ostream &outfile, const int seq) const
{
    if (seq < 0 || seq >= int(sequences.size())) Error("Reference to invalid sequence.");
    Assert(sequences[seq].length() == mapping.size(), "Inconsistent lengths.");
    
    for (size_t i = 1; i < mapping.size(); i++)
        outfile << i << ' ' << sequences[seq][i] << ' ' << mapping[i] << std::endl;
}

//////////////////////////////////////////////////////////////////////
// SStruct::WriteParens()
//
// Write sequence in parenthesized format.  This routine assumes that
// the structure does not contain pseudoknots.  All sequences are
// printed.
//////////////////////////////////////////////////////////////////////

void SStruct::WriteParens(std::ostream &outfile) const
{
    if (ContainsPseudoknots()) Error("Cannot write structure containing pseudoknots using parenthesized format.");
    
    // print sequences
    for (size_t k = 0; k < sequences.size(); k++)
    {
        outfile << ">" << names[k] << std::endl;
        outfile << sequences[k].substr(1) << std::endl;
    }

    // print structure
    outfile << ">structure" << std::endl;
    outfile << ConvertMappingToParens(mapping).substr(1) << std::endl;
}

//////////////////////////////////////////////////////////////////////
// SStruct::ComputePercentIdentity()
//
// Compute average pairwise percent identity for the alignment.
// The pairwise PID = # of identities / MIN(len1, len2).
//////////////////////////////////////////////////////////////////////

double SStruct::ComputePercentIdentity() const
{
    double pid = 0.0;
    for (size_t i = 0; i < sequences.size(); i++)
    {
        for (size_t j = i+1; j < sequences.size(); j++)
        {
            int identities = 0;
            int len1 = 0;
            int len2 = 0;

            const std::string &s = sequences[i];
            const std::string &t = sequences[j];

            for (size_t k = 0; k < s.length(); k++)
            {
                if (isalpha(s[k])) len1++;
                if (isalpha(t[k])) len2++;
                if (isalpha(s[k]) && toupper(s[k]) == toupper(t[k])) identities++;
            }

            int den = std::min(len1, len2);
            double pairwise_pid = (den == 0 ? 0.0 : double(identities) / den);

            pid += pairwise_pid;
        }
    }

    return pid / (sequences.size() * (sequences.size() - 1) / 2);
}

//////////////////////////////////////////////////////////////////////
// SStruct::ComputePositionBasedSequenceWeights()
//
// Compute sequence weights according to:
//    Henikoff, S., and Henikoff J.  1994.  Position-based
//    sequence weights.  J Mol Biol 243(4):574-578.
//////////////////////////////////////////////////////////////////////

std::vector<double> SStruct::ComputePositionBasedSequenceWeights() const
{
    std::vector<double> weights(sequences.size(), 0.0);
    std::vector<int> counts(256);

    for (size_t i = 1; i < sequences[0].length(); i++)
    {
        int diversity = 0;
        std::fill(counts.begin(), counts.end(), 0);
        
        for (size_t j = 0; j < sequences.size(); j++)
        {
            if (counts[BYTE(sequences[j][i])] == 0) diversity++;
            ++(counts[BYTE(sequences[j][i])]);
        }
        
        for (size_t j = 0; j < sequences.size(); j++)
            weights[j] += 1.0 / (diversity * counts[BYTE(sequences[j][i])]);
    }

    weights /= Sum(weights);
    return weights;
}

//////////////////////////////////////////////////////////////////////
// SStruct::SetMapping()
//
// Set secondary structure mapping.
//////////////////////////////////////////////////////////////////////

void SStruct::SetMapping(const std::vector<int> &mapping)
{
    this->mapping = mapping;
    ValidateMapping(mapping);
}
