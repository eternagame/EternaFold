////////////////////////////////////////////////////////////
// ScorePrediction.cpp
//
// Score a test prediction file against a reference.
////////////////////////////////////////////////////////////

#include "SStruct.hpp"
#include "Utilities.hpp"

///////////////////////////////////////////////////////////////////////////
// ComputeIntersection()
//
// Compute intersection of two sets.
///////////////////////////////////////////////////////////////////////////

std::set<std::vector<int> > ComputeIntersection (const std::set<std::vector<int> > &set1,
                                                 const std::set<std::vector<int> > &set2)
{
    std::set<std::vector<int> > ret;
    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                          std::insert_iterator<std::set<std::vector<int> > >(ret, ret.begin()));
    return ret;
}

///////////////////////////////////////////////////////////////////////////
// Match()
//
// Check if sequence matches.
///////////////////////////////////////////////////////////////////////////

bool Match(const std::pair<std::string, std::string> &s,
           const std::pair<std::string, std::string> &t)
{
    for (size_t i = 0; i < s.first.length(); i++)
    {
        if (toupper(s.first[i]) == 'N' || toupper(t.first[i]) == 'N') continue;
        if (toupper(s.first[i]) != toupper(t.first[i])) return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////
// Canonicalize()
//
// Convert sequence into standard unique form.
///////////////////////////////////////////////////////////////////////////

std::vector<std::pair<std::pair<std::string, std::string>, int> > Canonicalize(const SStruct &ref, const bool use_protein)
{
    std::vector<std::pair<std::pair<std::string, std::string>, int> > ret;
    for (int i = 0; i < ref.GetNumSequences(); i++)
    {
        const std::string &sequence = ref.GetSequences()[i];
        const std::string &name = ref.GetNames()[i];

        std::string projected_sequence;
        for (size_t j = 1; j < sequence.length(); j++)
        {
            const bool upper = isupper(sequence[j]);
            char ch = toupper(sequence[j]);
            if (isalpha(ch))
            {
                if (use_protein)
                {
                    if (std::string("ACDEFGHIKLMNPQRSTVWY").find(ch) == std::string::npos)
                        ch = 'N';
                }
                else
                {
                    if (ch == 'T')
                        ch = 'U';
                    else if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'U')
                        ch = 'N';
                }
                if (!upper) ch = tolower(ch);
                projected_sequence.push_back(ch);
            }
        }

        ret.push_back(std::make_pair(std::make_pair(projected_sequence, name), i));
    }
    return ret;
}

///////////////////////////////////////////////////////////////////////////
// GetSequenceOrdering()
//
// Compute mapping from test sequence indices to reference sequence
// indices.
///////////////////////////////////////////////////////////////////////////

std::vector<int> GetSequenceOrdering(const SStruct &ref,
                                     const SStruct &test,
                                     const std::string &ref_filename,
                                     const std::string &test_filename,
                                     const bool use_protein)
{
    std::vector<std::pair<std::pair<std::string, std::string>, int> > ref_strings = Canonicalize(ref, use_protein);
    std::vector<std::pair<std::pair<std::string, std::string>, int> > test_strings = Canonicalize(test, use_protein);

    Assert(ref_strings.size() == test_strings.size(), "Dimension mismatch.");

    std::vector<int> mapping(ref_strings.size());
    std::vector<int> used(test_strings.size());

    // assign each reference sequence to a test sequence
    
    for (size_t i = 0; i < ref_strings.size(); i++)
    {
        bool found = false;
        for (size_t j = 0; !found && j < test_strings.size(); j++)
        {
            if (used[j]) continue;
            if (Match(ref_strings[i].first, test_strings[j].first))
            {
                used[j] = 1;
                mapping[i] = j;
                found = true;
            }
        }
        if (!found)
        {
            Error("Unable to find matching sequence (%s vs %s)\n>%s\n%s",
                  ref_filename.c_str(), test_filename.c_str(),
                  ref_strings[i].first.second.c_str(),
                  ref_strings[i].first.first.c_str());
        }
    }

    return mapping;
}

///////////////////////////////////////////////////////////////////////////
// AddPairings()
// 
// Add all base-pairing positions in a sequence.
///////////////////////////////////////////////////////////////////////////

void AddPairings(std::set<std::vector<int> > &pairings, int ii, const std::string &s, const std::vector<int> &mapping)
{
    const std::vector<int> s_mapping = GetSequenceMapping(s);
    
    std::vector<int> pairing(3);
    pairing[0] = ii;

    if (s.length() != mapping.size()) Error("Dimension mismatch for alignment.");
    
    for (int i = 1; i < int(s.length()); i++)
    {
        if (mapping[i] > i && isalpha(s[i]) && isalpha(s[mapping[i]]))
        {
            pairing[1] = s_mapping[i];
            pairing[2] = s_mapping[mapping[i]];
            pairings.insert(pairing);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// AddMatches()
// 
// Add all matched positions in an alignment.
///////////////////////////////////////////////////////////////////////////

void AddMatches(std::set<std::vector<int> > &matches, int ii, int jj, const std::string &s, const std::string &t, const bool core_blocks_only, const bool use_protein)
{
    const std::vector<int> s_mapping = GetSequenceMapping(s);
    const std::vector<int> t_mapping = GetSequenceMapping(t);

    std::vector<int> match(4);
    match[0] = ii;
    match[1] = jj;

    if (s.length() != t.length()) Error("Dimension mismatch for alignment.");
    
    for (size_t i = 1; i < s.length(); i++)
    {
        if (isalpha(s[i]) && isalpha(t[i]))
        {
            if (use_protein && core_blocks_only && (!isupper(s[i]) || !isupper(t[i])))
                continue;
            match[2] = s_mapping[i];
            match[3] = t_mapping[i];
            matches.insert(match);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// ComputeScores()
//
// Compute sensitivity and specificity.
///////////////////////////////////////////////////////////////////////////

void ComputeScores(const std::string &ref_filename,
                   const std::string &test_filename,
                   const bool use_core_blocks,
                   const bool use_protein)
{
    SStruct ref(ref_filename);
    SStruct test(test_filename);

    if (ref.GetNumSequences() != test.GetNumSequences())
        Error("%s (%d) and %s (%d) have different numbers of sequences.",
              ref_filename.c_str(), ref.GetNumSequences(),
              test_filename.c_str(), test.GetNumSequences());

    std::vector<int> test_ordering = GetSequenceOrdering(ref, test, ref_filename, test_filename, use_protein);

    std::set<std::vector<int> > reference_pairings;
    std::set<std::vector<int> > test_pairings;
    std::set<std::vector<int> > reference_matches;
    std::set<std::vector<int> > test_matches;
    
    for (int i = 0; i < ref.GetNumSequences(); i++)
    {
        AddPairings(reference_pairings, i, ref.GetSequences()[i], ref.GetMapping());
        AddPairings(test_pairings, i, test.GetSequences()[test_ordering[i]], test.GetMapping());
        
        for (int j = i+1; j < ref.GetNumSequences(); j++)
        {
            AddMatches(reference_matches, i, j, ref.GetSequences()[i], ref.GetSequences()[j], use_core_blocks, use_protein);
            AddMatches(test_matches, i, j, test.GetSequences()[test_ordering[i]], test.GetSequences()[test_ordering[j]], false, use_protein);
        }
    }

    std::set<std::vector<int> > correct_pairings = ComputeIntersection(reference_pairings, test_pairings);
    std::set<std::vector<int> > correct_matches = ComputeIntersection(reference_matches, test_matches);

    double Qscore = (reference_matches.size() == 0) ? 1.0 : double(correct_matches.size()) / reference_matches.size();
    double fM = (test_matches.size() == 0) ? 1.0 : double(correct_matches.size()) / test_matches.size();

    double sensitivity = (reference_pairings.size() == 0) ? 1.0 : double(correct_pairings.size()) / reference_pairings.size();
    double ppv = (test_pairings.size() == 0) ? 1.0 : double(correct_pairings.size()) / test_pairings.size();

    std::cout << "ref=" << ref_filename << "; test=" << test_filename << "; N=" << ref.GetNumSequences()
              << "; ref_len=" << ref.GetLength() << "; ref_pid=" << ref.ComputePercentIdentity()
              << "; test_len=" << test.GetLength() << "; test_pid=" << test.ComputePercentIdentity()
              << "; Q=" << Qscore << "; fM=" << fM
              << "; sens=" << sensitivity << "; ppv=" << ppv << std::endl;
}

///////////////////////////////////////////////////////////////////////////
// main()
//
// Main program.
///////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        std::cerr << std::endl
                  << "Usage: " << argv[0] << " [protein|rna] REF TEST [--core]" << std::endl
                  << std::endl
                  << "       where REF    is the name of the reference file (in BPSEQ or FASTA format)" << std::endl
                  << "             TEST   is the name of the test file (in BPSEQ or FASTA format)" << std::endl
                  << std::endl;
        exit(1);
    }

    std::vector<std::string> filenames;
    bool use_protein = false;
    bool use_core_blocks = false;
    
    if (std::string(argv[1]) != "protein" &&
        std::string(argv[1]) != "rna")
    {
        Error("First parameter must be either \"protein\" or \"rna\".");
    }
    else
    {
        use_protein = (std::string(argv[1]) == "protein");
    }
       
    for (int i = 2; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (std::string(argv[i]) == "--core")
            {
                use_core_blocks = true;
            }
            else
            {
                Error("Unknown argument: %s", argv[i]);
            }               
        }
        else
        {
            filenames.push_back(argv[i]);
        }
    }

    if (filenames.size() != 2) Error("Incorrect number of filenames specified.");
    
    ComputeScores(filenames[0], filenames[1], use_core_blocks, use_protein);
}
