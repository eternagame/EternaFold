//////////////////////////////////////////////////////////////////////
// InferenceEngine.ipp
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Wrapper macros for certain model features.
//////////////////////////////////////////////////////////////////////


#if defined(PARAMS_EVIDENCE)
template<class RealT>
inline RealT InferenceEngine<RealT>::ScorePairedUnpositionEvidenceRaw(int which_data,int i) const
{
    return score_paired_position_raw[which_data][i-1];
}
#else
#define ScorePairedPositionEvidenceRaw(i,j) ( RealT(0) )
#endif
#define CountPairedPositionEvidenceRaw(i,v)


#if defined(PARAMS_EVIDENCE)
template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreUnpairedPositionEvidenceRaw(int which_data,int i) const
{
    return score_unpaired_position_raw[which_data][i-1];
}
#else
#define ScoreUnpairedPositionEvidenceRaw(i,j) ( RealT(0) )
#endif
#define CountUnpairedPositionEvidenceRaw(i,v)


// PARAMS_EVIDENCE

// score for leaving s[i] (nucleotide i) unpaired based on evidence
// (note: score_unpaired_position[k] refers to nucleotide k+1)

#if defined(PARAMS_EVIDENCE)
template<class RealT>
RealT InferenceEngine<RealT>::ScoreUnpairedPositionEvidence(int i) const
{
   RealT sum = 0;
   for (int n = 0; n < num_data_sources; n++)
       sum = sum + kappa*score_unpaired_position[n][i-1];
   return sum;
}
#else
#define ScoreUnpairedPositionEvidence(i) ( RealT(0) )
#endif
#define CountUnpairedPositionEvidence(i,v)

// score for making s[i] (nucleotide i) paired based on evidence
// (note: score_paired_position[k] refers to nucleotide k+1)

#if defined(PARAMS_EVIDENCE)
template<class RealT>
RealT InferenceEngine<RealT>::ScorePairedPositionEvidence(int i) const
{
   RealT sum = 0;
   for (int n = 0; n < num_data_sources; n++)
       sum = sum + kappa * score_paired_position[n][i-1];
   return sum;
}
#else
#define ScorePairedPositionEvidence(i) ( RealT(0) )
#endif
#define CountPairedPositionEvidence(i,v)

// score for leaving s[i+1...j] (nucleotides i+1,...j) unpaired based on evidence
// (note: score_unpaired_position[k] refers to nucleotide k+1)

#if defined(PARAMS_EVIDENCE)
template<class RealT>
RealT InferenceEngine<RealT>::ScoreUnpairedEvidence(int i, int j) const
{
    RealT sum = 0; 
    for (int k = i+1; k <= j; k++) { sum += ScoreUnpairedPositionEvidence(k); } 
    return sum;  
}
#else
#define ScoreUnpairedEvidence(i,j) ( RealT(0) )
#endif
#define CountUnpairedEvidence(i,j,v)

// score for leaving s[i] unpaired

#if defined(HAMMING_LOSS)
#define ScoreUnpairedPosition(i) (loss_unpaired_position[i])
#else
#define ScoreUnpairedPosition(i) (RealT(0))
#endif
#define CountUnpairedPosition(i,v)

// score for leaving s[i+1...j] unpaired

#if defined(HAMMING_LOSS)
#define ScoreUnpaired(i,j) (loss_unpaired[offset[i]+j])
#else
#define ScoreUnpaired(i,j) (RealT(0))
#endif
#define CountUnpaired(i,j,v)

// score for a base pair which is not part of any helix

#if PARAMS_ISOLATED_BASE_PAIR
#define ScoreIsolated() score_isolated_base_pair.first
#define CountIsolated(v) { score_isolated_base_pair.second += (v); }
#else
#define ScoreIsolated() RealT(0)
#define CountIsolated(v)
#endif

// base score for a multi-branch loop

#if PARAMS_MULTI_LENGTH
#define ScoreMultiBase() score_multi_base.first
#define CountMultiBase(v) { score_multi_base.second += (v); }
#else
#define ScoreMultiBase() RealT(0)
#define CountMultiBase(v)
#endif

// score for a base-pair adjacent to a multi-branch loop

#if PARAMS_MULTI_LENGTH
#define ScoreMultiPaired() score_multi_paired.first
#define CountMultiPaired(v) { score_multi_paired.second += (v); }
#else
#define ScoreMultiPaired() RealT(0)
#define CountMultiPaired(v)
#endif

// score for each unpaired position in a multi-branch loop based on evidence

#if PARAMS_MULTI_LENGTH
#define ScoreMultiUnpairedEvidence(i) (score_multi_unpaired.first + ScoreUnpairedPosition(i) + ScoreUnpairedPositionEvidence(i))
#define CountMultiUnpairedEvidence(i,v) { score_multi_unpaired.second += (v); CountUnpairedPosition(i,v); CountUnpairedPositionEvidence(i,v); }
#else
#define ScoreMultiUnpairedEvidence(i) (ScoreUnpairedPosition(i) + ScoreUnpairedPositionEvidence(i))
#define CountMultiUnpairedEvidence(i,v) { CountUnpairedPosition(i,v); CountUnpairedPositionEvidence(i,v); }
#endif

// score for each unpaired position in a multi-branch loop

#if PARAMS_MULTI_LENGTH
#define ScoreMultiUnpaired(i) (score_multi_unpaired.first + ScoreUnpairedPosition(i))
#define CountMultiUnpaired(i,v) { score_multi_unpaired.second += (v); CountUnpairedPosition(i,v); }
#else
#define ScoreMultiUnpaired(i) (ScoreUnpairedPosition(i))
#define CountMultiUnpaired(i,v) { CountUnpairedPosition(i,v); }
#endif

// score for each base-pair adjacent to an external loop

#if PARAMS_EXTERNAL_LENGTH
#define ScoreExternalPaired() score_external_paired.first
#define CountExternalPaired(v) { score_external_paired.second += (v); }
#else
#define ScoreExternalPaired() RealT(0)
#define CountExternalPaired(v)
#endif

// score for each unpaired position in an external loop based on evidence

#if PARAMS_EXTERNAL_LENGTH
#define ScoreExternalUnpairedEvidence(i) (score_external_unpaired.first + ScoreUnpairedPosition(i) + ScoreUnpairedPositionEvidence(i))
#define CountExternalUnpairedEvidence(i,v) { score_external_unpaired.second += (v); CountUnpairedPosition(i,v); CountUnpairedPositionEvidence(i,v); }
#else
#define ScoreExternalUnpairedEvidence(i) (ScoreUnpairedPosition(i) + ScoreUnpairedPositionEvidence(i))
#define CountExternalUnpairedEvidence(i,v) { CountUnpairedPosition(i,v); CountUnpairedPositionEvidence(i,v); }
#endif

// score for each unpaired position in an external loop

#if PARAMS_EXTERNAL_LENGTH
#define ScoreExternalUnpaired(i) (score_external_unpaired.first + ScoreUnpairedPosition(i))
#define CountExternalUnpaired(i,v) { score_external_unpaired.second += (v); CountUnpairedPosition(i,v); }
#else
#define ScoreExternalUnpaired(i) (ScoreUnpairedPosition(i))
#define CountExternalUnpaired(i,v) { CountUnpairedPosition(i,v); }
#endif

// score for a helix stacking pair of the form:
//
//       |         |
//    s[i+1] == s[j-1]
//       |         |
//     s[i] ==== s[j]
//       |         |

#if PARAMS_HELIX_STACKING
#if PROFILE
#define ScoreHelixStacking(i,j) profile_score_helix_stacking[i*(L+1)+j].first
#define CountHelixStacking(i,j,v) { profile_score_helix_stacking[i*(L+1)+j].second += (v); }
#else
#define ScoreHelixStacking(i,j) score_helix_stacking[s[i]][s[j]][s[i+1]][s[j-1]].first
#define CountHelixStacking(i,j,v) { score_helix_stacking[s[i]][s[j]][s[i+1]][s[j-1]].second += (v); }
#endif
#else
#define ScoreHelixStacking(i,j) RealT(0)
#define CountHelixStacking(i,j,v)
#endif

//////////////////////////////////////////////////////////////////////
// UPDATE_MAX()
//
// Macro for updating a score/traceback pointer which does not
// evaluate t unless an update is needed.  Make sure that this is
// used as a stand-alone statement (i.e., not the "if" condition
// of an if-then-else statement.)
//////////////////////////////////////////////////////////////////////

#define UPDATE_MAX(bs,bt,s,t) { RealT work(s); if ((work)>(bs)) { (bs)=(work); (bt)=(t); } }

//////////////////////////////////////////////////////////////////////
// FillScores()
// FillCounts()
// 
// Routines for setting scores and counts quickly. 
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::FillScores(typename std::vector<std::pair<RealT, RealT> >::iterator begin, typename std::vector<std::pair<RealT, RealT> >::iterator end, RealT value)
{
    while (begin != end)
    {
        begin->first = value;
        ++begin;
    }
}

template<class RealT>
void InferenceEngine<RealT>::FillCounts(typename std::vector<std::pair<RealT,RealT> >::iterator begin, typename std::vector<std::pair<RealT,RealT> >::iterator end, RealT value)
{
    while (begin != end)
    {
        begin->second = value;
        ++begin;
    }
}

//////////////////////////////////////////////////////////////////////
// ComputeRowOffset()
//
// Consider an N x N upper triangular matrix whose elements are
// stored in a one-dimensional flat array using the following
// row-major indexing scheme:
//
//     0  1  2  3     <-- row 0
//        4  5  6     <-- row 1
//           7 [8]    <-- row 2
//              9     <-- row 3
//
// Assuming 0-based indexing, this function computes offset[i]
// for the ith row such that offset[i]+j is the index of the
// (i,j)th element of the upper triangular matrix in the flat
// array.
//
// For example, offset[2] = 5, so the (2,3)th element of the
// upper triangular matrix (marked in the picture above) can be 
// found at position offset[2]+3 = 5+3 = 8 in the flat array.
//////////////////////////////////////////////////////////////////////

template<class RealT>
int InferenceEngine<RealT>::ComputeRowOffset(int i, int N) const
{
    Assert(i >= 0 && i <= N, "Index out-of-bounds.");
    
    // equivalent to:
    //   return N*(N+1)/2 - (N-i)*(N-i+1)/2 - i;
    return i*(N+N-i-1)/2;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::IsComplementary()
//
// Determine if a pair of positions is considered "complementary."
//////////////////////////////////////////////////////////////////////

template<class RealT>
bool InferenceEngine<RealT>::IsComplementary(int i, int j) const
{
    Assert(1 <= i && i <= L, "Index out-of-bounds.");
    Assert(1 <= j && j <= L, "Index out-of-bounds.");

#if !PROFILE
    return is_complementary[s[i]][s[j]];
#else
    RealT complementary_weight = 0;
    RealT total_weight = 0;

    for (int k = 0; k < N; k++)
    {
        if (is_complementary[A[k*(L+1)+i]][A[k*(L+1)+j]]) complementary_weight += weights[k];
        total_weight += weights[k];
    }

    return complementary_weight / total_weight >= std::min(RealT(N-1) / RealT(N), RealT(0.5));
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::InferenceEngine()
//
// Constructor
//////////////////////////////////////////////////////////////////////

template<class RealT>
InferenceEngine<RealT>::InferenceEngine(bool allow_noncomplementary, const int num_data_sources, const double kappa) :
    allow_noncomplementary(allow_noncomplementary),
    cache_initialized(false),
    parameter_manager(NULL),
    num_data_sources(num_data_sources),
    kappa(kappa),
    L(0),
    SIZE(0)
#if PROFILE
    , N(0)
    , SIZE2(0)
#endif

{
    // precompute mapping from characters to index representation
    std::memset(char_mapping, BYTE(alphabet.size()), 256);
    for (size_t i = 0; i < alphabet.size(); i++)
    {
        char_mapping[BYTE(tolower(alphabet[i]))] = 
            char_mapping[BYTE(toupper(alphabet[i]))] = i;
    }
    
    // precompute complementary pairings
    for (int i = 0; i <= M; i++)
        for (int j = 0; j <= M; j++)
            is_complementary[i][j] = 0;
    
    is_complementary[char_mapping[BYTE('A')]][char_mapping[BYTE('U')]] = 
        is_complementary[char_mapping[BYTE('U')]][char_mapping[BYTE('A')]] = 
        is_complementary[char_mapping[BYTE('G')]][char_mapping[BYTE('U')]] = 
        is_complementary[char_mapping[BYTE('U')]][char_mapping[BYTE('G')]] = 
        is_complementary[char_mapping[BYTE('C')]][char_mapping[BYTE('G')]] = 
        is_complementary[char_mapping[BYTE('G')]][char_mapping[BYTE('C')]] = 1;

    log_score_evidence.resize(num_data_sources);

    score_unpaired_position.resize(num_data_sources);
    score_paired_position.resize(num_data_sources);
    score_unpaired_position_raw.resize(num_data_sources);
    score_paired_position_raw.resize(num_data_sources);

    for (int i = 0; i < num_data_sources; i++)
    {
        log_score_evidence[i].resize(2);
        for (int j = 0; j < 2; j++)
        {
            log_score_evidence[i][j].resize(M);
            for (int k = 0; k < M; k++)
            {
                log_score_evidence[i][j][k].resize(2);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::~InferenceEngine()
//
// Destructor.
//////////////////////////////////////////////////////////////////////

template<class RealT>
InferenceEngine<RealT>::~InferenceEngine()
{}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::RegisterParameters()
//
// Establish a mapping between parameters in the inference
// engine and parameters in the parameter manager.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::RegisterParameters(ParameterManager<RealT> &parameter_manager)
{
    char buffer[1000];
    char buffer2[1000];

    cache_initialized = false;
    this->parameter_manager = &parameter_manager;
    parameter_manager.ClearParameters();
    
#if SINGLE_HYPERPARAMETER
    parameter_manager.AddParameterGroup("all_params");
#endif
    
#if PARAMS_BASE_PAIR
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("base_pair");
#endif
    for (int i = 0; i <= M; i++)
    {
        for (int j = 0; j <= M; j++)
        {
            if (i == M || j == M)
            {
                score_base_pair[i][j] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "base_pair_%c%c", alphabet[i], alphabet[j]);
                sprintf(buffer2, "base_pair_%c%c", alphabet[j], alphabet[i]);
                if (strcmp(buffer, buffer2) < 0)
                    parameter_manager.AddParameterMapping(buffer, &score_base_pair[i][j]);
                else
                    parameter_manager.AddParameterMapping(buffer2, &score_base_pair[i][j]);
            }
        }
    }
#endif
    
#if PARAMS_BASE_PAIR_DIST
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("base_pair_dist_at_least");
#endif
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
    {
        sprintf(buffer, "base_pair_dist_at_least_%d", BP_DIST_THRESHOLDS[i]);
        parameter_manager.AddParameterMapping(buffer, &score_base_pair_dist_at_least[i]);
    }
#endif
    
#if PARAMS_TERMINAL_MISMATCH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("terminal_mismatch");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int i2 = 0; i2 <= M; i2++)
            {
                for (int j2 = 0; j2 <= M; j2++)
                {
                    if (i1 == M || j1 == M || i2 == M || j2 == M)
                    {
                        score_terminal_mismatch[i1][j1][i2][j2] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "terminal_mismatch_%c%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2], alphabet[j2]);
                        parameter_manager.AddParameterMapping(buffer, &score_terminal_mismatch[i1][j1][i2][j2]);
                    }
                }
            }
        }
    }
#endif
    
#if PARAMS_HAIRPIN_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("hairpin_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
    {
        sprintf(buffer, "hairpin_length_at_least_%d", i);
        parameter_manager.AddParameterMapping(buffer, &score_hairpin_length_at_least[i]);
    }
#endif
    
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("hairpin_3_nucleotides");  
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                if (i1 == M || i2 == M || i3 == M)
                {
                    score_hairpin_3_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "hairpin_3_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_hairpin_3_nucleotides[i1][i2][i3]);
                }
            }
        }
    }
#endif
    
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("hairpin_4_nucleotides");  
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                for (int i4 = 0; i4 <= M; i4++)
                {
                    if (i1 == M || i2 == M || i3 == M || i4 == M)
                    {
                        score_hairpin_4_nucleotides[i1][i2][i3][i4] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "hairpin_4_nucleotides_%c%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3], alphabet[i4]);
                        parameter_manager.AddParameterMapping(buffer, &score_hairpin_4_nucleotides[i1][i2][i3][i4]);
                    }
                }
            }
        }
    }
#endif
    
#if PARAMS_HELIX_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("helix_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
    {
        if (i < 3)
        {
            score_helix_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "helix_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_helix_length_at_least[i]);
        }
    }
#endif
    
#if PARAMS_ISOLATED_BASE_PAIR
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("isolated_base_pair");
#endif
    parameter_manager.AddParameterMapping("isolated_base_pair", &score_isolated_base_pair);
#endif
    
#if PARAMS_INTERNAL_EXPLICIT
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_explicit");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_EXPLICIT_LENGTH; i++)
    {
        for (int j = 0; j <= D_MAX_INTERNAL_EXPLICIT_LENGTH; j++)
        {
            if (i == 0 || j == 0)
            {
                score_internal_explicit[i][j] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "internal_explicit_%d_%d", std::min(i, j), std::max(i, j));
                parameter_manager.AddParameterMapping(buffer, &score_internal_explicit[i][j]);
            }
        }
    }
#endif

#if PARAMS_BULGE_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_BULGE_LENGTH; i++)
    {
        if (i == 0)
        {
            score_bulge_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "bulge_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_bulge_length_at_least[i]);
        }
    }
#endif
    
#if PARAMS_INTERNAL_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_LENGTH; i++)
    {
        if (i < 2)
        {
            score_internal_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "internal_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_internal_length_at_least[i]);
        }
    }
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_symmetric_length_at_least");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
    {
        if (i == 0)
        {
            score_internal_symmetric_length_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "internal_symmetric_length_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_internal_symmetric_length_at_least[i]);
        }
    }
#endif

#if PARAMS_INTERNAL_ASYMMETRY
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_asymmetry_at_least");
#endif
    for (int i = 0; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
    {
        if (i == 0)
        {
            score_internal_asymmetry_at_least[i] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "internal_asymmetry_at_least_%d", i);
            parameter_manager.AddParameterMapping(buffer, &score_internal_asymmetry_at_least[i]);
        }
    }
#endif

#if PARAMS_BULGE_0x1_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_0x1_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        if (i1 == M)
        {
            score_bulge_0x1_nucleotides[i1] = std::pair<RealT,RealT>(0, 0);
            score_bulge_1x0_nucleotides[i1] = std::pair<RealT,RealT>(0, 0);
        }
        else
        {
            sprintf(buffer, "bulge_0x1_nucleotides_%c", alphabet[i1]);
            parameter_manager.AddParameterMapping(buffer, &score_bulge_0x1_nucleotides[i1]);
            parameter_manager.AddParameterMapping(buffer, &score_bulge_1x0_nucleotides[i1]);
        }
    }
#endif

#if PARAMS_BULGE_0x2_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_0x2_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            if (i1 == M || i2 == M)
            {
                score_bulge_0x2_nucleotides[i1][i2] = std::pair<RealT,RealT>(0, 0);
                score_bulge_2x0_nucleotides[i1][i2] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "bulge_0x2_nucleotides_%c%c", alphabet[i1], alphabet[i2]);
                parameter_manager.AddParameterMapping(buffer, &score_bulge_0x2_nucleotides[i1][i2]);
                parameter_manager.AddParameterMapping(buffer, &score_bulge_2x0_nucleotides[i1][i2]);
            }
        }
    }
#endif
    
#if PARAMS_BULGE_0x3_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("bulge_0x3_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                if (i1 == M || i2 == M)
                {
                    score_bulge_0x3_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                    score_bulge_3x0_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "bulge_0x3_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_bulge_0x3_nucleotides[i1][i2][i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_bulge_3x0_nucleotides[i1][i2][i3]);
                }
            }
        }
    }
#endif

#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_1x1_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            if (i1 == M || i2 == M)
            {
                score_internal_1x1_nucleotides[i1][i2] = std::pair<RealT,RealT>(0, 0);
            }
            else 
            {          
                sprintf(buffer, "internal_1x1_nucleotides_%c%c", alphabet[i1], alphabet[i2]);
                sprintf(buffer2, "internal_1x1_nucleotides_%c%c", alphabet[i2], alphabet[i1]);
                if (strcmp(buffer, buffer2) < 0)
                    parameter_manager.AddParameterMapping(buffer, &score_internal_1x1_nucleotides[i1][i2]);
                else
                    parameter_manager.AddParameterMapping(buffer2, &score_internal_1x1_nucleotides[i1][i2]);
            }
        }
    }
#endif
  
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_1x2_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                if (i1 == M || i2 == M || i3 == M)
                {
                    score_internal_1x2_nucleotides[i1][i2][i3] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "internal_1x2_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_internal_1x2_nucleotides[i1][i2][i3]);
                    sprintf(buffer, "internal_2x1_nucleotides_%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3]);
                    parameter_manager.AddParameterMapping(buffer, &score_internal_2x1_nucleotides[i1][i2][i3]);
                }
            }
        }
    }
#endif
  
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("internal_2x2_nucleotides");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int i2 = 0; i2 <= M; i2++)
        {
            for (int i3 = 0; i3 <= M; i3++)
            {
                for (int i4 = 0; i4 <= M; i4++)
                {
                    if (i1 == M || i2 == M || i3 == M || i4 == M)
                    {
                        score_internal_2x2_nucleotides[i1][i2][i3][i4] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "internal_2x2_nucleotides_%c%c%c%c", alphabet[i1], alphabet[i2], alphabet[i3], alphabet[i4]);
                        sprintf(buffer2, "internal_2x2_nucleotides_%c%c%c%c", alphabet[i3], alphabet[i4], alphabet[i1], alphabet[i2]);
                        if (strcmp(buffer, buffer2) < 0)
                            parameter_manager.AddParameterMapping(buffer, &score_internal_2x2_nucleotides[i1][i2][i3][i4]);
                        else
                            parameter_manager.AddParameterMapping(buffer2, &score_internal_2x2_nucleotides[i1][i2][i3][i4]);
                    }
                }
            }
        }
    }
#endif

#if PARAMS_HELIX_STACKING
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("helix_stacking");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int i2 = 0; i2 <= M; i2++)
            {
                for (int j2 = 0; j2 <= M; j2++)
                {
                    if (i1 == M || j1 == M || i2 == M || j2 == M)
                    {
                        score_helix_stacking[i1][j1][i2][j2] = std::pair<RealT,RealT>(0, 0);
                    }
                    else
                    {
                        sprintf(buffer, "helix_stacking_%c%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2], alphabet[j2]);
                        sprintf(buffer2, "helix_stacking_%c%c%c%c", alphabet[j2], alphabet[i2], alphabet[j1], alphabet[i1]);
                        if (strcmp(buffer, buffer2) < 0)
                            parameter_manager.AddParameterMapping(buffer, &score_helix_stacking[i1][j1][i2][j2]);
                        else
                            parameter_manager.AddParameterMapping(buffer2, &score_helix_stacking[i1][j1][i2][j2]);
                    }
                }
            }
        }
    }
#endif

#if PARAMS_HELIX_CLOSING
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("helix_closing");
#endif
    for (int i = 0; i <= M; i++)
    {
        for (int j = 0; j <= M; j++)
        {
            if (i == M || j == M)
            {
                score_helix_closing[i][j] = std::pair<RealT,RealT>(0, 0);
            }
            else
            {
                sprintf(buffer, "helix_closing_%c%c", alphabet[i], alphabet[j]);
                parameter_manager.AddParameterMapping(buffer, &score_helix_closing[i][j]);
            }
        }
    }
#endif

#if PARAMS_MULTI_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("multi_length");
#endif
    parameter_manager.AddParameterMapping("multi_base", &score_multi_base);
    parameter_manager.AddParameterMapping("multi_unpaired", &score_multi_unpaired);
    parameter_manager.AddParameterMapping("multi_paired", &score_multi_paired);
#endif

#if PARAMS_DANGLE
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("dangle");
#endif
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int i2 = 0; i2 <= M; i2++)
            {
                if (i1 == M || j1 == M || i2 == M)
                {
                    score_dangle_left[i1][j1][i2] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "dangle_left_%c%c%c", alphabet[i1], alphabet[j1], alphabet[i2]);
                    parameter_manager.AddParameterMapping(buffer, &score_dangle_left[i1][j1][i2]);
                }
            }
        }
    }
  
    for (int i1 = 0; i1 <= M; i1++)
    {
        for (int j1 = 0; j1 <= M; j1++)
        {
            for (int j2 = 0; j2 <= M; j2++)
            {
                if (i1 == M || j1 == M || j2 == M)
                {
                    score_dangle_right[i1][j1][j2] = std::pair<RealT,RealT>(0, 0);
                }
                else
                {
                    sprintf(buffer, "dangle_right_%c%c%c", alphabet[i1], alphabet[j1], alphabet[j2]);
                    parameter_manager.AddParameterMapping(buffer, &score_dangle_right[i1][j1][j2]);
                }
            }
        }
    }
#endif

#if PARAMS_EXTERNAL_LENGTH
#if MULTIPLE_HYPERPARAMETERS
    parameter_manager.AddParameterGroup("external_length");
#endif
    parameter_manager.AddParameterMapping("external_unpaired", &score_external_unpaired);
    parameter_manager.AddParameterMapping("external_paired", &score_external_paired);
#endif

#if PARAMS_EVIDENCE

#if SINGLE_HYPERPARAMETER
    parameter_manager.AddParameterGroup("evidence_cpds");
#endif

    for (int num_data_sources_current = 0; num_data_sources_current < num_data_sources; num_data_sources_current++)
    {
#if MULTIPLE_HYPERPARAMETERS
        sprintf(buffer,"%s%d","log_score_evidence",num_data_sources_current);
        parameter_manager.AddParameterGroup(buffer); //HKWS: comment this out???
#endif
    for (int i0 = 0; i0 < 2; i0++)
    {
        for (int i1 = 0; i1 < M; i1++)
        {
            for (int i2 = 0; i2 < 2; i2++)
            {
                if (i0 == 0)
                {
                    sprintf(buffer, "log_score_evidence%d_k_%c%d", num_data_sources_current,alphabet[i1],i2);
                }
                else
                {
                    sprintf(buffer, "log_score_evidence%d_theta_%c%d", num_data_sources_current,alphabet[i1],i2);
                }

                parameter_manager.AddParameterMapping(buffer, &log_score_evidence[num_data_sources_current][i0][i1][i2]);
            }
        }
    }
    }
#endif

}

template<class RealT>
double InferenceEngine<RealT>::LogGammaProb(RealT data, int which_data, int isUnpaired, int seq)
{
    RealT k = exp(log_score_evidence[which_data][0][seq][isUnpaired].first);
    RealT theta = exp(log_score_evidence[which_data][1][seq][isUnpaired].first);

    if (data < DATA_LOW_THRESH)
        return log(DATA_LOW_THRESH);

    return (k-1)*log(data) - data/theta - lgamma(k) - k*log(theta);
}

template<class RealT>
void InferenceEngine<RealT>::UpdateREVIVec(std::vector<RealT> perturb_unpaired, std::vector<RealT> perturb_paired)
{
    int which_data = 0; // update if we ever want to do more than one

    for (int i = 0; i < L; i++)
    {
        score_unpaired_position[which_data][i] += perturb_unpaired[i];
        //score_paired_position[which_data][i] += perturb_paired[i];
    }

    std::cerr << "u_v " << score_unpaired_position[0] << std::endl;
    std::cerr << "p_v " << score_paired_position[0] << std::endl;

}

template<class RealT>
std::vector<std::vector<double> > InferenceEngine<RealT>::GetREVIvec_up()
{
    return score_unpaired_position;
}

template<class RealT>
std::vector<std::vector<double> > InferenceEngine<RealT>::GetREVIvec_pr()
{
    return score_paired_position;
}

template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::GetREVIError(std::vector<RealT> p_i)
{
    int which_data = 0;
    std::vector<RealT> error;

     for (int i = 0; i < L; i++)
     {

        // analytical posterior
    //     RealT pos_paired = 0;
    //     // get paired posterior
    //     for (int l = 1; l <= L; l++)
    //     {
    //             int offset1 = i <= l ? i : l;
    //             int offset2 = i > l ? i : l;
    //             if (i != l)
    //                 pos_paired += posterior[offset[offset1]+offset2];
    //     }
    //     // if (i==1){
    //     //     std::cerr << "pos_paired " << pos_paired << std::endl;
    //     // }

    // HWS: reintroduce when using reactivities instead of probabilities

    // RealT k_pr = exp(log_score_evidence[which_data][0][s[i+1]][0].first);
    // RealT theta_pr = exp(log_score_evidence[which_data][1][s[i+1]][0].first);
    // RealT k_unp = exp(log_score_evidence[which_data][0][s[i+1]][1].first);
    // RealT theta_unp = exp(log_score_evidence[which_data][1][s[i+1]][1].first);

        // predict reactivity under this model
        //RealT rhat = k_pr*theta_pr*(1-p_i[i]) + k_unp*theta_unp*p_i[i];
        //score_unpaired_position_raw: where raw reacitivty data is stored

        // fitting to raw diff in p(unp values)
        RealT err = (p_i[i] - score_unpaired_position_raw[which_data][i]);

        error.push_back(err);
}
    std::cerr << "sc_up_p_raw " << score_unpaired_position_raw[which_data] << std::endl;
    return error;
    
}

template<class RealT>
void InferenceEngine<RealT>::InitializeREVIVec()
{    
    score_unpaired_position[0].clear();
    score_paired_position[0].clear();
    for (int i = 0; i < num_data_sources; i++)
    {
        for (int j = 0; j < L; j++)
    {
    score_unpaired_position[i].push_back(RealT(0));
    score_paired_position[i].push_back(RealT(0));
    }
}
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::LoadSequence()
//
// Load an RNA sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::LoadSequence(const SStruct &sstruct)
{
    cache_initialized = false;
    
    // compute dimensions
    L = sstruct.GetLength();
    SIZE = (L+1)*(L+2) / 2;
#if PROFILE
    N = sstruct.GetNumSequences();
    SIZE2 = (L+1)*(L+1);
#endif

    
    // allocate memory
    s.resize(L+1);
#if PROFILE
    A.resize(N*(L+1));
    weights.resize(N);
#endif
    offset.resize(L+1);
    allow_unpaired_position.resize(L+1);
    allow_unpaired.resize(SIZE);
    allow_paired.resize(SIZE);
    loss_unpaired_position.resize(L+1);
    loss_unpaired.resize(SIZE);
    loss_paired.resize(SIZE);
        
#if PROFILE

#if PARAMS_BASE_PAIR
    profile_score_base_pair.clear();                 profile_score_base_pair.resize(SIZE2);
#endif
#if PARAMS_TERMINAL_MISMATCH
    profile_score_terminal_mismatch.clear();         profile_score_terminal_mismatch.resize(SIZE2);
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
    profile_score_hairpin_3_nucleotides.clear();     profile_score_hairpin_3_nucleotides.resize(L+1);
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
    profile_score_hairpin_4_nucleotides.clear();     profile_score_hairpin_4_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
    profile_score_bulge_0x1_nucleotides.clear();     profile_score_bulge_0x1_nucleotides.resize(L+1);
    profile_score_bulge_1x0_nucleotides.clear();     profile_score_bulge_1x0_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
    profile_score_bulge_0x2_nucleotides.clear();     profile_score_bulge_0x2_nucleotides.resize(L+1);
    profile_score_bulge_2x0_nucleotides.clear();     profile_score_bulge_2x0_nucleotides.resize(L+1);
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
    profile_score_bulge_0x3_nucleotides.clear();     profile_score_bulge_0x3_nucleotides.resize(L+1);
    profile_score_bulge_3x0_nucleotides.clear();     profile_score_bulge_3x0_nucleotides.resize(L+1);
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
    profile_score_internal_1x1_nucleotides.clear();  profile_score_internal_1x1_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
    profile_score_internal_1x2_nucleotides.clear();  profile_score_internal_1x2_nucleotides.resize(SIZE2);
    profile_score_internal_2x1_nucleotides.clear();  profile_score_internal_2x1_nucleotides.resize(SIZE2);
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
    profile_score_internal_2x2_nucleotides.clear();  profile_score_internal_2x2_nucleotides.resize(SIZE2);
#endif
#if PARAMS_HELIX_STACKING
    profile_score_helix_stacking.clear();            profile_score_helix_stacking.resize(SIZE2);
#endif
#if PARAMS_HELIX_CLOSING
    profile_score_helix_closing.clear();             profile_score_helix_closing.resize(SIZE2);
#endif
#if PARAMS_DANGLE
    profile_score_dangle_left.clear();               profile_score_dangle_left.resize(SIZE2);
    profile_score_dangle_right.clear();              profile_score_dangle_right.resize(SIZE2);
#endif

#endif

#if FAST_HELIX_LENGTHS
    cache_score_helix_sums.clear();                  cache_score_helix_sums.resize((2*L+1)*L);
#endif

    // convert sequences to index representation
    const std::string &sequence = sstruct.GetSequences()[0];
    s[0] = BYTE(alphabet.size());
    for (int i = 1; i <= L; i++)
    {
        s[i] = char_mapping[BYTE(sequence[i])];
    }

#if PROFILE
    const std::vector<std::string> &alignment = sstruct.GetSequences();
    for (int k = 0; k < N; k++)
    {
        A[k*(L+1)+0] = BYTE(alphabet.size());
        for (int i = 1; i <= L; i++)
        {
            A[k*(L+1)+i] = char_mapping[BYTE(alignment[k][i])];
        }
    }

    weights = ConvertVector<RealT>(sstruct.ComputePositionBasedSequenceWeights());
#endif
    
    // compute indexing scheme for upper triangular arrays;
    // also allow each position to be unpaired by default, and
    // set the loss for each unpaired position to zero
    for (int i = 0; i <= L; i++)
    {
        offset[i] = ComputeRowOffset(i,L+1);
        allow_unpaired_position[i] = 1;
        loss_unpaired_position[i] = RealT(0);
    }

    // allow all ranges to be unpaired, and all pairs of letters
    // to be paired; set the respective losses to zero    
    for (int i = 0; i < SIZE; i++)
    {
        allow_unpaired[i] = 1;
        allow_paired[i] = 1;
        loss_unpaired[i] = RealT(0);
        loss_paired[i] = RealT(0);
    }

    // prevent the non-letter before each sequence from pairing with anything;
    // also prevent each letter from pairing with itself
    for (int i = 0; i <= L; i++)
    {
        allow_paired[offset[0]+i] = 0;
        allow_paired[offset[i]+i] = 0;
    }

    // enforce complementarity of base-pairings
    if (!allow_noncomplementary)
    {
        // for each pair of non-complementary letters in the sequence, disallow the pairing
        for (int i = 1; i <= L; i++)
        {
            for (int j = i+1; j <= L; j++)
            {
                if (!IsComplementary(i,j))
                    allow_paired[offset[i]+j] = 0;
            }
        }
    }

#if PARAMS_EVIDENCE
    for (int i = 0; i < num_data_sources; i++)
    {
      score_unpaired_position_raw[i].clear();
      score_unpaired_position[i].clear();
    
      score_paired_position_raw[i].clear();
      score_paired_position[i].clear();

      // load raw data
      score_unpaired_position_raw[i] = sstruct.GetUnpairedPotentials(i);
      score_paired_position_raw[i] = sstruct.GetPairedPotentials(i);
    }
#else
    // load unpaired potentials
    score_unpaired_position.clear();
    score_unpaired_position = sstruct.GetUnpairedPotentials();

    score_paired_position.clear();
    score_paired_position = 1-sstruct.GetPairedPotentials();  // probability of paired = 1 - probability of unpaired (NOTE: currently PairedPotentials returns the probability of unpaired, so need the "1-" in front
#endif

}


template<class RealT>
void InferenceEngine<RealT>::UpdateEvidenceStructures()
{
    for (int i = 0; i < num_data_sources; i++)
        UpdateEvidenceStructures(i);
}

template<class RealT>
void InferenceEngine<RealT>::UpdateEvidenceStructures(int which_data)
{

    // NOTE: this must be done after s and parameters have been updated
    score_unpaired_position[which_data].clear();
    score_paired_position[which_data].clear();

    for (int i = 0; i < L; i++)
    {
        double score_unpaired = score_unpaired_position_raw[which_data][i] == SStruct::UNKNOWN_POTENTIAL ? 0 : LogGammaProb(score_unpaired_position_raw[which_data][i],which_data,1,s[i+1]);
        score_unpaired_position[which_data].push_back(score_unpaired);

        double score_paired = score_paired_position_raw[which_data][i] == SStruct::UNKNOWN_POTENTIAL ? 0 : LogGammaProb(score_paired_position_raw[which_data][i],which_data,0,s[i+1]);
        score_paired_position[which_data].push_back(score_paired);
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::InitializeCache()
//
// Initialize scoring cache prior to inference.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::InitializeCache()
{
    if (cache_initialized) return;
    cache_initialized = true;

    // initialize length and distance scoring
#if PARAMS_BASE_PAIR_DIST
    for (int j = 0; j <= BP_DIST_LAST_THRESHOLD; j++)
        cache_score_base_pair_dist[j].first = RealT(0);
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            cache_score_base_pair_dist[j].first += score_base_pair_dist_at_least[i].first;
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    cache_score_hairpin_length[0].first = score_hairpin_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_HAIRPIN_LENGTH; i++)
        cache_score_hairpin_length[i].first = cache_score_hairpin_length[i-1].first + score_hairpin_length_at_least[i].first;
#endif

#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[0].first = score_helix_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].first = cache_score_helix_length[i-1].first + score_helix_length_at_least[i].first;
#endif

#if PARAMS_BULGE_LENGTH
    RealT temp_cache_score_bulge_length[D_MAX_BULGE_LENGTH+1];
    temp_cache_score_bulge_length[0] = score_bulge_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_BULGE_LENGTH; i++)
        temp_cache_score_bulge_length[i] = temp_cache_score_bulge_length[i-1] + score_bulge_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_score_internal_length[D_MAX_INTERNAL_LENGTH+1];
    temp_cache_score_internal_length[0] = score_internal_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_LENGTH; i++)
        temp_cache_score_internal_length[i] = temp_cache_score_internal_length[i-1] + score_internal_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_score_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    temp_cache_score_internal_symmetric_length[0] = score_internal_symmetric_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        temp_cache_score_internal_symmetric_length[i] = temp_cache_score_internal_symmetric_length[i-1] + score_internal_symmetric_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_score_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    temp_cache_score_internal_asymmetry[0] = score_internal_asymmetry_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        temp_cache_score_internal_asymmetry[i] = temp_cache_score_internal_asymmetry[i-1] + score_internal_asymmetry_at_least[i].first;
#endif
    
    // precompute score for single-branch loops of length l1 and l2
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            cache_score_single[l1][l2].first = RealT(0);

            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;

            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)];
#endif
            }

            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    cache_score_single[l1][l2].first += score_internal_explicit[l1][l2].first;
#endif
#if PARAMS_INTERNAL_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)];
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    cache_score_single[l1][l2].first += temp_cache_score_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)];
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                cache_score_single[l1][l2].first += temp_cache_score_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))];
#endif
            }
        }
    }
    
#if PROFILE
    // initialize counts for profile scoring
    for (int i = 0; i <= L; i++)
    {
        for (int j = 0; j <= L; j++)
        {
#if PARAMS_BASE_PAIR
            {
                const int pos[2] = {i, j};
                ComputeProfileScore(profile_score_base_pair[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_base_pair));
            }
#endif
#if PARAMS_TERMINAL_MISMATCH
            {
                const int pos[4] = {i, j+1, i+1, j};
                ComputeProfileScore(profile_score_terminal_mismatch[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_terminal_mismatch));
            }
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ComputeProfileScore(profile_score_hairpin_3_nucleotides[i].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_3_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ComputeProfileScore(profile_score_hairpin_4_nucleotides[i].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_4_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[1] = {j};
                ComputeProfileScore(profile_score_bulge_0x1_nucleotides[j].first, pos, 1, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x1_nucleotides));
            }
            if (j == 0)
            {
                const int pos[1] = {i+1};
                ComputeProfileScore(profile_score_bulge_1x0_nucleotides[i].first, pos, 1, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_1x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[2] = {j-1, j};
                ComputeProfileScore(profile_score_bulge_0x2_nucleotides[j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x2_nucleotides));
            }
            if (j == 0)
            {
                const int pos[2] = {i+1, i+2};
                ComputeProfileScore(profile_score_bulge_2x0_nucleotides[i].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_2x0_nucleotides));
            }
#endif            
#if PARAMS_BULGE_0x3_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[3] = {j-2, j-1, j};
                ComputeProfileScore(profile_score_bulge_0x3_nucleotides[j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x3_nucleotides));
            }
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ComputeProfileScore(profile_score_bulge_3x0_nucleotides[i].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_3x0_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
            {
                const int pos[2] = {i+1, j};
                ComputeProfileScore(profile_score_internal_1x1_nucleotides[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
            {
                const int pos[3] = {i+1, j-1, j};
                ComputeProfileScore(profile_score_internal_1x2_nucleotides[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x2_nucleotides));
            }
            {
                const int pos[3] = {i+1, i+2, j};
                ComputeProfileScore(profile_score_internal_2x1_nucleotides[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
            {
                const int pos[4] = {i+1, i+2, j-1, j};
                ComputeProfileScore(profile_score_internal_2x2_nucleotides[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x2_nucleotides));
            }
#endif     
#if PARAMS_HELIX_STACKING
            {
                const int pos[4] = {i, j, i+1, j-1};
                ComputeProfileScore(profile_score_helix_stacking[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_helix_stacking));
            }
#endif
#if PARAMS_HELIX_CLOSING
            {
                const int pos[2] = {i, j+1};
                ComputeProfileScore(profile_score_helix_closing[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_helix_closing));
            }
#endif
#if PARAMS_DANGLE
            {
                const int pos[3] = {i, j+1, i+1};
                ComputeProfileScore(profile_score_dangle_left[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_dangle_left));
            }
            {
                const int pos[3] = {i, j+1, j};
                ComputeProfileScore(profile_score_dangle_right[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_dangle_right));
            }
#endif
        }
    }

#endif

#if FAST_HELIX_LENGTHS
    // precompute helix partial sums
    FillScores(cache_score_helix_sums.begin(), cache_score_helix_sums.end(), RealT(0));
    for (int i = L; i >= 1; i--)
    {
        for (int j = i+3; j <= L; j++)
        {
            cache_score_helix_sums[(i+j)*L+j-i].first = cache_score_helix_sums[(i+j)*L+j-i-2].first;
            if (allow_paired[offset[i+1]+j-1])
            {
                cache_score_helix_sums[(i+j)*L+j-i].first += ScoreBasePair(i+1,j-1);
                if (allow_paired[offset[i]+j])
                    cache_score_helix_sums[(i+j)*L+j-i].first += ScoreHelixStacking(i,j);
            }
        }
    }
#endif

}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::InitializeCacheESS()
//
// This is similar to InitializeCache() but uses ScoreBasePairEvidence() instead.
// Initialize scoring cache prior to inference.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::InitializeCacheESS()
{
    if (cache_initialized) return;
    cache_initialized = true;
    
    // initialize length and distance scoring
#if PARAMS_BASE_PAIR_DIST
    for (int j = 0; j <= BP_DIST_LAST_THRESHOLD; j++)
        cache_score_base_pair_dist[j].first = RealT(0);
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            cache_score_base_pair_dist[j].first += score_base_pair_dist_at_least[i].first;
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    cache_score_hairpin_length[0].first = score_hairpin_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_HAIRPIN_LENGTH; i++)
        cache_score_hairpin_length[i].first = cache_score_hairpin_length[i-1].first + score_hairpin_length_at_least[i].first;
#endif
    
#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[0].first = score_helix_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].first = cache_score_helix_length[i-1].first + score_helix_length_at_least[i].first;
#endif
    
#if PARAMS_BULGE_LENGTH
    RealT temp_cache_score_bulge_length[D_MAX_BULGE_LENGTH+1];
    temp_cache_score_bulge_length[0] = score_bulge_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_BULGE_LENGTH; i++)
        temp_cache_score_bulge_length[i] = temp_cache_score_bulge_length[i-1] + score_bulge_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_score_internal_length[D_MAX_INTERNAL_LENGTH+1];
    temp_cache_score_internal_length[0] = score_internal_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_LENGTH; i++)
        temp_cache_score_internal_length[i] = temp_cache_score_internal_length[i-1] + score_internal_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_score_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    temp_cache_score_internal_symmetric_length[0] = score_internal_symmetric_length_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        temp_cache_score_internal_symmetric_length[i] = temp_cache_score_internal_symmetric_length[i-1] + score_internal_symmetric_length_at_least[i].first;
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_score_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    temp_cache_score_internal_asymmetry[0] = score_internal_asymmetry_at_least[0].first;
    for (int i = 1; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        temp_cache_score_internal_asymmetry[i] = temp_cache_score_internal_asymmetry[i-1] + score_internal_asymmetry_at_least[i].first;
#endif
    
    // precompute score for single-branch loops of length l1 and l2
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            cache_score_single[l1][l2].first = RealT(0);
            
            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;
            
            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)];
#endif
            }
            
            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    cache_score_single[l1][l2].first += score_internal_explicit[l1][l2].first;
#endif
#if PARAMS_INTERNAL_LENGTH
                cache_score_single[l1][l2].first += temp_cache_score_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)];
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    cache_score_single[l1][l2].first += temp_cache_score_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)];
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                cache_score_single[l1][l2].first += temp_cache_score_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))];
#endif
            }
        }
    }
    
#if PROFILE
    // initialize counts for profile scoring
    for (int i = 0; i <= L; i++)
    {
        for (int j = 0; j <= L; j++)
        {
#if PARAMS_BASE_PAIR
            {
                const int pos[2] = {i, j};
                ComputeProfileScore(profile_score_base_pair[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_base_pair));
            }
#endif
#if PARAMS_TERMINAL_MISMATCH
            {
                const int pos[4] = {i, j+1, i+1, j};
                ComputeProfileScore(profile_score_terminal_mismatch[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_terminal_mismatch));
            }
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ComputeProfileScore(profile_score_hairpin_3_nucleotides[i].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_3_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ComputeProfileScore(profile_score_hairpin_4_nucleotides[i].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_hairpin_4_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[1] = {j};
                ComputeProfileScore(profile_score_bulge_0x1_nucleotides[j].first, pos, 1, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x1_nucleotides));
            }
            if (j == 0)
            {
                const int pos[1] = {i+1};
                ComputeProfileScore(profile_score_bulge_1x0_nucleotides[i].first, pos, 1, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_1x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[2] = {j-1, j};
                ComputeProfileScore(profile_score_bulge_0x2_nucleotides[j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x2_nucleotides));
            }
            if (j == 0)
            {
                const int pos[2] = {i+1, i+2};
                ComputeProfileScore(profile_score_bulge_2x0_nucleotides[i].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_2x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[3] = {j-2, j-1, j};
                ComputeProfileScore(profile_score_bulge_0x3_nucleotides[j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_0x3_nucleotides));
            }
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ComputeProfileScore(profile_score_bulge_3x0_nucleotides[i].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_bulge_3x0_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
            {
                const int pos[2] = {i+1, j};
                ComputeProfileScore(profile_score_internal_1x1_nucleotides[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x1_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
            {
                const int pos[3] = {i+1, j-1, j};
                ComputeProfileScore(profile_score_internal_1x2_nucleotides[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_1x2_nucleotides));
            }
            {
                const int pos[3] = {i+1, i+2, j};
                ComputeProfileScore(profile_score_internal_2x1_nucleotides[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x1_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
            {
                const int pos[4] = {i+1, i+2, j-1, j};
                ComputeProfileScore(profile_score_internal_2x2_nucleotides[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_internal_2x2_nucleotides));
            }
#endif
#if PARAMS_HELIX_STACKING
            {
                const int pos[4] = {i, j, i+1, j-1};
                ComputeProfileScore(profile_score_helix_stacking[i*(L+1)+j].first, pos, 4, reinterpret_cast<std::pair<RealT,RealT> *>(score_helix_stacking));
            }
#endif
#if PARAMS_HELIX_CLOSING
            {
                const int pos[2] = {i, j+1};
                ComputeProfileScore(profile_score_helix_closing[i*(L+1)+j].first, pos, 2, reinterpret_cast<std::pair<RealT,RealT> *>(score_helix_closing));
            }
#endif
#if PARAMS_DANGLE
            {
                const int pos[3] = {i, j+1, i+1};
                ComputeProfileScore(profile_score_dangle_left[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_dangle_left));
            }
            {
                const int pos[3] = {i, j+1, j};
                ComputeProfileScore(profile_score_dangle_right[i*(L+1)+j].first, pos, 3, reinterpret_cast<std::pair<RealT,RealT> *>(score_dangle_right));
            }
#endif
        }
    }
    
#endif
    
#if FAST_HELIX_LENGTHS
    // precompute helix partial sums
    FillScores(cache_score_helix_sums.begin(), cache_score_helix_sums.end(), RealT(0));
    for (int i = L; i >= 1; i--)
    {
        for (int j = i+3; j <= L; j++)
        {
            cache_score_helix_sums[(i+j)*L+j-i].first = cache_score_helix_sums[(i+j)*L+j-i-2].first;
            if (allow_paired[offset[i+1]+j-1])
            {
                cache_score_helix_sums[(i+j)*L+j-i].first += ScoreBasePairEvidence(i+1,j-1);
                if (allow_paired[offset[i]+j])
                    cache_score_helix_sums[(i+j)*L+j-i].first += ScoreHelixStacking(i,j);
            }
        }
    }
#endif
    
}


#if PROFILE

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeProfileScore()
//
// Compute profile score for a single location.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeProfileScore(RealT &profile_score, const int *pos, int dimensions, std::pair<RealT,RealT> *table)
{
    profile_score = 0;

    // consider all sequences
    for (int k = 0; k < N; k++)
    {
        bool valid = true;
        int index = 0;
        int *seq = &A[k*(L+1)];
        
        // extract letters of the pattern for the current sequence
        for (int d = 0; valid && d < dimensions; d++)
        {
            if (pos[d] < 1 || pos[d] > L)
                valid = false;
            else
            {
                BYTE c = seq[pos[d]];
                if (c == BYTE(alphabet.size()))
                    valid = false;
                else
                    index = index * (M+1) + c;
            }
        }

        // add contribution of pattern to score
        if (valid) profile_score += weights[k] * table[index].first;
    }
}

#endif
            
//////////////////////////////////////////////////////////////////////
// InferenceEngine::LoadValues()
//
// Load parameter values.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::LoadValues(const std::vector<RealT> &values)
{

    if (values.size() != parameter_manager->GetNumLogicalParameters()) Error("Parameter size mismatch.");

    cache_initialized = false;
    for (size_t i = 0; i < values.size(); i++)
    {
        std::vector<std::pair<RealT,RealT> *> physical_parameters = parameter_manager->GetPhysicalParameters(i);
        for (size_t j = 0; j < physical_parameters.size(); j++) {
            physical_parameters[j]->first = values[i];
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetCounts()
//
// Return counts.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::GetCounts()
{
    std::vector<RealT> counts(parameter_manager->GetNumLogicalParameters());
    
    // clear counts for physical parameters
    for (size_t i = 0; i < parameter_manager->GetNumLogicalParameters(); i++)
    {
        std::vector<std::pair<RealT,RealT> *> physical_parameters = parameter_manager->GetPhysicalParameters(i);
        for (size_t j = 0; j < physical_parameters.size(); j++)
            counts[i] += physical_parameters[j]->second;
    }

    return counts;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ClearCounts()
//
// Set all counts to zero.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ClearCounts()
{
    // clear counts for physical parameters
    for (size_t i = 0; i < parameter_manager->GetNumLogicalParameters(); i++)
    {
        std::vector<std::pair<RealT,RealT> *> physical_parameters = parameter_manager->GetPhysicalParameters(i);
        for (size_t j = 0; j < physical_parameters.size(); j++)
            physical_parameters[j]->second = RealT(0);
    }

    // clear counts for cache
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i <= BP_DIST_LAST_THRESHOLD; i++)
        cache_score_base_pair_dist[i].second = RealT(0);
#endif

#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        cache_score_hairpin_length[i].second = RealT(0);
#endif

#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        cache_score_helix_length[i].second = RealT(0);
#endif

    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
        for (int l2 = 0; l2 <= C_MAX_SINGLE_LENGTH; l2++)
            cache_score_single[l1][l2].second = RealT(0);
    
#if FAST_HELIX_LENGTHS
    FillCounts(cache_score_helix_sums.begin(), cache_score_helix_sums.end(), RealT(0));
#endif

    // clear counts for profiles
#if PROFILE

#if PARAMS_BASE_PAIR
    FillCounts(profile_score_base_pair.begin(), profile_score_base_pair.end(), RealT(0));
#endif
#if PARAMS_TERMINAL_MISMATCH
    FillCounts(profile_score_terminal_mismatch.begin(), profile_score_terminal_mismatch.end(), RealT(0));
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
    FillCounts(profile_score_hairpin_3_nucleotides.begin(), profile_score_hairpin_3_nucleotides.end(), RealT(0));
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
    FillCounts(profile_score_hairpin_4_nucleotides.begin(), profile_score_hairpin_4_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x1_nucleotides.begin(), profile_score_bulge_0x1_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_1x0_nucleotides.begin(), profile_score_bulge_1x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x2_nucleotides.begin(), profile_score_bulge_0x2_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_2x0_nucleotides.begin(), profile_score_bulge_2x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
    FillCounts(profile_score_bulge_0x3_nucleotides.begin(), profile_score_bulge_0x3_nucleotides.end(), RealT(0));
    FillCounts(profile_score_bulge_3x0_nucleotides.begin(), profile_score_bulge_3x0_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
    FillCounts(profile_score_internal_1x1_nucleotides.begin(), profile_score_internal_1x1_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
    FillCounts(profile_score_internal_1x2_nucleotides.begin(), profile_score_internal_1x2_nucleotides.end(), RealT(0));
    FillCounts(profile_score_internal_2x1_nucleotides.begin(), profile_score_internal_2x1_nucleotides.end(), RealT(0));
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
    FillCounts(profile_score_internal_2x2_nucleotides.begin(), profile_score_internal_2x2_nucleotides.end(), RealT(0));
#endif
#if PARAMS_HELIX_STACKING
    FillCounts(profile_score_helix_stacking.begin(), profile_score_helix_stacking.end(), RealT(0));
#endif
#if PARAMS_HELIX_CLOSING
    FillCounts(profile_score_helix_closing.begin(), profile_score_helix_closing.end(), RealT(0));
#endif
#if PARAMS_DANGLE
    FillCounts(profile_score_dangle_left.begin(), profile_score_dangle_left.end(), RealT(0));
    FillCounts(profile_score_dangle_right.begin(), profile_score_dangle_right.end(), RealT(0));
#endif

#endif

}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::FinalizeCounts()
//
// Apply any needed transformations to counts.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::FinalizeCounts()
{
#if FAST_HELIX_LENGTHS

    // reverse helix partial sums    
    std::vector<std::pair<RealT,RealT> > reverse_sums(cache_score_helix_sums);
    
    for (int i = 1; i <= L; i++)
    {
        for (int j = L; j >= i+3; j--)
        {
            // the "if" conditions here can be omitted
            
            if (allow_paired[offset[i+1]+j-1])
            {
                CountBasePair(i+1,j-1,reverse_sums[(i+j)*L+j-i].second);
                if (allow_paired[offset[i]+j])
                {
                    CountHelixStacking(i,j,reverse_sums[(i+j)*L+j-i].second);
                }
                else
                {
                    Assert(Abs(double(reverse_sums[(i+j)*L+j-i].second)) < 1e-8, "Should be zero.");
                }
            }
            else
            {
                Assert(Abs(double(reverse_sums[(i+j)*L+j-i-2].second)) < 1e-8, "Should be zero.");
            }
            
            reverse_sums[(i+j)*L+j-i-2].second += reverse_sums[(i+j)*L+j-i].second;
        }
    }
#endif

    // perform transformations
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            score_base_pair_dist_at_least[i].second += cache_score_base_pair_dist[j].second;
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        for (int j = i; j <= D_MAX_HAIRPIN_LENGTH; j++)
            score_hairpin_length_at_least[i].second += cache_score_hairpin_length[j].second;
#endif
    
#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        for (int j = i; j <= D_MAX_HELIX_LENGTH; j++)
            score_helix_length_at_least[i].second += cache_score_helix_length[j].second;
#endif

    // allocate temporary storage
#if PARAMS_BULGE_LENGTH
    RealT temp_cache_counts_bulge_length[D_MAX_BULGE_LENGTH+1];
    std::fill(temp_cache_counts_bulge_length, temp_cache_counts_bulge_length + D_MAX_BULGE_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_counts_internal_length[D_MAX_INTERNAL_LENGTH+1];
    std::fill(temp_cache_counts_internal_length, temp_cache_counts_internal_length + D_MAX_INTERNAL_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_counts_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    std::fill(temp_cache_counts_internal_symmetric_length, temp_cache_counts_internal_symmetric_length + D_MAX_INTERNAL_SYMMETRIC_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_counts_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    std::fill(temp_cache_counts_internal_asymmetry, temp_cache_counts_internal_asymmetry + D_MAX_INTERNAL_ASYMMETRY+1, RealT(0));
#endif

    // compute contributions
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;

            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                temp_cache_counts_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
            }

            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    score_internal_explicit[l1][l2].second += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_LENGTH
                temp_cache_counts_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    temp_cache_counts_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                temp_cache_counts_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))] += cache_score_single[l1][l2].second;
#endif
            }
        }
    }

#if PARAMS_BULGE_LENGTH
    for (int i = 0; i <= D_MAX_BULGE_LENGTH; i++)
        for (int j = i; j <= D_MAX_BULGE_LENGTH; j++)
            score_bulge_length_at_least[i].second += temp_cache_counts_bulge_length[j];
#endif
    
#if PARAMS_INTERNAL_LENGTH
    for (int i = 0; i <= D_MAX_INTERNAL_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_LENGTH; j++)
            score_internal_length_at_least[i].second += temp_cache_counts_internal_length[j];
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; j++)
            score_internal_symmetric_length_at_least[i].second += temp_cache_counts_internal_symmetric_length[j];
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        for (int j = i; j <= D_MAX_INTERNAL_ASYMMETRY; j++)
            score_internal_asymmetry_at_least[i].second += temp_cache_counts_internal_asymmetry[j];
#endif

    // finalize profile counts
#if PROFILE
    for (int i = 0; i <= L; i++)
    {
        for (int j = 0; j <= L; j++)
        {
#if PARAMS_BASE_PAIR
            {
                const int pos[2] = {i, j};
                ConvertProfileCount(profile_score_base_pair[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_base_pair));
            }
#endif
#if PARAMS_TERMINAL_MISMATCH
            {
                const int pos[4] = {i, j+1, i+1, j};
                ConvertProfileCount(profile_score_terminal_mismatch[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_terminal_mismatch));
            }
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ConvertProfileCount(profile_score_hairpin_3_nucleotides[i].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_3_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ConvertProfileCount(profile_score_hairpin_4_nucleotides[i].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_4_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[1] = {j};
                ConvertProfileCount(profile_score_bulge_0x1_nucleotides[j].second, pos, 1, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x1_nucleotides));
            }
            if (j == 0)
            {
                const int pos[1] = {i+1};
                ConvertProfileCount(profile_score_bulge_1x0_nucleotides[i].second, pos, 1, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_1x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[2] = {j-1, j};
                ConvertProfileCount(profile_score_bulge_0x2_nucleotides[j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x2_nucleotides));
            }
            if (j == 0)
            {
                const int pos[2] = {i+1, i+2};
                ConvertProfileCount(profile_score_bulge_2x0_nucleotides[i].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_2x0_nucleotides));
            }
#endif            
#if PARAMS_BULGE_0x3_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[3] = {j-2, j-1, j};
                ConvertProfileCount(profile_score_bulge_0x3_nucleotides[j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x3_nucleotides));
            }
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ConvertProfileCount(profile_score_bulge_3x0_nucleotides[i].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_3x0_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
            {
                const int pos[2] = {i+1, j};
                ConvertProfileCount(profile_score_internal_1x1_nucleotides[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
            {
                const int pos[3] = {i+1, j-1, j};
                ConvertProfileCount(profile_score_internal_1x2_nucleotides[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x2_nucleotides));
            }
            {
                const int pos[3] = {i+1, i+2, j};
                ConvertProfileCount(profile_score_internal_2x1_nucleotides[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x1_nucleotides));
            }
#endif            
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
            {
                const int pos[4] = {i+1, i+2, j-1, j};
                ConvertProfileCount(profile_score_internal_2x2_nucleotides[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x2_nucleotides));
            }
#endif     
#if PARAMS_HELIX_STACKING
            {
                const int pos[4] = {i, j, i+1, j-1};
                ConvertProfileCount(profile_score_helix_stacking[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_helix_stacking));
            }
#endif
#if PARAMS_HELIX_CLOSING
            {
                const int pos[2] = {i, j+1};
                ConvertProfileCount(profile_score_helix_closing[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_helix_closing));
            }
#endif
#if PARAMS_DANGLE
            {
                const int pos[3] = {i, j+1, i+1};
                ConvertProfileCount(profile_score_dangle_left[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_dangle_left));
            }
            {
                const int pos[3] = {i, j+1, j};
                ConvertProfileCount(profile_score_dangle_right[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_dangle_right));
            }
#endif
        }
    }
#endif

}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::FinalizeCountsESS()
//
// This is similar to FinalizeCounts() but uses CountBasePairEvidence instead.
// Apply any needed transformations to counts.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::FinalizeCountsESS()
{
#if FAST_HELIX_LENGTHS
    
    // reverse helix partial sums
    std::vector<std::pair<RealT,RealT> > reverse_sums(cache_score_helix_sums);
    
    for (int i = 1; i <= L; i++)
    {
        for (int j = L; j >= i+3; j--)
        {
            // the "if" conditions here can be omitted
            
            if (allow_paired[offset[i+1]+j-1])
            {
                CountBasePairEvidence(i+1,j-1,reverse_sums[(i+j)*L+j-i].second);
                if (allow_paired[offset[i]+j])
                {
                    CountHelixStacking(i,j,reverse_sums[(i+j)*L+j-i].second);
                }
                else
                {
                    Assert(Abs(double(reverse_sums[(i+j)*L+j-i].second)) < 1e-8, "Should be zero.");
                }
            }
            else
            {
                Assert(Abs(double(reverse_sums[(i+j)*L+j-i-2].second)) < 1e-8, "Should be zero.");
            }
            
            reverse_sums[(i+j)*L+j-i-2].second += reverse_sums[(i+j)*L+j-i].second;
        }
    }
#endif
    
    // perform transformations
#if PARAMS_BASE_PAIR_DIST
    for (int i = 0; i < D_MAX_BP_DIST_THRESHOLDS; i++)
        for (int j = BP_DIST_THRESHOLDS[i]; j <= BP_DIST_LAST_THRESHOLD; j++)
            score_base_pair_dist_at_least[i].second += cache_score_base_pair_dist[j].second;
#endif
    
#if PARAMS_HAIRPIN_LENGTH
    for (int i = 0; i <= D_MAX_HAIRPIN_LENGTH; i++)
        for (int j = i; j <= D_MAX_HAIRPIN_LENGTH; j++)
            score_hairpin_length_at_least[i].second += cache_score_hairpin_length[j].second;
#endif
    
#if PARAMS_HELIX_LENGTH
    for (int i = 0; i <= D_MAX_HELIX_LENGTH; i++)
        for (int j = i; j <= D_MAX_HELIX_LENGTH; j++)
            score_helix_length_at_least[i].second += cache_score_helix_length[j].second;
#endif
    
    // allocate temporary storage
#if PARAMS_BULGE_LENGTH
    RealT temp_cache_counts_bulge_length[D_MAX_BULGE_LENGTH+1];
    std::fill(temp_cache_counts_bulge_length, temp_cache_counts_bulge_length + D_MAX_BULGE_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_LENGTH
    RealT temp_cache_counts_internal_length[D_MAX_INTERNAL_LENGTH+1];
    std::fill(temp_cache_counts_internal_length, temp_cache_counts_internal_length + D_MAX_INTERNAL_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    RealT temp_cache_counts_internal_symmetric_length[D_MAX_INTERNAL_SYMMETRIC_LENGTH+1];
    std::fill(temp_cache_counts_internal_symmetric_length, temp_cache_counts_internal_symmetric_length + D_MAX_INTERNAL_SYMMETRIC_LENGTH+1, RealT(0));
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    RealT temp_cache_counts_internal_asymmetry[D_MAX_INTERNAL_ASYMMETRY+1];
    std::fill(temp_cache_counts_internal_asymmetry, temp_cache_counts_internal_asymmetry + D_MAX_INTERNAL_ASYMMETRY+1, RealT(0));
#endif
    
    // compute contributions
    for (int l1 = 0; l1 <= C_MAX_SINGLE_LENGTH; l1++)
    {
        for (int l2 = 0; l1+l2 <= C_MAX_SINGLE_LENGTH; l2++)
        {
            // skip over stacking pairs
            if (l1 == 0 && l2 == 0) continue;
            
            // consider bulge loops
            if (l1 == 0 || l2 == 0)
            {
#if PARAMS_BULGE_LENGTH
                temp_cache_counts_bulge_length[std::min(D_MAX_BULGE_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
            }
            
            // consider internal loops
            else
            {
#if PARAMS_INTERNAL_EXPLICIT
                if (l1 <= D_MAX_INTERNAL_EXPLICIT_LENGTH && l2 <= D_MAX_INTERNAL_EXPLICIT_LENGTH)
                    score_internal_explicit[l1][l2].second += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_LENGTH
                temp_cache_counts_internal_length[std::min(D_MAX_INTERNAL_LENGTH, l1+l2)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_SYMMETRY
                if (l1 == l2)
                    temp_cache_counts_internal_symmetric_length[std::min(D_MAX_INTERNAL_SYMMETRIC_LENGTH, l1)] += cache_score_single[l1][l2].second;
#endif
#if PARAMS_INTERNAL_ASYMMETRY
                temp_cache_counts_internal_asymmetry[std::min(D_MAX_INTERNAL_ASYMMETRY, Abs(l1-l2))] += cache_score_single[l1][l2].second;
#endif
            }
        }
    }
    
#if PARAMS_BULGE_LENGTH
    for (int i = 0; i <= D_MAX_BULGE_LENGTH; i++)
        for (int j = i; j <= D_MAX_BULGE_LENGTH; j++)
            score_bulge_length_at_least[i].second += temp_cache_counts_bulge_length[j];
#endif
    
#if PARAMS_INTERNAL_LENGTH
    for (int i = 0; i <= D_MAX_INTERNAL_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_LENGTH; j++)
            score_internal_length_at_least[i].second += temp_cache_counts_internal_length[j];
#endif
    
#if PARAMS_INTERNAL_SYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; i++)
        for (int j = i; j <= D_MAX_INTERNAL_SYMMETRIC_LENGTH; j++)
            score_internal_symmetric_length_at_least[i].second += temp_cache_counts_internal_symmetric_length[j];
#endif
    
#if PARAMS_INTERNAL_ASYMMETRY
    for (int i = 0; i <= D_MAX_INTERNAL_ASYMMETRY; i++)
        for (int j = i; j <= D_MAX_INTERNAL_ASYMMETRY; j++)
            score_internal_asymmetry_at_least[i].second += temp_cache_counts_internal_asymmetry[j];
#endif
    
    // finalize profile counts
#if PROFILE
    for (int i = 0; i <= L; i++)
    {
        for (int j = 0; j <= L; j++)
        {
#if PARAMS_BASE_PAIR
            {
                const int pos[2] = {i, j};
                ConvertProfileCount(profile_score_base_pair[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_base_pair));
            }
#endif
#if PARAMS_TERMINAL_MISMATCH
            {
                const int pos[4] = {i, j+1, i+1, j};
                ConvertProfileCount(profile_score_terminal_mismatch[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_terminal_mismatch));
            }
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ConvertProfileCount(profile_score_hairpin_3_nucleotides[i].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_3_nucleotides));
            }
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
            if (j == 0)
            {
                const int pos[4] = {i+1, i+2, i+3, i+4};
                ConvertProfileCount(profile_score_hairpin_4_nucleotides[i].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_hairpin_4_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x1_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[1] = {j};
                ConvertProfileCount(profile_score_bulge_0x1_nucleotides[j].second, pos, 1, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x1_nucleotides));
            }
            if (j == 0)
            {
                const int pos[1] = {i+1};
                ConvertProfileCount(profile_score_bulge_1x0_nucleotides[i].second, pos, 1, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_1x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[2] = {j-1, j};
                ConvertProfileCount(profile_score_bulge_0x2_nucleotides[j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x2_nucleotides));
            }
            if (j == 0)
            {
                const int pos[2] = {i+1, i+2};
                ConvertProfileCount(profile_score_bulge_2x0_nucleotides[i].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_2x0_nucleotides));
            }
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
            if (i == 0)
            {
                const int pos[3] = {j-2, j-1, j};
                ConvertProfileCount(profile_score_bulge_0x3_nucleotides[j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_0x3_nucleotides));
            }
            if (j == 0)
            {
                const int pos[3] = {i+1, i+2, i+3};
                ConvertProfileCount(profile_score_bulge_3x0_nucleotides[i].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_bulge_3x0_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
            {
                const int pos[2] = {i+1, j};
                ConvertProfileCount(profile_score_internal_1x1_nucleotides[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x1_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
            {
                const int pos[3] = {i+1, j-1, j};
                ConvertProfileCount(profile_score_internal_1x2_nucleotides[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_1x2_nucleotides));
            }
            {
                const int pos[3] = {i+1, i+2, j};
                ConvertProfileCount(profile_score_internal_2x1_nucleotides[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x1_nucleotides));
            }
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
            {
                const int pos[4] = {i+1, i+2, j-1, j};
                ConvertProfileCount(profile_score_internal_2x2_nucleotides[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_internal_2x2_nucleotides));
            }
#endif
#if PARAMS_HELIX_STACKING
            {
                const int pos[4] = {i, j, i+1, j-1};
                ConvertProfileCount(profile_score_helix_stacking[i*(L+1)+j].second, pos, 4, reinterpret_cast<std::pair<RealT, RealT> *>(score_helix_stacking));
            }
#endif
#if PARAMS_HELIX_CLOSING
            {
                const int pos[2] = {i, j+1};
                ConvertProfileCount(profile_score_helix_closing[i*(L+1)+j].second, pos, 2, reinterpret_cast<std::pair<RealT, RealT> *>(score_helix_closing));
            }
#endif
#if PARAMS_DANGLE
            {
                const int pos[3] = {i, j+1, i+1};
                ConvertProfileCount(profile_score_dangle_left[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_dangle_left));
            }
            {
                const int pos[3] = {i, j+1, j};
                ConvertProfileCount(profile_score_dangle_right[i*(L+1)+j].second, pos, 3, reinterpret_cast<std::pair<RealT, RealT> *>(score_dangle_right));
            }
#endif
        }
    }
#endif
    
}


#if PROFILE

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ConvertProfileCount()
//
// Convert profile count for a single location.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ConvertProfileCount(const RealT &profile_score, const int *pos, int dimensions, std::pair<RealT,RealT> *table)
{
    // consider all sequences
    for (int k = 0; k < N; k++)
    {
        bool valid = true;
        int index = 0;
        int *seq = &A[k*(L+1)];
        
        // extract letters of the pattern for the current sequence
        for (int d = 0; valid && d < dimensions; d++)
        {
            if (pos[d] < 1 || pos[d] > L)
                valid = false;
            else
            {
                BYTE c = seq[pos[d]];
                if (c == BYTE(alphabet.size()))
                    valid = false;
                else
                    index = index * (M+1) + c;
            }
        }

        // add contribution of pattern to score
        if (valid) table[index].second += weights[k] * profile_score;
    }
}
            
#endif

//////////////////////////////////////////////////////////////////////
// InferenceEngine::UseLoss()
//
// Use per-position loss.  A loss is incurred if true_mapping[i] !=
// UNKNOWN && solution[i] != true_mapping[i].
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::UseLoss(const std::vector<int> &true_mapping, RealT example_loss)
{
    Assert(int(true_mapping.size()) == L+1, "Mapping of incorrect length!");
    cache_initialized = false;
    
    // compute number of pairings
    int num_pairings = 0;
    for (int i = 1; i <= L; i++)
        if (true_mapping[i] != SStruct::UNKNOWN && true_mapping[i] != SStruct::UNPAIRED)
            ++num_pairings;

    RealT per_position_loss = example_loss / RealT(num_pairings);
        
    // compute the penalty for each position that we declare to be unpaired
    for (int i = 0; i <= L; i++)
    {
        loss_unpaired_position[i] =
            ((i == 0 || true_mapping[i] == SStruct::UNKNOWN || true_mapping[i] == SStruct::UNPAIRED) ? RealT(0) : per_position_loss);
    }

    // now, compute the penalty for declaring ranges of positions to be unpaired;
    // also, compute the penalty for matching positions s[i] and s[j].
    for (int i = 0; i <= L; i++)
    {
        loss_unpaired[offset[i]+i] = RealT(0);
        loss_paired[offset[i]+i] = RealT(NEG_INF);
        for (int j = i+1; j <= L; j++)
        {
            loss_unpaired[offset[i]+j] = 
                loss_unpaired[offset[i]+j-1] +
                loss_unpaired_position[j];
            loss_paired[offset[i]+j] = 
                ((i == 0 || true_mapping[i] == SStruct::UNKNOWN || true_mapping[i] == SStruct::UNPAIRED || true_mapping[i] == j) ? RealT(0) : per_position_loss) +
                ((i == 0 || true_mapping[j] == SStruct::UNKNOWN || true_mapping[j] == SStruct::UNPAIRED || true_mapping[j] == i) ? RealT(0) : per_position_loss);
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::UseConstraints()
//
// Use known secondary structure mapping.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::UseConstraints(const std::vector<int> &true_mapping)
{
    Assert(int(true_mapping.size()) == L+1, "Supplied mapping of incorrect length!");
    cache_initialized = false;
    
    // determine whether we allow each position to be unpaired
    for (int i = 1; i <= L; i++)
    {
        allow_unpaired_position[i] =
            (true_mapping[i] == SStruct::UNKNOWN || 
             true_mapping[i] == SStruct::UNPAIRED);
    }

    // determine whether we allow ranges of positions to be unpaired;
    // also determine which base-pairings we allow
    for (int i = 0; i <= L; i++)
    {
        allow_unpaired[offset[i]+i] = 1;
        allow_paired[offset[i]+i] = 0;
        for (int j = i+1; j <= L; j++)
        {
            allow_unpaired[offset[i]+j] = 
                allow_unpaired[offset[i]+j-1] && 
                allow_unpaired_position[j];
            allow_paired[offset[i]+j] =
                (i > 0 &&
                 (true_mapping[i] == SStruct::UNKNOWN || true_mapping[i] == j) &&
                 (true_mapping[j] == SStruct::UNKNOWN || true_mapping[j] == i) &&
                 (allow_noncomplementary || IsComplementary(i,j)));
        }
    }
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreJunctionA()
// InferenceEngine::CountJunctionA()
//
// Returns the score for an asymmetric junction at positions i
// and j such that (i,j+1) are base-paired and (i+1,j) are free.
//
// In an RNA structure, this would look like
//
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
// Note that the difference between ScoreJunctionA() and
// ScoreJunctionB() is that the former applies to multi-branch
// loops whereas the latter is used for hairpin loops and
// single-branch loops.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreJunctionA(int i, int j) const
{
    // i and j must be bounded away from the edges so that s[i] and s[j+1]
    // refer to actual nucleotides.  To allow us to use this macro when
    // scoring the asymmetric junction for an exterior loop whose closing
    // base pair include the first and last nucleotides of the sequence,
    // we allow i to be as large as L and j to be as small as 0.
    
    Assert(0 < i && i <= L && 0 <= j && j < L, "Invalid indices.");

    return
        RealT(0)
#if PARAMS_HELIX_CLOSING
#if PROFILE
        + profile_score_helix_closing[i*(L+1)+j].first
#else                                          
        + score_helix_closing[s[i]][s[j+1]].first
#endif
#endif
#if PARAMS_DANGLE
#if PROFILE
        + (i < L ? profile_score_dangle_left[i*(L+1)+j].first : RealT(0))
        + (j > 0 ? profile_score_dangle_right[i*(L+1)+j].first : RealT(0))
#else
        + (i < L ? score_dangle_left[s[i]][s[j+1]][s[i+1]].first : RealT(0))
        + (j > 0 ? score_dangle_right[s[i]][s[j+1]][s[j]].first : RealT(0))
#endif
#endif
        ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountJunctionA(int i, int j, RealT value)
{
    Assert(0 < i && i <= L && 0 <= j && j < L, "Invalid indices.");
    
#if PARAMS_HELIX_CLOSING
#if PROFILE
    profile_score_helix_closing[i*(L+1)+j].second += value;
#else
    score_helix_closing[s[i]][s[j+1]].second += value;
#endif
#endif
#if PARAMS_DANGLE
#if PROFILE
    if (i < L) profile_score_dangle_left[i*(L+1)+j].second += value;
    if (j > 0) profile_score_dangle_right[i*(L+1)+j].second += value;
#else                                                               
    if (i < L) score_dangle_left[s[i]][s[j+1]][s[i+1]].second += value;
    if (j > 0) score_dangle_right[s[i]][s[j+1]][s[j]].second += value;
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreJunctionB()
// InferenceEngine::CountJunctionB()
//
// Returns the score for a symmetric junction at positions i
// and j such that (i,j+1) are base-paired and (i+1,j) are free.
//
// In an RNA structure, this would look like
//
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
// Note that the difference between ScoreJunctionA() and
// ScoreJunctionB() is that the former applies to multi-branch
// loops whereas the latter is used for hairpin loops and
// single-branch loops.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreJunctionB(int i, int j) const
{
    // The bounds here are similar to the asymmetric junction case, with
    // the main difference being that symmetric junctions are not allowed
    // for the exterior loop.  For this reason, i and j are bounded away
    // from the edges of the sequence (i.e., i < L && j > 0).
    
    Assert(0 < i && i < L && 0 < j && j < L, "Invalid indices.");
    
    return RealT(0)
#if PARAMS_HELIX_CLOSING
#if PROFILE
        + profile_score_helix_closing[i*(L+1)+j].first
#else
        + score_helix_closing[s[i]][s[j+1]].first
#endif
#endif
#if PARAMS_TERMINAL_MISMATCH
#if PROFILE
        + profile_score_terminal_mismatch[i*(L+1)+j].first
#else                                           
        + score_terminal_mismatch[s[i]][s[j+1]][s[i+1]][s[j]].first
#endif
#endif
        ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountJunctionB(int i, int j, RealT value)
{
    Assert(0 < i && i < L && 0 < j && j < L, "Invalid indices.");
    
#if PARAMS_HELIX_CLOSING
#if PROFILE
    profile_score_helix_closing[i*(L+1)+j].second += value;
#else
    score_helix_closing[s[i]][s[j+1]].second += value;
#endif
#endif
#if PARAMS_TERMINAL_MISMATCH
#if PROFILE
    profile_score_terminal_mismatch[i*(L+1)+j].second += value;
#else
    score_terminal_mismatch[s[i]][s[j+1]][s[i+1]][s[j]].second += value;
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreBasePair()
// InferenceEngine::CountBasePair()
//
// Returns the score for a base-pairing between letters i and j.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreBasePair(int i, int j) const
{

    // Clearly, i and j must refer to actual letters of the sequence,
    // and no letter may base-pair to itself.
    
    Assert(0 < i && i <= L && 0 < j && j <= L && i != j, "Invalid base-pair");
    
    return RealT(0)
#if defined(HAMMING_LOSS)
        + loss_paired[offset[i]+j]
#endif
#if PARAMS_BASE_PAIR
#if PROFILE
        + profile_score_base_pair[i*(L+1)+j].first
#else
        + score_base_pair[s[i]][s[j]].first
#endif
#endif
#if PARAMS_BASE_PAIR_DIST
        + cache_score_base_pair_dist[std::min(Abs(j - i), BP_DIST_LAST_THRESHOLD)].first
#endif
    ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountBasePair(int i, int j, RealT value)
{
    Assert(0 < i && i <= L && 0 < j && j <= L && i != j, "Invalid base-pair");
    
#if PARAMS_BASE_PAIR
#if PROFILE
    profile_score_base_pair[i*(L+1)+j].second += value;
#else
    score_base_pair[s[i]][s[j]].second += value;
#endif
#endif
#if PARAMS_BASE_PAIR_DIST
    cache_score_base_pair_dist[std::min(Abs(j - i), BP_DIST_LAST_THRESHOLD)].second += value;
#endif
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreBasePairEvidence()
// InferenceEngine::CountBasePairEvidence()
//
// Returns the score for a base-pairing between letters i and j.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreBasePairEvidence(int i, int j) const
{
    return ScoreBasePair(i,j) + ScorePairedPositionEvidence(i) + ScorePairedPositionEvidence(j);  
}

template<class RealT>
inline void InferenceEngine<RealT>::CountBasePairEvidence(int i, int j, RealT value)
{
    Assert(0 < i && i <= L && 0 < j && j <= L && i != j, "Invalid base-pair");

    CountBasePair(i,j,value);    
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreHairpinEvidence()
// InferenceEngine::CountHairpinEvidence()
//
// Similar to ScoreHairpin but based on evidence.
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHairpinEvidence(int i, int j) const
{
	return ScoreHairpin(i,j) + ScoreUnpairedEvidence(i,j);
}
template<class RealT>
inline void InferenceEngine<RealT>::CountHairpinEvidence(int i, int j, RealT value)
{  
	CountHairpin(i,j,value);
    CountUnpairedEvidence(i,j,value);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreHairpin()
// InferenceEngine::CountHairpin()
//
// Returns the score for a hairpin spanning positions i to j.
//
// In an RNA structure, this would look like
//
//                           ...
//                       /         \. 
//                   x[i+2]       x[j-1]
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHairpin(int i, int j) const
{
    // The constraints i > 0 && j < L ensure that s[i] and s[j+1] refer to
    // nucleotides which could base-pair.  The remaining constraint ensures
    // that only valid hairpins are considered.
    
    Assert(0 < i && i + C_MIN_HAIRPIN_LENGTH <= j && j < L, "Hairpin boundaries invalid.");
    
    return 
        ScoreUnpaired(i,j)
        + ScoreJunctionB(i,j)
#if PARAMS_HAIRPIN_LENGTH
        + cache_score_hairpin_length[std::min(j - i, D_MAX_HAIRPIN_LENGTH)].first
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
#if PROFILE
        + (j - i == 3 ? profile_score_hairpin_3_nucleotides[i].first : RealT(0))
#else
        + (j - i == 3 ? score_hairpin_3_nucleotides[s[i+1]][s[i+2]][s[i+3]].first : RealT(0))
#endif                                          
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
#if PROFILE
        + (j - i == 4 ? profile_score_hairpin_4_nucleotides[i].first : RealT(0))
#else
        + (j - i == 4 ? score_hairpin_4_nucleotides[s[i+1]][s[i+2]][s[i+3]][s[i+4]].first : RealT(0))
#endif
#endif
        ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountHairpin(int i, int j, RealT value)
{
    Assert(0 < i && i + C_MIN_HAIRPIN_LENGTH <= j && j < L, "Hairpin boundaries invalid.");
    
    CountUnpaired(i,j,value);
    CountJunctionB(i,j,value);
#if PARAMS_HAIRPIN_LENGTH
    cache_score_hairpin_length[std::min(j - i, D_MAX_HAIRPIN_LENGTH)].second += value;
#endif
#if PARAMS_HAIRPIN_3_NUCLEOTIDES
#if PROFILE
    if (j - i == 3) profile_score_hairpin_3_nucleotides[i].second += value;
#else
    if (j - i == 3) score_hairpin_3_nucleotides[s[i+1]][s[i+2]][s[i+3]].second += value;
#endif
#endif
#if PARAMS_HAIRPIN_4_NUCLEOTIDES
#if PROFILE
    if (j - i == 4) profile_score_hairpin_4_nucleotides[i].second += value;
#else
    if (j - i == 4) score_hairpin_4_nucleotides[s[i+1]][s[i+2]][s[i+3]][s[i+4]].second += value;
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreHelix()
// InferenceEngine::CountHelix()
//
// Returns the score for a helix of length m starting at positions
// i and j.  All base-pairs except for x[i+1]-x[j] are scored.
//
// In an RNA structure, this would look like
//
//                           ...
//                       \          /
// position i+m ------->  o        o  <----- position j-m
//                     x[i+3] -- x[j-2]
//                        |        |
//                     x[i+2] -- x[j-1]
//                        |        |
//                     x[i+1] -- x[j]
// position i --------->  o        o  <----- position j
//                       /          \.
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHelix(int i, int j, int m) const
{
    // First, i >= 0 && j <= L are obvious sanity-checks to make sure that
    // things are within range.  The check that i+2*m <= j ensures that there
    // are enough nucleotides to allow a helix of length m.
    
    Assert(0 <= i && i + 2 * m <= j && j <= L, "Helix boundaries invalid.");
    Assert(2 <= m && m <= D_MAX_HELIX_LENGTH, "Helix length invalid.");
    
#if FAST_HELIX_LENGTHS
    
    return
        cache_score_helix_sums[(i+j+1)*L+j-i-1].first - cache_score_helix_sums[(i+j+1)*L+j-i-m-m+1].first
#if PARAMS_HELIX_LENGTH
        + cache_score_helix_length[m].first
#endif
        ;
    
#else 
    
    RealT ret = RealT(0);
    for (int k = 1; k < m; k++)
        ret += ScoreHelixStacking(i+k,j-k+1) + ScoreBasePair(i+k+1,j-k);
    
#if PARAMS_HELIX_LENGTH
    ret += cache_score_helix_length[m].first;
#endif

    return ret;

#endif
  
}

template<class RealT>
inline void InferenceEngine<RealT>::CountHelix(int i, int j, int m, RealT value)
{
    Assert(0 <= i && i + 2 * m <= j && j <= L, "Helix boundaries invalid.");
    Assert(2 <= m && m <= D_MAX_HELIX_LENGTH, "Helix length invalid.");
    
#if FAST_HELIX_LENGTHS
    
    cache_score_helix_sums[(i+j+1)*L+j-i-1].second += value;
    cache_score_helix_sums[(i+j+1)*L+j-i-m-m+1].second -= value;
    
#else
    
    for (int k = 1; k < m; k++)
    {
        CountHelixStacking(i+k,j-k+1,value);
        CountBasePair(i+k+1,j-k,value);
    }
    
#endif
    
#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[m].second += value;
#endif
    
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreHelix()
// InferenceEngine::CountHelix()
//
// Returns the score for a helix of length m starting at positions
// i and j.  All base-pairs except for x[i+1]-x[j] are scored.
//
// In an RNA structure, this would look like
//
//                           ...
//                       \          /
// position i+m ------->  o        o  <----- position j-m
//                     x[i+3] -- x[j-2]
//                        |        |
//                     x[i+2] -- x[j-1]
//                        |        |
//                     x[i+1] -- x[j]
// position i --------->  o        o  <----- position j
//                       /          \.
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreHelixEvidence(int i, int j, int m) const
{
    // First, i >= 0 && j <= L are obvious sanity-checks to make sure that
    // things are within range.  The check that i+2*m <= j ensures that there
    // are enough nucleotides to allow a helix of length m.
    
    Assert(0 <= i && i + 2 * m <= j && j <= L, "Helix boundaries invalid.");
    Assert(2 <= m && m <= D_MAX_HELIX_LENGTH, "Helix length invalid.");
    
#if FAST_HELIX_LENGTHS
    
    return
    cache_score_helix_sums[(i+j+1)*L+j-i-1].first - cache_score_helix_sums[(i+j+1)*L+j-i-m-m+1].first
#if PARAMS_HELIX_LENGTH
    + cache_score_helix_length[m].first
#endif
    ;
    
#else
    
    RealT ret = RealT(0);
    for (int k = 1; k < m; k++)
        ret += ScoreHelixStacking(i+k,j-k+1) + ScoreBasePairEvidence(i+k+1,j-k);
    
#if PARAMS_HELIX_LENGTH
    ret += cache_score_helix_length[m].first;
#endif
    
    return ret;
    
#endif
    
}

template<class RealT>
inline void InferenceEngine<RealT>::CountHelixEvidence(int i, int j, int m, RealT value)
{
    Assert(0 <= i && i + 2 * m <= j && j <= L, "Helix boundaries invalid.");
    Assert(2 <= m && m <= D_MAX_HELIX_LENGTH, "Helix length invalid.");
    
#if FAST_HELIX_LENGTHS
    
    cache_score_helix_sums[(i+j+1)*L+j-i-1].second += value;
    cache_score_helix_sums[(i+j+1)*L+j-i-m-m+1].second -= value;
    
#else
    
    for (int k = 1; k < m; k++)
    {
        CountHelixStacking(i+k,j-k+1,value);
        CountBasePairEvidence(i+k+1,j-k,value);
    }
    
#endif
    
#if PARAMS_HELIX_LENGTH
    cache_score_helix_length[m].second += value;
#endif
    
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreSingleNucleotidesEvidence()
// InferenceEngine::CountSingleNucleotidesEvidence()
//
// Similar to ScoreSingleNucleotide but based on evidence
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreSingleNucleotidesEvidence(int i, int j, int p, int q) const
{
	return ScoreSingleNucleotides(i,j,p,q) + ScoreUnpairedEvidence(i,p) + ScoreUnpairedEvidence(q,j);
}
template<class RealT>
inline void InferenceEngine<RealT>::CountSingleNucleotidesEvidence(int i, int j, int p, int q, RealT value)
{
	CountSingleNucleotides(i,j,p,q,value);
    CountUnpairedEvidence(i,p,value);
    CountUnpairedEvidence(q,j,value);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreSingleNucleotides()
// InferenceEngine::CountSingleNucleotides()
//
// Returns the score for nucleotides in a single-branch loop 
// spanning i to j and p to q.
//
// In an RNA structure, this would look like
//
//                       ...      ...
//                        |        |
//                      x[p+1] -- x[q]
// position p -------->  o          o  <----- position q
//                    x[p]        x[q+1]
//                      |            |
//                     ...          ...
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreSingleNucleotides(int i, int j, int p, int q) const
{
    // Nucleotides s[i] and s[j+1] must exist, hence the conditions i > 0 and j < L.
    // the condition p+2 <= q comes from the fact that there must be enough room for
    // at least one nucleotide on the other side of the single-branch loop.  This
    // loop should only be used for dealing with single-branch loops, not stacking pairs.
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");

#if (!defined(NDEBUG) || PARAMS_BULGE_0x1_NUCLEOTIDES || PARAMS_BULGE_0x2_NUCLEOTIDES || PARAMS_BULGE_0x3_NUCLEOTIDES || PARAMS_INTERNAL_1x1_NUCLEOTIDES || PARAMS_INTERNAL_1x2_NUCLEOTIDES || PARAMS_INTERNAL_2x2_NUCLEOTIDES)
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
#endif
    
    return 
        ScoreUnpaired(i,p)
        + ScoreUnpaired(q,j)
#if PARAMS_BULGE_0x1_NUCLEOTIDES
#if PROFILE
        + (l1 == 0 && l2 == 1 ? profile_score_bulge_0x1_nucleotides[j].first : RealT(0))
        + (l1 == 1 && l2 == 0 ? profile_score_bulge_1x0_nucleotides[i].first : RealT(0))
#else
        + (l1 == 0 && l2 == 1 ? score_bulge_0x1_nucleotides[s[j]].first : RealT(0))
        + (l1 == 1 && l2 == 0 ? score_bulge_1x0_nucleotides[s[i+1]].first : RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
#if PROFILE
        + (l1 == 0 && l2 == 2 ? profile_score_bulge_0x2_nucleotides[j].first : RealT(0))
        + (l1 == 2 && l2 == 0 ? profile_score_bulge_2x0_nucleotides[i].first : RealT(0))
#else
        + (l1 == 0 && l2 == 2 ? score_bulge_0x2_nucleotides[s[j-1]][s[j]].first : RealT(0))
        + (l1 == 2 && l2 == 0 ? score_bulge_2x0_nucleotides[s[i+1]][s[i+2]].first : RealT(0))
#endif
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
#if PROFILE
        + (l1 == 0 && l2 == 3 ? profile_score_bulge_0x3_nucleotides[j].first : RealT(0))
        + (l1 == 3 && l2 == 0 ? profile_score_bulge_3x0_nucleotides[i].first : RealT(0))
#else
        + (l1 == 0 && l2 == 3 ? score_bulge_0x3_nucleotides[s[j-2]][s[j-1]][s[j]].first : RealT(0))
        + (l1 == 3 && l2 == 0 ? score_bulge_3x0_nucleotides[s[i+1]][s[i+2]][s[i+3]].first : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
#if PROFILE
        + (l1 == 1 && l2 == 1 ? profile_score_internal_1x1_nucleotides[i*(L+1)+j].first : RealT(0))
#else
        + (l1 == 1 && l2 == 1 ? score_internal_1x1_nucleotides[s[i+1]][s[j]].first : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
#if PROFILE
        + (l1 == 1 && l2 == 2 ? profile_score_internal_1x2_nucleotides[i*(L+1)+j].first : RealT(0))
        + (l1 == 2 && l2 == 1 ? profile_score_internal_2x1_nucleotides[i*(L+1)+j].first : RealT(0))
#else
        + (l1 == 1 && l2 == 2 ? score_internal_1x2_nucleotides[s[i+1]][s[j-1]][s[j]].first : RealT(0))
        + (l1 == 2 && l2 == 1 ? score_internal_2x1_nucleotides[s[i+1]][s[i+2]][s[j]].first : RealT(0))
#endif
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
#if PROFILE
        + (l1 == 2 && l2 == 2 ? profile_score_internal_2x2_nucleotides[i*(L+1)+j].first : RealT(0))
#else                                                                   
        + (l1 == 2 && l2 == 2 ? score_internal_2x2_nucleotides[s[i+1]][s[i+2]][s[j-1]][s[j]].first : RealT(0))
#endif
#endif
        ;
}

template<class RealT>
inline void InferenceEngine<RealT>::CountSingleNucleotides(int i, int j, int p, int q, RealT value)
{
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");

#if (!defined(NDEBUG) || PARAMS_BULGE_0x1_NUCLEOTIDES || PARAMS_BULGE_0x2_NUCLEOTIDES || PARAMS_BULGE_0x3_NUCLEOTIDES || PARAMS_INTERNAL_1x1_NUCLEOTIDES || PARAMS_INTERNAL_1x2_NUCLEOTIDES || PARAMS_INTERNAL_2x2_NUCLEOTIDES)
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
#endif
    
    CountUnpaired(i,p,value);
    CountUnpaired(q,j,value);
#if PARAMS_BULGE_0x1_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 1) profile_score_bulge_0x1_nucleotides[j].second += value;
    if (l1 == 1 && l2 == 0) profile_score_bulge_1x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 1) score_bulge_0x1_nucleotides[s[j]].second += value;
    if (l1 == 1 && l2 == 0) score_bulge_1x0_nucleotides[s[i+1]].second += value;
#endif
#endif
#if PARAMS_BULGE_0x2_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 2) profile_score_bulge_0x2_nucleotides[j].second += value;
    if (l1 == 2 && l2 == 0) profile_score_bulge_2x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 2) score_bulge_0x2_nucleotides[s[j-1]][s[j]].second += value;
    if (l1 == 2 && l2 == 0) score_bulge_2x0_nucleotides[s[i+1]][s[i+2]].second += value;
#endif
#endif
#if PARAMS_BULGE_0x3_NUCLEOTIDES
#if PROFILE
    if (l1 == 0 && l2 == 3) profile_score_bulge_0x3_nucleotides[j].second += value;
    if (l1 == 3 && l2 == 0) profile_score_bulge_3x0_nucleotides[i].second += value;
#else
    if (l1 == 0 && l2 == 3) score_bulge_0x3_nucleotides[s[j-2]][s[j-1]][s[j]].second += value;
    if (l1 == 3 && l2 == 0) score_bulge_3x0_nucleotides[s[i+1]][s[i+2]][s[i+3]].second += value;
#endif
#endif
#if PARAMS_INTERNAL_1x1_NUCLEOTIDES
#if PROFILE
    if (l1 == 1 && l2 == 1) profile_score_internal_1x1_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 1 && l2 == 1) score_internal_1x1_nucleotides[s[i+1]][s[j]].second += value;
#endif
#endif
#if PARAMS_INTERNAL_1x2_NUCLEOTIDES
#if PROFILE
    if (l1 == 1 && l2 == 2) profile_score_internal_1x2_nucleotides[i*(L+1)+j].second += value;
    if (l1 == 2 && l2 == 1) profile_score_internal_2x1_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 1 && l2 == 2) score_internal_1x2_nucleotides[s[i+1]][s[j-1]][s[j]].second += value;
    if (l1 == 2 && l2 == 1) score_internal_2x1_nucleotides[s[i+1]][s[i+2]][s[j]].second += value;
#endif    
#endif
#if PARAMS_INTERNAL_2x2_NUCLEOTIDES
#if PROFILE
    if (l1 == 2 && l2 == 2) profile_score_internal_2x2_nucleotides[i*(L+1)+j].second += value;
#else
    if (l1 == 2 && l2 == 2) score_internal_2x2_nucleotides[s[i+1]][s[i+2]][s[j-1]][s[j]].second += value;
#endif
#endif
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreSingleEvidence()
// InferenceEngine::CountSingleEvidence()
//
// Similar to ScoreSingle but based on evidence
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreSingleEvidence(int i, int j, int p, int q) const
{
    const int l1 = p - i;
    const int l2 = j - q;
    
    // Nucleotides s[i] and s[j+1] must exist, hence the conditions i > 0 and j < L.
    // the condition p+2 <= q comes from the fact that there must be enough room for
    // at least one nucleotide on the other side of the single-branch loop.  This
    // loop should only be used for dealing with single-branch loops, not stacking pairs.
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
    
    return 
        cache_score_single[l1][l2].first
        + ScoreBasePairEvidence(p+1,q)
        + ScoreJunctionB(i,j) 
        + ScoreJunctionB(q,p)
        + ScoreSingleNucleotidesEvidence(i,j,p,q);
}

template<class RealT>
inline void InferenceEngine<RealT>::CountSingleEvidence(int i, int j, int p, int q, RealT value)
{
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
    
    cache_score_single[l1][l2].second += value;
    CountBasePairEvidence(p+1,q,value);
    CountJunctionB(i,j,value);
    CountJunctionB(q,p,value);
    CountSingleNucleotidesEvidence(i,j,p,q,value);
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ScoreSingle()
// InferenceEngine::CountSingle()
//
// Returns the score for a single-branch loop spanning i to j and
// p to q.
//
// In an RNA structure, this would look like
//
//                       ...      ...
//                        |        |
//                      x[p+1] -- x[q]
// position p -------->  o          o  <----- position q
//                    x[p]        x[q+1]
//                      |            |
//                     ...          ...
//                      |            |
//                   x[i+1]        x[j]
// position i -------->  o          o  <----- position j
//                      x[i] -- x[j+1]
//                        |        |
//                     x[i-1] -- x[j+2]
//
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ScoreSingle(int i, int j, int p, int q) const
{
    const int l1 = p - i;
    const int l2 = j - q;
    
    // Nucleotides s[i] and s[j+1] must exist, hence the conditions i > 0 and j < L.
    // the condition p+2 <= q comes from the fact that there must be enough room for
    // at least one nucleotide on the other side of the single-branch loop.  This
    // loop should only be used for dealing with single-branch loops, not stacking pairs.
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
    
    return 
        cache_score_single[l1][l2].first
        + ScoreBasePair(p+1,q)
        + ScoreJunctionB(i,j) 
        + ScoreJunctionB(q,p)
        + ScoreSingleNucleotides(i,j,p,q);
}

template<class RealT>
inline void InferenceEngine<RealT>::CountSingle(int i, int j, int p, int q, RealT value)
{
    const int l1 = p - i;
    const int l2 = j - q;
    
    Assert(0 < i && i <= p && p + 2 <= q && q <= j && j < L, "Single-branch loop boundaries invalid.");
    Assert(l1 + l2 > 0 && l1 >= 0 && l2 >= 0 && l1 + l2 <= C_MAX_SINGLE_LENGTH, "Invalid single-branch loop size.");
    
    cache_score_single[l1][l2].second += value;
    CountBasePair(p+1,q,value);
    CountJunctionB(i,j,value);
    CountJunctionB(q,p,value);
    CountSingleNucleotides(i,j,p,q,value);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::EncodeTraceback()
// InferenceEngine::DecodeTraceback()
//
// Encode and decode traceback as an integer.  Here, i encodes
// a traceback type, and j encodes a length.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline int InferenceEngine<RealT>::EncodeTraceback(int i, int j) const
{
    Assert(0 <= i && i < NUM_TRACEBACK_TYPES && j >= 0, "Invalid values to encode as traceback.");
    return (j * NUM_TRACEBACK_TYPES) + i;
}

template<class RealT>
inline std::pair<int,int> InferenceEngine<RealT>::DecodeTraceback(int s) const
{
    return std::make_pair (s % NUM_TRACEBACK_TYPES, s / NUM_TRACEBACK_TYPES);
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeViterbi()
//
// Run Viterbi algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeViterbi()
{
    InitializeCache();
   
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif

#if CANDIDATE_LIST
    std::vector<int> candidates;
    candidates.reserve(L+1);
    long long int candidates_seen = 0;
    long long int candidates_possible = 0;
#endif
    
    // initialization

    F5t.clear(); F5t.resize(L+1, -1);
    FCt.clear(); FCt.resize(SIZE, -1);
    FMt.clear(); FMt.resize(SIZE, -1);
    FM1t.clear(); FM1t.resize(SIZE, -1);

    F5v.clear(); F5v.resize(L+1, RealT(NEG_INF));
    FCv.clear(); FCv.resize(SIZE, RealT(NEG_INF));
    FMv.clear(); FMv.resize(SIZE, RealT(NEG_INF));
    FM1v.clear(); FM1v.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEt.clear(); FEt.resize(SIZE, -1);
    FNt.clear(); FNt.resize(SIZE, -1);
    FEv.clear(); FEv.resize(SIZE, RealT(NEG_INF));
    FNv.clear(); FNv.resize(SIZE, RealT(NEG_INF));
#endif
    
    for (int i = L; i >= 0; i--)
    {
        
#if CANDIDATE_LIST
        candidates.clear();
#endif
        
        for (int j = i; j <= L; j++)
        {
            // FM2[i,j] = MAX (i<k<j : FM1[i,k] + FM[k,j])

            RealT FM2v = RealT(NEG_INF);
            int FM2t = -1;
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                UPDATE_MAX(FM2v, FM2t, FM1v[offset[i]+k] + FMv[offset[k]+j], k);
            
#else
            
#if !CANDIDATE_LIST
            
            if (i+2 <= j)
            {
                RealT *p1 = &(FM1v[offset[i]+i+1]);
                RealT *p2 = &(FMv[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    UPDATE_MAX(FM2v, FM2t, (*p1) + (*p2), k);
                    ++p1;
                    p2 += L-k;
                }
            }
            
#else
            
            for (register size_t kp = 0; kp < candidates.size(); kp++)
            {
                register const int k = candidates[kp];
                UPDATE_MAX(FM2v, FM2t, FM1v[offset[i]+k] + FMv[offset[k]+j], k);
            }
            
            candidates_seen += (long long int) candidates.size();
            candidates_possible += (long long int) std::max(j-i-1,0);
            
#endif
            
#endif

#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
      
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = MAX [ScoreHairpin(i,j),
            //                MAX (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + MAX (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    UPDATE_MAX(best_v, best_t, ScoreHairpin(i,j), EncodeTraceback(TB_FN_HAIRPIN,0));
                
                // compute MAX (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        UPDATE_MAX(best_v, best_t,
                                   ScoreSingle(i,j,p,q) + FCv[offset[p+1]+q-1],
                                   EncodeTraceback(TB_FN_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q));
                    }
                }
                
#else
                
                {
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    int bp = -1, bq = -1;
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCv[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT score = (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                                           ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            if (score > best_v)
                            {
                                best_v = score;
                                bp = p;
                                bq = q;
                            }
                        }
                    }
                    
                    if (bp != -1 && bq != -1)
                        best_t = EncodeTraceback(TB_FN_SINGLE,(bp-i)*(C_MAX_SINGLE_LENGTH+1)+j-bq);
                }
#endif
                
                // compute MAX (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                UPDATE_MAX(best_v, best_t,
                           FM2v + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase(), 
                           EncodeTraceback(TB_FN_BIFURCATION,FM2t));
                
                FNv[offset[i]+j] = best_v;
                FNt[offset[i]+j] = best_t;
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = MAX [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    UPDATE_MAX(best_v, best_t, 
                               ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEv[offset[i+1]+j-1],
                               EncodeTraceback(TB_FE_STACKING,0));
                }
                
                // compute FN(i,j)
                
                UPDATE_MAX(best_v, best_t, FNv[offset[i]+j], EncodeTraceback(TB_FE_FN,0));
                
                FEv[offset[i]+j] = best_v;
                FEt[offset[i]+j] = best_t;
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = MAX [ScoreIsolated() + FN(i,j),
            //                MAX (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreIsolated() + FN(i,j)
                
                UPDATE_MAX(best_v, best_t, ScoreIsolated() + FNv[offset[i]+j], EncodeTraceback(TB_FC_FN,0));
                
                // compute MAX (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    UPDATE_MAX(best_v, best_t, ScoreHelix(i-1,j+1,k) + FNv[offset[i+k-1]+j-k+1], EncodeTraceback(TB_FC_HELIX,k));
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        UPDATE_MAX(best_v, best_t, ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) +
                                   FEv[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1],
                                   EncodeTraceback(TB_FC_FE,0));
                }
                FCv[offset[i]+j] = best_v;
                FCt[offset[i]+j] = best_t;
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = MAX [ScoreHairpin(i,j),
            //                MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + MAX (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    UPDATE_MAX(best_v, best_t, ScoreHairpin(i,j), EncodeTraceback(TB_FC_HAIRPIN,0));
                
                // compute MAX (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])

#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        
                        UPDATE_MAX(best_v, best_t,
                                   FCv[offset[p+1]+q-1] +
                                   (p == i && q == j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingle(i,j,p,q)),
                                   EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q));
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    int bp = -1, bq = -1;
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCv[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            RealT score = (p == i && q == j) ?
                                (score_helix + FCptr[q]) :
                                (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            if (score > best_v)
                            {
                                best_v = score;
                                bp = p;
                                bq = q;
                            }
                        }
                    }
                    
                    if (bp != -1 && bq != -1)
                        best_t = EncodeTraceback(TB_FC_SINGLE,(bp-i)*(C_MAX_SINGLE_LENGTH+1)+j-bq);
                }
#endif
                
                // compute MAX (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                UPDATE_MAX(best_v, best_t,
                           FM2v + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase(), 
                           EncodeTraceback(TB_FC_BIFURCATION,FM2t));
                
                FCv[offset[i]+j] = best_v;
                FCt[offset[i]+j] = best_t;
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = MAX [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                {
                    UPDATE_MAX(best_v, best_t, 
                               FCv[offset[i+1]+j-1] + ScoreJunctionA(j,i) +
                               ScoreMultiPaired() + ScoreBasePair(i+1,j), 
                               EncodeTraceback(TB_FM1_PAIRED,0));
                }
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                {
                    UPDATE_MAX(best_v, best_t,
                               FM1v[offset[i+1]+j] + ScoreMultiUnpaired(i+1),
                               EncodeTraceback(TB_FM1_UNPAIRED,0));
                }
                
                FM1v[offset[i]+j] = best_v;
                FM1t[offset[i]+j] = best_t;
            }
            
#if CANDIDATE_LIST
            
            // If there exists some i <= k < j for which
            //   FM1[i,k] + FM[k,j] >= FM1[i,j]
            // then for all j' > j, we know that
            //   FM1[i,k] + FM[k,j'] >= FM1[i,j] + FM[j,j'].
            // since 
            //   FM[k,j'] >= FM[k,j] + FM[j,j'].
            //
            // From this, it follows that we only need to consider
            // j as a candidate partition point for future j' values
            // only if FM1[i,j] > FM1[i,k] + FM[k,j] for all k.
            
            if (FM1v[offset[i]+j] > FM2v)
                candidates.push_back(j);
#endif
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = MAX [MAX (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                RealT best_v = RealT(NEG_INF);
                int best_t = -1;
                
                // compute MAX (i<k<j : FM1[i,k] + FM[k,j])
                
                UPDATE_MAX(best_v, best_t, FM2v, EncodeTraceback(TB_FM_BIFURCATION,FM2t));
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                {
                    UPDATE_MAX(best_v, best_t,
                               FMv[offset[i]+j-1] + ScoreMultiUnpaired(j), 
                               EncodeTraceback(TB_FM_UNPAIRED,0));
                }
                
                // compute FM1[i,j]
                
                UPDATE_MAX(best_v, best_t, FM1v[offset[i]+j], EncodeTraceback(TB_FM_FM1,0));
                
                FMv[offset[i]+j] = best_v;
                FMt[offset[i]+j] = best_t;
            }
        }
    }
    
    F5v[0] = RealT(0);
    F5t[0] = EncodeTraceback(TB_F5_ZERO,0);
    for (int j = 1; j <= L; j++)
    {
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = MAX [F5[j-1] + ScoreExternalUnpaired(),
        //              MAX (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT best_v = RealT(NEG_INF);
        int best_t = -1;
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
        {
            UPDATE_MAX(best_v, best_t, 
                       F5v[j-1] + ScoreExternalUnpaired(j),
                       EncodeTraceback(TB_F5_UNPAIRED,0));
        }
        
        // compute MAX (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
            {
                UPDATE_MAX(best_v, best_t,
                           F5v[k] + FCv[offset[k+1]+j-1] + ScoreExternalPaired() +
                           ScoreBasePair(k+1,j) + ScoreJunctionA(j,k),
                           EncodeTraceback(TB_F5_BIFURCATION,k));
            }
        }
        
        F5v[j] = best_v;
        F5t[j] = best_t;
    }

#if SHOW_TIMINGS
    std::cerr << "Viterbi score: " << F5v[L] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
    
#if CANDIDATE_LIST
    //std::cerr << "Candidates: " << candidates_seen << "/" << candidates_possible << " = " << double(candidates_seen)/candidates_possible << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetViterbiScore()
//
// Return Viterbi score for a sequence.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::GetViterbiScore() const
{
    return F5v[L];
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::PredictPairingsViterbi()
// 
// Use Viterbi decoding to predict pairings.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsViterbi() const
{
    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;
    //return solution;
    
    std::queue<triple<const int *,int,int> > traceback_queue;
    traceback_queue.push(make_triple(&F5t[0], 0, L));
    
    while (!traceback_queue.empty())
    {
        triple<const int *,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int *V = t.first;
        const int i = t.second;
        const int j = t.third;
        
        std::pair<int,int> traceback = DecodeTraceback(V == &F5t[0] ? V[j] : V[offset[i]+j]);
        
        //std::cerr << (V == FCt ? "FC " : V == FMt ? "FM " : V == FM1t ? "FM1 " : "F5 ");
        //std::cerr << i << " " << j << ": " << traceback.first << " " << traceback.second << std::endl;
        
        switch (traceback.first)
        {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case TB_FN_HAIRPIN: 
                break;
            case TB_FN_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                solution[p+1] = q;
                solution[q] = p+1;
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FN_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FE_STACKING: 
            {
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(make_triple(&FEt[0], i+1, j-1));
            }
            break;
            case TB_FE_FN: 
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_FN:
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_HELIX:
            {
                const int m = traceback.second;
                for (int k = 2; k <= m; k++)
                {
                    solution[i+k-1] = j-k+2;
                    solution[j-k+2] = i+k-1;
                }
                traceback_queue.push(make_triple(&FNt[0], i+m-1, j-m+1));
            }
            break;
            case TB_FC_FE:
            {
                const int m = D_MAX_HELIX_LENGTH;
                for (int k = 2; k <= m; k++)
                {
                    solution[i+k-1] = j-k+2;
                    solution[j-k+2] = i+k-1;
                }
                traceback_queue.push(make_triple(&FEt[0], i+m-1, j-m+1));
            }
            break;
#else
            case TB_FC_HAIRPIN: 
                break;
            case TB_FC_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                solution[p+1] = q;
                solution[q] = p+1;
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FC_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
#endif
            case TB_FM1_PAIRED:
            {
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(make_triple(&FCt[0], i+1, j-1));
            }
      break;
            case TB_FM1_UNPAIRED:
            {
                traceback_queue.push(make_triple(&FM1t[0], i+1, j));
            }
            break;
            case TB_FM_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FM_UNPAIRED:
            {
                traceback_queue.push(make_triple(&FMt[0], i, j-1));
            }
            break;
            case TB_FM_FM1: 
            {
                traceback_queue.push(make_triple(&FM1t[0], i, j));
            }
            break;
            case TB_F5_ZERO:
                break;
            case TB_F5_UNPAIRED:
            {
                traceback_queue.push(make_triple(&F5t[0], 0, j-1));
            }
            break;
            case TB_F5_BIFURCATION:
            {
                const int k = traceback.second;
                solution[k+1] = j;
                solution[j] = k+1;
                traceback_queue.push(make_triple(&F5t[0], 0, k));
                traceback_queue.push(make_triple(&FCt[0], k+1, j-1));
            }
            break;
            default:
                Assert(false, "Bad traceback.");
        }
    }
    
    return solution;
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeViterbiFeatureCounts()
// 
// Use feature counts from Viterbi decoding.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::ComputeViterbiFeatureCounts()
{
    std::queue<triple<int *,int,int> > traceback_queue;
    traceback_queue.push(make_triple(&F5t[0], 0, L));

    ClearCounts();
    
    while (!traceback_queue.empty())
    {
        triple<int *,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int *V = t.first;
        const int i = t.second;
        const int j = t.third;
        
        std::pair<int,int> traceback = DecodeTraceback (V == &F5t[0] ? V[j] : V[offset[i]+j]);
        
        //std::cout << (V == FCt ? "FC " : V == FMt ? "FM " : V == FM1t ? "FM1 " : "F5 ");
        //std::cout << i << " " << j << ": " << traceback.first << " " << traceback.second << std::endl;
        
        switch (traceback.first)
        {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case TB_FN_HAIRPIN: 
                CountHairpin(i,j,RealT(1));
                break;
            case TB_FN_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                CountSingle(i,j,p,q,RealT(1));
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FN_BIFURCATION:
            {
                const int k = traceback.second;
                CountJunctionA(i,j,RealT(1));
                CountMultiPaired(RealT(1));
                CountMultiBase(RealT(1));
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FE_STACKING: 
            {
                CountBasePair(i+1,j,RealT(1));
                CountHelixStacking(i,j+1,RealT(1));
                traceback_queue.push(make_triple(&FEt[0], i+1, j-1));
            }
            break;
            case TB_FE_FN: 
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_FN:
            {
                CountIsolated(RealT(1));
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_HELIX:
            {
                const int m = traceback.second;
                CountHelix(i-1,j+1,m,RealT(1));
                traceback_queue.push(make_triple(&FNt[0], i+m-1, j-m+1));
            }
            break;
            case TB_FC_FE:
            {
                const int m = D_MAX_HELIX_LENGTH;
                CountHelix(i-1,j+1,m,RealT(1));
                traceback_queue.push(make_triple(&FEt[0], i+m-1, j-m+1));
            }
            break;
#else
            case TB_FC_HAIRPIN: 
                CountHairpin(i,j,RealT(1));
                break;
            case TB_FC_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);

                if (p == i && q == j)
                {
                    CountBasePair(i+1,j,RealT(1));
                    CountHelixStacking(i,j+1,RealT(1));
                }
                else
                {
                    CountSingle(i,j,p,q,RealT(1));
                }
                
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FC_BIFURCATION:
            {
                const int k = traceback.second;
                CountJunctionA(i,j,RealT(1));
                CountMultiPaired(RealT(1));
                CountMultiBase(RealT(1));
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
#endif
            case TB_FM1_PAIRED:
            {
                CountJunctionA(j,i,RealT(1));
                CountMultiPaired(RealT(1));
                CountBasePair(i+1,j,RealT(1));
                traceback_queue.push(make_triple(&FCt[0], i+1, j-1));
            }
            break;
            case TB_FM1_UNPAIRED:
            {
                CountMultiUnpaired(i+1,RealT(1));
                traceback_queue.push(make_triple(&FM1t[0], i+1, j));
            }
            break;
            case TB_FM_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FM_UNPAIRED:
            {
                CountMultiUnpaired(j,RealT(1));
                traceback_queue.push(make_triple(&FMt[0], i, j-1));
            }
            break;
            case TB_FM_FM1: 
                traceback_queue.push(make_triple(&FM1t[0], i, j));
                break;
            case TB_F5_ZERO:
                break;
            case TB_F5_UNPAIRED:
                CountExternalUnpaired(j,RealT(1));
                traceback_queue.push(make_triple(&F5t[0], 0, j-1));
                break;
            case TB_F5_BIFURCATION:
            {
                const int k = traceback.second;
                CountExternalPaired(RealT(1));
                CountBasePair(k+1,j,RealT(1));
                CountJunctionA(j,k,RealT(1));
                traceback_queue.push(make_triple(&F5t[0], 0, k));
                traceback_queue.push(make_triple(&FCt[0], k+1, j-1));
            }
            break;
            default:
                Assert(false, "Bad traceback.");
        }
    }

    FinalizeCounts();
    return GetCounts();
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetViterbiFeatures()
// 
// Use feature counts from Viterbi decoding.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::GetViterbiFeatures()
{
    std::queue<triple<int *,int,int> > traceback_queue;
    traceback_queue.push(make_triple(&F5t[0], 0, L));

    std::vector<triple<int, int, int> > multiloops;
    ClearCounts();

    std::map<int, RealT> eos_cb_map;
    for (int k = -1; k <= L; k++){
        eos_cb_map[k] = 0;
    }

    while (!traceback_queue.empty())
    {
        triple<int *,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int *V = t.first;
        const int i = t.second;
        const int j = t.third;
        
        std::pair<int,int> traceback = DecodeTraceback (V == &F5t[0] ? V[j] : V[offset[i]+j]);
        
        //std::cout << (V == FCt ? "FC " : V == FMt ? "FM " : V == FM1t ? "FM1 " : "F5 ");
        //std::cout << i << " " << j << ": " << traceback.first << " " << traceback.second << std::endl;
        
        switch (traceback.first)
        {

    // HWS: we're not using these params
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case TB_FN_HAIRPIN: 
                CountHairpin(i,j,RealT(1));
                break;
            case TB_FN_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                CountSingle(i,j,p,q,RealT(1));

                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FN_BIFURCATION:
            {
                const int k = traceback.second;
                CountJunctionA(i,j,RealT(1));
                CountMultiPaired(RealT(1));
                CountMultiBase(RealT(1));

                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FE_STACKING: 
            {
                CountBasePair(i+1,j,RealT(1));
                CountHelixStacking(i,j+1,RealT(1));

                traceback_queue.push(make_triple(&FEt[0], i+1, j-1));
            }
            break;
            case TB_FE_FN: 
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_FN:
            {
                CountIsolated(RealT(1));
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_HELIX:
            {
                const int m = traceback.second;
                CountHelix(i-1,j+1,m,RealT(1));
                traceback_queue.push(make_triple(&FNt[0], i+m-1, j-m+1));
            }
            break;
            case TB_FC_FE:
            {
                const int m = D_MAX_HELIX_LENGTH;
                CountHelix(i-1,j+1,m,RealT(1));
                traceback_queue.push(make_triple(&FEt[0], i+m-1, j-m+1));
            }
            break;

    //HWS: our param sections start here
#else
            case TB_FC_HAIRPIN: 
                //CountHairpin(i,j,RealT(1));
                //////////////////////////////////////////////////////////////////////
                // InferenceEngine::ScoreHairpin()
                // InferenceEngine::CountHairpin()
                //
                // Returns the score for a hairpin spanning positions i to j.
                //
                // In an RNA structure, this would look like
                //
                //                           ...
                //                       /         \. 
                //                   x[i+2]       x[j-1]
                //                      |            |
                //                   x[i+1]        x[j]
                // position i -------->  o          o  <----- position j
                //                      x[i] -- x[j+1]
                //                        |        |
                //                     x[i-1] -- x[j+2]
                //
                //////////////////////////////////////////////////////////////////////

                //HWS: subtracting off 1 to get zero-indexing for eos_cb

                //std::cerr << "Hairpin " << i-1 << " " << ScoreHairpin(i,j) << std::endl;
                eos_cb_map[i-1] += ScoreHairpin(i,j);

                break;
            case TB_FC_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);

                if (p == i && q == j)
                {
                    // CountBasePair(i+1,j,RealT(1));
                    // CountHelixStacking(i,j+1,RealT(1));

                //HWS: subtracting off 1 (from i+1) to get zero-indexing for eos_cb


                //std::cerr << "BasePair " << i << " " << ScoreBasePair(i+1,j) << std::endl;
                eos_cb_map[i] += ScoreBasePair(i+1,j);


                // ScoreHelixStacking(i,j): score for a helix stacking pair of the form:
                //
                //       |         |
                //    s[i+1] == s[j-1]
                //       |         |
                //     s[i] ==== s[j]
                //       |         |
                //HWS: subtracting off 1 (from i) to get zero-indexing for eos_cb

                //std::cerr << "HelixStacking " << i-1 << " " << ScoreHelixStacking(i,j+1) << std::endl;
                eos_cb_map[i-1] += ScoreHelixStacking(i,j+1);

                }
                else
                {
                    //CountSingle(i,j,p,q,RealT(1));
                    //here we need to break this down because they sum a bunch of things to make ScoreSingle
                    // here we need to subtract off the second base pair from ScoreSingle
                // std::cerr << "Single " << i-1 << " " << ScoreSingle(i,j,p,q)-ScoreBasePair(p+1,q) << std::endl;
                // std::cerr << "ScoreBasePair " << p << " " << ScoreBasePair(p+1,q) << std::endl;

                eos_cb_map[i-1] += ScoreSingle(i,j,p,q)-ScoreBasePair(p+1,q);
                eos_cb_map[p] += ScoreBasePair(p+1, q);

                }
                
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FC_BIFURCATION:
            {
                const int k = traceback.second;
                // CountJunctionA(i,j,RealT(1));
                // CountMultiPaired(RealT(1));
                // CountMultiBase(RealT(1));

                // This is where we first recognize we're in a new multiloop. Store the things at i
                multiloops.push_back(make_triple(i,j,0));
                // std::cerr << "New multiloop! " << i << ":" << j << std::endl;

                // std::cerr << "FCBifurc:JunctionA " << i << " " << ScoreJunctionA(i,j) << std::endl;
                // std::cerr << "FCBifurc:MultiPaired " << i << " " << ScoreMultiPaired() << std::endl;
                // std::cerr << "FCBifurc:MultiBase " << i << " " << ScoreMultiBase() << std::endl;

                eos_cb_map[i] += ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase();

                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
#endif
            case TB_FM1_PAIRED:
            {
                // CountJunctionA(j,i,RealT(1));
                // CountMultiPaired(RealT(1));
                // CountBasePair(i+1,j,RealT(1));

                // HWS: We're in a multiloop -- this gets the index of the 5' base 
                // of the current multiloop and puts in curr_multiloop ind.

                int min_dist = 100000;
                int curr_multiloop_ind = -1;

                for(int ind=0; ind < multiloops.size(); ind++){
                    if (multiloops[ind].first < i < multiloops[ind].second){
                        if (i - multiloops[ind].first < min_dist){
                            min_dist = i - multiloops[ind].first;
                            curr_multiloop_ind = multiloops[ind].first;
                        }
                    }
                }

                // These all go in the multiloop ind
                //std::cerr << "Paired:JunctionA " << curr_multiloop_ind << " " << ScoreJunctionA(j,i) << std::endl;
                //std::cerr << "Paired:MultiPaired " << curr_multiloop_ind << " " << ScoreMultiPaired() << std::endl;
                
                eos_cb_map[curr_multiloop_ind] += ScoreJunctionA(j,i) + ScoreMultiPaired();

                //This just goes with the next base pair
                //std::cerr << "Paired:BasePair " << i << " " << ScoreBasePair(i+1, j) << std::endl;
                eos_cb_map[i] += ScoreBasePair(i+1, j);

                traceback_queue.push(make_triple(&FCt[0], i+1, j-1));
            }
            break;
            case TB_FM1_UNPAIRED:
            {
                // CountMultiUnpaired(i+1,RealT(1));

                // HWS: We're in a multiloop -- this gets the index of the 5' base 
                // of the current multiloop and puts in curr_multiloop ind.

                int min_dist = 100000;
                int curr_multiloop_ind = -1;

                for(int ind=0; ind < multiloops.size(); ind++){
                    if (multiloops[ind].first < i+1 < multiloops[ind].second){
                        if (i+1 - multiloops[ind].first < min_dist){
                            min_dist = i+1 - multiloops[ind].first;
                            curr_multiloop_ind = multiloops[ind].first;
                        }
                    }
                }
                // HWS: a score for unpaired bases in multiloops, put in curr_multiloop ind
                //std::cerr << "MultiUnpaired " << curr_multiloop_ind << " " << ScoreMultiUnpaired(i+1) << std::endl;
                eos_cb_map[curr_multiloop_ind] += ScoreMultiUnpaired(i+1);

                traceback_queue.push(make_triple(&FM1t[0], i+1, j));
            }
            break;
            case TB_FM_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FM_UNPAIRED:
            {
                // CountMultiUnpaired(j,RealT(1));

                // HWS: We're in a multiloop -- this gets the index of the 5' base 
                // of the current multiloop and puts in curr_multiloop ind.

                int min_dist = 100000;
                int curr_multiloop_ind = -1;

                for(int ind=0; ind < multiloops.size(); ind++){
                    if (multiloops[ind].first < j < multiloops[ind].second){
                        if (j - multiloops[ind].first < min_dist){
                            min_dist = j - multiloops[ind].first;
                            curr_multiloop_ind = multiloops[ind].first;
                        }
                    }
                }
                // HWS: a score for unpaired bases in multiloops, put in curr_multiloop ind

                //std::cerr << "MultiUnpaired " << curr_multiloop_ind << " " << ScoreMultiUnpaired(j) << std::endl;
                eos_cb_map[curr_multiloop_ind] += ScoreMultiUnpaired(j);

                traceback_queue.push(make_triple(&FMt[0], i, j-1));
            }
            break;
            case TB_FM_FM1: 
                traceback_queue.push(make_triple(&FM1t[0], i, j));
                break;
            case TB_F5_ZERO:
                break;
            case TB_F5_UNPAIRED:
                CountExternalUnpaired(j,RealT(1));

                // HWS: a score for unpaired things in the external loop, dump into the -1 eos_cb

                //std::cerr << "ExternalUnpaired -1 " << ScoreExternalUnpaired(j) << std::endl;
                eos_cb_map[-1] += ScoreExternalUnpaired(j);

                traceback_queue.push(make_triple(&F5t[0], 0, j-1));
                break;
            case TB_F5_BIFURCATION:
            {
                const int k = traceback.second;
                // CountExternalPaired(RealT(1));
                // CountBasePair(k+1,j,RealT(1));
                // CountJunctionA(j,k,RealT(1));
                traceback_queue.push(make_triple(&F5t[0], 0, k));
                traceback_queue.push(make_triple(&FCt[0], k+1, j-1));

                //HWS: this is base pairs facing external loop. Hard-coding to print -1 as thinking these 
                //should get summed in to the -1 "external loop " eos_cb category
                
                //std::cerr << "F5Bifurc:ExternalPaired -1 " << ScoreExternalPaired() << std::endl;
                //std::cerr << "F5Bifurc:JunctionA -1 " << ScoreJunctionA(j,k) << std::endl;
                eos_cb_map[-1] += ScoreExternalPaired()+ScoreJunctionA(j,k);

                // HWS: another base pair, print (k+1) - 1 to zero-index
                //std::cerr << "F5Bifurc:BasePair " << k << " " << ScoreBasePair(k+1, j) << std::endl;
                eos_cb_map[k] += ScoreBasePair(k+1, j);

            }
            break;
            default:
                Assert(false, "Bad traceback.");
        }
    }
    RealT local_sum = 0;
    std::cerr << "EOS_CB: ";
    for (std::pair<int, RealT> element : eos_cb_map) {
        if (element.first == -1 || element.second != 0){
        std::cerr << element.first << ", " << element.second << ", ";
        local_sum += element.second;
    }
    }
    std::cerr << std::endl;
    std::cerr << "local sum " << local_sum << std::endl;

}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetViterbiFeaturesESS()
// 
// Use feature counts from Viterbi decoding, using evidence-based energies.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::GetViterbiFeaturesESS()
{
    std::queue<triple<int *,int,int> > traceback_queue;
    traceback_queue.push(make_triple(&F5t[0], 0, L));

    std::vector<triple<int, int, int> > multiloops;
    ClearCounts();

    std::map<int, RealT> eos_cb_map;
    for (int k = -1; k <= L; k++){
        eos_cb_map[k] = 0;
    }

    while (!traceback_queue.empty())
    {
        triple<int *,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int *V = t.first;
        const int i = t.second;
        const int j = t.third;
        
        std::pair<int,int> traceback = DecodeTraceback (V == &F5t[0] ? V[j] : V[offset[i]+j]);
        
        //std::cout << (V == FCt ? "FC " : V == FMt ? "FM " : V == FM1t ? "FM1 " : "F5 ");
        //std::cout << i << " " << j << ": " << traceback.first << " " << traceback.second << std::endl;
        
        switch (traceback.first)
        {

    // HWS: we're not using these params
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case TB_FN_HAIRPIN: 
                CountHairpin(i,j,RealT(1));
                break;
            case TB_FN_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                CountSingle(i,j,p,q,RealT(1));

                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FN_BIFURCATION:
            {
                const int k = traceback.second;
                CountJunctionA(i,j,RealT(1));
                CountMultiPaired(RealT(1));
                CountMultiBase(RealT(1));

                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FE_STACKING: 
            {
                CountBasePair(i+1,j,RealT(1));
                CountHelixStacking(i,j+1,RealT(1));

                traceback_queue.push(make_triple(&FEt[0], i+1, j-1));
            }
            break;
            case TB_FE_FN: 
            {
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_FN:
            {
                CountIsolated(RealT(1));
                traceback_queue.push(make_triple(&FNt[0], i, j));
            }
            break;
            case TB_FC_HELIX:
            {
                const int m = traceback.second;
                CountHelix(i-1,j+1,m,RealT(1));
                traceback_queue.push(make_triple(&FNt[0], i+m-1, j-m+1));
            }
            break;
            case TB_FC_FE:
            {
                const int m = D_MAX_HELIX_LENGTH;
                CountHelix(i-1,j+1,m,RealT(1));
                traceback_queue.push(make_triple(&FEt[0], i+m-1, j-m+1));
            }
            break;

    //HWS: our param sections start here
#else
            case TB_FC_HAIRPIN: 
                //CountHairpin(i,j,RealT(1));
                //////////////////////////////////////////////////////////////////////
                // InferenceEngine::ScoreHairpin()
                // InferenceEngine::CountHairpin()
                //
                // Returns the score for a hairpin spanning positions i to j.
                //
                // In an RNA structure, this would look like
                //
                //                           ...
                //                       /         \. 
                //                   x[i+2]       x[j-1]
                //                      |            |
                //                   x[i+1]        x[j]
                // position i -------->  o          o  <----- position j
                //                      x[i] -- x[j+1]
                //                        |        |
                //                     x[i-1] -- x[j+2]
                //
                //////////////////////////////////////////////////////////////////////

                //HWS: subtracting off 1 to get zero-indexing for eos_cb

                //std::cerr << "Hairpin " << i-1 << " " << ScoreHairpin(i,j) << std::endl;
                eos_cb_map[i-1] += ScoreHairpinEvidence(i,j);

                break;
            case TB_FC_SINGLE: 
            {
                const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);

                if (p == i && q == j)
                {
                    // CountBasePair(i+1,j,RealT(1));
                    // CountHelixStacking(i,j+1,RealT(1));

                //HWS: subtracting off 1 (from i+1) to get zero-indexing for eos_cb


                //std::cerr << "BasePair " << i << " " << ScoreBasePair(i+1,j) << std::endl;
                eos_cb_map[i] += ScoreBasePairEvidence(i+1,j);


                // ScoreHelixStacking(i,j): score for a helix stacking pair of the form:
                //
                //       |         |
                //    s[i+1] == s[j-1]
                //       |         |
                //     s[i] ==== s[j]
                //       |         |
                //HWS: subtracting off 1 (from i) to get zero-indexing for eos_cb

                //std::cerr << "HelixStacking " << i-1 << " " << ScoreHelixStacking(i,j+1) << std::endl;
                eos_cb_map[i-1] += ScoreHelixStacking(i,j+1);

                }
                else
                {
                    //CountSingle(i,j,p,q,RealT(1));
                    //here we need to break this down because they sum a bunch of things to make ScoreSingle
                    // here we need to subtract off the second base pair from ScoreSingle
                // std::cerr << "Single " << i-1 << " " << ScoreSingle(i,j,p,q)-ScoreBasePair(p+1,q) << std::endl;
                // std::cerr << "ScoreBasePair " << p << " " << ScoreBasePair(p+1,q) << std::endl;

                eos_cb_map[i-1] += ScoreSingleEvidence(i,j,p,q)-ScoreBasePairEvidence(p+1,q);
                eos_cb_map[p] += ScoreBasePairEvidence(p+1, q);

                }
                
                traceback_queue.push(make_triple(&FCt[0], p+1, q-1));
            }
            break;
            case TB_FC_BIFURCATION:
            {
                const int k = traceback.second;
                // CountJunctionA(i,j,RealT(1));
                // CountMultiPaired(RealT(1));
                // CountMultiBase(RealT(1));

                // This is where we first recognize we're in a new multiloop. Store the things at i
                multiloops.push_back(make_triple(i,j,0));
                // std::cerr << "New multiloop! " << i << ":" << j << std::endl;

                // std::cerr << "FCBifurc:JunctionA " << i << " " << ScoreJunctionA(i,j) << std::endl;
                // std::cerr << "FCBifurc:MultiPaired " << i << " " << ScoreMultiPaired() << std::endl;
                // std::cerr << "FCBifurc:MultiBase " << i << " " << ScoreMultiBase() << std::endl;

                eos_cb_map[i] += ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase();

                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
#endif
            case TB_FM1_PAIRED:
            {
                // CountJunctionA(j,i,RealT(1));
                // CountMultiPaired(RealT(1));
                // CountBasePair(i+1,j,RealT(1));

                // HWS: We're in a multiloop -- this gets the index of the 5' base 
                // of the current multiloop and puts in curr_multiloop ind.

                int min_dist = 100000;
                int curr_multiloop_ind = -1;

                for(int ind=0; ind < multiloops.size(); ind++){
                    if (multiloops[ind].first < i < multiloops[ind].second){
                        if (i - multiloops[ind].first < min_dist){
                            min_dist = i - multiloops[ind].first;
                            curr_multiloop_ind = multiloops[ind].first;
                        }
                    }
                }

                // These all go in the multiloop ind
                //std::cerr << "Paired:JunctionA " << curr_multiloop_ind << " " << ScoreJunctionA(j,i) << std::endl;
                //std::cerr << "Paired:MultiPaired " << curr_multiloop_ind << " " << ScoreMultiPaired() << std::endl;
                
                eos_cb_map[curr_multiloop_ind] += ScoreJunctionA(j,i) + ScoreMultiPaired();

                //This just goes with the next base pair
                //std::cerr << "Paired:BasePair " << i << " " << ScoreBasePair(i+1, j) << std::endl;
                eos_cb_map[i] += ScoreBasePairEvidence(i+1, j);

                traceback_queue.push(make_triple(&FCt[0], i+1, j-1));
            }
            break;
            case TB_FM1_UNPAIRED:
            {
                // CountMultiUnpaired(i+1,RealT(1));

                // HWS: We're in a multiloop -- this gets the index of the 5' base 
                // of the current multiloop and puts in curr_multiloop ind.

                int min_dist = 100000;
                int curr_multiloop_ind = -1;

                for(int ind=0; ind < multiloops.size(); ind++){
                    if (multiloops[ind].first < i+1 < multiloops[ind].second){
                        if (i+1 - multiloops[ind].first < min_dist){
                            min_dist = i+1 - multiloops[ind].first;
                            curr_multiloop_ind = multiloops[ind].first;
                        }
                    }
                }
                // HWS: a score for unpaired bases in multiloops, put in curr_multiloop ind
                //std::cerr << "MultiUnpaired " << curr_multiloop_ind << " " << ScoreMultiUnpaired(i+1) << std::endl;
                eos_cb_map[curr_multiloop_ind] += ScoreMultiUnpairedEvidence(i+1);

                traceback_queue.push(make_triple(&FM1t[0], i+1, j));
            }
            break;
            case TB_FM_BIFURCATION:
            {
                const int k = traceback.second;
                traceback_queue.push(make_triple(&FM1t[0], i, k));
                traceback_queue.push(make_triple(&FMt[0], k, j));
            }
            break;
            case TB_FM_UNPAIRED:
            {
                // CountMultiUnpaired(j,RealT(1));

                // HWS: We're in a multiloop -- this gets the index of the 5' base 
                // of the current multiloop and puts in curr_multiloop ind.

                int min_dist = 100000;
                int curr_multiloop_ind = -1;

                for(int ind=0; ind < multiloops.size(); ind++){
                    if (multiloops[ind].first < j < multiloops[ind].second){
                        if (j - multiloops[ind].first < min_dist){
                            min_dist = j - multiloops[ind].first;
                            curr_multiloop_ind = multiloops[ind].first;
                        }
                    }
                }
                // HWS: a score for unpaired bases in multiloops, put in curr_multiloop ind

                //std::cerr << "MultiUnpaired " << curr_multiloop_ind << " " << ScoreMultiUnpaired(j) << std::endl;
                eos_cb_map[curr_multiloop_ind] += ScoreMultiUnpairedEvidence(j);

                traceback_queue.push(make_triple(&FMt[0], i, j-1));
            }
            break;
            case TB_FM_FM1: 
                traceback_queue.push(make_triple(&FM1t[0], i, j));
                break;
            case TB_F5_ZERO:
                break;
            case TB_F5_UNPAIRED:
                CountExternalUnpaired(j,RealT(1));

                // HWS: a score for unpaired things in the external loop, dump into the -1 eos_cb

                //std::cerr << "ExternalUnpaired -1 " << ScoreExternalUnpaired(j) << std::endl;
                eos_cb_map[-1] += ScoreExternalUnpairedEvidence(j);


                traceback_queue.push(make_triple(&F5t[0], 0, j-1));
                break;
            case TB_F5_BIFURCATION:
            {
                const int k = traceback.second;
                CountExternalPaired(RealT(1));
                CountBasePair(k+1,j,RealT(1));
                CountJunctionA(j,k,RealT(1));
                traceback_queue.push(make_triple(&F5t[0], 0, k));
                traceback_queue.push(make_triple(&FCt[0], k+1, j-1));
                //HWS: this is base pairs facing external loop. Hard-coding to print -1 as thinking these 
                //should get summed in to the -1 "external loop " eos_cb category
                
                //std::cerr << "F5Bifurc:ExternalPaired -1 " << ScoreExternalPaired() << std::endl;
                //std::cerr << "F5Bifurc:JunctionA -1 " << ScoreJunctionA(j,k) << std::endl;
                eos_cb_map[-1] += ScoreExternalPaired()+ScoreJunctionA(j,k); //TODO: ScoreExternalPairedEvidence
 
                // HWS: another base pair, print (k+1) - 1 to zero-index
                //std::cerr << "F5Bifurc:BasePair " << k << " " << ScoreBasePair(k+1, j) << std::endl;
                eos_cb_map[k] += ScoreBasePairEvidence(k+1, j);

            }
            break;
            default:
                Assert(false, "Bad traceback.");
        }
    }

    RealT local_sum = 0;
    std::cerr << "EOS_CB: ";
    for (std::pair<int, RealT> element : eos_cb_map) {
        if (element.first == -1 || element.second != 0){
        std::cerr << element.first << ", " << element.second << ", ";
        local_sum += element.second;
    }
    }
    std::cerr << std::endl;
    std::cerr << "local sum " << local_sum << std::endl;

}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeInside()
//
// Run inside algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeInside()
{
    InitializeCache();
        
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif
    
    // initialization

    F5i.clear(); F5i.resize(L+1, RealT(NEG_INF));
    FCi.clear(); FCi.resize(SIZE, RealT(NEG_INF));
    FMi.clear(); FMi.resize(SIZE, RealT(NEG_INF));
    FM1i.clear(); FM1i.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEi.clear(); FEi.resize(SIZE, RealT(NEG_INF));
    FNi.clear(); FNi.resize(SIZE, RealT(NEG_INF));
#endif

    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i, FM1i[offset[i]+k] + FMi[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i[offset[i]+i+1]);
                const RealT *p2 = &(FMi[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR

            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    Fast_LogPlusEquals(sum_i, ScoreHairpin(i,j));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        Fast_LogPlusEquals(sum_i, ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                    }
                }
                
#else

                if (i+2 <= j)                
                {
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT score = (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                                           ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            Fast_LogPlusEquals(sum_i, score);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(sum_i, FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                FNi[offset[i]+j] = sum_i;
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    Fast_LogPlusEquals(sum_i, ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEi[offset[i+1]+j-1]);
                }
                
                // compute FN(i,j)

                Fast_LogPlusEquals(sum_i, FNi[offset[i]+j]);
                
                FEi[offset[i]+j] = sum_i;
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreIsolated() + FN(i,j)
                
                Fast_LogPlusEquals(sum_i, ScoreIsolated() + FNi[offset[i]+j]);
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    Fast_LogPlusEquals(sum_i, ScoreHelix(i-1,j+1,k) + FNi[offset[i+k-1]+j-k+1]);
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        Fast_LogPlusEquals(sum_i, ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]);
                }
                
                FCi[offset[i]+j] = sum_i;
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    Fast_LogPlusEquals(sum_i, ScoreHairpin(i,j));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        Fast_LogPlusEquals(sum_i,
                                           FCi[offset[p+1]+q-1] +
                                           (p == i && q == j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingle(i,j,p,q)));
                    }
                }
                
#else

                {
                    RealT score_helix = (i+2 <= j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            RealT score = (p == i && q == j) ?
                                (score_helix + FCptr[q]) :
                                (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) +
                                 ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            
                            Fast_LogPlusEquals(sum_i, score);
                        }
                    }
                }
#endif

                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(sum_i, FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                FCi[offset[i]+j] = sum_i;
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    Fast_LogPlusEquals(sum_i, FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j));
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                    Fast_LogPlusEquals(sum_i, FM1i[offset[i+1]+j] + ScoreMultiUnpaired(i+1));
                
                FM1i[offset[i]+j] = sum_i;
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                
                Fast_LogPlusEquals(sum_i, FM2i);
                
                // compute FM[i,j-1] + b
                //hzhang: debug
                //if (allow_unpaired_position[j])
                //    Fast_LogPlusEquals(sum_i, FMi[offset[i]+j-1] + ScoreMultiUnpaired(j));
                //
                float temp = 0.0;
                for (int zh=j; zh>=i; zh--)
                {
                   if (allow_unpaired_position[zh])
                     {
                        temp += ScoreMultiUnpaired(zh);
                         Fast_LogPlusEquals(sum_i, FM1i[offset[i]+zh-1] + temp);
			}
		else break;
		}
                // compute FM1[i,j]
                
                Fast_LogPlusEquals(sum_i, FM1i[offset[i]+j]);
                
                FMi[offset[i]+j] = sum_i;
            }
        }
    }
    
    F5i[0] = RealT(0);
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT sum_i = RealT(NEG_INF);
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            Fast_LogPlusEquals(sum_i, F5i[j-1] + ScoreExternalUnpaired(j));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
            if (allow_paired[offset[k+1]+j])
                Fast_LogPlusEquals(sum_i, F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
        
        F5i[j] = sum_i;
    }

#if SHOW_TIMINGS
    std::cerr << "Inside score: " << F5i[L] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeOutside()
//
// Run outside algorithm.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputeOutside()
{
    InitializeCache();
    
#if SHOW_TIMINGS    
    double starting_time = GetSystemTime();
#endif
    
    // initialization
    
    F5o.clear(); F5o.resize(L+1, RealT(NEG_INF));
    FCo.clear(); FCo.resize(SIZE, RealT(NEG_INF));
    FMo.clear(); FMo.resize(SIZE, RealT(NEG_INF));
    FM1o.clear(); FM1o.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEo.clear(); FEo.resize(SIZE, RealT(NEG_INF));
    FNo.clear(); FNo.resize(SIZE, RealT(NEG_INF));
#endif
    
    F5o[L] = RealT(0);  
    for (int j = L; j >= 1; j--)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            Fast_LogPlusEquals(F5o[j-1], F5o[j] + ScoreExternalUnpaired(j));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        {
            for (int k = 0; k < j; k++)
            {
                if (allow_paired[offset[k+1]+j])
                {
                    RealT temp = F5o[j] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k);
                    Fast_LogPlusEquals(F5o[k], temp + FCi[offset[k+1]+j-1]);
                    Fast_LogPlusEquals(FCo[offset[k+1]+j-1], temp + F5i[k]);
                }
            }
        }
    }
    
    for (int i = 0; i <= L; i++)
    {
        for (int j = L; j >= i; j--)
        {
            RealT FM2o = RealT(NEG_INF);
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                
                Fast_LogPlusEquals(FM2o, FMo[offset[i]+j]);
                
                // compute FM[i,j-1] + b
			
		//hzhang: debug                
                //if (allow_unpaired_position[j])
                //    Fast_LogPlusEquals(FMo[offset[i]+j-1], FMo[offset[i]+j] + ScoreMultiUnpaired(j));
                float temp = 0.0; 
		for (int zh=j; zh>i+1;zh--)
		{
		    if (allow_unpaired_position[zh])
		    {
                    temp += ScoreMultiUnpaired(zh);
                    Fast_LogPlusEquals(FM1o[offset[i]+zh-1], FMo[offset[i]+j] + temp);
		    }
		    else break;
		}

                // compute FM1[i,j]
                
                Fast_LogPlusEquals(FM1o[offset[i]+j], FMo[offset[i]+j]);
            }
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    Fast_LogPlusEquals(FCo[offset[i+1]+j-1], FM1o[offset[i]+j] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j));
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                    Fast_LogPlusEquals(FM1o[offset[i+1]+j], FM1o[offset[i]+j] + ScoreMultiUnpaired(i+1));
                
            }
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreIsolated() + FN(i,j)
                
                Fast_LogPlusEquals(FNo[offset[i]+j], ScoreIsolated() + FCo[offset[i]+j]);
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    Fast_LogPlusEquals(FNo[offset[i+k-1]+j-k+1], ScoreHelix(i-1,j+1,k) + FCo[offset[i]+j]);
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        Fast_LogPlusEquals(FEo[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1],
                                           ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FCo[offset[i]+j]);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    Fast_LogPlusEquals(FEo[offset[i+1]+j-1], FEo[offset[i]+j] + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1));
                }
                
                // compute FN(i,j)
                
                Fast_LogPlusEquals(FNo[offset[i]+j], FEo[offset[i]+j]);
            }
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreHairpin(i,j) -- do nothing
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                {
                    RealT temp = FNo[offset[i]+j];
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            Fast_LogPlusEquals(FCo[offset[p+1]+q-1], temp + ScoreSingle(i,j,p,q));
                        }
                    }
                }
#else
                
                {
                    RealT score_other = FNo[offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        RealT *FCptr = &(FCo[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            Fast_LogPlusEquals(FCptr[q], score_other + cache_score_single[p-i][j-q].first + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(FM2o, FNo[offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                // compute ScoreHairpin(i,j) -- do nothing
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])

#if !FAST_SINGLE_BRANCH_LOOPS
                {
                    RealT temp = FCo[offset[i]+j];
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            Fast_LogPlusEquals(FCo[offset[p+1]+q-1],
                                               temp + (p == i && q == j ? ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingle(i,j,p,q)));
                        }
                    }
                }
#else
                
                {
                    RealT score_helix = (i+2 <= j ? FCo[offset[i]+j] + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = FCo[offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        RealT *FCptr = &(FCo[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            Fast_LogPlusEquals(FCptr[q], 
                                               (p == i && q == j) ?
                                               score_helix :
                                               score_other + cache_score_single[p-i][j-q].first + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(FM2o, FCo[offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
            }
            
#endif
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
            {
                Fast_LogPlusEquals(FM1o[offset[i]+k], FM2o + FMi[offset[k]+j]);
                Fast_LogPlusEquals(FMo[offset[k]+j], FM2o + FM1i[offset[i]+k]);
            }

#else
            if (i+2 <= j)
            {
                RealT *p1i = &(FM1i[offset[i]+i+1]);
                RealT *p2i = &(FMi[offset[i+1]+j]);
                RealT *p1o = &(FM1o[offset[i]+i+1]);
                RealT *p2o = &(FMo[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(*p1o, FM2o + *p2i);
                    Fast_LogPlusEquals(*p2o, FM2o + *p1i);
                    ++p1i;
                    ++p1o;
                    p2i += L-k;
                    p2o += L-k;
                }
            }
            
#endif
        }
    }
    
#if SHOW_TIMINGS
    std::cerr << "Outside score: " << F5o[0] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeLogPartitionCoefficient()
//
// Return partition coefficient.
//////////////////////////////////////////////////////////////////////

template<class RealT>
inline RealT InferenceEngine<RealT>::ComputeLogPartitionCoefficient() const
{
    // NOTE: This should be equal to F5o[0]. 
    
    return F5i[L];
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeFeatureCountExpectations()
// 
// Combine the results of the inside and outside algorithms
// in order to compute feature count expectations.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::ComputeFeatureCountExpectations()
{
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif

    //std::cerr << "Inside score: " << F5i[L].GetLogRepresentation() << std::endl;
    //std::cerr << "Outside score: " << F5o[0].GetLogRepresentation() << std::endl;
    
    const RealT Z = ComputeLogPartitionCoefficient();

    ClearCounts();
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {

            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i, FM1i[offset[i]+k] + FMi[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i[offset[i]+i+1]);
                const RealT *p2 = &(FMi[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FNo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
		//hzhang
               CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        CountSingle(i,j,p,q,Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]));
                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                            cache_score_single[p-i][j-q].second += value;
                            CountBasePair(p+1,q,value);
                            CountJunctionB(i,j,value);
                            CountJunctionB(q,p,value);
                            CountSingleNucleotides(i,j,p,q,value);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FEo[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEi[offset[i+1]+j-1]);
                    CountBasePair(i+1,j,value);
                    CountHelixStacking(i,j+1,value);
                }
                
                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j)
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    CountHelix(i-1,j+1,k,Fast_Exp(outside + ScoreHelix(i-1,j+1,k) + FNi[offset[i+k-1]+j-k+1]));
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        CountHelix(i-1,j+1,D_MAX_HELIX_LENGTH,
                                   Fast_Exp(outside + ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]));
                }
            }

#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            RealT value = Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FCi[offset[p+1]+q-1]);
                            CountBasePair(i+1,j,value);
                            CountHelixStacking(i,j+1,value);
                        }
                        else
                        {
                            CountSingle(i,j,p,q,Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]));
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            if (p == i && q == j)
                            {
                                RealT value = Fast_Exp(score_helix + FCptr[q]);
                                cache_score_single[0][0].second += value;
                                CountBasePair(i+1,j,value);
                                CountHelixStacking(i,j+1,value);
                            }
                            else
                            {
                                RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + 
                                                       ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                                cache_score_single[p-i][j-q].second += value;
                                CountBasePair(p+1,q,value);
                                CountJunctionB(i,j,value);
                                CountJunctionB(q,p,value);
                                CountSingleNucleotides(i,j,p,q,value);
                            }
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(FM1o[offset[i]+j] + FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j) - Z);
                    CountJunctionA(j,i,value);
                    CountMultiPaired(value);
                    CountBasePair(i+1,j,value);
                }
                
                // compute FM1[i+1,j] + b
               // HWS -- change here? 
                if (allow_unpaired_position[i+1])
                {
                    CountMultiUnpaired(i+1,Fast_Exp(FM1o[offset[i]+j] + FM1i[offset[i+1]+j] + ScoreMultiUnpaired(i+1) - Z));
                }
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
                
                // compute FM[i,j-1] + b
               //hzhang,hws 
               // if (allow_unpaired_position[j])
               //     CountMultiUnpaired(j,Fast_Exp(FMo[offset[i]+j] + FMi[offset[i]+j-1] + ScoreMultiUnpaired(j) - Z));
               float temp = 0.0;
		for (int zh=j; zh>=i; zh--)
                 { if (allow_unpaired_position[zh]) { 
                  temp += ScoreMultiUnpaired(zh);
                  CountMultiUnpaired(j,Fast_Exp(FMo[offset[i]+j] + FMi[offset[i]+zh-1] + temp - Z));;
                    }
                 else break;
                 }
                // compute FM1[i,j] -- do nothing
            }
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            CountExternalUnpaired(j,Fast_Exp(outside + F5i[j-1] + ScoreExternalUnpaired(j)));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
            {
                RealT value = Fast_Exp(outside + F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
                CountExternalPaired(value);
                CountBasePair(k+1,j,value);
                CountJunctionA(j,k,value);
            }      
        }
    }
    
    FinalizeCounts();

#if SHOW_TIMINGS
    std::cerr << "Feature expectations (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif

    return GetCounts();
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeFeatureCountExpectationsESS()
// 
// Combine the results of the inside and outside algorithms
// in order to compute feature count expectations.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::ComputeFeatureCountExpectationsESS()
{
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif

    //std::cerr << "Inside score: " << F5i[L].GetLogRepresentation() << std::endl;
    //std::cerr << "Outside score: " << F5o[0].GetLogRepresentation() << std::endl;
    
    const RealT Z = ComputeLogPartitionCoefficientESS();

    ClearCounts();
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {

            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i_ess = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i_ess, FM1i_ess[offset[i]+k] + FMi_ess[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i_ess[offset[i]+i+1]);
                const RealT *p2 = &(FMi_ess[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i_ess, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FNo_ess[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpinEvidence(i,j,Fast_Exp(outside + ScoreHairpinEvidence(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        CountSingleEvidence(i,j,p,q,Fast_Exp(outside + ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]));
                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                            cache_score_single[p-i][j-q].second += value;
                            CountBasePairEvidence(p+1,q,value);
                            CountJunctionB(i,j,value);
                            CountJunctionB(q,p,value);
                            CountSingleNucleotidesEvidence(i,j,p,q,value);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i_ess + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FEo_ess[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FEi_ess[offset[i+1]+j-1]);
                    CountBasePairEvidence(i+1,j,value);
                    CountHelixStacking(i,j+1,value);
                }
                
                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo_ess[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j)
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi_ess[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    CountHelix(i-1,j+1,k,Fast_Exp(outside + ScoreHelixEvidence(i-1,j+1,k) + FNi_ess[offset[i+k-1]+j-k+1]));
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        CountHelix(i-1,j+1,D_MAX_HELIX_LENGTH,
                                   Fast_Exp(outside + ScoreHelixEvidence(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi_ess[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]));
                }
            }

#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo_ess[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpinEvidence(i,j,Fast_Exp(outside + ScoreHairpinEvidence(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            RealT value = Fast_Exp(outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FCi_ess[offset[p+1]+q-1]);
                            CountBasePairEvidence(i+1,j,value);
                            CountHelixStacking(i,j+1,value);
                        }
                        else
                        {
                            CountSingleEvidence(i,j,p,q,Fast_Exp(outside + ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]));
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            if (p == i && q == j)
                            {
                                RealT value = Fast_Exp(score_helix + FCptr[q]);
                                cache_score_single[0][0].second += value;
                                CountBasePairEvidence(i+1,j,value);
                                CountHelixStacking(i,j+1,value);
                            }
                            else
                            {
                                RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) +
                                                       ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                                cache_score_single[p-i][j-q].second += value;
                                CountBasePairEvidence(p+1,q,value);
                                CountJunctionB(i,j,value);
                                CountJunctionB(q,p,value);
                                CountSingleNucleotidesEvidence(i,j,p,q,value);
                            }
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i_ess + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(FM1o_ess[offset[i]+j] + FCi_ess[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePairEvidence(i+1,j) - Z);
                    CountJunctionA(j,i,value);
                    CountMultiPaired(value);
                    CountBasePairEvidence(i+1,j,value);
                }
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                {
                    CountMultiUnpairedEvidence(i+1,Fast_Exp(FM1o_ess[offset[i]+j] + FM1i_ess[offset[i+1]+j] + ScoreMultiUnpairedEvidence(i+1) - Z));
                }
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    CountMultiUnpairedEvidence(j,Fast_Exp(FMo_ess[offset[i]+j] + FMi_ess[offset[i]+j-1] + ScoreMultiUnpairedEvidence(j) - Z));
                
                // compute FM1[i,j] -- do nothing
            }
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o_ess[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            CountExternalUnpairedEvidence(j,Fast_Exp(outside + F5i_ess[j-1] + ScoreExternalUnpairedEvidence(j)));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
            {
                RealT value = Fast_Exp(outside + F5i_ess[k] + FCi_ess[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePairEvidence(k+1,j) + ScoreJunctionA(j,k));
                CountExternalPaired(value);
                CountBasePairEvidence(k+1,j,value);
                CountJunctionA(j,k,value);
            }      
        }
    }
    
    FinalizeCountsESS();

#if SHOW_TIMINGS
    std::cerr << "Feature expectations (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif

    return GetCounts();
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputePosterior()
// 
// Combine the results of the inside and outside algorithms
// in order to compute posterior probabilities of base pairing.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputePosterior()
{ 
    posterior.clear();
    posterior.resize(SIZE, RealT(0));
    
    //double starting_time = GetSystemTime();

    const RealT Z = ComputeLogPartitionCoefficient();
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i = RealT(NEG_INF);
            
#if SIMPLE_FM2
      
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i, FM1i[offset[i]+k] + FMi[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i[offset[i]+i+1]);
                const RealT *p2 = &(FMi[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
      
#endif

#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FNo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            posterior[offset[p+1]+q] += 
                                Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
#endif

                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c) -- do nothing
                
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FEo[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                    posterior[offset[i]+j] += Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FEi[offset[i+1]+j-1]);
                
                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j) -- do nothing
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    RealT value = Fast_Exp(outside + ScoreHelix(i-1,j+1,k) + FNi[offset[i+k-1]+j-k+1]);
                    for (int p = 1; p < k; p++)
                        posterior[offset[i+p]+j-p+1] += value;
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2]) {
                        RealT value = Fast_Exp(outside + ScoreHelix(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]);
                        
                        for (int k = 1; k < D_MAX_HELIX_LENGTH; k++)
                            posterior[offset[i+k]+j-k+1] += value;
                    }
                }
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FCo[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpin(i,j,Fast_Exp(outside + ScoreHairpin(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) + FCi[offset[p+1]+q-1]);
                        }
                        else
                        {
                            posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            posterior[offset[p+1]+q] +=
                                Fast_Exp(p == i && q == j ?
                                         score_helix + FCptr[q] :
                                         score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePair(p+1,q) + 
                                         ScoreJunctionB(q,p) + ScoreSingleNucleotides(i,j,p,q));
                        }
                    }
                }
                
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c) -- do nothing
                
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // Compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    posterior[offset[i+1]+j] += Fast_Exp(FM1o[offset[i]+j] + FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePair(i+1,j) - Z);
                
                // Compute FM1[i+1,j] + b -- do nothing
                
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            // Compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
            
            // Compute FM[i,j-1] + b -- do nothing
            
            // Compute FM1[i,j] -- do nothing
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired() -- do nothing
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
                posterior[offset[k+1]+j] += Fast_Exp(outside + F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePair(k+1,j) + ScoreJunctionA(j,k));
        }
    }

    for (int i = 1; i <= L; i++)
    {
     	for (int j = i+1; j <= L; j++)
        {
            posterior[offset[i]+j] = Clip(posterior[offset[i]+j], RealT(0), RealT(1));
        }
    }

}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputePosteriorESS()
// 
// Combine the results of the inside and outside algorithms
// in order to compute posterior probabilities of base pairing.
//////////////////////////////////////////////////////////////////////

template<class RealT>
void InferenceEngine<RealT>::ComputePosteriorESS()
{ 
    posterior.clear();
    posterior.resize(SIZE, RealT(0));
    
    //double starting_time = GetSystemTime();

    const RealT Z = ComputeLogPartitionCoefficientESS();

    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i_ess = RealT(NEG_INF);
            
#if SIMPLE_FM2
      
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i_ess, FM1i_ess[offset[i]+k] + FMi_ess[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i_ess[offset[i]+i+1]);
                const RealT *p2 = &(FMi_ess[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i_ess, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
      
#endif

#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FNo_ess[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpinEvidence(i,j,Fast_Exp(outside + ScoreHairpinEvidence(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]);

                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            posterior[offset[p+1]+q] += 
                                Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));

                        }
                    }
                }
#endif

                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c) -- do nothing
                
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FEo_ess[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j]) {
                    posterior[offset[i]+j] += Fast_Exp(outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FEi_ess[offset[i+1]+j-1]);

                }

                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FCo_ess[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j) -- do nothing
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi_ess[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    RealT value = Fast_Exp(outside + ScoreHelixEvidence(i-1,j+1,k) + FNi_ess[offset[i+k-1]+j-k+1]);
                    for (int p = 1; p < k; p++)
                        posterior[offset[i+p]+j-p+1] += value;
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2]) {
                        RealT value = Fast_Exp(outside + ScoreHelixEvidence(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi_ess[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]);
                        
                        for (int k = 1; k < D_MAX_HELIX_LENGTH; k++)
                            posterior[offset[i+k]+j-k+1] += value;
                    }
                }
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FCo_ess[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpinEvidence(i,j,Fast_Exp(outside + ScoreHairpinEvidence(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FCi_ess[offset[p+1]+q-1]);
                        }
                        else
                        {
                            posterior[offset[p+1]+q] += Fast_Exp(outside + ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]);
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            posterior[offset[p+1]+q] +=
                                Fast_Exp(p == i && q == j ?
                                         score_helix + FCptr[q] :
                                         score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) +
                                         ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                        }
                    }
                }
                
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c) -- do nothing
                
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // Compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    posterior[offset[i+1]+j] += Fast_Exp(FM1o_ess[offset[i]+j] + FCi_ess[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePairEvidence(i+1,j) - Z);
                
                // Compute FM1[i+1,j] + b -- do nothing
                
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            // Compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
            
            // Compute FM[i,j-1] + b -- do nothing
            
            // Compute FM1[i,j] -- do nothing
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o_ess[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired() -- do nothing
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j]) {
                posterior[offset[k+1]+j] += Fast_Exp(outside + F5i_ess[k] + FCi_ess[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePairEvidence(k+1,j) + ScoreJunctionA(j,k));
            }
        }
    }


    for (int i = 1; i <= L; i++)
    {
        for (int j = i+1; j <= L; j++)
        {
            posterior[offset[i]+j] = Clip(posterior[offset[i]+j], RealT(0), RealT(1));
        }
    }

}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::PredictPairingsPosterior()
//
// Use posterior decoding to predict pairings.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsPosterior(const RealT gamma) const
{
    Assert(gamma > 0, "Non-negative gamma expected.");
    
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif
    RealT* unpaired_posterior  = new RealT[L+1];
    RealT* score               = new RealT[SIZE];
    int* traceback             = new int[SIZE];
    
    // compute the scores for unpaired nucleotides
    
    for (int i = 1; i <= L; i++)
    {
        unpaired_posterior[i] = RealT(1);
        for (int j = 1; j < i; j++) unpaired_posterior[i] -= posterior[offset[j]+i];
        for (int j = i+1; j <= L; j++) unpaired_posterior[i] -= posterior[offset[i]+j];
    }
    
    for (int i = 1; i <= L; i++) unpaired_posterior[i] /= 2 * gamma;
    
    // initialize matrices
    
    std::fill(score, score+SIZE, RealT(-1.0));
    std::fill(traceback, traceback+SIZE, -1);
    
    // dynamic programming
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            RealT &this_score = score[offset[i]+j];
            int &this_traceback = traceback[offset[i]+j];
            
            if (i == j)
            {
                UPDATE_MAX(this_score, this_traceback, RealT(0), 0);
            }
            else
            {
                if (allow_unpaired_position[i+1])
                    UPDATE_MAX(this_score, this_traceback, unpaired_posterior[i+1] + score[offset[i+1]+j], 1);
                if (allow_unpaired_position[j])
                    UPDATE_MAX(this_score, this_traceback, unpaired_posterior[j] + score[offset[i]+j-1], 2);
                if (i+2 <= j)
                { 
                    if (allow_paired[offset[i+1]+j])
                        UPDATE_MAX(this_score, this_traceback, posterior[offset[i+1]+j] + score[offset[i+1]+j-1], 3);
                    
#if SIMPLE_FM2
                    
                    for (int k = i+1; k < j; k++)
                        UPDATE_MAX(this_score, this_traceback, score[offset[i]+k] + score[offset[k]+j], k+4);	
                    
#else
                    
                    RealT *p1 = &(score[offset[i]+i+1]);
                    RealT *p2 = &(score[offset[i+1]+j]);
                    for (register int k = i+1; k < j; k++)
                    {
                        UPDATE_MAX(this_score, this_traceback, (*p1) + (*p2), k+4);
                        ++p1;
                        p2 += L-k;
                    }
                    
#endif
                }
            }
        }
    }

#if SHOW_TIMINGS
    std::cerr << "Time: " << GetSystemTime() - starting_time << std::endl;
#endif
    
    // perform traceback
    
    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;
    
    std::queue<std::pair<int,int> > traceback_queue;
    traceback_queue.push(std::make_pair(0, L));
    
    while (!traceback_queue.empty())
    {
        std::pair<int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int i = t.first;
        const int j = t.second;
        
        switch (traceback[offset[i]+j])
        {
            case -1:
                Assert(false, "Should not get here.");
                break;
            case 0: 
                break;
            case 1: 
                traceback_queue.push(std::make_pair(i+1,j));
                break;
            case 2: 
                traceback_queue.push(std::make_pair(i,j-1));
                break;
            case 3:
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(std::make_pair(i+1,j-1));
                break;
            default:
            {
                const int k = traceback[offset[i]+j] - 4;
                traceback_queue.push(std::make_pair(i,k));
                traceback_queue.push(std::make_pair(k,j));
            }
            break;       
        }
    }
    
    return solution;
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::PredictPairingsPosteriorCentroid()
//
// Use Centroid estimator (CentroidFold, Hamada et. al.,  2009)
// to predict pairings.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsPosteriorCentroid(const RealT gamma) const
{
    Assert(gamma > 0, "Non-negative gamma expected.");

#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif
    RealT* score               = new RealT[SIZE];
    int* traceback             = new int[SIZE];

    // initialize matrices

    std::fill(score, score+SIZE, RealT(-1.0));
    std::fill(traceback, traceback+SIZE, -1);

    // dynamic programming

    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            RealT &this_score = score[offset[i]+j];
            int &this_traceback = traceback[offset[i]+j];

            if (i == j)
            {
                UPDATE_MAX(this_score, this_traceback, RealT(0), 0);
            }
            else
            {
                if (allow_unpaired_position[i+1])
                    UPDATE_MAX(this_score, this_traceback, score[offset[i+1]+j], 1);
                if (allow_unpaired_position[j])
                    UPDATE_MAX(this_score, this_traceback, score[offset[i]+j-1], 2);
                if (i+2 <= j)
                { 
                    if (allow_paired[offset[i+1]+j])
                        UPDATE_MAX(this_score, this_traceback,
                                   (gamma + 1)*posterior[offset[i+1]+j] - 1 + score[offset[i+1]+j-1], 3);

#if SIMPLE_FM2
                    for (int k = i+1; k < j; k++)
                        UPDATE_MAX(this_score, this_traceback, score[offset[i]+k] + score[offset[k]+j], k+4);
#else
                    RealT *p1 = &(score[offset[i]+i+1]);
                    RealT *p2 = &(score[offset[i+1]+j]);
                    for (register int k = i+1; k < j; k++)
                    {
                        UPDATE_MAX(this_score, this_traceback, (*p1) + (*p2), k+4);
                        ++p1;
                        p2 += L-k;
                    }
#endif
                }
            }
        }
    }

#if SHOW_TIMINGS
    std::cerr << "Time: " << GetSystemTime() - starting_time << std::endl;
#endif

    // perform traceback

    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;

    std::queue<std::pair<int,int> > traceback_queue;
    traceback_queue.push(std::make_pair(0, L));

    while (!traceback_queue.empty())
    {
        std::pair<int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int i = t.first;
        const int j = t.second;

        switch (traceback[offset[i]+j])
        {
            case -1:
                Assert(false, "Should not get here.");
                break;
            case 0:
                break;
            case 1:
                traceback_queue.push(std::make_pair(i+1,j));
                break;
            case 2:
                traceback_queue.push(std::make_pair(i,j-1));
                break;
            case 3:
                solution[i+1] = j;
                solution[j] = i+1;
                traceback_queue.push(std::make_pair(i+1,j-1));
                break;
            default:
            {
                const int k = traceback[offset[i]+j] - 4;
                traceback_queue.push(std::make_pair(i,k));
                traceback_queue.push(std::make_pair(k,j));
            }
            break;
        }
    }

    return solution;
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::GetPosterior()
//
// Return posterior probability matrix, thresholded.
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT *InferenceEngine<RealT>::GetPosterior(const RealT posterior_cutoff) const
{
    RealT *ret = new RealT[SIZE];
    for (int i = 0; i < SIZE; i++)
        ret[i] = (posterior[i] >= posterior_cutoff ? posterior[i] : RealT(0));
    return ret;
}

template<class RealT>
void InferenceEngine<RealT>::ComputeInsideESS() 
{
    InitializeCacheESS();
        
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif
    
    // initialization

    F5i_ess.clear(); F5i_ess.resize(L+1, RealT(NEG_INF));
    FCi_ess.clear(); FCi_ess.resize(SIZE, RealT(NEG_INF));
    FMi_ess.clear(); FMi_ess.resize(SIZE, RealT(NEG_INF));
    FM1i_ess.clear(); FM1i_ess.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEi_ess.clear(); FEi_ess.resize(SIZE, RealT(NEG_INF));
    FNi_ess.clear(); FNi_ess.resize(SIZE, RealT(NEG_INF));
#endif

    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i_ess = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i_ess, FM1i_ess[offset[i]+k] + FMi_ess[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i_ess[offset[i]+i+1]);
                const RealT *p2 = &(FMi_ess[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i_ess, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR

            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    Fast_LogPlusEquals(sum_i, ScoreHairpinEvidence(i,j));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        Fast_LogPlusEquals(sum_i, ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]);
                    }
                }
                
#else

                if (i+2 <= j)                
                {
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT score = (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) +
                                           ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                            
                            Fast_LogPlusEquals(sum_i, score);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(sum_i, FM2i_ess + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                FNi_ess[offset[i]+j] = sum_i;
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    Fast_LogPlusEquals(sum_i, ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FEi_ess[offset[i+1]+j-1]);
                }
                
                // compute FN(i,j)

                Fast_LogPlusEquals(sum_i, FNi_ess[offset[i]+j]);
                
                FEi_ess[offset[i]+j] = sum_i;
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreIsolated() + FN(i,j)
                
                Fast_LogPlusEquals(sum_i, ScoreIsolated() + FNi_ess[offset[i]+j]);
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    Fast_LogPlusEquals(sum_i, ScoreHelixEvidence(i-1,j+1,k) + FNi_ess[offset[i+k-1]+j-k+1]);
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        Fast_LogPlusEquals(sum_i, ScoreHelixEvidence(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi_ess[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]);
                }
                
                FCi_ess[offset[i]+j] = sum_i;
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT sum_i = RealT(NEG_INF);
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    Fast_LogPlusEquals(sum_i, ScoreHairpinEvidence(i,j));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        Fast_LogPlusEquals(sum_i,
                                           FCi_ess[offset[p+1]+q-1] +
                                           (p == i && q == j ? ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingleEvidence(i,j,p,q)));
                    }
                }
                
#else

                {
                    RealT score_helix = (i+2 <= j ? ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            RealT score = (p == i && q == j) ?
                                (score_helix + FCptr[q]) :
                                (score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) +
                                 ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                            
                            Fast_LogPlusEquals(sum_i, score);
                        }
                    }
                }
#endif

                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(sum_i, FM2i_ess + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
                FCi_ess[offset[i]+j] = sum_i;
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    Fast_LogPlusEquals(sum_i, FCi_ess[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePairEvidence(i+1,j));
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                    Fast_LogPlusEquals(sum_i, FM1i_ess[offset[i+1]+j] + ScoreMultiUnpairedEvidence(i+1));
                
                FM1i_ess[offset[i]+j] = sum_i;
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                RealT sum_i = RealT(NEG_INF);
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                
                Fast_LogPlusEquals(sum_i, FM2i_ess);
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    Fast_LogPlusEquals(sum_i, FMi_ess[offset[i]+j-1] + ScoreMultiUnpairedEvidence(j));
                
                // compute FM1[i,j]
                
                Fast_LogPlusEquals(sum_i, FM1i_ess[offset[i]+j]);
                
                FMi_ess[offset[i]+j] = sum_i;
            }
        }
    }
    
    F5i_ess[0] = RealT(0);
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT sum_i = RealT(NEG_INF);
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            Fast_LogPlusEquals(sum_i, F5i_ess[j-1] + ScoreExternalUnpairedEvidence(j));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
            if (allow_paired[offset[k+1]+j])
                Fast_LogPlusEquals(sum_i, F5i_ess[k] + FCi_ess[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePairEvidence(k+1,j) + ScoreJunctionA(j,k));
        
        F5i_ess[j] = sum_i;
    }

#if SHOW_TIMINGS
    std::cerr << "Inside score: " << F5i_ess[L] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
}


//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeOutsideESS()
//
// Run outside algorithm for ESS (includes experimental data).
//////////////////////////////////////////////////////////////////////
template<class RealT>
void InferenceEngine<RealT>::ComputeOutsideESS() 
{
    InitializeCacheESS();
    
#if SHOW_TIMINGS    
    double starting_time = GetSystemTime();
#endif
    
    // initialization
    
    F5o_ess.clear(); F5o_ess.resize(L+1, RealT(NEG_INF));
    FCo_ess.clear(); FCo_ess.resize(SIZE, RealT(NEG_INF));
    FMo_ess.clear(); FMo_ess.resize(SIZE, RealT(NEG_INF));
    FM1o_ess.clear(); FM1o_ess.resize(SIZE, RealT(NEG_INF));
    
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
    FEo_ess.clear(); FEo_ess.resize(SIZE, RealT(NEG_INF));
    FNo_ess.clear(); FNo_ess.resize(SIZE, RealT(NEG_INF));
#endif
    
    F5o_ess[L] = RealT(0);  
    for (int j = L; j >= 1; j--)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            Fast_LogPlusEquals(F5o_ess[j-1], F5o_ess[j] + ScoreExternalUnpairedEvidence(j));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        {
            for (int k = 0; k < j; k++)
            {
                if (allow_paired[offset[k+1]+j])
                {
                    RealT temp = F5o_ess[j] + ScoreExternalPaired() + ScoreBasePairEvidence(k+1,j) + ScoreJunctionA(j,k);
                    Fast_LogPlusEquals(F5o_ess[k], temp + FCi_ess[offset[k+1]+j-1]);
                    Fast_LogPlusEquals(FCo_ess[offset[k+1]+j-1], temp + F5i_ess[k]);
                }
            }
        }
    }
    
    for (int i = 0; i <= L; i++)
    {
        for (int j = L; j >= i; j--)
        {
            RealT FM2o_ess = RealT(NEG_INF);
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j])
                
                Fast_LogPlusEquals(FM2o_ess, FMo_ess[offset[i]+j]);
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    Fast_LogPlusEquals(FMo_ess[offset[i]+j-1], FMo_ess[offset[i]+j] + ScoreMultiUnpairedEvidence(j));
                
                // compute FM1[i,j]
                
                Fast_LogPlusEquals(FM1o_ess[offset[i]+j], FMo_ess[offset[i]+j]);
            }
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                    Fast_LogPlusEquals(FCo_ess[offset[i+1]+j-1], FM1o_ess[offset[i]+j] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePairEvidence(i+1,j));
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                    Fast_LogPlusEquals(FM1o_ess[offset[i+1]+j], FM1o_ess[offset[i]+j] + ScoreMultiUnpairedEvidence(i+1));
                
            }
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreIsolated() + FN(i,j)
                
                Fast_LogPlusEquals(FNo_ess[offset[i]+j], ScoreIsolated() + FCo_ess[offset[i]+j]);
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    Fast_LogPlusEquals(FNo_ess[offset[i+k-1]+j-k+1], ScoreHelixEvidence(i-1,j+1,k) + FCo_ess[offset[i]+j]);
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        Fast_LogPlusEquals(FEo_ess[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1],
                                           ScoreHelixEvidence(i-1,j+1,D_MAX_HELIX_LENGTH) + FCo_ess[offset[i]+j]);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    Fast_LogPlusEquals(FEo_ess[offset[i+1]+j-1], FEo_ess[offset[i]+j] + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1));
                }
                
                // compute FN(i,j)
                
                Fast_LogPlusEquals(FNo_ess[offset[i]+j], FEo_ess[offset[i]+j]);
            }
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                // compute ScoreHairpin(i,j) -- do nothing
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                {
                    RealT temp = FNo_ess[offset[i]+j];
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            Fast_LogPlusEquals(FCo_ess[offset[p+1]+q-1], temp + ScoreSingleEvidence(i,j,p,q));
                        }
                    }
                }
#else
                
                {
                    RealT score_other = FNo_ess[offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        RealT *FCptr = &(FCo_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            Fast_LogPlusEquals(FCptr[q], score_other + cache_score_single[p-i][j-q].first + ScoreBasePairEvidence(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(FM2o_ess, FNo_ess[offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
            }
            
#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                // compute ScoreHairpin(i,j) -- do nothing
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])

#if !FAST_SINGLE_BRANCH_LOOPS
                {
                    RealT temp = FCo_ess[offset[i]+j];
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            Fast_LogPlusEquals(FCo_ess[offset[p+1]+q-1],
                                               temp + (p == i && q == j ? ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : ScoreSingleEvidence(i,j,p,q)));
                        }
                    }
                }
#else
                
                {
                    RealT score_helix = (i+2 <= j ? FCo_ess[offset[i]+j] + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = FCo_ess[offset[i]+j] + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        RealT *FCptr = &(FCo_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            Fast_LogPlusEquals(FCptr[q], 
                                               (p == i && q == j) ?
                                               score_helix :
                                               score_other + cache_score_single[p-i][j-q].first + ScoreBasePairEvidence(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                Fast_LogPlusEquals(FM2o_ess, FCo_ess[offset[i]+j] + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                
            }
            
#endif
            
            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
            {
                Fast_LogPlusEquals(FM1o_ess[offset[i]+k], FM2o_ess + FMi_ess[offset[k]+j]);
                Fast_LogPlusEquals(FMo_ess[offset[k]+j], FM2o_ess + FM1i_ess[offset[i]+k]);
            }

#else
            if (i+2 <= j)
            {
                RealT *p1i = &(FM1i_ess[offset[i]+i+1]);
                RealT *p2i = &(FMi_ess[offset[i+1]+j]);
                RealT *p1o = &(FM1o_ess[offset[i]+i+1]);
                RealT *p2o = &(FMo_ess[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(*p1o, FM2o_ess + *p2i);
                    Fast_LogPlusEquals(*p2o, FM2o_ess + *p1i);
                    ++p1i;
                    ++p1o;
                    p2i += L-k;
                    p2o += L-k;
                }
            }
            
#endif
        }
    }
    
#if SHOW_TIMINGS
    std::cerr << "Outside score: " << F5o_ess[0] << " (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeLogPartitionCoefficientESS()
//
// Return partition coefficient (probability of evidence).
//////////////////////////////////////////////////////////////////////

template<class RealT>
RealT InferenceEngine<RealT>::ComputeLogPartitionCoefficientESS() const
{
    // NOTE: This should be equal to F5o_ess[0]. 
    
    return F5i_ess[L];
}

//////////////////////////////////////////////////////////////////////
// InferenceEngine::ComputeESS()
// 
// Combine the results of the inside and outside algorithms
// in order to compute feature count expectations (ESS) for EM.
//////////////////////////////////////////////////////////////////////

template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::ComputeESS()
{
#if SHOW_TIMINGS
    double starting_time = GetSystemTime();
#endif

    //std::cerr << "Inside score: " << F5i_ess[L].GetLogRepresentation() << std::endl;
    //std::cerr << "Outside score: " << F5o_ess[0].GetLogRepresentation() << std::endl;
    
    const RealT Z = ComputeLogPartitionCoefficientESS();

    ClearCounts();
    
    for (int i = L; i >= 0; i--)
    {
        for (int j = i; j <= L; j++)
        {

            // FM2[i,j] = SUM (i<k<j : FM1[i,k] + FM[k,j])
            
            RealT FM2i_ess = RealT(NEG_INF);
            
#if SIMPLE_FM2
            
            for (int k = i+1; k < j; k++)
                Fast_LogPlusEquals(FM2i_ess, FM1i_ess[offset[i]+k] + FMi_ess[offset[k]+j]);
            
#else
            
            if (i+2 <= j)
            {
                const RealT *p1 = &(FM1i_ess[offset[i]+i+1]);
                const RealT *p2 = &(FMi_ess[offset[i+1]+j]);
                for (register int k = i+1; k < j; k++)
                {
                    Fast_LogPlusEquals(FM2i_ess, (*p1) + (*p2));
                    ++p1;
                    p2 += L-k;
                }
            }
            
#endif
            
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            
            // FN[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           and the next interaction is not a stacking pair
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FNo_ess[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpinEvidence(i,j,Fast_Exp(outside + ScoreHairpinEvidence(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j, p-i+j-q>0 : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;
                        if (i == p && j == q) continue;
                        
                        CountSingleEvidence(i,j,p,q,Fast_Exp(outside + ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]));
                    }
                }
                
#else
                
                {
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            if (i == p && j == q) continue;
                            
                            RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) + ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                            cache_score_single[p-i][j-q].second += value;
                            CountBasePairEvidence(p+1,q,value);
                            CountJunctionB(i,j,value);
                            CountJunctionB(q,p,value);
                            CountSingleNucleotidesEvidence(i,j,p,q,value);
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i_ess + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
            // FE[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) ... (i-D+1,j+D) are 
            //           already base-paired
            //
            //         = SUM [ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]   if i+2<=j,
            //                FN(i,j)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                
                RealT outside = FEo_ess[offset[i]+j] - Z;
                
                // compute ScoreBP(i+1,j) + ScoreHelixStacking(i,j+1) + FE[i+1,j-1]
                
                if (i+2 <= j && allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FEi_ess[offset[i+1]+j-1]);
                    CountBasePairEvidence(i+1,j,value);
                    CountHelixStacking(i,j+1,value);
                }
                
                // compute FN(i,j) -- do nothing
                
            }
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //           but (i-1,j+2) are not
            //
            //         = SUM [ScoreIsolated() + FN(i,j),
            //                SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k)),
            //                FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
            //
            //           (assuming 0 < i <= j < L)
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo_ess[offset[i]+j] - Z;
                
                // compute ScoreIsolated() + FN(i,j)
                
                CountIsolated(Fast_Exp(outside + ScoreIsolated() + FNi_ess[offset[i]+j]));
                
                // compute SUM (2<=k<D : FN(i+k-1,j-k+1) + ScoreHelix(i-1,j+1,k))
                
                bool allowed = true;
                for (int k = 2; k < D_MAX_HELIX_LENGTH; k++)
                {
                    if (i + 2*k - 2 > j) break;
                    if (!allow_paired[offset[i+k-1]+j-k+2]) { allowed = false; break; }
                    CountHelix(i-1,j+1,k,Fast_Exp(outside + ScoreHelixEvidence(i-1,j+1,k) + FNi_ess[offset[i+k-1]+j-k+1]));
                }
                
                // compute FE(i+D-1,j-D+1) + ScoreHelix(i-1,j+1,D)]
                
                if (i + 2*D_MAX_HELIX_LENGTH-2 <= j)
                {
                    if (allowed && allow_paired[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+2])
                        CountHelix(i-1,j+1,D_MAX_HELIX_LENGTH,
                                   Fast_Exp(outside + ScoreHelixEvidence(i-1,j+1,D_MAX_HELIX_LENGTH) + FEi_ess[offset[i+D_MAX_HELIX_LENGTH-1]+j-D_MAX_HELIX_LENGTH+1]));
                }
            }

#else
            
            // FC[i,j] = optimal energy for substructure between positions
            //           i and j such that letters (i,j+1) are base-paired
            //
            //         = SUM [ScoreHairpin(i,j),
            //                SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1]),
            //                ScoreJunctionA(i,j) + a + c + SUM (i<k<j : FM1[i,k] + FM[k,j])]
            //
            //           (assuming 0 < i <= j < L)
            //
            // Multi-branch loops are scored as [a + b * (# unpaired) + c * (# branches)]
            
            if (0 < i && j < L && allow_paired[offset[i]+j+1])
            {
                RealT outside = FCo_ess[offset[i]+j] - Z;
                
                // compute ScoreHairpin(i,j)
                
                if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                    CountHairpinEvidence(i,j,Fast_Exp(outside + ScoreHairpinEvidence(i,j)));
                
                // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                
#if !FAST_SINGLE_BRANCH_LOOPS
                for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                {
                    if (p > i && !allow_unpaired_position[p]) break;
                    int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                    for (int q = j; q >= q_min; q--)
                    {
                        if (q < j && !allow_unpaired_position[q+1]) break;
                        if (!allow_paired[offset[p+1]+q]) continue;

                        if (p == i && q == j)
                        {
                            RealT value = Fast_Exp(outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) + FCi_ess[offset[p+1]+q-1]);
                            CountBasePairEvidence(i+1,j,value);
                            CountHelixStacking(i,j+1,value);
                        }
                        else
                        {
                            CountSingleEvidence(i,j,p,q,Fast_Exp(outside + ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]));
                        }
                    }
                }
                
#else
                
                {
                    RealT score_helix = (i+2 <= j ? outside + ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) : 0);
                    RealT score_other = outside + ScoreJunctionB(i,j);
                    
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        const RealT *FCptr = &(FCi_ess[offset[p+1]-1]);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;
                            
                            if (p == i && q == j)
                            {
                                RealT value = Fast_Exp(score_helix + FCptr[q]);
                                cache_score_single[0][0].second += value;
                                CountBasePairEvidence(i+1,j,value);
                                CountHelixStacking(i,j+1,value);
                            }
                            else
                            {
                                RealT value = Fast_Exp(score_other + cache_score_single[p-i][j-q].first + FCptr[q] + ScoreBasePairEvidence(p+1,q) +
                                                       ScoreJunctionB(q,p) + ScoreSingleNucleotidesEvidence(i,j,p,q));
                                cache_score_single[p-i][j-q].second += value;
                                CountBasePairEvidence(p+1,q,value);
                                CountJunctionB(i,j,value);
                                CountJunctionB(q,p,value);
                                CountSingleNucleotidesEvidence(i,j,p,q,value);
                            }
                        }
                    }
                }
#endif
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                
                {
                    RealT value = Fast_Exp(outside + FM2i_ess + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase());
                    CountJunctionA(i,j,value);
                    CountMultiPaired(value);
                    CountMultiBase(value);
                }
            }
            
#endif
            
            // FM1[i,j] = optimal energy for substructure belonging to a
            //            multibranch loop containing a (k+1,j) base pair
            //            preceded by 5' unpaired nucleotides from i to k
            //            for some i <= k <= j-2
            //
            //          = SUM [FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)  if i+2<=j,
            //                 FM1[i+1,j] + b                                          if i+2<=j]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                
                if (allow_paired[offset[i+1]+j])
                {
                    RealT value = Fast_Exp(FM1o_ess[offset[i]+j] + FCi_ess[offset[i+1]+j-1] + ScoreJunctionA(j,i) + ScoreMultiPaired() + ScoreBasePairEvidence(i+1,j) - Z);
                    CountJunctionA(j,i,value);
                    CountMultiPaired(value);
                    CountBasePairEvidence(i+1,j,value);
                }
                
                // compute FM1[i+1,j] + b
                
                if (allow_unpaired_position[i+1])
                {
                    CountMultiUnpairedEvidence(i+1,Fast_Exp(FM1o_ess[offset[i]+j] + FM1i_ess[offset[i+1]+j] + ScoreMultiUnpairedEvidence(i+1) - Z));
                }
            }
            
            // FM[i,j] = optimal energy for substructure belonging to a
            //           multibranch loop which contains at least one 
            //           helix
            //
            //         = SUM [SUM (i<k<j : FM1[i,k] + FM[k,j]),
            //                FM[i,j-1] + b,
            //                FM1[i,j]]
            //
            //            (assuming 0 < i < i+2 <= j < L)
            
            if (0 < i && i+2 <= j && j < L)
            {
                
                // compute SUM (i<k<j : FM1[i,k] + FM[k,j]) -- do nothing
                
                // compute FM[i,j-1] + b
                
                if (allow_unpaired_position[j])
                    CountMultiUnpairedEvidence(j,Fast_Exp(FMo_ess[offset[i]+j] + FMi_ess[offset[i]+j-1] + ScoreMultiUnpairedEvidence(j) - Z));
                
                // compute FM1[i,j] -- do nothing
            }
        }
    }
    
    for (int j = 1; j <= L; j++)
    {
        
        // F5[j] = optimal energy for substructure between positions 0 and j
        //         (or 0 if j = 0)
        //
        //       = SUM [F5[j-1] + ScoreExternalUnpaired(),
        //              SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))]
        
        RealT outside = F5o_ess[j] - Z;
        
        // compute F5[j-1] + ScoreExternalUnpaired()
        
        if (allow_unpaired_position[j])
            CountExternalUnpairedEvidence(j,Fast_Exp(outside + F5i_ess[j-1] + ScoreExternalUnpairedEvidence(j)));
        
        // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
        
        for (int k = 0; k < j; k++)
        {
            if (allow_paired[offset[k+1]+j])
            {
                RealT value = Fast_Exp(outside + F5i_ess[k] + FCi_ess[offset[k+1]+j-1] + ScoreExternalPaired() + ScoreBasePairEvidence(k+1,j) + ScoreJunctionA(j,k));
                CountExternalPaired(value);
                CountBasePairEvidence(k+1,j,value);
                CountJunctionA(j,k,value);
            }      
        }
    }
    
    FinalizeCountsESS();

#if SHOW_TIMINGS
    std::cerr << "Feature expectations (" << GetSystemTime() - starting_time << " seconds)" << std::endl;
#endif

    return GetCounts();
}


template<class RealT>
RealT InferenceEngine<RealT>::ComputeGammaMLESum(std::vector<int> ev_cpd_id, bool ignorePairing, bool usePosterior, int which_data) 
{
    RealT sum1 = 0; 
    int i = ev_cpd_id[0];
    int j = ev_cpd_id[1];

    double score_up = 0;

    if (ignorePairing)
    {
        for (int k = 1; k <= L; k++) {

            score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

            if (s[k] == i) 
            {
                if (usePosterior)
                {
                RealT pos = 0;
                for (int l = 1; l <= L; l++)
                {
                    int offset1 = k <= l ? k : l;
                    int offset2 = k > l ? k : l;
                    if (k != l)
                        pos += posterior[offset[offset1]+offset2];
                }
                if (j == 1)
                {
                   pos = 1-pos;
                }
                sum1 += score_up*pos;
                }
                else
                sum1 += score_up;
            }
        }
    }
    else
    {
        for (int k = 1; k <= L; k++) {

            score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);


            if (s[k] == i && allow_unpaired_position[k] == j)
            {
                if (usePosterior)
                {
                RealT pos = 0;
                for (int l = 1; l <= L; l++)
                {
                    int offset1 = k <= l ? k : l;
                    int offset2 = k > l ? k : l;
                    if (k != l)
                        pos += posterior[offset[offset1]+offset2];
                }
                if (j == 1)
                {
                   pos = 1-pos;
                }
                sum1 += score_up*pos;
                }
                else
                sum1 += score_up;
            }
        }
    }
    return sum1;
}


template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::ComputeGammaMLESS(std::vector<int> ev_cpd_id, bool ignoreZeros, bool useMLE, int which_data) 
{
    // Needs to be called after clamping to true structures so that
    // allowed_unpaired_position is properly set.

    std::vector<RealT> result;

    // If there are zeros in the data, the SS for the method of moments estimator is used (sum d^2)
    RealT sum1 = 0;
    RealT sum2 = 0;
    int N = 0;
    int i = ev_cpd_id[0]; // nucleotide type
    int j = ev_cpd_id[1]; // paired or unpaired

    double score_up = 0;

    for (int k = 1; k <= L; k++)
    {
        score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

        if (s[k] == i && allow_unpaired_position[k] == j)  // 0 = paired, 1 = unpaired
        {
            if (useMLE)  // if using MLE, return log-score unless we need to ignore the zeros and score is 0
            {
                    if (ignoreZeros && score_up == 0)
                        continue;

                    sum1 += score_up;
                    sum2 += log(score_up);
            }
            else
            {
                sum1 += score_up;
                sum2 += score_up*score_up;
            }
            N++;
        }
    }

    result.push_back(sum1);
    result.push_back(sum2);
    result.push_back(N);
    return result;
}


template<class RealT>
std::vector<RealT> InferenceEngine<RealT>::ComputeGammaMLEESS(std::vector<int> ev_cpd_id, bool ignoreZeros, bool useMLE, int which_data) 
{
    std::vector<RealT> result;

    // If there are zeros in the data, the SS for the method of moments estimator is used (sum d^2)
    RealT sum1 = 0;
    RealT sum2 = 0;
    RealT sum3 = 0;
    int i = ev_cpd_id[0]; // base identity
    int j = ev_cpd_id[1]; // pairedness

    double score_up = 0;

    for (int k = 1; k <= L; k++)
    {
        score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

        if (s[k] == i) 
        {

            RealT pos = 0;
            for (int l = 1; l <= L; l++)
            {
                    int offset1 = k <= l ? k : l;
                    int offset2 = k > l ? k : l;
                    if (k != l)
                        pos += posterior[offset[offset1]+offset2];
            }
            if (j == 1)
            {
               pos = 1-pos;
            }

            if (useMLE)  // if using MLE, return log-score unless we need to ignore the zeros and score is 0
            {
                    if (ignoreZeros && score_up <= 0)
                        continue;

                    sum1 += score_up*pos;
                    sum2 += log(score_up)*pos;
            }
            else
            {
                sum1 += score_up*pos;
                sum2 += score_up*score_up*pos;
            }
            sum3 += pos;
        }
    }

    result.push_back(sum1);
    result.push_back(sum2);
    result.push_back(sum3);

    return result;
}


template<class RealT>
RealT InferenceEngine<RealT>::GetNumExamplesSeqPairing(std::vector<int> ev_cpd_id, bool ignoreZeros, int which_data) 
{
    RealT sum = 0; 
    int i = ev_cpd_id[0];
    int j = ev_cpd_id[1];
    double score_up = 0;
    for (int k = 1; k <= L; k++)
    {
        score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

        if (s[k] == i && allow_unpaired_position[k] == j)  // 0 = paired, 1 = unpaired
        {
            if (ignoreZeros && score_up == 0)
                continue;
            sum += 1;
        } 
    }
    return sum;
}

template<class RealT>
RealT InferenceEngine<RealT>::GetNumExamplesSeq(std::vector<int> ev_cpd_id, bool ignoreZeros, int which_data) 
{
    RealT sum = 0; 
    int i = ev_cpd_id[0];
    int j = ev_cpd_id[1];

    double score_up = 0;
    for (int k = 1; k <= L; k++)
    {
        score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

        if (s[k] == i)  // 0 = paired, 1 = unpaired
        {
            if (ignoreZeros && score_up == 0)
                continue;

            RealT pos = 0;
            for (int l = 1; l <= L; l++)
            {
                    int offset1 = k <= l ? k : l;
                    int offset2 = k > l ? k : l;
                    if (k != l)
                        pos += posterior[offset[offset1]+offset2];
            }
            if (j == 1)
            {
               pos = 1-pos;
            }

            sum += pos;
        }
    }
    return sum;

}

template<class RealT>
int InferenceEngine<RealT>::AreZerosInSeqPairing(int id_base, int id_pairing, int which_data) 
{

    int i = id_base;
    int j = id_pairing;

    int areZeros = 0;

    double score_up = 0;
    for (int k = 1; k <= L; k++)
    {
        score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

        if (s[k] == i && allow_unpaired_position[k] == j && score_up == 0)  // 0 = paired, 1 = unpaired
        {
            areZeros += 1;
        }
    }
    return areZeros;
}

template<class RealT>
int InferenceEngine<RealT>::AreZerosInSeq(int id_base, int which_data)  
{

    int i = id_base; 

    int areZeros = 0;

    double score_up = 0;
    for (int k = 1; k <= L; k++)
    {
        score_up = ScoreUnpairedPositionEvidenceRaw(which_data,k);

        if (s[k] == i && score_up == 0)  // 0 = paired, 1 = unpaired
        {
            areZeros += 1;
        }
    }
    return areZeros;
}

template<class RealT>
void InferenceEngine<RealT>::InitRand(unsigned int seed)
{
    if (die) delete die;
    die = new Die(seed);
}

// stochastic traceback algorithm
template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsStochasticTraceback() const
{
    enum { ST_FC, ST_F5, ST_FM, ST_FM1, ST_FE, ST_FN };

    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;

    std::queue<triple<int,int,int> > traceback_queue;
    traceback_queue.push(make_triple(int(ST_F5), 0, L));

    while (!traceback_queue.empty())
    {
        triple<int,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int i = t.second;
        const int j = t.third;
        int L2 = L; // ed to remove max bp dist

        switch (t.first)
        {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case ST_FC:
                break;
            case ST_FE:
                break;
            case ST_FN:
                break;
#else
            case ST_FC:
            {
                if (0 < i && j < L2 && allow_paired[offset[i]+j+1]) // ???
                {
                    Roulette<int,RealT> roulette(*die);
                
                    // compute ScoreHairpin(i,j)
                    if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                        roulette.add(EncodeTraceback(TB_FC_HAIRPIN,0), ScoreHairpin(i,j));
                
                    // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;

                            if (p == i && q == j)
                            {
                                roulette.add(EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q),
                                             ScoreBasePair(i+1,j) + ScoreHelixStacking(i,j+1) +
                                             FCi[offset[p+1]+q-1]);
                            }
                            else
                            {
                                roulette.add(EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q),
                                             ScoreSingle(i,j,p,q) + FCi[offset[p+1]+q-1]);
                            }
                        }
                    }

                    // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                    for (int k=i+1; k < j; k++)
                    {
                        RealT FM2i = FM1i[offset[i]+k] + FMi[offset[k]+j];
                        RealT val = FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase();
                        roulette.add(EncodeTraceback(TB_FC_BIFURCATION, k), val);
                    }

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_FC_HAIRPIN: 
                            break;
                        case TB_FC_SINGLE: 
                        {
                            const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                            const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                            solution[p+1] = q;
                            solution[q] = p+1;
                            traceback_queue.push(make_triple(int(ST_FC), p+1, q-1));
                        }
                        break;
                        case TB_FC_BIFURCATION:
                        {
                            const int k = traceback.second;
                            traceback_queue.push(make_triple(int(ST_FM1), i, k));
                            traceback_queue.push(make_triple(int(ST_FM), k, j));
                        }
                        break;
                    }
                } //else { Assert(!, "unreachable"); }
            } 
            break;
#endif

            case ST_FM:
                if (0 < i && i+2 <= j && j < L2) // ???
                {
                    Roulette<int,RealT> roulette(*die);

                    // compute SUM (i<k<j : FM1[i,k] + FM[k,j]) 
                    for (int k=i+1; k < j; k++)
                    {
                        RealT FM2i = FM1i[offset[i]+k] + FMi[offset[k]+j];
                        roulette.add(EncodeTraceback(TB_FM_BIFURCATION, k), FM2i);
                    }

                    // compute FM[i,j-1] + b
                    if (allow_unpaired_position[j])
                    {
                        roulette.add(EncodeTraceback(TB_FM_UNPAIRED,0),
                                     FMi[offset[i]+j-1] + ScoreMultiUnpaired(j));
                    }

                    // compute FM1[i,j]
                    roulette.add(EncodeTraceback(TB_FM_FM1,0), FM1i[offset[i]+j]);

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_FM_BIFURCATION:
                        {
                            const int k = traceback.second;
                            traceback_queue.push(make_triple(int(ST_FM1), i, k));
                            traceback_queue.push(make_triple(int(ST_FM), k, j));
                        }
                        break;
                        case TB_FM_UNPAIRED:
                        {
                            traceback_queue.push(make_triple(int(ST_FM), i, j-1));
                        }
                        break;
                        case TB_FM_FM1: 
                        {
                            traceback_queue.push(make_triple(int(ST_FM1), i, j));
                        }
                        break;
                    }
                
                } //else { Assert(!, "unreachable"); }
                break;

            case ST_FM1:
                if (0 < i && i+2 <= j && j < L2) // ???
                {
                    Roulette<int,RealT> roulette(*die);

                    // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                    if (allow_paired[offset[i+1]+j])
                    {
                        RealT value = FCi[offset[i+1]+j-1] + ScoreJunctionA(j,i) +
                            ScoreMultiPaired() + ScoreBasePair(i+1,j);
                        roulette.add(EncodeTraceback(TB_FM1_PAIRED, 0), value);
                    }
                
                    // compute FM1[i+1,j] + b
                    if (allow_unpaired_position[i+1])
                    {
                        roulette.add(EncodeTraceback(TB_FM1_UNPAIRED,0),
                                     FM1i[offset[i+1]+j] + ScoreMultiUnpaired(i+1));
                    }

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_FM1_PAIRED:
                        {
                            solution[i+1] = j;
                            solution[j] = i+1;
                            traceback_queue.push(make_triple(int(ST_FC), i+1, j-1));
                        }
                        break;
                        case TB_FM1_UNPAIRED:
                        {
                            traceback_queue.push(make_triple(int(ST_FM1), i+1, j));
                        }
                        break;
                    }
                } //else { Assert(!, "unreachable"); }
                break;

            case ST_F5:
                if (j!=0)
                {
                    Roulette<int,RealT> roulette(*die);

                    // compute F5[j-1] + ScoreExternalUnpaired()
                    if (allow_unpaired_position[j])
                    {
                        roulette.add(EncodeTraceback(TB_F5_UNPAIRED,0),
                                     F5i[j-1] + ScoreExternalUnpaired(j));
                    }
        
                    // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
                    int l = 0; // ed remove max bp dist
                    for (int k = l; k < j; k++)
                    {
                        if (allow_paired[offset[k+1]+j])
                        {
                            RealT value = F5i[k] + FCi[offset[k+1]+j-1] + ScoreExternalPaired() +
                                ScoreBasePair(k+1,j) + ScoreJunctionA(j,k);
                            roulette.add(EncodeTraceback(TB_F5_BIFURCATION,k), value);
                        }      
                    }

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_F5_ZERO:
                            break;
                        case TB_F5_UNPAIRED:
                        {
                            traceback_queue.push(make_triple(int(ST_F5), 0, j-1));
                        }
                        break;
                        case TB_F5_BIFURCATION:
                        {
                            const int k = traceback.second;
                            solution[k+1] = j;
                            solution[j] = k+1;
                            traceback_queue.push(make_triple(int(ST_F5), 0, k));
                            traceback_queue.push(make_triple(int(ST_FC), k+1, j-1));
                        }
                        break;
                    }
                }
                break;

            default:
                break;
        }
    }

    return solution;
}

// stochastic traceback algorithm with chemical mapping data
template<class RealT>
std::vector<int> InferenceEngine<RealT>::PredictPairingsStochasticTracebackESS() const
{
    enum { ST_FC, ST_F5, ST_FM, ST_FM1, ST_FE, ST_FN };

    std::vector<int> solution(L+1,SStruct::UNPAIRED);
    solution[0] = SStruct::UNKNOWN;

    std::queue<triple<int,int,int> > traceback_queue;
    traceback_queue.push(make_triple(int(ST_F5), 0, L));

    while (!traceback_queue.empty())
    {
        triple<int,int,int> t = traceback_queue.front();
        traceback_queue.pop();
        const int i = t.second;
        const int j = t.third;
        int L2 = L; // ed to remove max bp dist

        switch (t.first)
        {
#if PARAMS_HELIX_LENGTH || PARAMS_ISOLATED_BASE_PAIR
            case ST_FC:
                break;
            case ST_FE:
                break;
            case ST_FN:
                break;
#else
            case ST_FC:
            {
                if (0 < i && j < L2 && allow_paired[offset[i]+j+1]) // ???
                {
                    Roulette<int,RealT> roulette(*die);
                
                    // compute ScoreHairpin(i,j)
                    if (allow_unpaired[offset[i]+j] && j-i >= C_MIN_HAIRPIN_LENGTH)
                        roulette.add(EncodeTraceback(TB_FC_HAIRPIN,0), ScoreHairpinEvidence(i,j));
                
                    // compute SUM (i<=p<p+2<=q<=j : ScoreSingle(i,j,p,q) + FC[p+1,q-1])
                    for (int p = i; p <= std::min(i+C_MAX_SINGLE_LENGTH,j); p++)
                    {
                        if (p > i && !allow_unpaired_position[p]) break;
                        int q_min = std::max(p+2,p-i+j-C_MAX_SINGLE_LENGTH);
                        for (int q = j; q >= q_min; q--)
                        {
                            if (q < j && !allow_unpaired_position[q+1]) break;
                            if (!allow_paired[offset[p+1]+q]) continue;

                            if (p == i && q == j)
                            {
                                roulette.add(EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q),
                                             ScoreBasePairEvidence(i+1,j) + ScoreHelixStacking(i,j+1) +
                                             FCi_ess[offset[p+1]+q-1]);
                            }
                            else
                            {
                                roulette.add(EncodeTraceback(TB_FC_SINGLE,(p-i)*(C_MAX_SINGLE_LENGTH+1)+j-q),
                                             ScoreSingleEvidence(i,j,p,q) + FCi_ess[offset[p+1]+q-1]);
                            }
                        }
                    }

                    // compute SUM (i<k<j : FM1[i,k] + FM[k,j] + ScoreJunctionA(i,j) + a + c)
                    for (int k=i+1; k < j; k++)
                    {
                        RealT FM2i = FM1i_ess[offset[i]+k] + FMi_ess[offset[k]+j];
                        RealT val = FM2i + ScoreJunctionA(i,j) + ScoreMultiPaired() + ScoreMultiBase();
                        roulette.add(EncodeTraceback(TB_FC_BIFURCATION, k), val);
                    }

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_FC_HAIRPIN: 
                            break;
                        case TB_FC_SINGLE: 
                        {
                            const int p = i + traceback.second / (C_MAX_SINGLE_LENGTH+1);
                            const int q = j - traceback.second % (C_MAX_SINGLE_LENGTH+1);
                            solution[p+1] = q;
                            solution[q] = p+1;
                            traceback_queue.push(make_triple(int(ST_FC), p+1, q-1));
                        }
                        break;
                        case TB_FC_BIFURCATION:
                        {
                            const int k = traceback.second;
                            traceback_queue.push(make_triple(int(ST_FM1), i, k));
                            traceback_queue.push(make_triple(int(ST_FM), k, j));
                        }
                        break;
                    }
                } //else { Assert(!, "unreachable"); }
            } 
            break;
#endif

            case ST_FM:
                if (0 < i && i+2 <= j && j < L2) // ???
                {
                    Roulette<int,RealT> roulette(*die);

                    // compute SUM (i<k<j : FM1[i,k] + FM[k,j]) 
                    for (int k=i+1; k < j; k++)
                    {
                        RealT FM2i = FM1i_ess[offset[i]+k] + FMi_ess[offset[k]+j];
                        roulette.add(EncodeTraceback(TB_FM_BIFURCATION, k), FM2i);
                    }

                    // compute FM[i,j-1] + b
                    if (allow_unpaired_position[j])
                    {
                        roulette.add(EncodeTraceback(TB_FM_UNPAIRED,0),
                                     FMi_ess[offset[i]+j-1] + ScoreMultiUnpairedEvidence(j));
                    }

                    // compute FM1[i,j]
                    roulette.add(EncodeTraceback(TB_FM_FM1,0), FM1i_ess[offset[i]+j]);

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_FM_BIFURCATION:
                        {
                            const int k = traceback.second;
                            traceback_queue.push(make_triple(int(ST_FM1), i, k));
                            traceback_queue.push(make_triple(int(ST_FM), k, j));
                        }
                        break;
                        case TB_FM_UNPAIRED:
                        {
                            traceback_queue.push(make_triple(int(ST_FM), i, j-1));
                        }
                        break;
                        case TB_FM_FM1: 
                        {
                            traceback_queue.push(make_triple(int(ST_FM1), i, j));
                        }
                        break;
                    }
                
                } //else { Assert(!, "unreachable"); }
                break;

            case ST_FM1:
                if (0 < i && i+2 <= j && j < L2) // ???
                {
                    Roulette<int,RealT> roulette(*die);

                    // compute FC[i+1,j-1] + ScoreJunctionA(j,i) + c + ScoreBP(i+1,j)
                    if (allow_paired[offset[i+1]+j])
                    {
                        RealT value = FCi_ess[offset[i+1]+j-1] + ScoreJunctionA(j,i) +
                            ScoreMultiPaired() + ScoreBasePairEvidence(i+1,j);
                        roulette.add(EncodeTraceback(TB_FM1_PAIRED, 0), value);
                    }
                
                    // compute FM1[i+1,j] + b
                    if (allow_unpaired_position[i+1])
                    {
                        roulette.add(EncodeTraceback(TB_FM1_UNPAIRED,0),
                                     FM1i_ess[offset[i+1]+j] + ScoreMultiUnpaired(i+1));
                    }

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_FM1_PAIRED:
                        {
                            solution[i+1] = j;
                            solution[j] = i+1;
                            traceback_queue.push(make_triple(int(ST_FC), i+1, j-1));
                        }
                        break;
                        case TB_FM1_UNPAIRED:
                        {
                            traceback_queue.push(make_triple(int(ST_FM1), i+1, j));
                        }
                        break;
                    }
                } //else { Assert(!, "unreachable"); }
                break;

            case ST_F5:
                if (j!=0)
                {
                    Roulette<int,RealT> roulette(*die);

                    // compute F5[j-1] + ScoreExternalUnpaired()
                    if (allow_unpaired_position[j])
                    {
                        roulette.add(EncodeTraceback(TB_F5_UNPAIRED,0),
                                     F5i_ess[j-1] + ScoreExternalUnpairedEvidence(j));
                    }
        
                    // compute SUM (0<=k<j : F5[k] + FC[k+1,j-1] + ScoreExternalPaired() + ScoreBP(k+1,j) + ScoreJunctionA(j,k))
                    int l = 0; // ed remove max bp dist
                    for (int k = l; k < j; k++)
                    {
                        if (allow_paired[offset[k+1]+j])
                        {
                            RealT value = F5i_ess[k] + FCi_ess[offset[k+1]+j-1] + ScoreExternalPaired() +
                                ScoreBasePairEvidence(k+1,j) + ScoreJunctionA(j,k);
                            roulette.add(EncodeTraceback(TB_F5_BIFURCATION,k), value);
                        }      
                    }

                    // choose
                    std::pair<int,int> traceback = DecodeTraceback(roulette.choose());
                    switch (traceback.first)
                    {
                        case TB_F5_ZERO:
                            break;
                        case TB_F5_UNPAIRED:
                        {
                            traceback_queue.push(make_triple(int(ST_F5), 0, j-1));
                        }
                        break;
                        case TB_F5_BIFURCATION:
                        {
                            const int k = traceback.second;
                            solution[k+1] = j;
                            solution[j] = k+1;
                            traceback_queue.push(make_triple(int(ST_F5), 0, k));
                            traceback_queue.push(make_triple(int(ST_FC), k+1, j-1));
                        }
                        break;
                    }
                }
                break;

            default:
                break;
        }
    }

    return solution;
}


