#include <iostream>
#include <set>
#include <map>
#include <cassert>
#include <cctype>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <algorithm>

using namespace std;

struct Projection {
  vector<pair<char,int> > mapping;
  double quality;
};

struct Family {
  string name;
  string annot;
  bool predicted;
  map<string,string> data;
  map<string,Projection> full_projections;
  map<string,Projection> partial_projections;
  
  Family() : predicted (false) {}
  Family (const Family &rhs) : 
    name(rhs.name), annot(rhs.annot), predicted(rhs.predicted), data(rhs.data){}
};

vector<Family> rfam;
char group = 'A';
int numA = 0, numB = 0;

////////////////////////////////////////////////////////////
// NonCanonical()
////////////////////////////////////////////////////////////

int NonCanonical (char c, char d){
  if (c == 'A' && d == 'U') return false;
  if (c == 'U' && d == 'A') return false;
  if (c == 'G' && d == 'U') return false;
  if (c == 'U' && d == 'G') return false;
  if (c == 'C' && d == 'G') return false;
  if (c == 'G' && d == 'C') return false;
  return true;
}

////////////////////////////////////////////////////////////
// Trim()
////////////////////////////////////////////////////////////

string Trim (const string &s){
  int left = 0, right = s.length();
  while (left < right && isspace(s[left])) left++;
  while (left < right && isspace(s[right-1])) right--;
  return s.substr(left,right-left);
}

////////////////////////////////////////////////////////////
// ReadData()
////////////////////////////////////////////////////////////

void ReadData (const string &filename){
  ifstream data (filename.c_str());

  // read entries

  while (true){
    string s, t;
    getline (data, s);
    s = Trim(s);
    if (s != "# STOCKHOLM 1.0") break;

    // read entry

    Family family;

    while (getline (data, s)){
      s = Trim(s);
      if (s == "") continue;
      if (s == "//") break;
      if (s[0] == '#'){
        if (!strncmp ("#=GF AC", s.c_str(), 7)){
          family.name = s.substr(10);
	} else if (!strncmp ("#=GF SS   Pred", s.c_str(), 14)){
	  family.predicted = true;
        } else if (!strncmp ("#=GC SS", s.c_str(), 7)){
          istringstream iss (s);
          iss >> s >> s >> t;
          family.annot += t;
        }
      } else {
        istringstream iss (s);
        iss >> s >> t;
        family.data[s] += t;        
      }
    }

    // upper-case letters in sequence and convert unknowns to Ns
    
    for (map<string,string>::iterator iter = family.data.begin();
         iter != family.data.end(); ++iter){
      string &sequence = iter->second;
      for (int j = 0; j < sequence.length(); j++){
        sequence[j] = toupper(sequence[j]);
        if (isalpha(sequence[j]) && 
            sequence[j] != 'A' && 
            sequence[j] != 'G' && 
            sequence[j] != 'U' && 
            sequence[j] != 'C') 
          sequence[j] = 'N';
      }
    }
    
    rfam.push_back (family);
  }
  data.close();
}

////////////////////////////////////////////////////////////
// IsLeftPairing()
// IsRightPairing()
// PairIndex()
////////////////////////////////////////////////////////////

bool IsLeftPairing (char c, bool allow_pseudoknot){
  if (!allow_pseudoknot) return (c == '<');
  return (c == '<' || isupper(c));
}

bool IsRightPairing (char c, bool allow_pseudoknot){
  if (!allow_pseudoknot) return (c == '>');
  return (c == '>' || islower(c));
}

unsigned int PairIndex (char c){
  if (c == '>') c = '<';
  return (unsigned int) toupper(c);
}

////////////////////////////////////////////////////////////
// ComputeMapping()
////////////////////////////////////////////////////////////

vector<int> ComputeMapping (const string &sstruct, bool allow_pseudoknot){
  vector<vector<int> > stack (256);
  vector<int> mapping;
  
  for (int i = 0; i < sstruct.length(); i++){
    char c = sstruct[i];

    if (IsRightPairing(c, allow_pseudoknot)){
      assert (stack[PairIndex(c)].size() > 0);
      mapping[stack[PairIndex(c)].back()] = mapping.size();
      mapping.push_back (stack[PairIndex(c)].back());
      stack[PairIndex(c)].pop_back();
    } else {
      mapping.push_back (0);
      if (IsLeftPairing(c, allow_pseudoknot))
        stack[PairIndex(c)].push_back(i);
    }
  }

  for (int i = 0; i < 256; i++)
    assert (stack[i].size() == 0);
  
  return mapping;  
}

////////////////////////////////////////////////////////////
// ComputePositionNumbers()
////////////////////////////////////////////////////////////

vector<int> ComputePositionNumbers (const string &sequence){
  vector<int> numbers (1);
  for (int i = 1; i < sequence.length(); i++){
    numbers.push_back (numbers.back() + (isalpha(sequence[i]) != 0));
  }    
  return numbers;
}

////////////////////////////////////////////////////////////
// ProjectMapping()
////////////////////////////////////////////////////////////

vector<pair<char,int> > ProjectMapping (const vector<int> &mapping, const string &sequence, const vector<int> &numbers){
  vector<pair<char,int> > ret (1, make_pair('@',0));
  for (int i = 1; i < sequence.length(); i++){
    if (isalpha(sequence[i])){
      ret.push_back (make_pair(sequence[i], isalpha(sequence[mapping[i]]) ? numbers[mapping[i]] : 0));
    }
  }  
  return ret;
}

////////////////////////////////////////////////////////////
// ComputeMappingScore()
////////////////////////////////////////////////////////////

double ComputeMappingScore(const vector<pair<char,int> > &mapping){
  int numNs = 0, numNCs = 0;
  for (int i = 1; i < mapping.size(); i++){
    numNs += (mapping[i].first == 'N');
    numNCs += (mapping[i].second > i && NonCanonical(mapping[i].first, mapping[mapping[i].second].first));
  }
  return (double)(numNs + numNCs) / mapping.size();
}

////////////////////////////////////////////////////////////
// Pad()
////////////////////////////////////////////////////////////

string Pad (int i, int width){
  ostringstream oss;
  oss << i;
  string s = oss.str();
  while (s.length() < width) s = "0" + s;
  return s;
}

////////////////////////////////////////////////////////////
// MakeBPSEQ()
////////////////////////////////////////////////////////////

void MakeBPSEQ (Family &family){

  if (family.predicted){ 
    cerr << "Skipping predicted family: " << family.name << endl;
    return;
  }
  
  // convert annotation to mapping
  
  string sstruct = "@" + family.annot;
  vector<int> full_mapping = ComputeMapping (sstruct, true);
  vector<int> partial_mapping = ComputeMapping (sstruct, false);

  vector<pair<double, string> > scores;
  
  // compute all mappings
  
  for (map<string,string>::const_iterator iter = family.data.begin();
       iter != family.data.end(); ++iter){
    
    string sequence = "@" + iter->second;
    vector<int> numbers = ComputePositionNumbers(sequence);
    
    vector<pair<char,int> > full_projected_mapping = ProjectMapping (full_mapping, sequence, numbers);
    vector<pair<char,int> > partial_projected_mapping = ProjectMapping (partial_mapping, sequence, numbers);
    
    double full_score = ComputeMappingScore (full_projected_mapping);
    double partial_score = ComputeMappingScore (partial_projected_mapping);
    
    family.full_projections[iter->first].quality = full_score;
    family.full_projections[iter->first].mapping = full_projected_mapping;
    family.partial_projections[iter->first].quality = partial_score;
    family.partial_projections[iter->first].mapping = partial_projected_mapping;

    scores.push_back (make_pair(full_score,iter->first));
  }

  // pick out lowest scores subject to increasing percent identity cutoff

  sort (scores.begin(), scores.end());
  string selected = scores[0].second;

  // build sequences for family
  
  cerr << "Building equences for family: " << family.name << endl;
  
  ostringstream oss;
  oss << "conv/full_bpseq/" << family.name << "_" << group << ".bpseq";
  ofstream outfile (oss.str().c_str());
  const vector<pair<char,int> > &mapping = family.full_projections[selected].mapping;
  for (int j = 1; j < mapping.size(); j++)
    outfile << j << " " << mapping[j].first << " " << mapping[j].second << endl;
  outfile.close();
  outfile.clear();
  
  oss.clear(); oss.str("");
  oss << "conv/partial_bpseq/" << family.name << "_" << group << ".bpseq";
  outfile.open (oss.str().c_str());
  const vector<pair<char,int> > &mapping2 = family.partial_projections[selected].mapping;    
  for (int j = 1; j < mapping2.size(); j++)
    outfile << j << " " << mapping2[j].first << " " << mapping2[j].second << endl;
  outfile.close();
  
  // build mfa

  oss.str ("");
  oss.clear();
  oss << "conv/full_mfa/" << family.name << "_" << group << ".mfa";
  outfile.clear();
  outfile.open (oss.str().c_str());
  for (map<string,string>::const_iterator iter = family.data.begin(); iter != family.data.end(); ++iter)
    outfile << ">" << iter->first << endl << iter->second << endl;
  outfile << ">secondary structure" << endl  << family.annot << endl;
  outfile.close();
  outfile.clear();
  
  oss.clear(); oss.str("");
  oss << "conv/partial_mfa/" << family.name << "_" << group << ".mfa";
  outfile.open (oss.str().c_str());
  for (map<string,string>::const_iterator iter = family.data.begin(); iter != family.data.end(); ++iter)
    outfile << ">" << iter->first << endl << iter->second << endl;
  outfile << ">secondary structure" << endl;
  for (int i = 0; i < family.annot.length(); i++){
    if (family.annot[i] == '<' || family.annot[i] == '>') 
      outfile << family.annot[i];
    else
      outfile << '.';    
  }
  outfile << endl;
  outfile.close();
  
  // switch groups
  
  if (group == 'A') 
    numA++;
  else
    numB++;
  
  if (group == 'A' && numA > numB)
    group = 'B';
  else if (group == 'B' && numB > numA)
    group = 'A';
}

////////////////////////////////////////////////////////////
// main()
////////////////////////////////////////////////////////////

int main (int argc, char **argv){
  
  if (argc != 2){
    cerr << "Usage: convert-rfam RFAM_SEED_FILE" << endl;
    exit(1);
  }

  ReadData (argv[1]);

  system ("mkdir conv");
  system ("mkdir conv/full_bpseq");
  system ("mkdir conv/partial_bpseq");
  system ("mkdir conv/full_mfa");
  system ("mkdir conv/partial_mfa");
  
  for (int i = 0; i < rfam.size(); i++){
    cerr << "Converting " << rfam[i].name << "..." << endl;
    MakeBPSEQ (rfam[i]);
  }
}
