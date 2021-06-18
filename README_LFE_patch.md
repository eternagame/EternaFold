# Setup for LinearFold-E patch

The EternaFold parameters have also been adapted for use with the LinearFold algorithm, described in

Liang Huang, He Zhang, Dezhong Deng, Kai Zhao, Kaibo Liu,
David Hendrix, and David Mathews (2019). LinearFold: Linear-Time
Approximate RNA Folding by 5'-to-3' Dynamic Programming and Beam Search.
Bioinformatics, Vol. 35, July 2019, Special Issue of ISMB 2019 Proceedings.

LinearFold has a different license than EternaFold, please read the LinearFold [license](https://github.com/LinearFold/LinearFold/blob/master/LICENSE) before proceeding.

1. Clone the LinearFold repository at [https://github.com/LinearFold/LinearFold](https://github.com/LinearFold/LinearFold) :

```
git clone https://github.com/LinearFold/LinearFold.git
```

2. To apply the patch and compile:

```
cd LinearFold
git apply --whitespace=fix /path/to/EternaFold/LinearFold-E.patch
make
```

3. Add the following line to your Arnie file:

`linearfold: /path/to/LinearFold/bin`

4. To check LinearFold is correctly linked to Arnie:

```python
from arnie.mfe import mfe
ref_seq = 'AUGGCCGUUUACCCAUACGAUGUUCCUGACUAUGCGGGCUAUCCCUAUGACGUCCCGGACUAUGCAGGCUCCUAUCCAUAUGACGUUCCAGAUUACGCUGGAUCUGGCGUCUUCACACUCGAAGAUUUCGUUGGGGACUGGCGACAGACAGCCGGCUACAACCUGGACCAAGUCCUUGAACAGGGAGGUGUGUCCAGUUUGUUUCAGAAUCUCGGGGUGUCCGUAACUCCGAUCCAAAGGAUUGUCCUGAGCGGUGAAAAUGGGCUGAAGAUCGACAUCCAUGUCAUCAUCCCGUAUGAAGGUCUGAGCGGCGACCAAAUGGGCCAGAUCGAAAAAAUUUUUAAGGUGGUGUACCCUGUGGAUGAUCAUCACUUUAAGGUGAUCCUGCACUAUGGCACACUGGUAAUCGACGGGGUUACGCCGAACAUGAUCGACUAUUUCGGACGGCCGUAUGAAGGCAUCGCCGUGUUCGACGGCAAAAAGAUCACUGUAACAGGGACCCUGUGGAACGGCAACAAAAUUAUCGACGAGCGCCUGAUCAACCCCGACGGCUCCCUGCUGUUCCGAGUAACCAUCAACGGAGUGACCGGCUGGCGGCUGUGCGAACGCAUUCUGGCGUAA'

print('Vienna RNAfold:')
print(mfe(ref_seq))

print('\nLinearFold-V:')
print(mfe(ref_seq, linear=True))

print('\nEternaFold:')
print(mfe(ref_seq, package='eternafold'))

print('\nLinearFold-E:')
print(mfe(ref_seq, linear=True, package='eternafold'))
```

Expected output:

```
Vienna RNAfold:
...((((((...(((((....(((((((....(((((((.(((((...((((((..(((...............)))....))))))(((((((......))))))).((((((......))))))...((((.((.(((((........))))))).))))((((((((..(((((....)))))..)).))))))................))))))))))))....((((((...))))))......(((((((..(((((.(((.(((..........))).))).)))))..((..((((..(((((..((((((((....(((((...............(((((((((((((.((((((((((.......))))))).(((......)))........))).)))))).)))))))....)))))..)))))).))...)))))....)))).))(((((.....))))).......)))))))..)))))))..))))))))))).................(((((.((.......((((((((.((.(((((((((............)))))....)))).)).)))))).)).......)).)))))..

LinearFold-V:
...((((((...(((((....(((((((....(((((((.(((((...((((((..(((...............)))....))))))(((((((......))))))).((((((......))))))...((((.((.(((((........))))))).))))((((((((..(((((....)))))..)).))))))................))))))))))))....((((((...))))))......(((((((..(((((.(((.....((((....)))).))).)))))..((..((((..(((((..((((((((....(((((...............(((((((((((((.((((((((((.......))))))).(((......)))........))).))))))).))))))....)))))..)))))).))...)))))....)))).))(((((.....))))).......)))))))..)))))))..))))))))))).................(((((.((.......((((((((.((.(((((((((............)))))....)))).)).)))))).)).......)).)))))..

EternaFold:
((((((((...((((.((((....((((......(((((..((.....)).).)))).......))))...............((..(((((((......)))))))))(((((......)))))..))))))))..(((....))).(((((((..((((.((((((((..(((((....)))))..)).)))))).))))........(((((((((((.....(((((.((((.(((...((((((.(((((((..(((((.(...(((.((((....)))))))).)))))..(...(((.(((((((((((..((.((((.(.((((((..........(.(((((((((((((.((((((((((.......)))))))..((......)).........))).))))))).))))))..)...........)))))).)))))..))...))..)))))))..)))))..).......)))))))..))).)))))).)))).))).................)).....)))..))))))).)((.....)).((((((.(.......)..))))))..)))))))))))))))....((((......))))..

LinearFold-E:
.(((((.....((((.((((....((.((......(((....)))....(((((..(((....(......)...)))....)))))(((((......))))))).)).((((((......)))))).))))))))..(((((........)))))..((((.((((((((..(((((....)))))..)).)))))).))))((((......((((((((((((....(((.((...(((.....))))).))).....))))))....(((.((((....)))))))))))))..))))((((........)))).....))))).....................(((((.((..(((...(((((((.......)))))))...)))...)).))))).....((.((((((((.(.....((((((((..........((..(((.......))).))(((((.....))))).....))))).)))....).)))))))).))......................(((((.((.......((((((((.((.((((.((((............)))).....)))).)).)))))).)).......)).)))))..
```
(NB: The EternaFold and LinearFold-E structures differ more than the Vienna/LinearFold-V structures because EternaFold default structure is the MEA structure, whereas the other 3 are MFE (Viterbi) structures).
