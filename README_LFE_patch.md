# Setup for LinearFold-E patch

1. Clone the LinearFold repository at [https://github.com/LinearFold/LinearFold](https://github.com/LinearFold/LinearFold) :

```python
git clone [https://github.com/LinearFold/LinearFold.git](https://github.com/LinearFold/LinearFold.git)
```

2. To apply the patch and compile:

```python
cd LinearFold
git apply --whitespace=fix /path/to/EternaFold/LinearFold-E.patch
make
```

3. Add the following line to your arnie file:

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

```python
Vienna RNAfold:
...((((((...(((((....(((((((....(((((((.(((((...((((((..(((...............)))....))))))(((((((......))))))).((((((......))))))...((((.((.(((((........))))))).))))((((((((..(((((....)))))..)).))))))................))))))))))))....((((((...))))))......(((((((..(((((.(((.(((..........))).))).)))))..((..((((..(((((..((((((((....(((((...............(((((((((((((.((((((((((.......))))))).(((......)))........))).)))))).)))))))....)))))..)))))).))...)))))....)))).))(((((.....))))).......)))))))..)))))))..))))))))))).................(((((.((.......((((((((.((.(((((((((............)))))....)))).)).)))))).)).......)).)))))..

LinearFold-V:
...((((((...(((((....(((((((....(((((((.(((((...((((((..(((...............)))....))))))(((((((......))))))).((((((......))))))...((((.((.(((((........))))))).))))((((((((..(((((....)))))..)).))))))................))))))))))))....((((((...))))))......(((((((..(((((.(((.....((((....)))).))).)))))..((..((((..(((((..((((((((....(((((...............(((((((((((((.((((((((((.......))))))).(((......)))........))).))))))).))))))....)))))..)))))).))...)))))....)))).))(((((.....))))).......)))))))..)))))))..))))))))))).................(((((.((.......((((((((.((.(((((((((............)))))....)))).)).)))))).)).......)).)))))..

EternaFold:
((((((((...((((.((((....((((......(((((..((.....)).).)))).......))))...............((..(((((((......)))))))))(((((......)))))..))))))))..(((....))).(((((((..((((.((((((((..(((((....)))))..)).)))))).))))........(((((((((((.....(((((.((((.(((...((((((.(((((((..(((((.(...(((.((((....)))))))).)))))..(...(((.(((((((((((..((.((((.(.((((((..........(.(((((((((((((.((((((((((.......)))))))..((......)).........))).))))))).))))))..)...........)))))).)))))..))...))..)))))))..)))))..).......)))))))..))).)))))).)))).))).................)).....)))..))))))).)((.....)).((((((.(.......)..))))))..)))))))))))))))....((((......))))..

LinearFold-E:
.(((((.....((((.((((....((.((......(((....)))....(((((..(((....(......)...)))....)))))(((((......))))))).)).((((((......)))))).))))))))..(((((........)))))..((((.((((((((..(((((....)))))..)).)))))).))))((((......((((((((((((....(((.((...(((.....))))).))).....))))))....(((.((((....)))))))))))))..))))((((........)))).....))))).....................(((((.((..(((...(((((((.......)))))))...)))...)).))))).....((.((((((((.(.....((((((((..........((..(((.......))).))(((((.....))))).....))))).)))....).)))))))).))......................(((((.((.......((((((((.((.((((.((((............)))).....)))).)).)))))).)).......)).)))))..
```