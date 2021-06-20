# Setting up LinearFold/LinearPartition patches

The EternaFold parameters have also been adapted for use with the LinearFold and LinearPartition algorithms, described in


LinearFold: Linear-Time Approximate RNA Folding by 5’-to-3’ Dynamic Programming and Beam Search. Bioinformatics, Volume 35, Issue 14, July 2019, Pages i295–i304. ISMB 2019

Liang Huang, He Zhang, Dezhong Deng, Kai Zhao, Kaibo Liu, David Hendrix, David Mathews

LinearPartition: linear-time approximation of RNA folding partition function and base-pairing probabilities. Bioinformatics, Volume 36, Issue Supplement_1, July 2020, Pages i258–i267. ISMB 2020

He Zhang, Liang Zhang, David Mathews, Liang Huang

LinearFold and LinearPartition have different licenses than EternaFold, please read the LinearFold [license](https://github.com/LinearFold/LinearFold/blob/master/LICENSE) and LinearPartition [license](https://github.com/LinearFold/LinearPartition/blob/master/LICENSE) before proceeding.

1. Clone the LinearFold repository at [https://github.com/LinearFold/LinearFold](https://github.com/LinearFold/LinearFold). The most recently-tested working commit is 260c6bbb9bf8cc84b807fa7633b9cb731e639884 (June 06 2021). You can get this commit with

```
git clone https://github.com/LinearFold/LinearFold.git
cd LinearFold
git reset --hard 260c6bbb9bf8cc84b807fa7633b9cb731e639884
```

2. Clone the LinearPartition repository at [https://github.com/LinearFold/LinearPartition](https://github.com/LinearFold/LinearPartition). The most recently-tested working commit is ae6507f3053573decd2e4bdae60d5a96eac87783 (May 19 2021). You can get this commit with

```
git clone https://github.com/LinearFold/LinearPartition.git
cd LinearPartition
git reset --hard ae6507f3053573decd2e4bdae60d5a96eac87783
```

3. To apply the LinearFold patch and compile:

```
cd LinearFold
git apply --whitespace=fix /PATH/TO/ETERNAFOLD/LinearFold-E.patch
make
cd ..
```

4. To apply the LinearPartition patch and compile:

```
cd LinearPartition
git apply --whitespace=fix /PATH/TO/ETERNAFOLD/LinearPartition-E.patch
make
cd ..
```

3. Add the following lines to your Arnie file:

```
linearfold: /path/to/LinearFold/bin
linearpartition: /path/to/LinearPartition/bin
```

4. To check LinearFold executables are properly built:

```
$ echo CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG | ./bin/linearfold_v

CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG
((((((.((((((......)))))).......((((.....))))...)))))) (-15.20)

$ echo CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG | ./bin/linearfold_c

CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG
((((((.((((((......)))))).......((((.....))))...)))))) (4.27)

$ echo CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG | ./bin/linearfold_e

CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG
((((((.((((((......)))))).......((((.....))))...)))))) (10.94)
```

5. To check LinearPartition executables are properly built:

```
$ echo CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG | ./bin/linearpartition_v
Free Energy of Ensemble: -15.92 kcal/mol

$ echo CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG | ./bin/linearpartition_c
Log Partition Coefficient: 6.77346

$ echo CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG | ./bin/linearpartition_e
Log Partition Coefficient: 13.28986
```

Interpreting this output:

LinearFold returns the minimum free energy structure (MFE) and its free energy (dG(MFE)) for that energy model. For LinearFold-C and LinearFold-E, these are positive numbers because this is technically the log partition coefficient. Negative of this is the free energy. (Arnie automatically converts this.)

LinearPartition calculates dG(ensemble), again positive numbers from LinearPartition-C and LinearPartition-E.

6. To check LinearFold and LinearPartition are correctly linked to Arnie:

```python
from arnie.mfe import mfe
from arnie.free_energy import free_energy

seq = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'

# Arnie mfe utility calls LinearFold for structure. Arnie free_energy utility calls LinearPartition for free energy.

print('Vienna RNAfold:')
print(mfe(seq))
print(free_energy(seq))

print('\nLinearFold-V:')
print(mfe(seq, linear=True))
print(free_energy(seq, linear=True))

print('\nCONTRAfold (v2_02):')
print(mfe(seq, package='contrafold', viterbi=True)) # setting viterbi=True because CONTRAfold default is MEA structure, not MFE structure
print(free_energy(seq, package='contrafold'))

print('\nLinearFold-C:')
print(mfe(seq, linear=True, package='contrafold'))
print(free_energy(seq, linear=True, package='contrafold',beam_size=100000000))
```

Expected output:

```
Vienna RNAfold:
((((((.((((((......)))))).......((((.....))))...))))))
-15.92

LinearFold-V:
((((((.((((((......)))))).......((((.....))))...))))))
-15.92

CONTRAfold (v2_02):
((((((.((((((......)))))).......((((.....))))...))))))
-6.87394

LinearFold-C:
((((((.((((((......)))))).......((((.....))))...))))))
-6.77342
```

CONTRAfold (v2_02) and LinearPartition-C free energies are known to differ. Fixing multiloop bug in CONTRAfold v2_02 brings free energy to -6.83585, but that still does not match LinearPartition-C free energy.

