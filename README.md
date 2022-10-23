## EternaFold - Improving RNA structure prediction through multitask learning

An EternaFold server is available at [eternafold.eternagame.org](https://eternafold.eternagame.org/)!

EternaFold performs multitask learning to improve RNA structure prediction. Its training tasks include 1) predicting single structures, 2) maximizing the likelihood of structure probing data, and 3) predicting experimentally-measured affinities of RNA molecules to proteins and small molecules.

Its training data comes from diverse high-throughput experimental crowdsourced data from the [Eterna](www.eternagame.org) project.

EternaFold is possible thanks to [CONTRAfold-SE](https://github.com/csfoo/contrafold-se) (C.-S. Foo, C. Pop).

For scripts and datasets pertaining to benchmarking EternaFold on secondary structure prediction tasks, see the [EternaBench](https://github.com/eternagame/EternaBench) repo.

If you use EternaFold in your research, please cite [the paper](https://www.nature.com/articles/s41592-022-01605-0):

```
H.K. Wayment-Steele, W. Kladwang, A.I. Strom, J. Lee, A. Treuille, A. Becka,
Eterna Participants, R. Das. (2022). RNA secondary structure packages ranked
and improved by high-throughput experiments. Nature Methods 19, 1234–1242.
```

### Installation

Clone the repository and run `make` in `src` to compile.
Multithreaded version: run `make multi` in `src`.

Compiled with gcc 4.8.5 and openmpi 2.0.2. 

See instructions in [README_LinearFold-E_patch.md](README_LinearFold-E_patch.md) for using EternaFold parameters with LinearFold and LinearPartition algorithms.

### Prediction

#### Single-structure prediction
Predict the MEA structure of example test sequence (Hammerhead ribozyme), using the EternaFold parameters:

`./src/contrafold predict test.seq --params parameters/EternaFoldParams.v1`

Output:
```
Training mode: 
Use constraints: 0
Use evidence: 0
Predicting using MEA estimator.
>test.seq
CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG
>structure
(((((((((((((......))))))..)....((((.....))))...))))))
```

#### SHAPE-directed secondary structure prediction

Predict the maximum-likelihood structure of the given sequence, using SHAPE likelihood potentials learned from Cloud Lab SHAPE MAP-seq experiments (Wayment-Steele et. al, 2022).

Predicted structure of example construct without incorporating SHAPE data:
```
./src/contrafold predict test_SHAPE.seq --params parameters/EternaFoldParams.v1
```
Output:
``` 
Training mode: 
Use constraints: 0
Use evidence: 0
Predicting using MEA estimator.
>test_SHAPE.seq
UGUACCGGAAGGUGCGAAUCUUCCG
>structure
.....((((((((....))))))))
```

Alternate structure is predicted upon incorporating SHAPE data in `test_SHAPE.bpseq`:

```
./src/contrafold predict test_SHAPE.bpseq --evidence --numdatasources 1 --kappa 0.1 --params parameters/EternaFoldParams_PLUS_POTENTIALS.v1 
```

Output:
```
Training mode: 
Use constraints: 0
Use evidence: 1
Predicting using MEA estimator.
>test_SHAPE.bpseq
UGUACCGGAAGGUGCGAAUCUUCCG
>structure
((((((....)))))).........
``` 



#### Ensemble free energy prediction

```
$ ./src/contrafold predict test.seq --params parameters/EternaFoldParams.v1 --partition
```

Output (log partition coefficient)
```
Training mode: 
Use constraints: 0
Use evidence: 0
Log partition coefficient for "test.seq": 13.7489
```
#### Base-pairing probability prediction

```
./src/contrafold predict test.seq --params parameters/EternaFoldParams.v1 --posteriors 0.00001 bps.txt
```

Base-pairing probabilities are output to `bps.txt`:
```
1 C 9:3.55095e-05 28:0.000274751 31:0.0050855 33:0.000420935 46:0.00100593 52:0.000674974 54:0.815493
2 G 7:0.000290278 10:0.000150796 16:6.48946e-05 22:0.000711706 24:6.9622e-05 26:0.000379153 27:0.000149917 30:0.005751
06 34:0.00134091 35:0.00017805 45:0.000245854 50:0.000512436 53:0.913047
3 C 9:0.000150353 15:6.90445e-05 21:0.000968743 28:0.00245682 31:0.00417261 33:0.00229046 46:0.000703465 52:0.91348 54
:0.000561778
4 U 20:0.00104566 25:0.00225947 28:0.000790812 29:0.000620939 31:0.0282994 32:0.00417285 41:2.64421e-05 46:0.000350788
 47:0.00012951 48:0.000177715 49:0.000466242 51:0.825096 52:8.55103e-05 54:0.000171393
5 G 12:0.000356127 19:0.00100131 24:0.00327509 26:0.00645875 27:0.000742048 30:0.189649 45:0.00133716 50:0.74545 53:0.
000533715
6 U 11:0.00030157 17:0.000271142 23:0.00321799 25:0.0050542 28:0.0026176 29:0.230025 31:0.000505041 32:0.0156585 43:0.
000190485 44:0.00204912 46:0.00040755 47:0.0307307 48:0.0032191 49:0.561613 52:0.00014478
...
```

#### Sample structures

Stochastically samples structures from the underlying distribution. 
```
./src/contrafold sample test.seq --params parameters/EternaFoldParams.v1 --nsamples 10
```

Output
```
Training mode: 
Use constraints: 0
Use evidence: 0
(((((((..((((......)))).......))((((.....))))....)))))
..(.(((.((....(((....))))).)))).((((.....)))).........
................................((((.....)))).........
........(((((......)))))........((((.....)))).........
.(((((.((((((......)))))).......((((.....))))...))))).
.(((((..((((........))))........((((.....))))...))))).
.((((((.(((((......)))))........((((.....)))))..))))).
.(((((.((((((......)))))).......((((.....))))...))))).
....(((.(((((......)))))...)))..((((.....)))).........
....(((((((((......))))))..)))..((((.....)))).........
```

`sample` can be used in conjunction with SHAPE data to sample SHAPE-reweighted distribution:

```
./src/contrafold sample test_SHAPE.bpseq --params parameters/EternaFoldParams_PLUS_POTENTIALS.v1 --nsamples 10 --evidence --numdatasources 1 --kappa 0.1
```

Output:
```
Training mode: 
Use constraints: 0
Use evidence: 1
.(((((....)))))..........
((((((....)))))).........
((((((....)))))).........
((((((....)))))).........
((((((....)))))).........
.(((((....)))))..........
.(.(((....))).)..........
.(((((....)))))..........
...(((....)))............
.(((((....)))))..........
```



Please see the documentation of [CONTRAfold](http://contra.stanford.edu/contrafold/manual_v2_02.pdf) for further information on parameters and usage. See below for documented discrepancies (besides parameters) from CONTRAfold codebases.


### Training

Training data is in `input_data` (unzip first).

Text files containing the lists used for training, test, and holdout models for the EternaFold models reported in Wayment-Steele et al. (2020) are found in the `datalists` repo.

From CONTRAfold-SE:

"Learn parameters based on a set of sequences, in which sequences with associated probing data have data from 2 sources, and with a relative weight (specified by `hyperparam_data`) of 0.1.

Assumes that folder "trainset" has a set of sequences of type ".bpseq" in evidence format for the ones with data.

`contrafold train --regularize 1 --numdatasources 2 --maxiter 1000 --hyperparam_data 0.1 --initweights contrafold.params.complementary_data2 trainset/*.bpseq`

If there are a large number of input files used (> 1000 files; e.g. for training on RMDB data), provide a text file containing the list of example files instead with the `--examplefile` option.

`contrafold train --regularize 1 --numdatasources 1 --maxiter 500 --examplefile examples.txt`"

#### Training options for riboswitch data

`contrafold train --examplefile ../production_struct_riboswitches.txt --regularize 32  --kd_hyperparam_data 30 --ligand --ligand_bonus 90 --lig_hyperparam_data 30`

`kd_hyperparam_data`: weight placed on no-ligand kd values.
`lig_hyperparam_data`: weight placed on ligand kd values.
`ligand_bonus`: ligand bonus used.


### Input file formats

#### Chemical mapping

From CONTRAfold-SE:

"To support structure probing data we adapt the BPSEQ format in two ways to support sequences with only probing data (BPP2SEQ), and sequences with both probing data and known structure (base-pairings) (BPP2TSEQ).

The original BPSEQ format consists of a row for each base in the RNA molecule, describing the index (1-based), the actual base present, and the index of the pairing partner (0 if unpaired). 

In the BPP2SEQ format, there is an evidence string `e<N>` following the base, where `<N>` is an integer denoting how many probing data sources there are, followed by `<N>` (positive) values of the unpairedness potential (derived from probing data) for that base. For instance, `e2` denotes 2 probing sources, and should be followed by a two positive real numbers; all entries should have the same evidence string. An example is shown below:

```
1 G e2 7.070000e+00 -1.000000e+00
2 A e2 7.570000e+00 3.333333e-02
3 A e2 6.500000e+00 -1.000000e+00
4 A e2 5.310000e+00 4.444444e-02
```

In the BPP2TSEQ format, the evidence string is `t<N>` instead, and following the unpairedness potential values is the index of the pairing partner."

Values that are less than 1e-5 are ignored (treated as not present) for numerical stability.


#### Riboswitch format

To account for three different constrained structures, riboswitch molecules are input using a BPSEQ file that is modified to have three columns correpsonding to three different constrained structures:

```
k1.0 2.0 99
1 G -1 -1 -1
2 A -1 -1 -1
3 U -1 19 19
4 C -1 18 18
```

#### ❗️ Discrepancies from CONTRAfold-SE code

This code has been modified in two ways that means its output, even using the CONTRAfold parameters, will differ from the CONTRAfold codebase here and the CONTRAfold-SE codebase here.

1. A bug was fixed in the multiloop traceback `InferenceEngine.ipp` which was first identified by He Zhang (Oregon State). 

2. The minimum allowable hairpin size was increased from `0` to `3` to prevent structure predictions with `(())` hairpins. To revert back to the original CONTRAfold behavior, set `C_MIN_HP_LENGTH=0` in `Config.hpp` before compiling.

Predictions for Hammerhead Ribozyme sequence, using default CONTRAfold parameters: `CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG`

contrafold predict hhr.bpseq --partition 

| Version | hhr.bpseq Log Partition Coefficient |
| --- | ----------- |
|CONTRAfold v2.02| 6.87394|
|CONTRAfold-SE| 6.87394|
|EternaFold code, no ML fix and C_MIN_HP_LENGTH=0| 6.87394|
|EternaFold code, C_MIN_HP_LENGTH=0| 6.83585|
|EternaFold code | 6.77285 |
