## EternaFold - Improving RNA structure prediction through multitask learning

EternaFold performs multitask learning to improve RNA structure prediction. Its training tasks include 1) predicting single structures, 2) maximizing the likelihood of structure probing data, and 3) predicting experimentally-measured affinities of RNA molecules to proteins and small molecules.

Its training data comes from diverse high-throughput experimental crowdsourced data from the [Eterna](www.eternagame.org) project.

EternaFold is possible thanks to [CONTRAfold-SE](https://github.com/csfoo/contrafold-se) (C.-S. Foo, C. Pop).

### Installation

Clone the repository and run `make` in `src` to compile.
Multithreaded version: run `make multi` in `src`.

See instructions in [README_LinearFold-E_patch.md](README_LinearFold-E_patch.md) for using EternaFold parameters with LinearFold and LinearPartition algorithms.

### Prediction

#### Single-structure prediction
Predict the MEA structure of sequence "test", using the EternaFold parameters:

`contrafold predict test.bpseq --params parameters/EternaFoldParams.v1`

Please see the documentation of [CONTRAfold](http://contra.stanford.edu/contrafold/manual_v2_02.pdf) for further information on parameters and usage.

#### Fold change prediction
Predict log K_MS2 values for riboswitch molecules to MS2 in the presence and absence of small molecule aptamers.

`contrafold predict-foldchange test_ms2.bpseq`

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
