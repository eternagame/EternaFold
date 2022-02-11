# Datasets in fasta form

EternaBench_datasets/EB_ChemMapping_unfiltered.fasta # Original set of xx,xxx constructs from Cloud Labs
EternaBench_datasets/EB_ChemMapping_Filtered.fasta # EB_ChemMapping_unfiltered, filtered for <80% sequence similarity within each experiment
EternaBench_datasets/EB_Riboswitch_unfiltered.fasta # Original set of xx,xxx riboswitch constructs from Andreasson, 2019 and Wu, 2019
EternaBench_datasets/EB_Riboswitch_Filtered.fasta # EB_Riboswitch_unfiltered, filtered for <80% sequence similarity within each experiment

train_datasets/ChemMappingData_train.fasta # From Cloud labs data that was originally used in CONTRAfold-SE
train_datasets/RiboswitchData_train.fasta # training split of RiboLogic riboswitches
train_datasets/S-Processed-TRA.fasta # STRAND S-Processed training set
train_datasets/All_EternaFold_training_data.fasta # All of the above data, plus holdout datasets, concatenated


holdout_datasets/ChemMappingData_holdout.fasta # From Cloud labs data that was originally used in CONTRAfold-SE
holdout_datasets/RiboswitchData_holdout.fasta # holdout split of RiboLogic riboswitches
holdout_datasets/S-Processed-VAL.fasta # STRAND S-Processed validation set

test_datasets/ExternalData.fasta # Sequences from independent sources (no duplicates)
test_datasets/ExternalData_window300_uniq.fasta # ExternalData, windowed at size 300 (no duplicates)
test_datasets/ExternalData_window600_uniq.fasta # ExternalData, windowed at size 600 (no duplicates)
test_datasets/ExternalData_window900_uniq.fasta # ExternalData, windowed at size 900 (no duplicates)
test_datasets/ExternalData_window1200_uniq.fasta # ExternalData, windowed at size 1200 (no duplicates)
test_datasets/S-Processed-TES.fasta # STRAND S-Processed test set, demonstrated ~50% sequence redundancy with S-Processed train and test
test_datasets/archiveII_original.fasta # original ArchiveII (not used for benchmarking)
test_datasets/archiveII_filt_80_NR.fasta # ArchiveII, filtered to maintain <80% sequence similarity to all EternaFold training data
