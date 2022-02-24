# Script to evaluate fold change prediction.
# Usage: first arg is name of parameter file to evaluate.

../src/contrafold predict-foldchange /home/users/hannahw1/secstruct_software/contrafold-se/input_data/RiboswitchData/test/* --bpseq tmp --params $1
cat tmp/*.bpseq > fold_change_data.txt
python get_rmse.py fold_change_Data.txt > fold_change_stats.txt
