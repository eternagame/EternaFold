#!/bin/bash

#Usage: ./batch_score_prediction.sh /path/to/ref_files/*.bpseq /path/to/predicted/files

echo "$1"
echo $2

# /home/users/hannahw1/secstruct_software/contrafold-se/src/contrafold predict $@ --bpseq $pred_file_dir

# for arg in $@; do
# 	if [ "$arg" != "$pred_file_dir" ]; then
# 		/home/users/hannahw1/secstruct_software/contrafold-se/src/score_prediction rna $arg ${pred_file_dir}/${arg##*/} >> ${pred_file_dir}/output
# 	fi
# done