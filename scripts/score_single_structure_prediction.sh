# Evaluate single-state structure prediction.
# Usage: first argument is name of parameter file to evaluate.

../src/contrafold predict ../input_data/StructureData/test/*.bpseq  --params $1 --bpseq curr_outdir;
cd curr_outdir

for j in *.bpseq; do
../src/score_prediction rna ../input_data/StructureData/test/${j} $j; done > scoring
cut -f 10 -d';' scoring | cut -f2 -d'=' > sensitivity
cut -f 11 -d';' scoring | cut -f2 -d'=' > ppv
paste sensitivity ppv | awk '{print sqrt($1 * $2)}'  > fscore
awk '{sum += $1; n++} END {if (n>0) print sum/n;}' fscore > stats

cd ..
