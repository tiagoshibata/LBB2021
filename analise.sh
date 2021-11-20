#!/bin/bash
REF=Data/ref.fa.gz
ANC=Data/anc.fa.gz
ANGSD=~/src/angsd
NGSTOOLS=~/src/ngsTools

$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL.qc -r 11\
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 \
        -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500

ls Results
echo Quality scores
less -S Results/ALL.qc.qs
echo Per sample depth
less -S Results/ALL.qc.depthSample 
wc -l Results/ALL.qc.depthSample # 30 Results/ALL.qc.depthSample
echo Global depth
less -S Results/ALL.qc.depthGlobal

Rscript $NGSTOOLS/Scripts/plotQC.R Results/ALL.qc
less -S Results/ALL.qc.info
xdg-open Results/ALL.qc.pdf
