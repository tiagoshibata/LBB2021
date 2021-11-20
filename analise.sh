#!/bin/bash
set -e
REF=Data/ref.fa.gz
ANC=Data/anc.fa.gz
ANGSD=~/src/angsd
NGSTOOLS=~/src/ngsTools
FASTME=~/src/FastME/src/fastme

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

$ANGSD/angsd -P 4 -b ALL.bamlist -ref $REF -out Results/ALL -r 11 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 15 -setMinDepth 60 -setMaxDepth 400 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
        -SNP_pval 1e-3 \
        -doGeno 8 -doPost 1
NSITES=`zcat Results/ALL.mafs.gz | tail -n+2 | wc -l`
echo $NSITES

Rscript -e 'cat(paste(rep(c("LWK","TSI","PEL"),each=10), rep(1:10, 3), sep="_"), sep="\n", file="Data/pops.label")'
cat Data/pops.label

$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 30 -n_sites $NSITES -labels Data/pops.label -o Results/ALL.dist -n_threads 4
less -S Results/ALL.dist

$FASTME -D 1 -i Results/ALL.dist -o Results/ALL.tree -m b -n b
cat Results/ALL.tree

Rscript $NGSTOOLS/Scripts/plotTree.R Results/ALL.tree
xdg-open Results/ALL.tree.pdf

$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/ALL.geno.gz -probs -n_ind 30 -n_sites $NSITES -labels Data/pops.label -o Results/ALL.boot.dist -n_threads 4 -n_boot_rep 20 -boot_block_size 20

$FASTME -D 21 -i Results/ALL.boot.dist -o Results/ALL.boot.tree -m b -n b
Rscript $NGSTOOLS/Scripts/plotTreeBoots.R Results/ALL.boot.tree
xdg-open Results/ALL.boot.tree.pdf

NSAMPLES=30
tail -n +3 Results/ALL.dist | head -n $NSAMPLES | Rscript --vanilla --slave $NGSTOOLS/Scripts/getMDS.R --no_header --data_symm -n 4 -m "mds" -o Results/ALL.mds
less -S Results/ALL.mds

Rscript $NGSTOOLS/Scripts/plotMDS.R -i Results/ALL.mds -c 1-2 -a Results/ALL.clst -o Results/ALL.mds.pdf
xdg-open Results/ALL.mds.pdf

for POP in LWK TSI PEL
do
	echo $POP
	$ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $ANC -out Results/$POP -r 11 \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
		-GL 1 -doSaf 1 &> /dev/null
done

$ANGSD/misc/realSFS print Results/PEL.saf.idx | less -S

for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/misc/realSFS Results/$POP.saf.idx -P 4 2> /dev/null > Results/$POP.sfs
done

Rscript $NGSTOOLS/Scripts/plotSFS.R Results/LWK.sfs-Results/TSI.sfs-Results/PEL.sfs LWK-TSI-PEL 0 Results/ALL.sfs.pdf
xdg-open Results/ALL.sfs.pdf

for POP in LWK TSI PEL
do
        echo $POP
        $ANGSD/angsd -P 4 -b $POP.bamlist -ref $REF -anc $REF -out Results/${POP}.ref -r 11 \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 10 -setMinDepth 20 -setMaxDepth 200 -doCounts 1 \
                -GL 1 -doSaf 1 &> /dev/null
	$ANGSD/misc/realSFS Results/$POP.ref.saf.idx -P 4 2> /dev/null > Results/$POP.ref.sfs
	
done
Rscript $NGSTOOLS/Scripts/plotSFS.R Results/LWK.ref.sfs-Results/TSI.ref.sfs-Results/PEL.ref.sfs LWK-TSI-PEL 1 Results/ALL.ref.sfs.pdf
xdg-open Results/ALL.ref.sfs.pdf

$ANGSD/misc/realSFS Results/PEL.saf.idx -bootstrap 10  2> /dev/null > Results/PEL.boots.sfs
cat Results/PEL.boots.sfs

for POP in LWK TSI
do
	$ANGSD/misc/realSFS -P 4 Results/$POP.saf.idx Results/PEL.saf.idx > Results/$POP.PEL.sfs
done
