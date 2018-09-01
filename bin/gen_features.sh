#!/bin/bash

target=$1
if [ $# -eq 1 ]
then
    ncpu=1
else
    ncpu=$2
fi

blastbin=/usr/bin
dbname=/home1/shade/Database/uniref/20180131/uniref90.fasta

psipred=~/Programs/psipred
psipredbin=$psipred/bin
psipreddata=$psipred/data

condodir=~/ConDo
condobin=$condodir/bin
condodata=$condodir/data

sann=~/Programs/sann
#dir of database.rsa, database.vec
export NNDB_HOME=$sann/nndb

log="gen_features.log"

#blast
$blastbin/blastpgp -b 0 -v 5000 -j 3 -h 0.001 -a $ncpu -d $dbname -i $target.fasta -C $target.chk >& $target.blast

if [ ! -e $target.chk ]; then
    echo ">>$target.chk is not exist has problem in a3m, hhr, ..." 
    echo "check it"
    exit
fi

#$psipred/runpsipred_single $target

tmproot=psitmp$$$hostid
cp -f $target.fasta $tmproot.fasta
$psipredbin/seq2mtx $target.fasta > $tmproot.mtx
$psipredbin/psipred $tmproot.mtx $psipreddata/weights.dat $psipreddata/weights.dat2 $psipreddata/weights.dat3 > $target.ss
$psipredbin/psipass2 $psipreddata/weights_p2.dat 1 1.0 1.0 $target.ss2 $target.ss > $target.horiz

rm -f $tmproot* 

# sann2 (using new database)
$condobin/mkchk2.py $target $condodata/qij  | tee -a $log
$sann/bin/sann -i $target -np $ncpu
ln -s $target.a3 $target.sa2


