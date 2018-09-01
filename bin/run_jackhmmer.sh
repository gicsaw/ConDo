#!/bin/bash 

# set directory and file path
hhpath=/usr/local/hh-suite
condodir=~/ConDo
condobin=$condodir/bin
database=/home1/shade/Database/uniref/20180131/uniref90.fasta
jackhmmerbin=~/Programs/hmmer/bin

export HHLIB=$hhpath

target=${1%.*}

if [ $# -eq 1 ]
then
    nprocessor=1
else
    nprocessor=$2
fi

$jackhmmerbin/jackhmmer -N 4 --cpu $nprocessor -o $target.dat -A $target.align $target.fasta $database 
$jackhmmerbin/esl-reformat -o $target.a2m a2m $target.align
$hhpath/scripts/reformat.pl -r $target.a2m $target.hmm.fas
$condobin/jackhammer_aln.py $target
$condobin/jackhammer_si.py $target

rm -f $target.align $target.a2m $target.hmm.fas 

