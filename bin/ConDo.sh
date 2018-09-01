#!/bin/bash 

#condo target.fasta ncpu

condodir=$HOME/ConDo
condobin=$condodir/bin
weight_file=$condodir/data/weight.h5
target=${1%.*}

if [ $# -eq 1 ]
then
    nprocessor=1
else
    nprocessor=$2
fi

if [ ! -e $target.fasta ]; then
    echo ">>$target.fasta is not exist " 
    exit
fi

$condobin/run_jackhmmer.sh $target $nprocessor

if [ ! -e $target.aln ]; then
    echo ">>$target.aln is not exist " 
    exit
fi
$condobin/run_ccmpred.sh $target $nprocessor
$condobin/gen_features.sh $target $nprocessor
if [ ! -e $target.ss2 ]; then
    echo ">>$target.ss2 is not exist " 
    echo ">>check PSIPRED"
    exit
fi
if [ ! -e $target.a22 ]; then
    echo ">>$target.a22 is not exist "
    echo ">>check SANN"
    exit
fi
if [ ! -e $target.a3 ]; then
    echo ">>$target.a3 is not exist "
    echo ">>check SANN"
    exit
fi
if [ ! -e $target.msa ]; then
    echo ">>$target.msa is not exist "
    echo ">>check jackhammer"
    exit
fi
if [ ! -e $target.ccmpred ]; then
    echo ">>$target.ccmpred is not exist "
    echo ">>check ccmpred"
    exit
fi

$condobin/feature $target $nprocessor
if [ ! -e $target"_feature.txt" ]; then
    echo ">>$target"_feature.txt" is not exist "
    echo ">>check feature "
    exit
fi

$condobin/gather_input.py $target
#keras with TF or theano (CPU)
$condobin/prediction.py data_feature.dat.npz y_pred.dat.npz $weight_file
# GPU for theano
#THEANO_FLAGS=device=cuda,floatX=float32, python $condobin/prediction.py data_feature.dat.npz y_pred.dat.npz $weight_file

conf_cut=1.4
$condobin/gen_results.py $target $conf_cut

