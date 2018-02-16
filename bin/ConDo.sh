#!/bin/bash 

#ConDo target.fasta ncpu

ConDodir=$HOME/Programs/ConDo
ConDobin=$ConDodir/bin
weight_file=$ConDodir/data/weight.h5
target=${1%.*}

if [ $# -eq 1 ]
then
    nprocessor=1
else
    nprocessor=$2
fi

$ConDobin/run_jackhmmer.sh $target $nprocessor
$ConDobin/run_ccmpred.sh $target $nprocessor
$ConDobin/gen_features.sh $target $nprocessor
$ConDobin/feature $target $nprocessor
$ConDobin/gather_input.py $target
#keras with TF or theano (CPU)
$ConDobin/prediction.py data_feature.dat.npz y_pred.dat.npz $weight_file
# GPU for theano
#THEANO_FLAGS=device=cuda,floatX=float32, python $ConDobin/prediction.py data_feature.dat.npz y_pred.dat.npz $weight_file

$ConDobin/gen_results.py $target

