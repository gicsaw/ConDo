#!/usr/bin/env python

import numpy as np
import sys

USAGE="""
gen_input.py target
"""

def get_seq(seq_file) :
    fp=open(seq_file,"r")
    lines=fp.readlines()
    fp.close()
    seq=""
    for line in lines:
        if line[0]==">":
            continue
        else :
            seq+=line.strip()

    return seq

def main():
    if len(sys.argv)<2:
        print USAGE
        sys.exit()

    data_feature=[]
    k=0

    target=sys.argv[1]
    seq_file=target+".fasta"
    seq=get_seq(seq_file)
    Nseq=len(seq)
 
    feature_file=target+"_feature.txt"
    fp=open(feature_file)
    feature_lines=fp.readlines()
    fp.close()

    for feature_line in feature_lines:
        feature=feature_line.split()
        data_feature+=[np.array(feature[1:], dtype=np.float32)]

    data_feature=np.array(data_feature, dtype=np.float32)
    outfile_feature = "data_feature.dat"
    np.savez(outfile_feature, feature=data_feature)


if __name__ == '__main__':
    main()


