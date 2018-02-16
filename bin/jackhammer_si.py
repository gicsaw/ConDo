#!/usr/bin/env python

import numpy as np
import sys

USAGE="""
mi.py target
"""

seqcode1="ARNDCQEGHILKMFPSTWYV-"

def get_seq(target_f):
    fp = open(target_f+".fasta", "r")
    line = fp.readlines()
    fp.close()
    
    return line[1].strip()


def si(target) :

    seq=get_seq(target)
    Nseq=len(seq)

    msafile = target+".msa"
    fp_msa=open(msafile)
    msa_lines=fp_msa.readlines()
    fp_msa.close()

    k=0 
    count=0
    msa=[]
    si_a=[]
    for line in msa_lines:
        if (line[0]==">"):
            k+=1
            title=line[:-1]
            continue
        seq2=line[:-1]
        Mcount=0
        Gcount=0
        for i in range(0,Nseq):
            if seq2[i]=="-" :
                Gcount+=1
            elif seq[i]==seq2[i]:
                Mcount+=1
        Acount=Nseq-Gcount
        SI=Mcount/float(Acount)
#        print title, SI, Gcount, Mcount, Acount
        msa+=[[title,seq2,SI,Gcount,Mcount,Acount]]
        si_a+=[SI]
    Nmsa=k

    msa_cut=10000
    if Nmsa>=msa_cut : 
        si_array=np.array(si_a)
        sorder=np.argsort(-si_array)
        ind=sorder[msa_cut]
        si_cut=si_a[ind]
#        print si_cut
    else :
        si_cut = 0.0

    count=0
    output_file=target+".aln"
    output_file2=target+".sim"

    fp_out=open(output_file,"w")
    fp_out2=open(output_file2,"w")

    for i in range(0,Nmsa):
        gg=msa[i]
        if gg[2]<=si_cut:
            continue
        count+=1
        title=gg[0]
        seq2=gg[1]
        SI=gg[2]
        Gcount=gg[3]
        Mcount=gg[4]
        Acount=gg[5]
#            print title, SI, Gcount, Mcount, Acount
        print >> fp_out2, title, SI, Gcount, Mcount, Acount
        print >> fp_out, seq2

    fp_out.close()
    fp_out2.close()

    print count,k, si_cut
 

def main():

    if len(sys.argv)<2:
        print USAGE
        sys.exit()

    target= sys.argv[1]
    si(target)



if __name__ == '__main__':
    main()


