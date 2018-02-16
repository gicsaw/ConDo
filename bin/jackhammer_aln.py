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


def main():

    if len(sys.argv)<2:
        print USAGE
        sys.exit()

    target= sys.argv[1]
    a2mfile=target+".hmm.fas"

    seq=get_seq(target)
    Nseq=len(seq)

    fp_a2m=open(a2mfile)
    a2m_lines=fp_a2m.readlines()
    fp_a2m.close()

    k=0 
    seq=""
    aln=[]
    for line in a2m_lines:
        if (line[0]==">"):
            k+=1
            seq=""
            arr=line.split()
            arr2=arr[0].split("/")
            title=arr2[0]

#            inifin2=arr2[1].split("-")
#            if len(inifin2)>1 :
#                ini2=int(inifin2[0])
#                fin2=int(inifin2[1])
#            else :
#                ini2=0
#                fin2=Nseq-1

            continue

        seq+=line.strip()
        if(Nseq==len(seq)):
            ini1=0
            fin1=Nseq-1
            for i in range(0,Nseq):
                if seq[i]!="-" :
                    ini1=i
                    break
            for i in range(Nseq-1,-1,-1):
                if seq[i]!="-" :
                    fin1=i
                    break
            aln+=[[title,seq,ini1,fin1]]

    Naln=len(aln)
#    msafile=target+".msa"
#    fp_msa=open(msafile,"w")
#    for i in range(0,Naln):
#        title=aln[i][0]
#        seq=aln[i][1]
#        ini1=aln[i][2]
#        fin1=aln[i][3]
#        print >>fp_msa, "%s/%d-%d" %(title, ini1+1, fin1+1)
#        print >>fp_msa, seq

#    fp_msa.close()

    aln2=[]
    k=0
    title_old=""
    check=0
    for i in range(0,Naln) :
        title=aln[i][0]
        seq=aln[i][1]
        ini1=aln[i][2]
        fin1=aln[i][3]

        if title!=title_old:
            if check!=0:
                aln2+=[con]
            check=1
            k+=1
            title_old=title
            con=[aln[i]]
        else:
            con+=[aln[i]]
    aln2+=[con]

    alnfile=target+".aln"
    fp_aln=open(alnfile,"w")

    msafile=target+".msa"
    fp_msa=open(msafile,"w")

    Naln2=len(aln2)
    for k in range(0,Naln2) :
        dom=aln2[k]
        Ndom=len(dom)
        title= dom[0][0]
        if Ndom==1 :
            seq= dom[0][1]
            ini1=dom[0][2]
            fin1=dom[0][3]

            print >>fp_msa, "%s/%d-%d" %(title, ini1+1, fin1+1)
            print >>fp_msa, seq
            print >>fp_aln, seq
            continue
        check_use=[]
        for m in range(0,Ndom):
            check_use+=[0]
        for m in range(0,Ndom):
            dm=dom[m]
            seq1=dm[1]
            ini1=dm[2]
            fin1=dm[3]
            for n in range(m+1,Ndom):
                dn=dom[n]
                seq2=dn[1]
                ini2=dn[2]
                fin2=dn[3]
                if (ini1<ini2 and fin1<ini2):
                    seq_new=seq1[0:ini2]+seq2[ini2:]
                    if(ini2-fin1>25):
                        print >>fp_msa, "%s/%d-%d,%d-%d" %(title, ini1+1, fin1+1,ini2+1,fin2+1)
                    else :
                        print >>fp_msa, "%s/%d-%d" %(title, ini1+1, fin2+1)
                    print >>fp_msa, seq_new
                    print >>fp_aln, seq_new
                    check_use[m]=1
                    check_use[n]=1

                elif (ini1>ini2 and fin2<ini1):
                    seq_new=seq2[0:ini1]+seq1[ini1:]
                    if(ini1-fin2>25):
                        print >>fp_msa, "%s/%d-%d,%d-%d" %(title, ini2+1, fin2+1,ini1+1,fin1+1)
                    else :
                        print >>fp_msa, "%s/%d-%d" %(title, ini2+1, fin1+1)
                    print >>fp_msa, seq_new
                    print >>fp_aln, seq_new
                    check_use[m]=1
                    check_use[n]=1
            if (check_use[m]==0):
                print >>fp_msa, "%s/%d-%d" %(title, ini1+1, fin1+1)
                print >>fp_msa, seq1
                print >>fp_aln, seq1

    fp_aln.close()
    fp_msa.close()

if __name__ == '__main__':
    main()


