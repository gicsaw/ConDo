#!/usr/bin/env python

import os, sys
import numpy as np

USAGE="""
gen_results.py target
"""
NC=40

def gen_results(target,conf_cut):
    pred_npz_data="y_pred.dat.npz" 
    pred_boundary=np.load(pred_npz_data)
    boundary_pred= pred_boundary['y_pred']

    Nseq=len(boundary_pred)
    predfile=target+".ConDo"
    
    fp_dom=open(predfile,"w")

    scores=[]
    for i in range(0,Nseq):
        score=0
        line_out=""
        for j in range(0,4):
            score+=boundary_pred[i][j]
            line_out+="%f " %(boundary_pred[i][j])
        scores+=[score]
        print >> fp_dom, i+1, score, line_out
    npscores=np.array(scores)
    arg_scores=np.argsort(-npscores)
    fp_dom.close()

    boundary=[]
    boundary2=[]
    bd2=[]
    for i in range(0,Nseq):
        k=arg_scores[i]
        score=scores[k]
        if score<conf_cut:
            break
        if (k<1 or k>=Nseq-1) :
            continue
        check=0
        for bb in boundary:
            if abs(bb-(k+1))<40:
                check+=1
        if check==0:
            boundary+=[k+1]
            if scores[k+1]<scores[k-1]:
                boundary2+=[[k,k+1]]
            else:
                boundary2+=[[k+1,k+2]]
            bd2+=[[k+1,score]]

    
    Ncount=0
    for k in range(0,len(bd2)):
        i=bd2[k][0]
        score=bd2[k][1]
        if (score>=conf_cut) and (i>NC and i< Nseq-NC-1):
            Ncount+=1
        print "boundary:", boundary2[k]
    
    if Ncount==0:
        print target+": single-domain"
        return

    print target+ ": multi-domain"

    return

def main ():

    if len(sys.argv)<2:
        print USAGE
        sys.exit()

    target= sys.argv[1]

    conf_cut=1.4
    gen_results(target,conf_cut)

if __name__ == '__main__':

    main ()

