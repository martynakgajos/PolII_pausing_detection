#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np
import scipy.stats as st

def calculatePercentile(peak,reads,positions,bootstrap_N):
    bootstrap_table=np.zeros((bootstrap_N,reads),dtype=np.int)
    for i in range(bootstrap_N):
        bootstrap_table[i]=np.random.randint(positions,size=reads)
    maxi=np.apply_along_axis(getMax,1,bootstrap_table)
    return st.percentileofscore(maxi,peak,'mean')

def getMax(x):
    return max(np.bincount(x))

strands={'pos':'+',
         'neg':'-'}

if os.stat(snakemake.input[0]).st_size == 0:
    with open(snakemake.output[0],'w') as f:
        exit() 
        
df=pd.read_csv(snakemake.input[0],delimiter='\t',header=None)
df['strand']=df.apply(lambda x:
    strands[snakemake.params[2]] if snakemake.params[2] in strands.keys() else '.',1)

if snakemake.params[3]=='poisson':
    df['score']=df.apply(lambda x: st.poisson(x[4]*1./x[6]).cdf(x[3])*100,1)
else:
    df['score']=df.apply(lambda x: calculatePercentile(x[3],x[4],x[6],snakemake.params[0]),1)
df=df[[0,1,2,3,'score','strand']]
df=df.sort_values(by=[0,1])
df.to_csv(snakemake.output[0],header=None,index=None,sep='\t')