#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd
import numpy as np

if os.stat(snakemake.input[0]).st_size == 0:
    with open(snakemake.output[0],'w') as f:
        exit()
        
strands={'pos':'+',
         'neg':'-'}
try:
    strand=strands[snakemake.params[4]]
except:
    strand=snakemake.params[4]

df=pd.read_csv(snakemake.input[0],sep='\t',header=None)
    
flanks=(snakemake.params[0]-1)//2
pos_0=df[1].min()  
signal_len=df[2].max()-pos_0

signal=np.zeros(signal_len,dtype=np.int)
for row in df.iterrows():
    start=row[1][1]-pos_0
    stop=row[1][2]-pos_0
    value=row[1][3]
    for i in range(start,stop):
        signal[i]=value

max_left=np.sign(signal-np.pad(signal,1,'constant')[2:])
max_right=np.sign(signal-np.pad(signal,1,'constant')[:-2])

threshold=signal-snakemake.params[1]

window=np.ones(flanks*2+1)
surroundings=np.convolve(signal,window,'same')
nonzero=np.convolve((signal>0)*1,window,'same')

potential_max=max_left*max_right*threshold
potential_max_positions=np.where(potential_max>0)[0]
with open(snakemake.output[0],'w') as f:  
    for i in potential_max_positions:
            if max_left[i]>0 and max_right[i]>0 and threshold[i]>0:
                f.write(snakemake.params[3]+'\t'+str(i+pos_0)+'\t'+str(i+pos_0+1)+'\t'+str(signal[i])+'\t'+str(int(surroundings[i]))+'\t'+strand+'\t'+str(int(nonzero[i]))+'\n')