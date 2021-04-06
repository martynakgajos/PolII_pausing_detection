#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import pandas as pd

dfs=[]
for fn in snakemake.input:
    if os.stat(fn).st_size == 0:
        continue
    df=pd.read_csv(fn,delimiter='\t',header=None)
    dfs.append(df)
df=pd.concat(dfs,sort=[0,1],ignore_index=True)
df.to_csv(snakemake.output[1],sep='\t',index=False,header=False)
df[df[4]>=snakemake.params[0]].to_csv(snakemake.output[0],sep='\t',index=False,header=False)
