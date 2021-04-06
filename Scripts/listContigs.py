#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

contigs=set()
for fn in snakemake.input:
    df=pd.read_csv(fn, delimiter='\t', usecols=[0],header=None)[0].values
    contigs=contigs.union(set(list(df)))
contigs=list(contigs)
contigs.sort()
with open(snakemake.output[0],'w') as fn:
    fn.writelines(["%s\n" % item  for item in contigs])
    