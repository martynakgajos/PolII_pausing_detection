#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import pandas as pd
import pybedtools

if os.stat(snakemake.input[0]).st_size == 0:
    with open(snakemake.output[0],'w') as f:
        exit()

ps = pybedtools.BedTool(snakemake.input[0])
if snakemake.params[0]=='':
    df=ps
else:   
    mf = pybedtools.BedTool(snakemake.params[0])
    df = ps.subtract(mf, s=True)
df.moveto(snakemake.output[0])

