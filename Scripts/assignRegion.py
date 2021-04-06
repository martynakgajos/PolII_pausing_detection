#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import pandas as pd

def map2Regions(row,df,bed_type=6):
    data=df[df[0]==row[0]]
    data=data[data[1]<=row[1]]
    data=data[data[2]>=row[2]]
    if bed_type==6:
        data=data[data[5]==row[5]]
    if data.shape[0]==0:
        return 'intergenic_intergenic'
    elif data.shape[0]==1:
        return data[3].iloc[0]
    else:
        if set(data[3].values)==1:
            return 'multiple_'+data[3].iloc[0]
        else:
            return 'multiple_multiple'

df=pd.read_csv(snakemake.input[0],sep='\t',header=None)
bed_type=df.shape[1]
df['seqname']=df[0]
df['source']='peakCalling'
df['feature']='peak'
df['start']=df[1]+1
df['end']=df[2]
df['frame']=0
df['score']=(df[4] if bed_type>=6 else '.')
df['strand']=(df[5] if bed_type>=6 else '.')

if snakemake.params[0]=='':
    df['attribute']='.'
else:
    afile=pd.read_csv(snakemake.params[0],sep='\t',header=None)
    df['attribute']=df.apply(lambda x: map2Regions(x,afile,bed_type),1)

bed=df[['seqname',1,'end','attribute','score','strand']]
bed.to_csv(snakemake.output[1],header=None,index=None,sep='\t')
bed['region']=df.apply(lambda x: x['attribute'].split('_')[1],1)
bed=bed[['seqname',1,'end','region','score','strand']]
bed.to_csv(snakemake.output[2],header=None,index=None,sep='\t')

df=df[['seqname','source','feature','start','end','score','strand','frame','attribute']]
df.to_csv(snakemake.output[0],header=None,index=None,sep='\t')
                
