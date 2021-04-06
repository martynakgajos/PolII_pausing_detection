#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

def getarg():
    parser = argparse.ArgumentParser(description = "To be added")
    parser.add_argument("-f", "--File", help = "GFF file with pausing sites. Example: /path/*.gff", required = True)
    parser.add_argument("-o", "--Output", help = "Path to output plot. Example: /path/*.pdf", required = False, default = 'PausingDistribution.pdf')
    parser.add_argument("-a", "--Annotation", help = "GFF file with gene annotation. Example: /path/*.gff", required = True)
    parser.add_argument("-i", "--Isoforms", help = "Path to RSEM isoform file results (directory). Example: /path/", required = False, default = False)
    arg = parser.parse_args()
    return arg.File, arg.Output, arg.Annotation , arg.Isoforms

def returnExons(afile, activeTs=[]):
    d={}
    with open(afile,'r') as af:
        for linia in af:
            if linia[0]!='#':
                ln=linia.strip().split()
                if ln[2]=='exon':
                    if (ln[8].split('ID=exon:')[1].split(':')[0] in activeTs) or (len(activeTs)==0):
                        gene_id=ln[8].split('gene_id=')[1].split(';')[0]
                        try:
                            d[gene_id].append((int(ln[3]),int(ln[4])))
                        except:
                            d[gene_id]=[(int(ln[3]),int(ln[4]))]       
    return d
    
def isExon(row,dic):
    reg="intron"
    if row["gene"] in set(dic.keys()):
        for tup in dic[row["gene"]]:
            if tup[0]<=row["start"]<=tup[1]:
                reg="exon"
                break
    return reg

#inputfile, outputfile, annotation, isoforms=getarg()
inputfile='/project/owlmayerTemporary/Martyna/NET_pro/Results/GRCh38p12_HeLaS3.none.Rep1.gff'
annotation='/project/PausingDynamics/GeneralResources/Genecode29/gencode.v29.annotation.sorted.gff3'
outputfile='/project/owlmayerTemporary/Martyna/NET_pro/Results/GRCh38p12_HeLaS3.none.Rep1.pdf'
isoforms=False
df=pd.read_csv(inputfile,names=['start','score','attribute'],usecols=[3,5,8],sep='\t')
df['gene']=df.apply(lambda x: x['attribute'].split('_')[0],1)
df['region']=df.apply(lambda x: x['attribute'].split('_')[1],1)
df=df[~(df["region"].isin(["OP","NC","RNA"]))]
df=df.replace({'GB':'gene\nbody','PP':'promoter\nproximal',
               'CA':'convergent\nantisense','DA':'divergent\nantisense','AS':'antisense',
               'multiple':'undetermined','TW':'termination\nwindow'})
order=['gene\nbody', 'intergenic', 'convergent\nantisense','promoter\nproximal',
       'divergent\nantisense', 'antisense', 'undetermined', 'termination\nwindow']
colors={}
for i in range(len(order)):
    colors[order[i]]=sns.husl_palette(len(order), s=0.7 )[i]

#intron/exon
if isoforms:
    pct=10
    l=[]
    for fn in os.listdir(isoforms+'data/'):
        if fn.endswith('.isoforms.results'):
            sample=fn.split('.')[0]
            data=pd.read_csv(isoforms+'data/'+sample+'.isoforms.results',sep='\t',usecols=[0,1,7])
            data=data.set_index(['transcript_id','gene_id'])
            data=data.rename({"IsoPct":"P_"+sample},axis='columns')
            l.append(data)
    data=l[0]
    for i in range(1,len(l)):
        data=data.join(l[i],how="outer")
    data=data.fillna(0)
    data=data[data.mean(axis=1)>pct]
    activeTranscripts=set([i[0] for i in data.index.values])
    del(data)
    exons=returnExons(annotation,activeTranscripts)
else:
    exons=returnExons(annotation)
    
gb=df[df["region"]=="gene\nbody"]
gb["part"]=gb.apply(lambda x: isExon(x,exons),1)

cap='number of pausing sites of the type'
fig, ax = plt.subplots(figsize=(5,6)) 
cur_axes = plt.gca()
ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 
percentage=(df['region'].value_counts()/df['region'].value_counts().sum()*100).to_frame().sort_values(by='region',ascending=False)
sns.countplot(y="region",data=df,order = percentage.index,palette = [colors[i] for i in percentage.index])
ax.set_xlabel(cap)
for i in range(0,percentage.shape[0]):
    p=ax.patches[i]
    h=percentage.iloc[i]['region']
    ax.annotate('{:.1f}%'.format(h), (p.get_width()+10, p.get_y()-0.5*p.get_height()-0.15+1))
plt.subplots_adjust(left=0.27, right=0.9, top=0.9, bottom=0.05)
gb_index=list(percentage.index).index('gene\nbody')
rect = plt.Rectangle((ax.patches[gb_index].get_x(),ax.patches[gb_index].get_y()), gb.groupby("part").count()["region"]["exon"], ax.patches[gb_index].get_height(), color='k', alpha=0.3)
ax.add_patch(rect)
ax.annotate('E', (gb.groupby("part").count()["region"]["exon"]*0.5, ax.patches[gb_index].get_y()-0.5*ax.patches[gb_index].get_height()-0.15+1),color='white')
ax.annotate('I', (gb.groupby("part").count()["region"]["exon"]+gb.groupby("part").count()["region"]["intron"]*0.5, ax.patches[0].get_y()-0.5*ax.patches[gb_index].get_height()-0.15+1),color='white')
plt.savefig(outputfile)
