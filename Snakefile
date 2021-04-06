import pandas as pd
import shutil

#GLOBAL VALUES
STRANDS = config["sample"]
SAMPLES = STRANDS.values()

CONTIGS=set()
for s in STRANDS:
    fn="Data/"+STRANDS[s]
    c=list(pd.read_csv(fn, delimiter='\t', usecols=[0],header=None)[0].values)
    CONTIGS=CONTIGS.union(set(c))

for key in ["mask_bed", "region_file"]:
    if key not in config.keys():
        config[key]=""

if "test" not in config.keys():
    config["test"]="bootstrap"

rule all:
    input:
        "Results/Significant_"+config["prefix"]+".bed",
        "Results/Regions_"+config["prefix"]+".bed",
        "Results/"+config["prefix"]+".gff"

rule ExtractChromosome:
    input:
        lambda wildcards: expand("Data/{sample}", sample=STRANDS[wildcards.strand])
    output:
        temp("Temporary/Bed/{strand}/{chr}.bedgraph")
    run:
        df=pd.read_csv(input[0], delimiter='\t',header=None)
        df[df[0]==wildcards.chr].to_csv(output[0],sep='\t',header=False,index=False)

rule GetPausing:
    input:
        "Temporary/Bed/{strand}/{chr}.bedgraph"
    output:
        temp("Temporary/Pausing/{strand}/{chr}.bed")
    params:
        window_size = config["window_size"],
        minimum_intensity = config["minimum_intensity"],
        include_spikes = config["include_spikes"],
        chromosome = "{chr}",
        strand = "{strand}"
    script:
        "Scripts/getPausing.py"

rule MaskRegions:
    input:
        "Temporary/Pausing/{strand}/{chr}.bed"
    output:
        temp("Temporary/Masked/{strand}/{chr}.bed")
    params:
        config["mask_bed"] if config["mask_bed"]!="" else ""
    script:
        "Scripts/maskRegion.py"

rule TestPausing:
    input:
        "Temporary/Masked/{strand}/{chr}.bed"
    output:
        temp("Temporary/Tested/{strand}/{chr}.bed")
    params:
        bootstrap_N = config["bootstrap_N"],
        window_size = config["window_size"],
        strand = "{strand}",
        test=config["test"]
    script:
        "Scripts/testPausing.py"

rule FilterPeaks:
    input:
        lambda wildcards: expand("Temporary/Tested/{strand}/{chr}.bed", strand=STRANDS, chr=CONTIGS)
    output:
        "Results/Significant_{prefix}.bed",
        "Results/All_{prefix}.bed"
    params:
        config["percentile"]
    script:
        "Scripts/filterPausing.py"

rule AssignRegion:
    input:
        "Results/Significant_{prefix}.bed"
    output:
        "Results/{prefix}.gff",
        "Results/Regions_{prefix}.bed",
        "Results/Region_{prefix}.bed"
    params:
        config["region_bed"] if config["region_bed"]!="" else ""
    script:
        "Scripts/assignRegion.py"
