SIF_EXEC = "/cluster/tools/software/centos7/singularity/3.5.2/bin/singularity exec config_slurm/ml_repeat_pipeline_v1.0.sif"

import pandas as pd
from snakemake.utils import validate, min_version
import os.path as path
import itertools

##### set minimum snakemake version #####
min_version("5.5.4")




SCRIPTS_DIR = "scripts"
RESULTS_DIR = "results"
MATRIX_DIR = "matrix"
CHROM_DIR = "chromvar"
PLOTS_DIR = "plots_chromvar"

configfile: "config.yaml"

######Groups##########
GROUPS=[]
groups=set(pd.read_table(config["metadatafile"],sep='\s+')["Group"].tolist())
groupcombinations=list(itertools.combinations(groups, 2))
for i in groupcombinations:
	v=sorted([i[0],i[1]])
	GROUPS.append(v[0]+"_vs_"+v[1])


# ==============================================================================
# Meta Rules
# ==============================================================================
rule all:
    input:
        ## STEP1_generate_binarymatrix
        path.join(RESULTS_DIR, MATRIX_DIR, config["opname"]+".consensus.Binarymat.rds"),
        path.join(RESULTS_DIR, MATRIX_DIR, config["opname"]+".consensus.Binarymat.repeats.rds"),
        
        ## STEP2_run_chromvar
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".anno_ix.Rdata"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".counts_filtered.Rdata"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".devobj.Rdata"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Zscore.txt"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Variability.txt"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Deviations.txt"),
        
        ## STEP3_chromvar_analysis
        expand(path.join(RESULTS_DIR, PLOTS_DIR, "{groups}.DiffDeviations.with.medians.txt"),groups=GROUPS),
        expand(path.join(RESULTS_DIR, PLOTS_DIR, "Heatmap.{groups}.DiffDeviations.qval"+str(config["qval"])+".pdf"),groups=GROUPS)
        
        

# ==============================================================================
# Rules
# ==============================================================================
rule generate_binarymatrix:
    input:
        script = path.join(SCRIPTS_DIR, "01_SLURM_generate_binarymat.sh"),
        consensusfile = config["consensusfile"],
        PEAKS_DIR = config["PEAKS_DIR"],
        hg38TEs_DIR = config["hg38TEs_DIR"],
    output:
        path.join(RESULTS_DIR, MATRIX_DIR, config["opname"]+".consensus.Binarymat.rds"),
        path.join(RESULTS_DIR, MATRIX_DIR, config["opname"]+".consensus.Binarymat.repeats.rds")
    params:
        config["peaks_extension"]+" "+config["hg38TEs_extension"]+" "+path.join(RESULTS_DIR, MATRIX_DIR)+" "+config["opname"]
    shell:
        "{SIF_EXEC} bash {input.script} {input.consensusfile} {input.PEAKS_DIR} {input.hg38TEs_DIR} {params}"


rule run_chromvar:
    input:
        script = path.join(SCRIPTS_DIR, "02_run_chromvar_slurm2.h4h.sh"),
        peaks_binarymat = path.join(RESULTS_DIR, MATRIX_DIR, config["opname"]+".consensus.Binarymat.rds"),
        peaks_binarymat_repeats = path.join(RESULTS_DIR, MATRIX_DIR, config["opname"]+".consensus.Binarymat.repeats.rds")
    output:
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".anno_ix.Rdata"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".counts_filtered.Rdata"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".devobj.Rdata"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Zscore.txt"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Variability.txt"),
        path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Deviations.txt")
    params:
        config["opname"]
    shell:
        "{SIF_EXEC} bash {input.script} {input.peaks_binarymat} {input.peaks_binarymat_repeats} {params} {RESULTS_DIR}"
        
rule chromvar_analysis:
    input:
        script = path.join(SCRIPTS_DIR, "03_Chromvar_analysis.R"),
        metadatafile = config["metadatafile"],
        devobjfile = path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".devobj.Rdata"),
        zscorefile = path.join(RESULTS_DIR, CHROM_DIR, config["opname"]+".Zscore.txt"),
    output:
    	expand(path.join(RESULTS_DIR, PLOTS_DIR, "{groups}.DiffDeviations.with.medians.txt"),groups=GROUPS),
        expand(path.join(RESULTS_DIR, PLOTS_DIR, "Heatmap.{groups}.DiffDeviations.qval"+str(config["qval"])+".pdf"),groups=GROUPS)
    params:
        config["opname"]+" "+str(config["qval"])
    shell:
        "{SIF_EXEC} Rscript {input.script} {input.metadatafile} {input.zscorefile} {input.devobjfile} {RESULTS_DIR} {params} {GROUPS}"





