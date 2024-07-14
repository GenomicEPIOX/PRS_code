import os
import pandas as pd
# define chromosomes
CHROM = list(range(1, 2))
#CHROM.append("X") 
# read in list of summary stats
sumstats_info = pd.read_csv("input.txt", sep="\t", header=None)
sumstats_info.columns = ["sample", "sumstat_path", "N_gwas"]
sumstats_info = sumstats_info.set_index('sample')

def path_from_sample(wildcards):
  return sumstats_info.loc[wildcards.sample, "sumstat_path"]

def N_gwas_from_sample(wildcards):
  return sumstats_info.loc[wildcards.sample, "N_gwas"]



# define the final rule that depends on both prs_cs_chr and prs_cs_concat
rule all:
    input:
        expand("output/{sample}/{sample}.panel",  sample=sumstats_info.index.tolist() )


# rule to create a separate profile for each chromosome



rule untar:
    input:
        path_from_sample
    resources:
        mem='32G',
        time='00:20:00',
        parition="short"
    output:
        touch("output/{sample}/.naughty")
       
    params:
        folder = "output/{sample}/"
    shell:
        """
        mkdir -p {params.folder}
        tar xvf {input} -C {params.folder} 
        mv output/{wildcards.sample}/*/*gz output/{wildcards.sample}/
        
       
        """

rule preprocess:
    resources:
        mem='32G',
        time='00:10:00',
        parition="short"
    input:
        "output/{sample}/.naughty"
    params:
        chrom = "{chrom}",
        A1 = "effectAllele",
        A2 = "otherAllele",
        SNP = "rsids",
        P = "Pval",
        BETA = "Beta",
        pre_script = "/well/travis-prostate/projects/PRS_CS/annotate_chromsomes.py",
        sample = "{sample}" 
    output:
        "output/{sample}/chr{chrom}.temp"
    shell:
        """
        mv output/{params.sample}/*_chr{params.chrom}_*.gz output/{params.sample}/chr{params.chrom}.gz
        python {params.pre_script} -i output/{params.sample}/chr{params.chrom}.gz -k /well/travis-prostate/projects/annotation/windows/UKBB_rsids_keys/ -o output/{params.sample}/chr{params.chrom}.temp 

        """


rule prs_cs_chr:
    input:
        "output/{sample}/chr{chrom}.temp"
    output:
        "output/{sample}/{sample}.chr{chrom}.txt"
    params:
        ref_dir="/well/travis-prostate/projects/PRS_CS/ldblk_ukbb_eur",
        bim_prefix="/well/travis-prostate/projects/PRS_CS/1000G/",
        gwas_size=N_gwas_from_sample ,
        chrom = "{chrom}",
        sample_name = "{sample}"
    threads: 8
    resources:
        mem='32G',
        time='04:00:00',
        partition="short"
    shell:
        """
        python /well/travis-prostate/projects/PRS_CS/PRScs/PRScs.py \
        --ref_dir={params.ref_dir} \
        --bim_prefix={params.bim_prefix}chr{params.chrom} \
        --sst_file={input} \
        --n_gwas={params.gwas_size} \
        --out_dir=output/{params.sample_name}/{params.sample_name} \
        --chrom={params.chrom}

        mv output/{params.sample_name}/{params.sample_name}*chr{params.chrom}.txt output/{params.sample_name}/{params.sample_name}.chr{params.chrom}.txt 
        """
       




# rule to concatenate individual chromosome profiles
rule prs_cs_concat:
    input:
        expand("output/{sample}/{sample}.chr{chrom}.txt", chrom=CHROM, sample=sumstats_info.index.tolist())
    threads: 1
    resources:
        mem='8G',
        time='04:00:00',
        partition="short"
    output:
        panel="output/{sample}/{sample}.panel"
    shell:
        """
        cat {input} > {output.panel}
        #rm {input}
        #rm output/{wildcards.sample}/*.temp
        #rm output/{wildcards.sample}/*.gz
        #rm -r output/{wildcards.sample}/*_{wildcards.sample}_* 
    
    
	"""

