# Attempt to make a local running RNAseq standard pipeline.  Starting with Manu's simpler format but modifying it to run Kat's pipeline

# 20220718 changed "metadata" to "config"
configfile: "config/config.yaml"

import pandas as pd
import os

# 20220718 changed "units_table" to "samples"
# 20220718 changed "samples_file" to "samples"
samples = pd.read_csv(config["samples"], sep="," , dtype = str).set_index("Run", drop=False)

groups=pd.read_csv(config["comparisons"], sep =",")
# control_groups=list(groups.Control_Group)
# treatment_groups=list(groups.Treatment_Group)


# 20220718 changed "files" to "results"
rule all:
	input:
		expand('results/deltapsi_results/{comparison}.deltapsi.tsv',comparison=groups.Comparison),
		expand("results/bam_results/{sample}_Aligned.sortedByCoord.out.bam", sample=samples.sample_name),
		expand("results/salmon/{run}_quant/quant.sf",run=samples.Run)

def get_raw_reads(wildcards):
	return config["samples"][wildcards.sample]

rule SRR_download:
	input:
		get_raw_reads
	output:
		expand("/Volumes/Bioinfo/raw_data/tmp/{{run}}_{num}.fastq.gz", num=[1,2])
	conda:
		"workflow/envs/SraTools.yaml"
	threads: 10
	params:
		path="results/fastq_raw/"
	shell:
		"fastq-dump --gzip --split-3 -O {params.path} {input}"

#rule SRR_convert:
#	input:
#		srr="results/fastq_results/{run}/{run}.sra"
#	params:
#		path="results/fastq_results/"
#	output:
#		expand("results/fastq_results/{{run}}_{num}.fastq.gz", num=[1,2])
#	conda:
#		"envs/SraTools.yaml"
#	resources:
#		cpu=1,
#		mem=lambda wildcards, attempt: attempt * 4,
#		time=360
#	shell:
#		"fastq-dump --gzip --split-3 -O {params.path} {input.srr}"



#rule STAR_align:
#	input:
#		genome=config["STAR_genome"],
#		results=lambda wildcards: expand("results/fastq_results/{run}_{num}.fastq.gz", run=samples.Run[samples.sample_name == wildcards.sample], num=[1,2])
#	params:
#		path="results/bam_results/{sample}_"
#	resources:
#		cpu=10,
#		mem=lambda wildcards, attempt: attempt * 120
#	output:
#		"results/bam_results/{sample}_Aligned.sortedByCoord.out.bam"
#	shell:
#		"STAR --twopassMode Basic --genomeDir {input.genome} --outTmpKeep None "
#		"--readresultsIn {input.results} --readresultsCommand zcat "
#		"--runThreadN {resources.cpu} --outSAMtype BAM SortedByCoordinate "
#		"--outFileNamePrefix {params.path} --alignSJoverhangMin 8 "
#		"--limitBAMsortRAM {resources.mem}000000000 --outSAMattributes All "
#		"--quantMode GeneCounts"

# rule samtools_merge_bam:
# 	input:
# 		lambda wildcards: expand('results/bam_results/{run}_Aligned.sortedByCoord.out.bam',
# 			run= samples.Run[samples.source_name == wildcards.sample])
# 	output:
# 		bam = 'results/bam_results/{sample}_merged.bam'
# 	resources:
# 		cpu=2,
# 		mem=lambda wildcards, attempt: attempt * 8
# 	shell:
# 		"""
# 		PS1=dummy
# 		. $(conda info --base)/etc/profile.d/conda.sh
# 		conda activate samtools
# 		samtools merge {output.bam} {input}
# 		"""


r#ule samtools_index:
#	input:
#		'results/bam_results/{sample}_Aligned.sortedByCoord.out.bam'
#	output:
#		'results/bam_results/{sample}_Aligned.sortedByCoord.out.bam.bai'
#	conda:
#		"envs/SamTools_env.yaml"
#	resources:
#		cpu=10,
#		mem=lambda wildcards, attempt: attempt * 10
#	shell:
#		"samtools index -@ {resources.cpu} {input} {output}"


#rule config_majiq:
#	input:
#		bam=expand('results/bam_results/{sample}_Aligned.sortedByCoord.out.bam',sample=samples.sample_name),
#		bai=expand('results/bam_results/{sample}_Aligned.sortedByCoord.out.bam.bai',sample=samples.sample_name)
#	output:
#		'config/majiq_config_file.txt'
#	resources:
#		cpu=1,
#		mem=lambda wildcards, attempt: attempt * 1
#	params:
#		samples=list(samples.sample_name),
#		genome="hg38"
#	script:
#		'scripts/config_creator_majiq.py'

#rule majiq_build:
#	input:
#		bam=expand('results/bam_results/{sample}_Aligned.sortedByCoord.out.bam',sample=samples.sample_name),
#		bai=expand('results/bam_results/{sample}_Aligned.sortedByCoord.out.bam.bai',sample=samples.sample_name),
#		genome=config["genome"],
#		majiq_config=config["build_config"]
#	params:
#		path='results/majiq_results/'
#	conda:
#		"envs/majiq_env.yaml"
#	resources:
#		cpu=8,
#		mem=lambda wildcards, attempt: attempt * 8
#	output:
#		majiq_results=expand('results/majiq_results/{sample}_Aligned.sortedByCoord.out.majiq',sample=samples.sample_name),
#		sj_results=expand('results/majiq_results/{sample}_Aligned.sortedByCoord.out.sj',sample=samples.sample_name),
#		splicegraph='results/majiq_results/splicegraph.sql'
#		log='results/majiq_results/majiq.log'
#	shell:
#		"majiq build --conf {input.majiq_config} --nproc {resources.cpu} "
#		"--disable-ir --simplify -o {params.path} {input.genome}"

# #rule to create new column with underscore between names and loop
# #this will in sense merge the two columns

#rule delta_psi:
#	input:
#		control_majiq_results=lambda wildcards: expand(
#			'results/majiq_results/{group_control}{rep}_Aligned.sortedByCoord.out.majiq',
#			group_control=groups.Control_Group[groups.Comparison == wildcards.comparison],rep=["_rep1","_rep2","_rep3"]),
#		treatment_majiq_results=lambda wildcards: expand(
#			'results/majiq_results/{group_treatment}{rep}_Aligned.sortedByCoord.out.majiq',
#			group_treatment=groups.Treatment_Group[groups.Comparison == wildcards.comparison],rep=["_rep1","_rep2","_rep3"])
#	conda:
#		"envs/majiq_env.yaml"
#	resources:
#		cpu=4,
#		mem=lambda wildcards, attempt: attempt * 64
#	params:
#		path="results/deltapsi_results",
#		control=lambda wildcards: expand(
#			'{group_control}',group_control=groups.Control_Group[groups.Comparison == wildcards.comparison]),
#		treated=lambda wildcards: expand(
#			'{group_treatment}',group_treatment=groups.Treatment_Group[groups.Comparison == wildcards.comparison])
#	output:
#		'results/deltapsi_results/{comparison}.deltapsi.voila'
#	shell:
#		"majiq deltapsi --output-type voila -grp1 {input.control_majiq_results} -grp2 {input.treatment_majiq_results} "
#		"-n {params.control} {params.treated} -j {resources.cpu} -o {params.path}"

#rule voila_tsv:
#	input:
#		voila='results/deltapsi_results/{comparison}.deltapsi.voila',
#		sql='results/majiq_results/splicegraph.sql'
#	resources:
#		cpu=4,
#		mem=lambda wildcards, attempt: attempt * 64
#	conda:
#		"envs/majiq_env.yaml"
#	output:
#		'results/deltapsi_results/{comparison}.deltapsi.tsv'
#	shell:
#		"voila tsv -f {output} --show-all -j {resources.cpu} {input.voila} {input.sql}"

#rule psi:
#	input:
#		majiq_results=expand('results/majiq_results/{sample}_Aligned.sortedByCoord.out.majiq',
#			sample=samples)
#	resources:
#		cpu=4,
#		mem="4G"
#	conda:
#		"envs/majiq_env.yaml"
#	params:
#		path="results/psi_results"
#	output:
#		tsv_results='results/psi_results/{sample}.psi.tsv',
#		voila_results='results/psi_results/{sample}.psi.voila'
#	shell:
#		"majiq psi --output-type all -j {resources.cpu} {input.majiq_results} -o {params.path}"



#rule salmon_quant:
#        input:
#                "results/fastq_results/{run}_1.fastq.gz", "results/fastq_results/{run}_2.fastq.gz"
#        resources:
#                cpu=8,
#                time="96:00:00",
#                mem=lambda wildcards, attempt: attempt * 25
#        params:
#                genome=config["salmon_genome"]
#        output:
#                "results/salmon/{run}_quant/quant.sf"
#        shell:
#                "/home/torresdizm/bin/salmon quant "
#                "-i {params.genome} "
#                "-l A "
#                "-1 {input[0]} -2 {input[1]} "
#                "-p 8 --validateMapping "
#                "-o ./results/salmon/{wildcards.run}_quant"
