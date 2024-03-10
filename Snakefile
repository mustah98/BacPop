import os
import sys
import importlib
import pandas as pd
from datetime import datetime



# Set config file and input and absolut path
path_abs = os.path.abspath(os.getcwd()) + "/"
configfile: "Pipeline/config/config.yaml"
input_path = config["data"]


# Chek for old results if option is on
if config["global"]["check_results"]:
	results_folder = sorted(os.listdir(path_abs + "/results"))
	if len(results_folder) > 0:
		print("A run folder already exists.")
		user_input = input("Do you want to continue with the existing run folder? (yes/no): ")
		if user_input.lower() == "yes":
			date_folder = results_folder[-1]
		else:
			print("Generating new results folder.")
			current_datetime = datetime.now()
			date_folder = current_datetime.strftime("%d-%m-%Y-%H-%M")
			os.system("mkdir -p " + path_abs + "/results/" + date_folder)
else:
	print("Generating new results folder.")
	current_datetime = datetime.now()
	date_folder = current_datetime.strftime("%d-%m-%Y-%H-%M")
	os.system("mkdir -p " + path_abs + "/results/" + date_folder)


# Set table for Snakemake file forwarding

if config["global"]["files"] == "fastq":

	files = [file for file in (sorted(list(os.listdir(input_path)))) if file.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
	sample_names = [name.split(".")[0] for name in files][::2]
	R1 = list([input_path + name for name in files])[0::2]
	R2 = list([input_path + name for name in files])[1::2]

	sample_table = {
		"index" :   sample_names,
		"R1"    :   R1,
		"R2"    :   R2,
	}

	sample_table = pd.DataFrame.from_dict(data = sample_table)
	sample_table = sample_table.set_index("index")

elif config["global"]["files"] == "fasta":

	files = [file for file in (sorted(list(os.listdir(input_path)))) if file.endswith((".fasta", ".fasta.gz", ".fa", ".fa.gz"))]
	sample_names = [name.split(".")[0] for name in files]
	assemblys = list([input_path + file for file in files])

	sample_table = {
		"index" : sample_names,
		"assembly" : assemblys,
	}

	sample_table = pd.DataFrame.from_dict(data = sample_table)
	sample_table = sample_table.set_index("index")

print(sample_table)

def rule_all():
	rule_all = []

	if config["global"]["files"] == "fastq":
		if config["global"]["qc"]:
			rule_all.append(expand(path_abs + "results/" + date_folder + "/QC/{sample}_trimm_1_fastqc.html", sample = sample_names))
			rule_all.append(expand(path_abs + "results/" + date_folder + "/QC/{sample}_trimm_2_fastqc.html", sample = sample_names))

		if config["global"]["annotation"] and os.path.isdir(config["annotation"]["db"]):
			rule_all.append(expand(path_abs + "results/" + date_folder + "/annotations/{sample}/{sample}.gff3", sample = sample_names))
			rule_all.append(expand(path_abs + "results/" + date_folder + "/all_gffs/{sample}.gff3", sample = sample_names))

		if config["global"]["amr"]:
			rule_all.append(expand(path_abs + "results/" + date_folder + "/AMR/{sample}.tsv", sample = sample_names))

		if config["global"]["mlst"]:
			rule_all.append(path_abs + "results/" + date_folder + "/mlst/mlst.csv")

		rule_all.append(expand(path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_1.fastq.gz", sample = sample_names))
		rule_all.append(expand(path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_2.fastq.gz", sample = sample_names))
		rule_all.append(expand(path_abs + "results/" + date_folder + "/assemblies/{sample}", sample = sample_names))
		rule_all.append(expand(path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta", sample = sample_names))
		rule_all.append(path_abs + "results/" + date_folder + "/cgMLST/")
		rule_all.append(path_abs + "results/" + date_folder + "/cgMSA/")
		rule_all.append(path_abs + "results/" + date_folder + "/phylogenetic_tree/")

	elif config["global"]["files"] == "fasta":
		if config["global"]["annotation"] and os.path.isdir(config["annotation"]["db"]):
			rule_all.append(expand(path_abs + "results/" + date_folder + "/annotations/{sample}/{sample}.gff3", sample = sample_names))
			rule_all.append(expand(path_abs + "results/" + date_folder + "/all_gffs/{sample}.gff3", sample = sample_names))

		if config["global"]["amr"]:
			rule_all.append(expand(path_abs + "results/" + date_folder + "/AMR/{sample}.tsv", sample = sample_names))

		if config["global"]["mlst"]:
			rule_all.append(path_abs + "results/" + date_folder + "/mlst/mlst.csv")

		rule_all.append(expand(path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta", sample = sample_names))
		rule_all.append(path_abs + "results/" + date_folder + "/cgMLST/")
		rule_all.append(path_abs + "results/" + date_folder + "/cgMSA/")
		rule_all.append(path_abs + "results/" + date_folder + "/phylogenetic_tree/")


	print(rule_all)
	return rule_all


rule all:
	input:
		rule_all()


###########
## RULES ##
###########


rule qc:
	input:
		R1 = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_1.fastq.gz",
		R2 = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_2.fastq.gz",
	output:
		OUT1 = path_abs + "results/" + date_folder + "/QC/{sample}_trimm_1_fastqc.html",
		OUT2 = path_abs + "results/" + date_folder + "/QC/{sample}_trimm_2_fastqc.html",
	conda: "Pipeline/workflow/envs/qc.yml"
	threads: int(config["global"]["threads"]),
	shell:
		"fastqc {input.R1} {input.R2} -o results/" + date_folder +  "/QC/"

#----------------------------------------------------------------------------------------------------------------------------

rule trimm:
	input:
		R1 = lambda wildcards: sample_table.at[wildcards.sample, "R1"] if wildcards.sample in sample_table.index else "",
		R2 = lambda wildcards: sample_table.at[wildcards.sample, "R2"] if wildcards.sample in sample_table.index else "",
	output:
		R1 = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_1.fastq.gz",
		R2 = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_2.fastq.gz",
		HTML = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}.html",
		JSON = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}.json",
	conda: "Pipeline/workflow/envs/trimm.yml"
	threads: int(config["global"]["threads"]),
	shell:
		"fastp -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}  -h results/" + date_folder +  "/trimmed_sequences/{wildcards.sample}.html -j  results/" + date_folder +  "/trimmed_sequences/{wildcards.sample}.json --qualified_quality_phred 20 --length_required 50"

#----------------------------------------------------------------------------------------------------------------------------

rule assembly:
	input:
		R1 = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_1.fastq.gz",
		R2 = path_abs + "results/" + date_folder + "/trimmed_sequences/{sample}_trimm_2.fastq.gz",
	output:
		OUT = directory(path_abs + "results/" + date_folder + "/assemblies/{sample}"),
	params:
		MEMORY  = config["global"]["memory"],
		PATH = path_abs + "results/" + date_folder + "assemblies/"
	conda: "Pipeline/workflow/envs/assembly.yml"
	threads: int(config["global"]["threads"]),
	shell:
		"spades.py -t {threads} -m {params.MEMORY} --isolate -1 {input.R1} -2 {input.R2} -o {output.OUT}"

#----------------------------------------------------------------------------------------------------------------------------

rule move_assemblies:
	input:
		ASSEMBLY = path_abs + "results/" + date_folder + "/assemblies/{sample}/",
	output:
		ASSEMBLY = path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta",
	shell:
		"cp {input.ASSEMBLY}/contigs.fasta {output.ASSEMBLY}"

#----------------------------------------------------------------------------------------------------------------------------

rule annotation:
	input:
		ASSEMBLY = path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta",
	output:
		OUT = path_abs + "results/" + date_folder + "/annotations/{sample}/{sample}.gff3",
	params:
		DB = config["annotation"]["db"],
		GENUS = config["annotation"]["genus"],
		SPECIES = config["annotation"]["species"],
	conda: "Pipeline/workflow/envs/annotation.yml"
	threads: int(config["global"]["threads"]),
	shell:
		"bakta {input.ASSEMBLY} --db {params.DB} -o {output.OUT} --threads {threads} --genus {params.GENUS} --species {params.SPECIES}"

#----------------------------------------------------------------------------------------------------------------------------

rule move_annotations:
	input:
		ANNOATION = path_abs + "results/" + date_folder + "/annotations/{sample}/{sample}.gff3"
	output:
		ALL_GFF = path_abs + "results/" + date_folder + "/all_gffs/{sample}.gff3"
	shell:
		"cp {input.ANNOATION} {output.ALL_GFF}"

#----------------------------------------------------------------------------------------------------------------------------

rule mlst:
	input:
		ASSEMBLY = expand(path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta", sample = sample_names),
	output:
		OUT = path_abs + "results/" + date_folder + "/mlst/mlst.csv",
	params:
		SED = path_abs + "results/" + date_folder + "/assemblies_contigs/",
		PATH = path_abs + "results/" + date_folder + "/mlst/",
	threads: int(config["global"]["threads"]),
	conda: "Pipeline/workflow/envs/mlst.yml"
	shell:
		"""
			mlst  {input.ASSEMBLY} --nopath --csv --minid 90 --mincov 90 --threads {threads} > {output.OUT}
			cd {params.PATH} && sed -i 's#.fasta##' mlst.csv
		"""

#----------------------------------------------------------------------------------------------------------------------------

rule roary:
	input:
		ANNOTATIONS = expand(path_abs + "results/" + date_folder + "/all_gffs/{sample}.gff3", sample = sample_names)
	output:
		OUT = directory(path_abs + "results/" + date_folder + "/roary/"),
	conda: "Pipeline/workflow/envs/roary.yml",
	params:
		ANNOTATIONS_DIR = path_abs + "results/" + date_folder + "/all_gffs/",
	shell:
		"roary -i 90 -f {output.OUT} -p 32 {params.ANNOTATIONS_DIR}"

#----------------------------------------------------------------------------------------------------------------------------


rule amr:
	input:
		ASSEMBLY = path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta",
	output:
		OUT = path_abs + "results/" + date_folder + "/AMR/{sample}.tsv",
	params:
		MINID = config["amr"]["minid"],
		MINCOV = config["amr"]["mincov"],
		DATABASE = config["amr"]["database"],
	conda: "Pipeline/workflow/envs/amr.yml"
	threads: int(config["global"]["threads"]),
	shell:
		"abricate {input.ASSEMBLY} --threads {threads}  --minid {params.MINID} --mincov {params.MINCOV} > {output.OUT}"

#----------------------------------------------------------------------------------------------------------------------------

rule cgMLST:
	input:
		CHECK = expand(path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta", sample = sample_names),
	output:
		CGMLST = directory(path_abs + "results/" + date_folder + "/cgMLST/"),
	params:
		THRESHOLD = config["cgmlst"]["threshold"],
		ASSEMBLYS = path_abs + "results/" + date_folder + "/assemblies_contigs/",
	threads: int(config["global"]["threads"]),
	conda: "Pipeline/workflow/envs/cgmlst.yml"
	shell:
		"""
			if [ -d "{output.CGMLST}" ]; then
				rm -rf cgMLST
			fi
			chewBBACA.py CreateSchema -i {params.ASSEMBLYS} --cpu {threads} -o {output.CGMLST}
			cd {output.CGMLST} && chewBBACA.py AlleleCall -i {params.ASSEMBLYS} -g schema_seed/ -o results_cg --cpu {threads}
			cd {output.CGMLST} && chewBBACA.py RemoveGenes -i results_cg/results_alleles.tsv -g results_cg/paralogous_counts.tsv -o alleleCallMatrix_cg.tsv
			cd {output.CGMLST} && chewBBACA.py ExtractCgMLST -i alleleCallMatrix_cg.tsv -o cgMLST_ref --t {params.THRESHOLD}
			cd {output.CGMLST}/cgMLST_ref/ && sed 's/.fasta//g' cgMLST*.tsv > cgMLST.tsv
		"""

#----------------------------------------------------------------------------------------------------------------------------

rule cg_alignment:
	input:
		CGMLST = path_abs + "results/" + date_folder + "/cgMLST/",
	output:
		MSA = directory(path_abs + "results/" + date_folder + "/cgMSA/"),
	conda: "Pipeline/workflow/envs/cgmlst.yml",
	threads: int(config["global"]["threads"]),
	shell:
		"python3 Pipeline/workflow/scripts/cgMSA.py -i {input.CGMLST}/ -o {output.MSA} -t {threads}"

#----------------------------------------------------------------------------------------------------------------------------

rule tree:
	input:
		MSA = path_abs + "results/" + date_folder + "/cgMSA/"
	output:
		TREE = directory(path_abs + "results/" + date_folder + "/phylogenetic_tree/")
	params:
		MODEL = config["tree"]["model"],
	conda: "Pipeline/workflow/envs/cgmlst.yml",
	threads: int(config["global"]["threads"]),
	shell:
		"""
			mkdir -p {output.TREE}
			raxml-ng --all --force perf_threads --prefix {output.TREE}/CoreGenomePhylogeny --msa {input.MSA}/concat_alignment/concat_file.afa --model {params.MODEL} --threads  {threads}
		"""

#----------------------------------------------------------------------------------------------------------------------------

rule set_assemblys_path:
	input:
		ASSEMBLYS = lambda wildcards: sample_table.at[wildcards.sample, "assembly"] if wildcards.sample in sample_table.index else "",
		# ASSEMBLYS = config["data"]
	output:
		ASSEMBLYS = path_abs + "results/" + date_folder + "/assemblies_contigs/{sample}.fasta"
	shell:
		"cp {input.ASSEMBLYS} {output.ASSEMBLYS}"
