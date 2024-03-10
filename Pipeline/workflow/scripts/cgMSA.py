#Author:  Mustafa Helal
#!/usr/bin/python
import os
import pandas as pd
import glob
import numpy as np
import multiprocessing
from multiprocessing import Pool
import argparse


#generate multifasta files
def generate_multifasta(prot_files, df, schema_seed, dir_name_1):
	read_files = df.index.tolist()
	file_header = df.columns.tolist()
	for i in range(len(read_files)):
		fasta_name = read_files[i] + ".fasta"
		with open(schema_seed+fasta_name, "r") as ID_file, open(dir_name_1 + "/"+fasta_name,"w+") as multifasta:
			Lines = ID_file.readlines()
			for header in file_header:
				index = df.at[read_files[i], header]
				line_num = 0
				index_count = 0
				for line in Lines:
					if line.startswith(">"):
						index_2 = line.rsplit("_")[-1]
						index_2 = index_2.split("\n")[0]
						if str(index_2) == str(index):
							multifasta.write(Lines[line_num].strip("\n") + ":ID-" + header + "\n")
							multifasta.write(Lines[line_num+1])
							break
						if index == 0:
							multifasta.write(Lines[line_num].strip("\n") + ":ID-" + header + "\n")
							multifasta.write("X"+ "\n")
							break
					line_num += 1

#alignment function for parallelization
def align(list):
	thread = list[0]
	prot = list[1]
	name = list[2]
	dir_name_2 = list[3]

	for i in range(len(prot)):
		os.system("muscle -in " + prot[i] + " -out " + dir_name_2 + "/" + name[i] + ".afa")


def main():

	parser = argparse.ArgumentParser(description="Build a core genome phylogeny using assembly or paird-end read files.")
	parser.add_argument("-i", "--input", required=True, help="Input path to files directory.")
	parser.add_argument("-t", "--threads", default=4, required=False, help="Amount of threads to use. Default is 4.")
	parser.add_argument("-o", "--output_dir", required=True, help="Output path for the generated files.")
	# parser.add_argument("-m", "--model", required=False, default="GTR+G", help="Substitution Model für RaxML-ng. Default is GTR+G")
	# parser.add_argument("-th", "--threshold", required=False, default=0.95, help="Clustering Threshold for ChewBBACA. Default is 0.95")
	args = parser.parse_args()

	#define input
	threads = args.threads #amount of threads

	schema_seed = args.input + "schema_seed/"	#path to gene allele_type files
	table = args.input + "cgMLST_ref/cgMLST.tsv" #path to mlst_table.tsv

	#import prot_files
	prot_files = [f for f in glob.glob(schema_seed + "*.fasta")]
	df = pd.read_csv(table, sep = "\t")
	df.set_index('FILE', inplace=True)

	#created dir for multifasta results
	dir_name_1 = args.output_dir + "/multifasta_alignments_results"
	if os.path.exists(dir_name_1):
		os.system("rm -rf " + dir_name_1)
		os.system("mkdir -p " + dir_name_1)
		print("Directory " , dir_name_1 ,  " Created ")
	else:
		os.system("mkdir -p " + dir_name_1)
		print("Directory " , dir_name_1 ,  " Created ")


	generate_multifasta(prot_files, df.T, schema_seed, dir_name_1)

	# second file import for alignments
	multi_fasta_files = [f for f in glob.glob(dir_name_1 +"/"+"*.fasta")]
	result_names = [name.split(".")[0] for name in os.listdir(schema_seed) if name.endswith(".fasta")]


	#create dir for alignment results
	dir_name_2 = args.output_dir + "/gene_alignment_results"
	if os.path.exists(dir_name_2):
		os.system("rm -rf " + dir_name_2)
		os.system("mkdir -p " + dir_name_2)
		print("Directory " , dir_name_2 ,  " Created ")
	else:
		os.system("mkdir -p " + dir_name_2)
		print("Directory " , dir_name_2 ,  " Created ")

	cpu_threads = multiprocessing.cpu_count()


	# set amount of thread and parallele functions here

	amount_threads = int(threads)
	if cpu_threads >= amount_threads:
		print("using " + str(amount_threads) + " threads:")
		#split files for parallel usage
		prot_list = np.array_split(multi_fasta_files, amount_threads)
		name_list = [elem.split("/")[-1] for elem in multi_fasta_files]
		name_list = [elem.split(".")[0] for elem in name_list]
		name_list = np.array_split(name_list, amount_threads)
	else:
		amount_threads = cpu_threads-3 # -3 to not occupy all the cores in the system
		print("passed more threads than available -- running with " + str(amount_threads) + " Threads")
		#split files for parallel usage
		prot_list = np.array_split(multi_fasta_files, amount_threads)
		name_list = [elem.split("/")[-1] for elem in multi_fasta_files]
		name_list = [elem.split(".")[0] for elem in name_list]
		name_list = np.array_split(name_list, amount_threads)

	#Parallel function
	if cpu_threads >= amount_threads:
		pool = Pool(amount_threads)
		params = [[i+1, prot_list[i].tolist(), name_list[i].tolist(), dir_name_2] for i in range(amount_threads)]
		results = pool.map(align, params)

	else:
		pool = Pool(cpu_threads-3)
		params = [[i+1, prot_list[i].tolist(), name_list[i].tolist(), dir_name_2] for i in range(cpu_threads-3)]
		results = pool.map(align, params)

	#create dir for concat results
	dir_name_3 = args.output_dir + "/concat_alignment"
	if os.path.exists(dir_name_3):
		os.system("rm -rf " + dir_name_3)
		os.system("mkdir -p " + dir_name_3)
		print("Directory " , dir_name_3 ,  " Created ")
	else:
		os.system("mkdir -p " + dir_name_3)
		print("Directory " , dir_name_3 ,  " Created ")


	alignment_files = os.listdir(dir_name_2)
	id_names = df.index.tolist()
	prot_names = df.columns.tolist()


	#concat all files
	with open(dir_name_3 + "/concat_file.afa", "w") as write_file:
		for id in id_names:
			ID = "ID-" + id
			write_file.write(">"+id+"\n")
			print(ID)
			for file in prot_names:
				with open(dir_name_2 + "/" + file + ".afa", "r") as read_file:
					lines = read_file.readlines()
					num_lines_between = 0
					j = 1
					while lines[j][0] != ">":
						num_lines_between += 1
						j += 1
					for i in range(len(lines)):
						ID_FILE = lines[i].split(":")[-1]
						ID_FILE = ID_FILE.replace("\n", "")
						if ID_FILE == ID:
							END = 1
							while END != num_lines_between+1:
								write_line = lines[i+1].replace("X", "-")
								write_file.write(write_line)
								END += 1
								i += 1
						else:
							continue


if __name__ == "__main__":
	try:
		main()
	except Exception as e:
		print("\033[91m An error occurred: " + str(e))
