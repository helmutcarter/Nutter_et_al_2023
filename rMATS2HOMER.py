#!/usr/bin/env python3
import re
from pyfaidx import Fasta
from Bio.Seq import Seq
import os

#settings
FDR_cutoff = 0.05
use_count_cutoff = True
count_cutoff=5
event_list = ["SE"]
SJ_flank = 100
rMATS_folder = "/path/to/directory/containing/rMATS/output"
path_to_fasta = "/path/to/appropriate/genome/assembly/Mus_musculus.GRCm38.dna.primary_assembly.fa"

output_folder = rMATS_folder + "_HOMER_output_" + str(SJ_flank) + "bp"
splicing_dict = {
				"SE":{
						"alt_exon":{"5P":"exonStart_0base","3P":"exonEnd"},
						"US_intron_5P":{"5P":"upstreamEE","3P":"exonStart_0base"},
						"US_intron_3P":{"5P":"upstreamEE","3P":"exonStart_0base"},
						"DS_intron_5P":{"5P":"exonEnd","3P":"downstreamES"},
						"DS_intron_3P":{"5P":"exonEnd","3P":"downstreamES"}}
}

def parse_rMATS(rMATS_path, genome_dict):
	AS_event_type = rMATS_path.split("/")[-1].split(".")[0]
	with open(rMATS_path, "r") as f:
		header_dict = {}
		rMATS_dict = {}
		found_list = []
		temp_dict = {}
		species = ""
		for line in f:
			splitline = line.strip().replace('"',"").split("\t")
			if splitline[0] == "ID":
				for n in range(len(splitline)):
					header_dict[n] = splitline[n]
				continue
			event_id = splitline[0]

			temp_dict[event_id] = {}
			for index in header_dict:
				if header_dict[index] == "chr" and splitline[index][0:3] == "chr":
					temp_dict[event_id][header_dict[index]] = splitline[index][3:]
				else:
					temp_dict[event_id][header_dict[index]] = splitline[index]


			#for exon_part in splicing_dict["SE"]:
			temp_dict[event_id]["alt_exon"] = str(genome_dict[temp_dict[event_id]["chr"]][int(temp_dict[event_id][splicing_dict["SE"]["alt_exon"]["5P"]]):int(temp_dict[event_id][splicing_dict["SE"]["alt_exon"]["3P"]])])
			US_exon_seq = str(genome_dict[temp_dict[event_id]["chr"]][int(temp_dict[event_id][splicing_dict["SE"]["US_intron_5P"]["5P"]]):int(temp_dict[event_id][splicing_dict["SE"]["US_intron_5P"]["3P"]])])
			temp_dict[event_id]["US_intron_5P"] = US_exon_seq[:250]
			temp_dict[event_id]["US_intron_3P"] = US_exon_seq[-250:]
			DS_exon_seq = str(genome_dict[temp_dict[event_id]["chr"]][int(temp_dict[event_id][splicing_dict["SE"]["DS_intron_5P"]["5P"]]):int(temp_dict[event_id][splicing_dict["SE"]["DS_intron_5P"]["3P"]])])
			temp_dict[event_id]["DS_intron_5P"] = DS_exon_seq[:250]
			temp_dict[event_id]["DS_intron_3P"] = DS_exon_seq[-250:]

			if temp_dict[event_id]["strand"] == "-":
				for region in ("alt_exon","US_intron_5P","US_intron_3P",
				"DS_intron_5P","DS_intron_3P"):

					temp_dict[event_id][region] = Seq(temp_dict[event_id][region]).reverse_complement()
					#makes sure the sequence is a string
					temp_dict[event_id][region] = str(temp_dict[event_id][region])

			ensembl_gene = splitline[1]
			#get the average number of counts
			count_string = temp_dict[event_id]["IJC_SAMPLE_1"] + "," + temp_dict[event_id]["IJC_SAMPLE_2"] + "," + temp_dict[event_id]["SJC_SAMPLE_1"] + "," + temp_dict[event_id]["SJC_SAMPLE_2"]
			total = 0
			if "" in count_string.split(","):
				temp_dict[event_id]["avg_count"] = "NA"

			else:
				for count in count_string.split(","):
					total += int(count)
				temp_dict[event_id]["avg_count"] = str(total / len(count_string.split(",")))


			if temp_dict[event_id]["IncLevelDifference"] != "NA":
				temp_dict[event_id]["abs_IncLevelDifference"] = str(abs(float(temp_dict[event_id]["IncLevelDifference"])))
			else:
				temp_dict[event_id]["abs_IncLevelDifference"] = "NA"

			#adds the temp dict to the list of events for a given gene
			if ensembl_gene not in rMATS_dict:
				event_count = 0
				rMATS_dict[ensembl_gene] = {}
			identifier = ensembl_gene + "_" + str(event_count)
			rMATS_dict[ensembl_gene][identifier] = temp_dict[event_id]
			event_count += 1
	return rMATS_dict

if __name__ == "__main__":
	genome_dict1 = Fasta(path_to_fasta)
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	for region in splicing_dict["SE"]:
		with open(f"{output_folder}/{region}_up.fasta", "w+") as up_file, open(f"{output_folder}/{region}_down.fasta", "w+") as down_file, open(f"{output_folder}/{region}_both.fasta", "w+") as both_file, open(f"{output_folder}/{region}_bg.fasta", "w+") as bg_file:
			rMATS_path1 = f"{rMATS_folder}/SE.MATS.JCEC.txt"
			AS_event_type = rMATS_path1.split("/")[-1].split(".")[0]

			rMATS_dict1 = parse_rMATS(rMATS_path1, genome_dict1)
			for gene in rMATS_dict1:
				for event in list(rMATS_dict1[gene].keys()):
					if use_count_cutoff and float(rMATS_dict1[gene][event]["avg_count"]) < count_cutoff:
						if float(rMATS_dict1[gene][event]["FDR"]) < 0.05:
							if abs(float(rMATS_dict1[gene][event]["IncLevelDifference"])) > 0.1:
								if float(rMATS_dict1[gene][event]["IncLevelDifference"]) > 0:
									up_file.write(f">{event}\n{rMATS_dict1[gene][event][region]}\n")
								else:
									down_file.write(f">{event}\n{rMATS_dict1[gene][event][region]}\n")
								both_file.write(f">{event}\n{rMATS_dict1[gene][event][region]}\n")
							else:
								bg_file.write(f">{event}\n{rMATS_dict1[gene][event][region]}\n")
						else:
							bg_file.write(f">{event}\n{rMATS_dict1[gene][event][region]}\n")
					else:
						bg_file.write(f">{event}\n{rMATS_dict1[gene][event][region]}\n")
