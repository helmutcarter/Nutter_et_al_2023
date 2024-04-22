#!/usr/bin/env python3

#settings
use_dPSI_cutoff = True
use_FDR_cutoff = True
use_count_cutoff = True
count_cutoff = 5
FDR_cutoff = 0.05
dPSI_cutoff = 0.1
rMATS_folder_1 = "/path/to/directory/containing/first/rMATS/output"
rMATS_folder_2 = "/path/to/directory/containing/second/rMATS/output"
event_list = ["SE","RI","MXE","A5SS","A3SS"]




event_window_dict = {"SE":["exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"],
				"A5SS":["longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE"],
				"A3SS":["longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE"],
				"MXE":["1stExonStart_0base","1stExonEnd","2ndExonStart_0base","2ndExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"],
				"RI":["riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"]}

def parse_rMATS(rMATS_path):
	f = open(rMATS_path, "r")
	header_dict = {}
	rMATS_dict = {}
	temp_dict = {}
	for line in f:
		splitline = line.strip().replace('"',"").split("\t")
		if splitline[0] == "ID":
			for n in range(len(splitline)):
				header_dict[n] = splitline[n]
			continue
		event_id = splitline[0]

		temp_dict[event_id] = {}
		for index in header_dict:
			temp_dict[event_id][header_dict[index]] = splitline[index]
			
		if use_FDR_cutoff and float(temp_dict[event_id]["FDR"]) > FDR_cutoff:
			del temp_dict[event_id]
			continue

		count_string = temp_dict[event_id]["IJC_SAMPLE_1"] + "," + temp_dict[event_id]["IJC_SAMPLE_2"] + "," + temp_dict[event_id]["SJC_SAMPLE_1"] + "," + temp_dict[event_id]["SJC_SAMPLE_2"]
		total = 0
		if "" in count_string.split(","):
			temp_dict[event_id]["avg_count"] = "NA"
			print(warning)
		else:
			for count in count_string.split(","):
				total += int(count)
			temp_dict[event_id]["avg_count"] = str(total / len(count_string.split(",")))
			if use_count_cutoff and float(temp_dict[event_id]["avg_count"]) < count_cutoff:
				del temp_dict[event_id]
				continue
		if use_dPSI_cutoff and abs(float(temp_dict[event_id]["IncLevelDifference"])) < dPSI_cutoff:
			del temp_dict[event_id]
			continue

		rMATS_dict[event_id] = temp_dict[event_id]
	return rMATS_dict

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

if __name__ == "__main__":
	total1 = 0
	total2 = 0
	intersection_total = 0
	intersection_list = []
	unique_list1 = []
	unique_list2 = []
	rMATS_list1 = []
	rMATS_list2 = []
	for event in event_list:
		rMATS_dict1 = parse_rMATS(f"{rMATS_folder_1}/{event}.MATS.JCEC.txt")
		rMATS_dict2 = parse_rMATS(f"{rMATS_folder_2}/{event}.MATS.JCEC.txt")
		for id in rMATS_dict1:
			outvalue = rMATS_dict1[id]["geneSymbol"] + "_" + event
			for value in event_window_dict[event]:
				outvalue += "_" + rMATS_dict1[id][value]
			rMATS_list1.append(outvalue)
			print(outvalue + ", " + rMATS_dict1[id]["IncLevelDifference"])

		for id in rMATS_dict2:
			outvalue = rMATS_dict2[id]["geneSymbol"] + "_" + event
			for value in event_window_dict[event]:
				outvalue += "_" + rMATS_dict2[id][value]
			rMATS_list2.append(outvalue)

	rMATS_list1 = list(dict.fromkeys(rMATS_list1))
	rMATS_list2 = list(dict.fromkeys(rMATS_list2))
	total1 = len(rMATS_list1)
	total2 = len(rMATS_list2)
	this_intersection = intersection(rMATS_list1,rMATS_list2)
	intersection_len = len(this_intersection)

	for id in rMATS_list1:
		if id not in this_intersection:
			unique_list1.append(id)


	for id in rMATS_list2:
		if id not in this_intersection:
			unique_list2.append(id)
	print("number of sig events unique to rMATS 1:")
	print(len(unique_list1))
	print("number of sig events unique to rMATS 2:")
	print(len(unique_list2))
	print("number of sig events in rMATS 1:")
	print(total1)
	print("number of sig events in rMATS 2:")
	print(total2)
	print("number of sig events overlap between rMATS 1 and 2:")
	print(intersection_len)

