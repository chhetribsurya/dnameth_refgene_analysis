import numpy as np
import pandas as pd
import pybedtools
import re
import os
import time
import shutil
from os.path import basename
from os.path import join
from os.path import splitext

start_time = time.time()


input_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed")
input_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/MACS_SL68321_SL68325_15_peaks.bed")
input_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/MACS_SL60583_SL60585_10_peaks.bed")
input_file_4 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/MACS_SL63653_SL63655_10_peaks.bed")

#refgene_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_head.txt")
refgene_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_hg19")
#fasta_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/hg19-male.fa")
#cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/wgbs_GATA2.cov")
#cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL60583_GATA2_rep1.cov")
#cpg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/gata/SL63653_GATA2_rep2.cov")

cpg_file_1 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/wgbs_Pol2.cov")
cpg_file_2 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL68321_Pol2_rep1.cov")
cpg_file_3 = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/meth_chip_bs_data/pol2/SL88935_Pol2_rep2.cov")


#shutil.rmtree(os.path.expanduser("~/Dropbox/local_miscellaneous_data/tss_tes_peak_methplot/pol2"))
output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/tss_tes_peak_methplot/pol2")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

final_output_file = "final_wgbs_tss_tes_intersect.bed"



def intersect_all(*beds, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	bedlist = list(beds)
	x = pybedtools.BedTool(bedlist[0])

	for bed in bedlist[1:]:
		x = x.intersect(bed, **kwargs)
	return(x)


#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)
#cleaned_peak_file_1.saveas('/home/surya/surya/for_gata_parse/pybed_overlap.bed')



#Merging all the cpg files with the common key (chrom, start, end):
def merge_all(*cpgs, **kwargs):

	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	cpg_file_list = list(cpgs)
	x = pd.read_csv(cpg_file_list[0], sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])

	for cpg_file in cpg_file_list[1:]:
		current_file = pd.read_csv(cpg_file, sep="\t", names=["chrom","start","end", "perc_meth", "meth", "unmeth"])
		x = pd.merge(x, current_file, **kwargs)
	return(x)


#f(x):merged_df = merge_all(cpg_file_1, cpg_file_2, cpg_file_3, how='inner', on=['chrom', 'start', 'end'] )
#merged_df.to_csv(output_file, sep="\t",header=True, index= False)



#Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
#gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
def chrom_analysis_range(range_of_chrom_list=22, chrom_XY=False):
	chrom_list = []
	for chrom_no in range(range_of_chrom_list):
		chrom_list.append("chr" + str(chrom_no + 1))
	if chrom_XY:
		chrom_list.append("chrX")
		chrom_list.append("chrY")

	return(chrom_list)



def single_chrom_based_analysis(which_chrom):
	chrom_list = ["chr" + str(which_chrom)]	
	return(chrom_list)



def parse_refgene_model(refgene_file):

	with open(refgene_file, "r") as file:
		sorted_refgene_list = []

		ref_file_read = file.readlines()
		header_splitted = ref_file_read[0].split()
		header_tuple = (header_splitted[2], header_splitted[4], header_splitted[5], header_splitted[3],header_splitted[12])
		header = "\t".join(header_tuple)
		#single_chrom_analysis = single_chrom_based_analysis(19)
		#autosomal_chrom_list = chrom_analysis_range()
		normal_chrom_list = chrom_analysis_range(22, True)

		for line in ref_file_read[1:]:
			splitted = line.strip().split()
			#if str(splitted[2]) in single_chrom_analysis:
			if str(splitted[2]) in normal_chrom_list:
				if splitted[3] == "+":
					sorted_refgene_list.append([str(splitted[2]), int(splitted[4]), int(splitted[5]), str(splitted[3]),str(splitted[12])])
				if splitted[3] == "-":
					sorted_refgene_list.append([str(splitted[2]), int(splitted[5]), int(splitted[4]), str(splitted[3]), str(splitted[12])])

		return(sorted_refgene_list,header)



def final_refgene_model(refgene_file):

	sorted_refgene_list, header = parse_refgene_model(refgene_file)
	df_refgene = pd.DataFrame(sorted_refgene_list, columns = header.split())
	grouped_refgene = df_refgene.groupby(["chrom","strand","name2"], sort = False, as_index=False)
	refgene_model = grouped_refgene.agg({"txStart": np.min,"txEnd": np.max})	
 	colsorted_final_model = refgene_model.reindex(columns=header.split())
 	print "\nCurrent dimension of the refgene model:\n", colsorted_final_model.shape
 	colsorted_final_model = colsorted_final_model.drop_duplicates()
 	print "Dropping duplicates, current dimension of the refgene model:\n", colsorted_final_model.shape
	colsorted_final_model.to_csv(join(output_dir, "filtered_refgene_model_head.txt"), sep="\t", header = True, index= False)

	return(colsorted_final_model)



# Dividing the gene or peak or binding sites into the quantiles:
#( set groups = 4 for quartiles, groups = 10 for deciles, 
# & groups = 100 for percentiles)
def generate_gene_quantiles_coords(total_groups, refgene_coord_info):
	#groups = 4
	#refgene_coord_info = refgene_coord_df
	groups = total_groups
	refgene_concat_df = pd.concat([refgene_coord_info]*groups, ignore_index=True)
	refgene_sorted_df = refgene_concat_df.sort_values(["chrom","strand","name2"])
	
	### Sort refgene_sorted_df and orginal_refgene_df with same values to make the order of gene same:
	orginal_refgene_df = refgene_coord_info.sort_values(["chrom","strand","name2"]) 
	start_list = list(orginal_refgene_df["txStart"])
	end_list = list(orginal_refgene_df["txEnd"])

	zipped_tuple = zip(start_list, end_list)

	bins_dict = {"chrom_start": [], "chrom_end": []}
	#extra_bins_dict = {"bin_start": [], "bin_end":[]}
	for each_start,each_end in zip(start_list, end_list):
	    #print np.linspace(each_start, each_end, (groups+1))
	    list1 = list(np.linspace(each_start, each_end, (groups+1)).round().astype(int)[:-1])
	    list2 = list(np.linspace(each_start, each_end, (groups+1)).round().astype(int)[1:])
	    bins_dict["chrom_start"].extend(list1)
	    bins_dict["chrom_end"].extend(list2)

	### For negative strands, make sure that start strand is lesser than the end strand
	### so, switch the positions else pybedtools would get upset:
	bining_df = pd.DataFrame(bins_dict).iloc[:,[1,0]]
	df_pos = bining_df[bining_df["chrom_start"] <= bining_df["chrom_end"]]
	df_neg = bining_df[bining_df["chrom_start"] > bining_df["chrom_end"]]
	df_neg_switched = df_neg.rename(columns={"chrom_start" : "chrom_end", "chrom_end" : "chrom_start"}).iloc[:,[1,0]]	
	final_bining_df = pd.concat([df_pos, df_neg_switched]).sort_index()

	group_dict = {"group":[]}
	for i in range(len(start_list)):
		for j in range(groups):
			group_dict["group"].append("group_" + str(j+1))

	group_df = pd.DataFrame(group_dict)
	coord_df = pd.concat([refgene_sorted_df.reset_index(), final_bining_df, group_df], axis = 1)
	
	quantiles_coord_df1 = coord_df.drop(["index"], axis=1)
	quantiles_coord_df1.to_csv(join(output_dir, "tss_tes_gene_quantile_data.bed"), sep="\t", index=False, header=True)
	
	return(quantiles_coord_df1)

#quantiles_coord_df = generate_gene_quatiles_coords(4, refgene_coord_info)


##################################
##################################
##################################


def generate_tss_tes_binned_coords(refgene_coordinates_info, upstream_range, downstream_range, bin_size, genequantiles_coordinates_info):
	#upstream = 1000
	#downstream = 1000
	#bin_size = 100
	refgene_df =  refgene_coordinates_info.sort_values(["chrom","txStart","txEnd"])
	quantiles_coord_df = genequantiles_coordinates_info
	upstream = upstream_range
	downstream = downstream_range
	bin_size = bin_size
	nrows =  refgene_df.shape[0]

	bins = range(-upstream, (downstream), bin_size)
	bin_len = len(bins)
	refgene_concat_df = pd.concat([refgene_df]*bin_len, ignore_index="TRUE")
	refgene_sorted_df = refgene_concat_df.sort_values(["chrom","txStart","txEnd"])

	### Copy the bin list that is deep copy:
	bin_start_list = bins[:]
	bin_end_list = []
	for each in bin_start_list:
		bin_end_list.append(each+bin_size)

	bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
	bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
	bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

	### Combine the refgene df and bin df by cbind or column wise:
	temp_refgene_df = pd.concat([refgene_sorted_df.reset_index(), bin_concat_df], axis = 1)
	temp_refgene_df["tss_midpoint"] = temp_refgene_df["txStart"]
	temp_refgene_df["tes_midpoint"] = temp_refgene_df["txEnd"]
	

	# TSS bining operation df:
	# Select all the rows with upstream properties i.e negative dist:
	tss_upstream_df = temp_refgene_df[temp_refgene_df["bin_start"] < 0]
	#final_refgene_df = temp_refgene_df.loc[:,["chrom", "tss_midpoint", "bin_start","bin_end","name2","strand"]]

	""" For positive strand, i.e if strand1 == "+": 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end)"""

	refgene_pos_df = tss_upstream_df.loc[tss_upstream_df["strand"] == "+",]
	refgene_pos_df["chrom_start"] = refgene_pos_df["tss_midpoint"] + refgene_pos_df["bin_start"]
	refgene_pos_df["chrom_end"] = refgene_pos_df["tss_midpoint"] + refgene_pos_df["bin_end"]
	#refgene_pos_df.iloc[:,range(1,8)]

	""" For negative strand, i.e if strand1 == "-": 
	chrom_start = tss_midpt - (bin_end ); 
	chrom_end = tss_midpt - (bin_start);
	Meth_r concept, start and end switched, so as to always maintain the higher coords for end site """

	refgene_neg_df = tss_upstream_df.loc[tss_upstream_df["strand"] == "-",]
	refgene_neg_df["chrom_start"] = refgene_neg_df["tss_midpoint"] - refgene_neg_df["bin_end"]
	refgene_neg_df["chrom_end"] = refgene_neg_df["tss_midpoint"] - refgene_neg_df["bin_start"]
	#refgene_pos_df.iloc[:,range(1,8)]

	### Combine the positive and negative stranded genes:
	refgene_model_df = pd.concat([refgene_pos_df, refgene_neg_df])

	# TES bining operation df:
	# Select all the rows with upstream properties i.e negative dist:
	tes_downstream_df = temp_refgene_df[temp_refgene_df["bin_start"] >= 0]
	#final_refgene_df = temp_refgene_df.loc[:,["chrom", "tss_midpoint", "bin_start","bin_end","name2","strand"]]

	""" For positive strand, i.e if strand1 == "+": 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end)"""

	refgene_pos_df_1 = tes_downstream_df.loc[tes_downstream_df["strand"] == "+",]
	refgene_pos_df_1["chrom_start"] = refgene_pos_df_1["tes_midpoint"] + refgene_pos_df_1["bin_start"]
	refgene_pos_df_1["chrom_end"] = refgene_pos_df_1["tes_midpoint"] + refgene_pos_df_1["bin_end"]
	#refgene_pos_df_1.iloc[:,range(1,8)]

	""" For negative strand, i.e if strand1 == "-": 
	chrom_start = tss_midpt - (bin_end ); 
	chrom_end = tss_midpt - (bin_start);
	Meth_r concept, start and end switched, so as to always maintain the higher coords for end site """

	refgene_neg_df_1 = tes_downstream_df.loc[tes_downstream_df["strand"] == "-",]
	refgene_neg_df_1["chrom_start"] = refgene_neg_df_1["tes_midpoint"] - refgene_neg_df_1["bin_end"]
	refgene_neg_df_1["chrom_end"] = refgene_neg_df_1["tes_midpoint"] - refgene_neg_df_1["bin_start"]

	### Combine the positive and negative stranded genes:
	refgene_model_df_1 = pd.concat([refgene_pos_df_1, refgene_neg_df_1])
	
	concat_df =  pd.concat([refgene_model_df, refgene_model_df_1], ignore_index=True)
	df_sorted  = concat_df.sort_values(["name2", "bin_start"])
	#df_sorted[df_sorted["strand"] == "-"]
	refgene_select_df = df_sorted.loc[:, [u'chrom', u'txStart', u'txEnd', u'strand', u'name2', u'chrom_start',
       u'chrom_end', u'bin_start', 'bin_end']]

	### The dataframe from breaking the genes in to quantiles:
	### Import the gene quantiles data from generate_gene_quantiles_coords() funtion:
	if "group" in quantiles_coord_df.keys():
		quantiles_coord_df = quantiles_coord_df.drop(["group"],axis=1)

	quantiles_coord_df["bin_start"] = 0
	quantiles_coord_df["bin_end"] = 0
	if "index" not in quantiles_coord_df.keys():
		quantiles_coord_df = quantiles_coord_df.reset_index()

	if "index" not in refgene_select_df.keys():	
		refgene_select_df = refgene_select_df.reset_index()

	final_refgene_df = pd.concat([quantiles_coord_df, refgene_select_df])
	refgene_sorted = final_refgene_df.sort_values(["name2", "bin_start"])

	### Extract the sample gene for counting the no. of repetition, and multiply later with group df:
	gene_name_list = refgene_sorted["name2"].unique()	
	sample_gene = refgene_sorted[refgene_sorted["name2"] == gene_name_list[0]]
	group_len = sample_gene.shape[0]
	total_gene_num = len(list(refgene_df["name2"]))

	group_dict = {"group":[]}
	for i in range(group_len):
		group_dict["group"].append("group_" + str(i+1))
	group_df = pd.DataFrame(group_dict)
	final_group_df = pd.concat([group_df]*total_gene_num, ignore_index=True)

	combined_df = pd.concat([refgene_sorted.reset_index(), final_group_df], axis=1)
	### tss_tes_coord_df1.columns ---> [u'index', u'level_0', u'chrom', u'txStart', u'txEnd', u'strand', u'name2', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'group']
	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'strand', u'group', u'name2', u'txStart', u'txEnd', u'index']
	tss_tes_coord_df1 = combined_df.loc[:,select_cols]
	tss_tes_coord_df1.to_csv(join(output_dir, "tss_tes_coordinate_info.bed"), sep="\t", index=False, header=False)

	return(tss_tes_coord_df1)


#tss_tes_coord_df = generate_tss_tes_binned_coords(refgene_coord_df, 1000, 1000, 100, quantiles_coord_df)


def generate_tss_binned_coords(coordinates_info, upstream_range, downstream_range, bin_size):
	#upstream = -1000
	#downstream = 1000
	#bin_size=100
	refgene_df =  coordinates_info.sort_values(["chrom","txStart","txEnd"])
	upstream = upstream_range
	downstream = downstream_range
	bin_size = bin_size
	nrows =  refgene_df.shape[0]

	bins = range(-upstream, (downstream), bin_size)
	bin_len = len(bins)
	refgene_concat_df = pd.concat([refgene_df]*bin_len, ignore_index="TRUE")
	refgene_sorted_df = refgene_concat_df.sort_values(["chrom","txStart","txEnd"])

	### Copy the bin list that is deep copy:
	bin_start_list = bins[:]
	bin_end_list = []
	for each in bin_start_list:
		bin_end_list.append(each+bin_size)

	bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
	bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
	bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

	### Combine the refgene df and bin df by cbind or column wise:
	temp_refgene_df = pd.concat([refgene_sorted_df.reset_index(), bin_concat_df], axis = 1)
	temp_refgene_df["tss_midpoint"] = temp_refgene_df["txStart"]
	final_refgene_df = temp_refgene_df.loc[:,["chrom", "tss_midpoint", "bin_start","bin_end","name2","strand"]]

	""" For positive strand, i.e if strand1 == "+": 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end)"""

	refgene_pos_df = final_refgene_df.loc[final_refgene_df["strand"] == "+",]
	refgene_pos_df["chrom_start"] = refgene_pos_df["tss_midpoint"] + refgene_pos_df["bin_start"]
	refgene_pos_df["chrom_end"] = refgene_pos_df["tss_midpoint"] + refgene_pos_df["bin_end"]
	#refgene_pos_df.iloc[:,range(1,8)]

	""" For negative strand, i.e if strand1 == "-": 
	chrom_start = tss_midpt - (bin_end ); 
	chrom_end = tss_midpt - (bin_start);
	Meth_r concept, start and end switched, so as to always maintain the higher coords for end site """

	refgene_neg_df = final_refgene_df.loc[final_refgene_df["strand"] == "-",]
	refgene_neg_df["chrom_start"] = refgene_neg_df["tss_midpoint"] - refgene_neg_df["bin_end"]
	refgene_neg_df["chrom_end"] = refgene_neg_df["tss_midpoint"] - refgene_neg_df["bin_start"]
	#refgene_pos_df.iloc[:,range(1,8)]

	### Combine the positive and negative stranded genes:
	refgene_model_df = pd.concat([refgene_pos_df, refgene_neg_df])
	refgene_model_df.to_csv(join(output_dir, "tss_coordinate_info.bed"), sep="\t", index=False, header=False)

	return(refgene_model_df)


# def load_pybedObject_for_refseq_model_df(final_dict):
# 	refseq_file <- pybedtools.BedTool((output_dir, "tss_coordinate_out_test.bed"))
# 	return(refseq_bed_file)


def load_tss_tes_coord_pybedtool_object(file_name_with_full_path): 
	each_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file

    #if not os.path.exists(join(output_dir, sorted_peak_file)):
    CMD = 'sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file'
    os.system(CMD)   
  
    print "Generating bedtools object..."
    peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))
 
    return(peak_bed)


def load_cpg_pybedtool_object(file_name_with_full_path):
	print " \nProcessing Cpg bed file\n "
	cpg_bed_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["cpg_bed_file"] = cpg_bed_file
	cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_bedEdited" + splitext(cpg_bed_file)[1]
	os.environ["cpg_bed_edited"] = cpg_bed_edited

	#if not os.path.exists(join(output_dir,cpg_bed_edited)):
	CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n | grep -v "^*" > $output_dir/$cpg_bed_edited'''
	os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
	#else:
	 #   print "CpG file already converted to bed format..."

	print "Generating bedtools object..."
	Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))

	return(Cpg_bed_file)


def generate_tss_tes_binned_perc_meth(refseq_final_bedfile, meth_file_list, **kwargs):
	#file_name =  kwargs["files_basename"]
	print "kwargs: ", kwargs
	refseq_bedfile = refseq_final_bedfile

	master_dict = {}

	for idx, each_file in enumerate(meth_file_list):
		file_basename = splitext(basename(each_file))[0]	
		cpg_bed_file = load_cpg_pybedtool_object(each_file)

	    """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
	    print " Processing the intersection for", file_basename
	    pybed_outfile = join(output_dir, (file_basename + "_cpg_tss_tes_intersect.bed"))
	    pybed_outfile_v = join(output_dir, (file_basename + "_cpg_tss_tes_peaks_not_intersect.bed"))

	    #if not os.path.exists(pybed_outfile):
	    cpg_bed_intersect = cpg_bed_file.intersect(refseq_bedfile, wa = True, wb = True, output=pybed_outfile)	    
	    print(cpg_bed_intersect.head())  
	       #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
	    #else:
	    #	print "Pybedtools object already present"

	    ### Working with the dataframes; reading the output file of pybedtool intersect:
	    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
	    df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 13]]
	    #print df_ordered.head()
	    df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "group", "gene_name" ]
	    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
	    df_grouped =  df_ordered.groupby(["group"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
	    print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

	    if file_basename not in master_dict:
	      master_dict[file_basename] = df_grouped

		print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

	### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
	cpg_intersect_combined_df = pd.concat(master_dict).reset_index()
	cpg_intersect_combined_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

	return(cpg_intersect_combined_df)

###generate_tss_tes_binned_perc_meth(refseq_tss_bed_file, cpg_bed_file, files_basename=file_basename)


### Pass on the quantile input file generated by generate_gene_quantiles_coords():
def rearrange_quantiles_coord_for_pybedtool_load():
	### Arrange the quantiles into bedtool format for pybedtools load:
	quantiles_df = pd.read_csv(join(output_dir, "tss_tes_gene_quantile_data.bed"), sep="\t")
	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'group', u'name2', u'strand']
	quantiles_df = quantiles_df.loc[:,select_cols]
	quantiles_df.to_csv(join(output_dir, "tss_tes_gene_quantile_data_colRearranged.bed"), sep="\t", index=False, header=False)
	
	return(quantiles_df)


def load_quantiles_coord_pybedtool_object(file_name_with_full_path): 
	each_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file

    if not os.path.exists(join(output_dir, sorted_peak_file)):
    	CMD = 'sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file'
    	os.system(CMD)   
  
    print "Generating bedtools object..."
    peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))
 
    return(peak_bed)


### Pass on the quantile bedfile obtained from prior generate_gene_quantiles_coords():
### Usage: quantiles_coord_df = generate_gene_quantiles_coords(4, refgene_coord_info)

def generate_quantiles_binned_perc_meth(quantiles_final_bedfile, meth_file_list, **kwargs):
	#file_name =  kwargs["files_basename"]
	print "kwargs: ", kwargs
	quantiles_bedfile = quantiles_final_bedfile

	master_dict = {}

	for idx, each_file in enumerate(meth_file_list):
		file_basename = splitext(basename(each_file))[0]	
		cpg_bed_file = load_cpg_pybedtool_object(each_file)

	    """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
	    print " Processing the intersection for", file_basename
	    pybed_outfile = join(output_dir, (file_basename + "_cpg_quantiles_intersect.bed"))
	    pybed_outfile_v = join(output_dir, (file_basename + "_cpg_quantiles_peaks_not_intersect.bed"))

	    #if not os.path.exists(pybed_outfile):
	    cpg_bed_intersect = cpg_bed_file.intersect(quantiles_bedfile, wa = True, wb = True, output=pybed_outfile)	    
	    print(cpg_bed_intersect.head())  
	       #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
	    #else:
	    #	print "Pybedtools object already present"

	    ### Working with the dataframes; reading the output file of pybedtool intersect:
	    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
	    df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11]]
	    #print df_ordered.head()
	    df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "group", "gene_name", "strand"]
	    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
	    df_grouped =  df_ordered.groupby(["group"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
	    print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

	    if file_basename not in master_dict:
	      master_dict[file_basename] = df_grouped

		print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

	### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
	cpg_intersect_combined_df = pd.concat(master_dict).reset_index()
	cpg_intersect_combined_df.to_csv(join(output_dir, "final_wgbs_quantiles_intersect.bed"), sep ="\t", header = True, index = True)

	return(cpg_intersect_combined_df)

#generate_tss_tes_binned_perc_meth(refseq_tss_bed_file, cpg_bed_file, files_basename=file_basename)


def capture_binned_meth_perc_list(meth_file_list, refseq_tss_bed_file):	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename1 = os.path.splitext(basename(each_file))[0]	
		#cpg_bed_file = load_stringf_for_cpg_after_edits(each_file)
		cpg_bed_file = load_pybedObject_for_cpg_after_edits(each_file)
		meth_data = generate_tss_binned_perc_meth(refseq_tss_bed_file, cpg_bed_file, file_basename=file_basename1)
		meth_data["file_name"] = os.path.splitext(basename(each_file))[0]
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	appended_data.to_csv(output_dir + "/" + ggplot_outfile_suffix, sep ="\t", header = True, index = False)

	return(appended_data)


def main():

	refgene_coord_df = final_refgene_model(refgene_file)
	quantiles_coord_df = generate_gene_quantiles_coords(100, refgene_coord_df)
	#quantiles_coord_df[quantiles_coord_df["chrom_start"].isnull()]
	tss_tes_coord_df = generate_tss_tes_binned_coords(refgene_coord_df, 2000, 2000, 100, quantiles_coord_df)

	#join(output_dir, "tss_tes_coordinate_info.bed") : this file is a output of generate_tss_tes_binned_coords(), last function:
	refseq_bed_file = load_tss_tes_coord_pybedtool_object(join(output_dir, "tss_tes_coordinate_info.bed"))

	list_of_files = [cpg_file_2, cpg_file_3, cpg_file_1]
	plot_data_set = generate_tss_tes_binned_perc_meth(refseq_bed_file, list_of_files)

### Only run the following functions, if you need quantiles methylation:
################################################################################

	### Rearrange the quantiles generated by:
	### generate_gene_quantiles_coords(4, refgene_coord_df):	
	final_quantiles_coord_df = rearrange_quantiles_coord_for_pybedtool_load()
    
    ### join(output_dir, "tss_tes_gene_quantile_data_colRearranged.bed") : this file is a output of 
    ### generate_gene_quantiles_coords() & rearrange_quantiles_coord_for_pybedtool_load() function:
	quantiles_bed_file = load_quantiles_coord_pybedtool_object(join(output_dir, "tss_tes_gene_quantile_data_colRearranged.bed"))
	
	list_of_files = [cpg_file_2, cpg_file_3, cpg_file_1]
	plot_quantiles_data_set = generate_quantiles_binned_perc_meth(quantiles_bed_file, list_of_files)

################################################################################


if __name__ == '__main__':
	main()

else:
	print "Functions Imported from other module"


print "Time for analysis = ", time.time()-start_time


##########################


