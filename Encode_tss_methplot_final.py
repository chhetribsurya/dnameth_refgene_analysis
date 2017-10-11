import numpy as np
import pandas as pd
import pybedtools
import re
import os
import time
from os.path import basename
from os.path import join
from os.path import splitext
import shutil

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


#shutil.rmtree(os.path.expanduser("~/Dropbox/local_miscellaneous_data/tss_peak_methplot/pol2"))
output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/tss_peak_methplot/pol2")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

final_output_file = "final_wgbs_tss_intersect.bed"



def intersect_all(*beds, **kwargs):
	#convert args tuple to a list, so as to have the flexibility, or mutability advantage of modifying the list later on, if needed:
	bedlist = list(beds)
	x = pybedtools.BedTool(bedlist[0])

	for bed in bedlist[1:]:
		x = x.intersect(bed, **kwargs)
	return(x)

#f(x):cleaned_peak_file_1 = intersect_all(input_file_1, input_file_2, input_file_3, input_file_4, wa = True, u = True)


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


def generate_tss_binned_coords(refgene_coordinates_info, upstream_range, downstream_range, bin_size):
	#upstream = 1000
	#downstream = 1000
	#bin_size=100
	refgene_df =  refgene_coordinates_info.sort_values(["chrom","txStart","txEnd"])
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
	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'strand', u'name2', u'tss_midpoint']
	tss_coord_df1 = refgene_model_df.loc[:,select_cols]
	tss_coord_df1.to_csv(join(output_dir, "tss_coordinate_info.bed"), sep="\t", index=False, header=False)

	return(tss_coord_df1)



def load_tss_coord_pybedtool_object(file_name_with_full_path): 
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
	#    print "CpG file already converted to bed format..."

	print "Generating bedtools object..."
	Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))

	return(Cpg_bed_file)


def generate_tss_binned_perc_meth(refseq_final_bedfile, meth_file_list, **kwargs):
	#file_name =  kwargs["files_basename"]
	print "kwargs: ", kwargs
	refseq_bedfile = refseq_final_bedfile

	master_dict = {}

	for idx, each_file in enumerate(meth_file_list):
		file_basename = splitext(basename(each_file))[0]	
		cpg_bed_file = load_cpg_pybedtool_object(each_file)

	    """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
	    print " Processing the intersection for", file_basename
	    pybed_outfile = join(output_dir, (file_basename + "_cpg_tss_intersect.bed"))
	    pybed_outfile_v = join(output_dir, (file_basename + "_cpg_tss_not_intersect.bed"))

	    #if not os.path.exists(pybed_outfile):
	    cpg_bed_intersect = cpg_bed_file.intersect(refseq_bedfile, wa = True, wb = True, output=pybed_outfile)	    
	    print(cpg_bed_intersect.head())  
	       #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
	    #else:
	    #	print "Pybedtools object already present"

	    ### Working with the dataframes; reading the output file of pybedtool intersect:
	    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
	    df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12]]
	    df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "gene_name" ]
	    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
	    df_grouped =  df_ordered.groupby(["bin_start", "bin_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
	    print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

	    if file_basename not in master_dict:
	      master_dict[file_basename] = df_grouped

		print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

	### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
	cpg_intersect_combined_df = pd.concat(master_dict).reset_index()
	cpg_intersect_combined_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

	return(cpg_intersect_combined_df)

#generate_tss_binned_perc_meth(refseq_tss_bed_file, cpg_bed_file, files_basename=file_basename)


def capture_binned_meth_perc_list(meth_file_list, refseq_tss_bed_file):	
	frames = []
	for idx, each_file in enumerate(meth_file_list):
		file_basename = splitext(basename(each_file))[0]	
		cpg_bed_file = load_cpg_pybedtool_object(each_file)
		meth_data = generate_tss_binned_perc_meth(refseq_tss_bed_file, cpg_bed_file, files_basename=file_basename)
		meth_data["file_name"] = file_basename
		frames.append(meth_data)
	appended_data = pd.concat(frames, ignore_index=True)
	appended_data.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = False)

	return(appended_data)


def main():

	refgene_coord_df = final_refgene_model(refgene_file)
	tss_coord_df = generate_tss_binned_coords(refgene_coord_df, 5000, 5000, 100 )

	#join(output_dir, "tss_coordinate_info.bed") : this file is a output of generate_tss_binned_coords(), last function:
	refseq_bed_file = load_tss_coord_pybedtool_object(join(output_dir, "tss_coordinate_info.bed"))

	list_of_files = [cpg_file_2, cpg_file_3, cpg_file_1]
	plot_data_set = generate_tss_binned_perc_meth(refseq_bed_file, list_of_files)

	#refseq_bed_file = None

if __name__ == '__main__':
	main()

else:
	print "Functions Imported from other module"


print "Time for analysis = ", time.time()-start_time


