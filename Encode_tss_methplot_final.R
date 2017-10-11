
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(gtools)

# tf_name = "GATA"

# ### For motif based methylation plots:
# input_file <- "~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/motifs_methplot/final_wgbs_motifs_intersect_2kb_GATA.bed"
# binned_perc_meth_table <- fread(input_file, sep="\t", header= TRUE)

# names(binned_perc_meth_table) <- c("index", "annotation", "bin_start", "bin_end", "meth_percent")
# binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)

# :waiver
# output_file_name <- paste0("~/Dropbox", "/", "motif_composite_plot.pdf")       
# pdf(output_file_name, width=6, height=7)

# this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, meth_percent)) +
# ggplot2::geom_line(aes(color=annotation)) +
# ggplot2::ylab("Percent Methylation") +
# ggplot2::xlab("Distance from Center (bp)") +
# ggplot2::ggtitle(paste(tf_name,"motifs methylation profile"))+
# #ggplot2::ggsave(output_file_name, width = 4, height = 10, dpi = 200)+
# ggplot2::scale_x_continuous() +
# ggplot2::scale_y_continuous(limits=c(0,1)) +
# ggplot2::theme( legend.key = element_blank(),
#       #legend.title=element_text(size=8, face="bold"),
#       #legend.text=element_text(size=5),
#       #legend.key.size = unit(0.2, "cm"),
#       axis.line = element_line(),
#       panel.background = element_blank(),
#       axis.text.x = element_text(color="black", size=12),
#       axis.title.x = element_text(color="black", size=14, vjust=-1.5),
#       axis.text.y = element_text(color="black", size=12),
#       axis.title.y = element_text(color="black", size=14, vjust=3),
#       plot.title=element_text(size=12, face="bold"),
#       plot.margin = grid::unit(c(1,1,1,1), "cm"),
#       panel.border=element_blank(),
#       axis.ticks=element_line(size=0.6, color="black")) + theme_bw()

# this_plot

# dev.off()




### For tss plots:
file <- "~/Dropbox/local_miscellaneous_data/tss_peak_methplot/pol2/final_wgbs_tss_intersect.bed"
binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)

names(binned_perc_meth_table) <- c("index", "annotation", "bin_start", "bin_end", "meth_percent")
binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)


output_file_name <- paste0("~/Dropbox", "/", "tss_composite_plot.pdf")       
pdf(output_file_name)

this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, meth_percent)) +
ggplot2::geom_line(aes(color=annotation)) +
ggplot2::ylab("Percent Methylation") +
ggplot2::xlab("Distance from Center (bp)") +
ggplot2::scale_x_continuous() +
ggplot2::scale_y_continuous(limits=c(0,1)) +
ggplot2::theme( legend.key = element_blank(),
      axis.line = element_line(),
      panel.background = element_blank(),
      axis.text.x = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=14, vjust=-1.5),
      axis.text.y = element_text(color="black", size=12),
      axis.title.y = element_text(color="black", size=14, vjust=3),
      plot.margin = grid::unit(c(1,1,1,1), "cm"),
      panel.border=element_blank(),
      axis.ticks=element_line(size=0.6, color="black")) + theme_bw()

this_plot

dev.off()



### For tss_tes plots:
tss_tes_file <- "/Users/suryachhetri/Dropbox/local_miscellaneous_data/tss_tes_peak_methplot/pol2/final_wgbs_tss_tes_intersect.bed"
binned_perc_meth_table <- fread(tss_tes_file, sep="\t", header= TRUE)

names(binned_perc_meth_table) <- c("index", "annotation", "group", "meth_percent")

output_file_name <- paste0("~/Dropbox", "/", "tss_tes_composite_plot.pdf")       
pdf(output_file_name)

this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(group, meth_percent, group="name")) +
ggplot2::geom_line(aes(color=annotation)) +
ggplot2::ylab("Percent Methylation") +
ggplot2::xlab("Distance from Center (bp)") +
ggplot2::scale_x_discrete(limits=paste0("group_", seq(0,140)), breaks = c("group_1", "group_20","group_120", "group_139"), labels=c("group_1"= "-2kb", "group_20" = "TSS", "group_120" = "TTS", "group_139" = "+2kb")) +
ggplot2::scale_y_continuous(limits=c(0,1))
ggplot2::theme( legend.key = element_blank(),
      axis.line = element_line(),
      panel.background = element_blank(),
      axis.text.x = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=14, vjust=-1.5),
      axis.text.y = element_text(color="black", size=12),
      axis.title.y = element_text(color="black", size=14, vjust=3),
      plot.margin = grid::unit(c(1,1,1,1), "cm"),
      panel.border=element_blank(),
      axis.ticks=element_line(size=0.6, color="black")) + theme_bw()

this_plot

dev.off()



### For peaks plot:
file <- "~/Dropbox/local_miscellaneous_data/pol2/peaks_methplot/final_wgbs_peaks_intersect.bed"
binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)

names(binned_perc_meth_table) <- c("index", "annotation", "bin_start", "bin_end", "meth_percent")
binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)


output_file_name <- paste0("~/Dropbox", "/", "peaks_composite_plot.pdf")       
pdf(output_file_name)

this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, meth_percent)) +
ggplot2::geom_line(aes(color=annotation)) +
ggplot2::ylab("Percent Methylation") +
ggplot2::xlab("Distance from Center (bp)") +
ggplot2::scale_x_continuous() +
ggplot2::scale_y_continuous(limits=c(0,1)) +
ggplot2::theme( legend.key = element_blank(),
      axis.line = element_line(),
      panel.background = element_blank(),
      axis.text.x = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=14, vjust=-1.5),
      axis.text.y = element_text(color="black", size=12),
      axis.title.y = element_text(color="black", size=14, vjust=3),
      plot.margin = grid::unit(c(1,1,1,1), "cm"),
      panel.border=element_blank(),
      axis.ticks=element_line(size=0.6, color="black")) + theme_bw()

this_plot

dev.off()



### For chromHMM plot:
file <- "~/Dropbox/local_miscellaneous_data/pol2/chromHMM_methplot/final_wgbs_chromHMM_intersect.bed"
binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)
bin_size <-  100

names(binned_perc_meth_table) <- c("index", "annotation", "bin_start", "bin_end", "state", "meth_percent")
binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)

### Maintain the order of the list as per the targetl list:
ordered_df <- binned_perc_meth_table[mixedorder(binned_perc_meth_table$state),]
binned_perc_meth_table$state <- factor(binned_perc_meth_table$state, levels= ordered_df$state)
factor(binned_perc_meth_table$state)

bin_size <-  100
min_size <- min(binned_perc_meth_table$bin_mid) #original upstream
max_size <- max(binned_perc_meth_table$bin_mid) #original downstream

label_min_size <- min(binned_perc_meth_table$bin_mid) + (-bin_size/2) #Just to label upstream (but not original upstream)
label_max_size <- max(binned_perc_meth_table$bin_mid) + (bin_size/2) #Just to label downstream (but not original downstream)

# second_max_size <- sort(unique(binned_perc_meth_table$bin_mid))[2]
# second_min_size <- sort(unique(binned_perc_meth_table$bin_mid), decreasing=TRUE)[2]

output_file_name <- paste0("~/Dropbox", "/", "chromHMM_composite_plot.pdf")       
pdf(output_file_name)

ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, meth_percent)) +
ggplot2::geom_line(aes(color=annotation)) +
ggplot2::ylab("Percent Methylation") +
ggplot2::xlab("Distance from Center (bp)") +
ggplot2::scale_x_continuous(breaks=c(min_size,0,max_size), labels=c(label_min_size/1000, "0", label_max_size/1000)) +
#ggplot2::scale_x_continuous(breaks=c(second_min_size,0,second_max_size), labels=c(paste0(second_min_size/1000), "0", paste0(second_max_size/1000))) +
ggplot2::scale_y_continuous(limits=c(0,1)) +
ggplot2::theme( legend.key = element_blank(),
      axis.line = element_line(),
      panel.background = element_blank(),
      axis.text.x = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=14, vjust=-1.5),
      axis.text.y = element_text(color="black", size=12),
      axis.title.y = element_text(color="black", size=14, vjust=3),
      plot.margin = grid::unit(c(1,1,1,1), "cm"),
      panel.border=element_blank(),
      axis.ticks=element_line(size=0.6, color="black")) + theme_bw() +
ggplot2::facet_wrap(~state) 

this_plot

dev.off()










plot_percent_meth <- function(file){

	binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)
	names(binned_perc_meth_table) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")
	binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
    this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, Percent_meth)) +
    ggplot2::geom_line(aes(color=annotation)) +
    ggplot2::ylab("Percent Methylation") +
    ggplot2::xlab("Distance from Center (bp)") +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) +
    ggplot2::theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          # legend.position="none",
          legend.key = element_blank(),
          axis.line = element_line(),
          panel.background = element_blank(),
          axis.text.x = element_text(color="black", size=12),
          axis.title.x = element_text(color="black", size=14, vjust=-1.5),
          axis.text.y = element_text(color="black", size=12),
          axis.title.y = element_text(color="black", size=14, vjust=3),
          plot.margin = grid::unit(c(1,1,1,1), "cm"),
          panel.border=element_blank(),
          axis.ticks=element_line(size=0.6, color="black"))

  this_plot %>% return
}


meth_plot <- plot_percent_meth(file)
meth_plot.save(dir + "")


#######################
#######################
#######################



read_meth_data <- fread("/home/surya/Desktop/scripts/data/groupby_peaks_final/group_data_motifs_concating.txt", sep="\t", header=TRUE)
read_meth_data %>% head
names(read_meth_data) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")

theme_1 <-theme(panel.grid.minor = element_blank(),
			panel.grid.major = element_blank(),
			# legend.position="none",
			legend.key = element_blank(),
			axis.line = element_line(),
			panel.background = element_blank(),
			axis.text.x = element_text(color="black", size=12),
			axis.title.x = element_text(color="black", size=14, vjust=-1.5),
			axis.text.y = element_text(color="black", size=12),
			axis.title.y = element_text(color="black", size=14, vjust=3),
			plot.margin = grid::unit(c(1,1,1,1), "cm"),
			panel.border=element_blank(),
			axis.ticks=element_line(size=0.6, color="black"))


read_meth_data$bin_mid <- (read_meth_data$bin_start + read_meth_data$bin_end)/2
meth_plot <- ggplot2::ggplot(read_meth_data, aes(x=bin_mid, y=Percent_meth)) + 
			ggplot2::geom_line(aes(color=annotation)) +
    		ggplot2::ylab("Percent Methylation") +
    		ggplot2::xlab("Distance from Center (bp)") +
    		ggplot2::scale_x_continuous(expand=c(0,0)) +
    		ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1))
    		
meth_plot


plot_percent_meth_with_depth <- function(binned_perc_meth_table){
  binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table,
                                         (bin_start + bin_end)/2)
  merged_df <- melt(binned_perc_meth_table, id.vars = c("bin_start", "bin_end", "meth", "unmeth", "group", "bin_mid"))
  this_plot <- ggplot2::ggplot(merged_df) +
    ggplot2::geom_line(aes(x = bin_mid, y = value, color = annotation)) +
    ggplot2::geom_line(aes(x = bin_mid, y = value, color = annotation)) +
    ggplot2::facet_grid(variable~., scales = "free_y")

  this_plot %>% return
}



#Loops can come in handy on numerous occasions. While loops are like repeated if statements; 
#the for loop is designed to iterate over all elements in a sequence.

#Also, whenever you're using a for loop, you might want to revise your code and see whether you can use the 
#lapply function instead. Learn all about this intuitive way of applying a function over a list or a 
#vector, and its variants sapply and vapply.


file1 <- "/home/surya/Desktop/scripts/data/groupby_peaks_final/group_data_motifs_concating.txt"
file_list <- c(file1)

plot_percent_meth <- function(file_list){

	for (file in file_list){
		#or use lapply
		print(file)
		binned_perc_meth_table <- fread(file, sep="\t", header= TRUE)
		names(binned_perc_meth_table) <- c("bin_start", "bin_end", "meth", "unmeth", "Percent_meth", "annotation")
		print(names(binned_perc_meth_table))
		binned_perc_meth_table$bin_mid <- with(binned_perc_meth_table, (bin_start + bin_end)/2)
	    this_plot <- ggplot2::ggplot(binned_perc_meth_table, aes(bin_mid, Percent_meth)) +
	    ggplot2::geom_line(aes(color=annotation)) +
	    ggplot2::ylab("Percent Methylation") +
	    ggplot2::xlab("Distance from Center (bp)") +
	    ggplot2::scale_x_continuous(expand=c(0,0)) +
	    ggplot2::scale_y_continuous(expand=c(0,0.0001), limits = c(0,1)) 
	    print("cool")

	    #this_plot.save(dir + "")
		return(this_plot)

  	}
}


meth_plot <- plot_percent_meth(file_list)
meth_plot
meth_plot.save(dir + "")


#Testing the loop:
for (num in c(1,2,3)){
	print(num)
}