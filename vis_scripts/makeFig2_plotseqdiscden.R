library("ggplot2")
library("cowplot")

# Args/Output name
args = commandArgs(trailingOnly=TRUE)
if (length(args)>1) {
  stop("Only one argument must be supplied (output file name).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  outname = args[1]
} else {
  print("Using default outputname 'seqdiscplots.eps', you can specify a name as the first argument")
  outname='seqdiscplots.eps'
}


#Function that controls the different aspects of the plot
mydenplotfunc<- function(blastdf){
  return(ggplot(blastdf, aes(x=PID, colour=SAG))+geom_density(alpha=0.5, size=1)+ theme(panel.background = element_blank()) + xlab('Percent Identity') + scale_fill_discrete(name="Metagenome") + ylab('Density') + theme(axis.text = element_text(size=12), legend.text=element_text(size=10), axis.title=element_text(size=15), legend.title=element_text(size=15)) + theme_bw())
}

# Forming each panel (all same using function above)
LD12_ME <- mydenplotfunc(read.delim('LD12/PTXW.len150-vs-LD12_norrna_short.blast.len200.bbh.onlyME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG')))+ggtitle("LD12-Mendota")

FW_ME <- mydenplotfunc(read.delim('FWset/PTXW.len150-vs-FWset_norrna_short.blast.len200.bbh.onlyME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG')))+ggtitle("other-SAGs-Mendota")

acI_ME <- mydenplotfunc(read.delim('acI/PTXW.len150-vs-acI_norrna.fna_short.blast.len200.bbh.onlyME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG')))+ggtitle("acI-Mendota")

LD12_nonME <- mydenplotfunc(read.delim('LD12/PTXW.len150-vs-LD12_norrna_short.blast.len200.bbh.nonME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG')))+ggtitle("LD12-non-Mendota")

acI_nonME <- mydenplotfunc(read.delim('acI/PTXW.len150-vs-acI_norrna.fna_short.blast.len200.bbh.nonME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG')))+ggtitle("acI-non-Mendota")

## Putting plots together and saving them
all_plots <- plot_grid(acI_ME,LD12_ME, FW_ME, acI_nonME, LD12_nonME, labels=c("A", "B", "C", "D", "E"), ncol = 3, nrow = 2)
save_plot(outname, all_plots, base_height = 8, base_width = 18)

