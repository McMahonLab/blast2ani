library("ggplot2")
library("cowplot")
library("dplyr")
library("stringr")

# Args/Output name
args = commandArgs(trailingOnly=TRUE)
if (length(args)>1) {
  stop("Only one argument must be supplied (output file name).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  outname = args[1]
} else {
  print("Using default outputname 'DATE_Fig3_seqdiscplots' both eps and pdf, you can specify a name as the first argument")
  outname=paste0(Sys.Date(),'_Fig3_seqdiscplots')
}

# Metadata for colors
metadata <- read.delim("../Table-1-allSAGs-Metadata.txt")
metadata$Color <- as.character(metadata$Color)

# Reformatting Metadata
metadata$Genome.name <- gsub("-", "", paste(metadata$Genome.name))
metadata$SAG_lake <- paste(str_sub(metadata$Genome.name,start = -3, end=-1),metadata$Lake,sep="-")
meta_sub <- metadata %>%
  select(Genome.name, Lake, Tribe, Color, SAG_lake)
## Getting colors
jColors <- metadata$Color
names(jColors) <- metadata$SAG_lake


#Function that controls the different aspects of the plot
mydenplotfunc<- function(blastdf){
  return(blastdf  %>% 
           ggplot(aes(x=PID, colour=SAG_lake))+geom_density(alpha=0.5, size=1)+
           scale_color_manual(values=jColors) + 
           theme_bw() +
           xlab('Percent Identity') + 
           ylab('Density') + 
           theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), legend.text= element_text(size=15),axis.text = element_text(size=15),axis.text.x  = element_text(angle=90, vjust=0.5)))
}

# Forming each panel (all same using function above)
LD12_ME_df <- read.delim('LD12/PTXW.len150-vs-LD12_norrna_short.blast.len200.bbh.onlyME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG'))
LD12_ME_df <- left_join(LD12_ME_df,meta_sub,by = c("SAG"="Genome.name"))
LD12_ME <- mydenplotfunc(LD12_ME_df)+ggtitle("LD12-Mendota")

FW_ME_df <- read.delim('FWset/PTXW.len150-vs-FWset_norrna_short.blast.len200.bbh.onlyME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG'))
FW_ME_df <- left_join(FW_ME_df,meta_sub,by = c("SAG"="Genome.name"))
FW_ME <- mydenplotfunc(FW_ME_df)+ggtitle("other-SAGs-Mendota")

acI_ME_df <- read.delim('acI/PTXW.len150-vs-acI_norrna.fna_short.blast.len200.bbh.onlyME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG'))
acI_ME_df <- left_join(acI_ME_df,meta_sub,by = c("SAG"="Genome.name"))
acI_ME <- mydenplotfunc(acI_ME_df)+ggtitle("acI-Mendota")

LD12_nonME_df <- read.delim('LD12/PTXW.len150-vs-LD12_norrna_short.blast.len200.bbh.nonME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG'))
LD12_nonME_df <- left_join(LD12_nonME_df,meta_sub,by = c("SAG"="Genome.name"))
LD12_nonME <- mydenplotfunc(LD12_nonME_df)+ggtitle("LD12-non-Mendota")

acI_nonME_df <- read.delim('acI/PTXW.len150-vs-acI_norrna.fna_short.blast.len200.bbh.nonME', header = FALSE, sep = '\t',col.names=c('Read', 'Contig', 'PID', 'aln_len','mismatches','gaps','q_start','q_end','s_start','s_end', 'evalue','bitscore','SAG'))
acI_nonME_df <- left_join(acI_nonME_df,meta_sub,by = c("SAG"="Genome.name"))
acI_nonME <- mydenplotfunc(acI_nonME_df)+ggtitle("acI-non-Mendota")

## Putting plots together and saving them
all_plots <- plot_grid(acI_ME,LD12_ME, FW_ME, acI_nonME, LD12_nonME, labels=c("A", "B", "C", "D", "E"), ncol = 3, nrow = 2, align='hv')
save_plot(paste0(outname,'.eps'), all_plots, base_height = 8, base_width = 18)
save_plot(paste0(outname,'.pdf'), all_plots, base_height = 8, base_width = 18)

