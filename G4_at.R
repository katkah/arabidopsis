if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("pqsfinder")
BiocManager::install("rtracklayer")
BiocManager::install("Gviz")
library("Biostrings")
library("pqsfinder")
library("rtracklayer")
library("Gviz")


#upload sequence
###############################################################################
s <- readDNAStringSet("WT_consensuses.fasta")
my_seq <- s[[20]]
s_string <- toString(my_seq)


#pqsfinder
###############################################################################
pqs <- pqsfinder(my_seq, min_score = 30) 

#export the result to .gff and .fa
#gr <- as(pqs, "GRanges")
#export(gr, "rDNA_pqs30_result.gff", version = "3")
#dss <- as(pqs, "DNAStringSet")
#writeXStringSet(dss, file = "rDNA_pqs30_result.fa", format = "fasta")


#rDNA sequence annotation
###############################################################################
ir <- IRanges(start = c(1,1133,1752,2835,4728,6800,7151), end = c(777,1193,1813,2896,6531,6963,10525), names = c("3ETS","spacer_promoter1","spacer_promoter2", "gene_promoter","18S","5.8S","25S") )
gr_annotation <- GRanges(seqnames = "chr1", strand = rep("+",7),
                         ranges = ir, type = c("3ETS","SP1","SP2", "GP","18S","5.8S","25S"))
#export(gr_annotation, "rDNA_annotation.gff", version = "3")


#plot the result
###############################################################################
gtrack <- AnnotationTrack(gr_annotation, 
                          name = "rDNA",
                          shape = "box",
                          group = c("3ETS","SP1","SP2", "GP","18S","5.8S","25S")
                          )
atrack <- GenomeAxisTrack()

pqs_df <- data.frame(chrom = "chr1", start = start(pqs), end = end(pqs), strand = strand(pqs), plus_strand = score(pqs), minus_strand = score(pqs))
pqs_df[which(pqs_df$strand == "+"),"minus_strand"] <- NA
pqs_df[which(pqs_df$strand == "-"),"plus_strand"] <- NA
pqs_df[,"strand"] <- "*"
pqs_df <- pqs_df[,c(1,2,3,4,6,5)] #switch the columns otherwise the legend in the plot has switched colors
  

#colnames(pqs_df)<-c("chrom","start","end","strand","minus_strand","plus_strand")
gr <- as(pqs_df, "GRanges")
pqs_track <- DataTrack(
  gr,
  name = "pqsfinder G4 score",
  col = c("cornflowerblue","red"),
  groups=c("plus_strand","minus_strand")
)



tiff("G4_rDNA_at.tiff", width = 6, height = 4, units = 'in', res = 300)
suppressWarnings(plotTracks(c(gtrack, pqs_track, atrack), groupAnnotation = "group", just.group = "below",type='h', labelPos="below", cex.group = 0.5, fontsize = 9))
dev.off()





