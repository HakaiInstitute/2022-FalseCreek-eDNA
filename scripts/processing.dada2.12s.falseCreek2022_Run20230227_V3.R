#notes for processing falseCreek2022_Run20230227_V3 run
#author: Evan Morien
#created: 2023/03/01
#last modified: 2023/03/01
#working directory: ~/projects/12S_runs/falseCreek2022_Run20230227_V3/

####Libraries####
library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)

####Environment Setup####
theme_set(theme_bw())
setwd("~/projects/12S_runs/falseCreek2022_Run20230227_V3/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "~/projects/12S_runs/falseCreek2022_Run20230227_V3/raw_data/" # CHANGE ME to the directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name

####fastq Quality Plots####
numplot <- 49 #CHANGE ME to the number of samples you want to include in the quality plot
a <- sample(fnFs, numplot)
b <- which(fnFs %in% a)
plotfnFs <- fnFs[b]
plotfnRs <- fnRs[b]
pdf("quality_plots.dada2.R1s.pdf", width = 16, height = 9) # define plot width and height. completely up to user.
plotQualityProfile(plotfnFs) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
dev.off()
pdf("quality_plots.dada2.R2s.pdf", width = 16, height = 9) # define plot width and height. completely up to user.
plotQualityProfile(plotfnRs)
dev.off()

#quality plots must be inspected at this point. if they look okay, we can proceed with the rest of the processing. if not, consult Rute & Matt, explain what looks off to them

####primer checking####
#these are the MiFish-U primers
FWD <- "GTCGGTAAAACTCGTGCCAGC"
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, trimRight = c(0,0), truncLen=c(180,150), trimLeft = c(0,0), maxN = 0, multithread = 16, compress = TRUE, matchIDs=TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
index <- 2 #CHANGE ME (if necessary) this is the index of the file we want to check for primers, within the lists "fn*s.filtN", it can be any number from 1 to N, where N is the number of samples you are processing. don't use extraction, pcr, or field blanks for this sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[index]]), #the index of the sample you'd like to use for this test is used here (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[index]]))

####OPTIONAL!!!! change primer orientation####
#REV <- REV.orients[["Reverse"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS_workflow guide section "Identify Primers" for details

#### primer removal ####
cutadapt <- "/usr/local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-j", 16,# -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
#sanity check, should report zero for all orientations and read sets
index <- 2 #this is the index of the file we want to check for primers, within the lists "fn*s.cut", it can be any number from 1 to N, where N is the number of samples you are processing
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[index]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[index]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[index]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[index]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE)) #remember to change this so it matches ALL your file names!
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE)) #remember to change this so it matches ALL your file names!

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#IMPORTANT: DO NOT CHANGE THESE PARAMETERS FOR ICO EXPERIMENTS, THEY ARE SET AND MUST BE RETAINED FOR DATA COMPARABILITY
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimLeft = c(0,0), trimRight = c(0,0), minLen = c(110,110),
                     maxN=c(0,0), maxEE=c(4,6), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=16)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
write.table(retained, "retained_reads.filterAndTrim_step.length_var.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=16)
errR <- learnErrors(filtRs, multithread=16)

pdf("error_rates.dada2.R1s.length_var.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
plotErrors(errF, nominalQ=TRUE) #assess this graph after processing. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered carefully before moving on.
dev.off()
pdf("error_rates.dada2.R2s.length_var.pdf", width = 10, height = 10) # define plot width and height. completely up to user.
plotErrors(errR, nominalQ=TRUE) #assess this graph after processing. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered carefully before moving on.
dev.off()

#sometimes, samples with low sequence count don't produce any filtered output, so we have to do something extra here to get rid of their entries in filtFs and filtRs
filtFs <- sort(list.files("~/projects/12S_runs/falseCreek2022_Run20230227_V3/raw_data/cutadapt/filtered/", pattern = "_R1_001.fastq.gz", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!
filtRs <- sort(list.files("~/projects/12S_runs/falseCreek2022_Run20230227_V3/raw_data/cutadapt/filtered/", pattern = "_R2_001.fastq.gz", full.names = TRUE)) #remember to change the pattern so it matches ALL your file names!

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=16)
dadaRs <- dada(derepRs, err=errR, multithread=16)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 50 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing
write.table(names(which(samples_to_keep == TRUE)), "samples_retained.txt", row.names=FALSE, quote=F, sep="\n")
write.table(setdiff(sample.names, names(which(samples_to_keep == TRUE))), "samples_removed.txt", row.names=FALSE, quote=F, sep="\n")


####merge paired reads####
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
pdf("length_histogram.merged_reads.length_var.pdf", width = 10, height = 8) # define plot width and height. completely up to user.
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot
dev.off()
#review after processing as sanity check. most ASVs should be w/in a few bp of the amplicon length

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
dim(otus) # (this should be the same as last command, but the dimensions reversed)

otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
removed <- length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

#Start ASV report
cat("dim(otus): ", dim(otus),file="ASV_report.txt",sep="\t",append=TRUE)
cat(" ", file="ASV_report.txt", sep="\n", append=TRUE)
cat("# ASVs Removed: ", length(intersect(a,b)),file="ASV_report.txt",sep="\t",append=TRUE)
cat(" ", file="ASV_report.txt", sep="\n", append=TRUE)
cat("dim(otus_filt): ", dim(otus_filt),file="ASV_report.txt",sep="\t",append=TRUE)
cat(" ", file="ASV_report.txt", sep="\n", append=TRUE)

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=16, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
cat("dim(seqtab.nosingletons.nochim): ", dim(seqtab.nosingletons.nochim),file="ASV_report.txt",sep="\t",append=TRUE)
cat(" ", file="ASV_report.txt", sep="\n", append=TRUE)
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low
cat("proprtion of chimeras: ", sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons),file="ASV_report.txt",sep="\t",append=TRUE)

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100, track[,7]/track[,1]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras", "percent_retained_of_total")

####save output from sequence table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.12S_merged.length_var.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.12S.merged.length_var.txt", row.names=FALSE, quote=F, sep="\t")


#OPTIONAL: read in sequence table on cluster, for taxonomy assignment
#code to read in before taxonomy assignment, if done on a separate machine
seqtab.nosingletons.nochim <- fread("sequence_table.12S.merged.length_var.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix
mode(seqtab.nosingletons.nochim) <- "numeric"


#### save sequences and do taxonomy assignment with blast ####
##### replace the long ASV names (the actual sequences) with human-readable names####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.length_var.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "12S_ASV_sequences.length_var.fasta") #save sequences with new names in fasta format
#IMPORTANT: sanity checks
colnames(seqtab.nosingletons.nochim) == ASV.seq #only proceed if this tests as true for all elements

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtab.nosingletons.nochim) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.12S.merged.w_ASV_names.length_var.txt", row.names=FALSE, quote=F, sep="\t")


#assign taxonomy with blast NT database at 96% similarity threshold
mkdir blast_96_sim_LCA_besthit
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db /data/taxonomyDBs/NCBI_NT/2023-02-07/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query 12S_ASV_sequences.length_var.fasta  -out blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out -t /data/taxonomyDBs/NCBI_taxonomy/2023-03-02/rankedlineage.dmp -m /data/taxonomyDBs/NCBI_taxonomy/2023-03-02/merged.dmp -o blast_96_sim_LCA_besthit/taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_12S_ASV_sequences.length_var.blast.out > tmp
python2 ~/programs/galaxy-tool-lca/lca.species.py -i tmp -o blast_96_sim_LCA_besthit/taxonomy_table.12S.NCBI_NT.96sim.LCA+besthit.txt -b 100 -id 96 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentified,kingdom -flh unclassified

mkdir blast_96_sim_LCA_only
python2 ~/programs/galaxy-tool-lca/lca.species.py -i tmp -o blast_96_sim_LCA_only/taxonomy_table.12S.NCBI_NT.96sim.LCA_only.txt -b 100 -id 96 -cov 50 -t only_lca -fh environmental,unidentified -flh unclassified

#cleanup
rm blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out #remove blast output without taxonomy
rm taxonomy_12S_ASV_sequences.length_var.blast.out #remove redundant file
mv tmp blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out #replace with taxonomy added blast output

