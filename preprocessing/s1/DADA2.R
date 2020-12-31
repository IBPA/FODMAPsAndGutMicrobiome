#!/usr/bin/env Rscript

# Run DADA2 analysis for 16s rRNA sequence files of single study. 
#   Results are saved in: 
#     (V1_seqtab.nochim.rds,	V2_taxa.rds,	V3_fittedPhylTree.rds	V4_ps.rds)
#   Command to run it hpc1 cluster (after running fetch.sh on dataset's directory):
#     #sbatch -c 32 -n 1 -N 1 -t 3600 DADA2_a6.R
#     

# Process 16s rRNA
library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")

src_dir <- "./microbiome_data/"
shared_dir <- "../shared/"

# Forward and reverse fastq filenames have format: SRR*_1.fastq.gz and SRR*_2.fastq.gz
fnFs <- sort(list.files(src_dir, pattern="SRR.*_1.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SRR*_X.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(src_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

names(filtFs) <- sample.names


# Filter and Trim
out <- filterAndTrim(fnFs, filtFs, trimLeft = 12,
                     maxN=0, truncQ=2, maxEE = 2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)

# Learn Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

# Core DADA2 Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

# Merge paired reads (no merging, because forward only)
mergers <- dadaFs

# Create Sequence Table
seqtab <- makeSequenceTable(mergers)
print(sprintf("Log| dim(seqtab):%s", toString(dim(seqtab))))

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sprintf("Log| chimeras: %f %%", 100*(sum(seqtab)-sum(seqtab.nochim))/sum(seqtab))
saveRDS(seqtab.nochim, file.path(src_dir, "V1_seqtab.nochim.rds"))
print("************* Saved V1_seqtab.nochim.rds *************")

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(derepFs, getN) , sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "drepF", "merged", "nonchim")
rownames(track) <- sample.names
print("Log| Read counts in the pipeline:")
track

# Assign Taxonomy
dna <- DNAStringSet(getSequences(seqtab.nochim))
load(file.path(shared_dir, "SILVA_SSU_r132_March2018.RData")) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxa <- taxid
saveRDS(taxa, file.path(src_dir, "V2_taxa.rds"))
print("************* Saved V2_taxa.rds *************")

# Construct Phylogenetic Tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
library(phangorn)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) # Long running
fittedPhylTree <- fitGTR$tree
saveRDS(fittedPhylTree, file.path(src_dir, "V3_fittedPhylTree.rds"))
print("************* Saved V3_fittedPhylTree.rds *************")
detach("package:phangorn", unload=TRUE)

# **********  polyseq **********
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

#Load sample metadata from file:
samdf <- data.frame("name" = rownames(seqtab.nochim))
rownames(samdf) <- samdf$name

#Phyloseq
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(taxa),
               phy_tree(fittedPhylTree))

saveRDS(ps, file.path(src_dir, "V4_ps.rds"))
print("************* Saved V4_ps.rds *************")
