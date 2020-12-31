#!/usr/bin/env Rscript
# Run DADA2 analysis for 16s rRNA sequence files of a single study. 
#   Results are saved in separate files:
#     (V1_seqtab.nochim.rds,	V2_taxa.rds,	V3_fittedPhylTree.rds	V4_ps.rds)
#   Command to run on SLURM cluster (after running fetch.sh on under "./data/study_dir/"):
#     #sbatch -c 32 -n 1 -N 1 -t 3600 ./s_DADA2.R ./data/study_dir/

library("testit")
library("optparse")

this_parse_args <- function() {
  # Parse command line arguments
  #
  # Args: None
  #
  # Returns:
  #   A list containing: shared_dir and src_dir
  
  parser <- OptionParser(usage = "%prog [options] ./data/study_dir/")
  parser <- add_option(parser, c("-s", "--shared_dir"),
                       type="character", 
                       default="./data/shared/",
                       help="Directory including necessary files for DADA2 processing (i.e. Silva database SILVA_SSU_r132_March2018.RData) (required). [default: %default]" , 
                       metavar="path")
  arguments <- parse_args(parser, positional_arguments = 1)
  shared_dir <- arguments$options$shared_dir
  src_dir <- arguments$args
  assert(sprintf("Directory %s does not exist!", shared_dir), dir.exists(shared_dir))
  assert(sprintf("Directory %s does not exist!", src_dir), dir.exists(src_dir))
  
    return(list(shared_dir = shared_dir, src_dir = src_dir))
}

this_get_sample_filenames <- function(src_dir) {
  # Get microbiome sample file names in src_dir (and name for filtered )
  #
  # Args: src_dir in which  *_[1-2].fastq.gz microbiome files reside.
  #
  # Returns:
  #   A list sample filenames
  
  fnFWs <- sort(list.files(src_dir, pattern=".*R1_001.fastq.gz", full.names = TRUE))
  fnREVs <- sort(list.files(src_dir, pattern=".*R2_001.fastq.gz", full.names = TRUE))
  
  # Extract sample names, assuming filenames have format: *_[1-2].fastq.gz
  sample.names <- sapply(strsplit(basename(fnFWs), "_"), `[`, 1)
  
  # Place filtered files in filtered/ subdirectory
  filtFWs <- file.path(src_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtREVs <- file.path(src_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFWs) <- sample.names
  names(filtREVs) <- sample.names
  
  return(list(FWs = fnFWs,
              REVs = fnREVs,
              filtFWs = filtFWs,
              filtREVs = filtREVs))
}

this_log <- function(message, obj = NULL){
  print(sprintf("** LOG ** %s", message))
  if (!is.null(obj)){
    print(obj)
  }
}

# Get sample filenames
#opts <- this_parse_args()

opts <- list(src_dir = "./microbiome_data/", shared_dir = "../shared/")
sfns <- this_get_sample_filenames(opts$src_dir)

# Start DADA2 calls
library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")

# Filter and Trim [ToDo: find the right parameters automatically!]
out <- filterAndTrim(sfns$FWs, sfns$filtFWs, sfns$REVs, sfns$filtREVs, 
                     truncLen=c(230,150), trimLeft = c(16, 16),
                     maxN=0, maxEE=c(4,7), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
this_log("filterAndTrim output", out)

# Learn Error Rates
errF <- learnErrors(sfns$filtFWs, multithread=TRUE)
errR <- learnErrors(sfns$filtREVs, multithread=TRUE)

# Dereplicate identical reads
derepFs <- derepFastq(sfns$filtFWs, verbose = TRUE)
derepRs <- derepFastq(sfns$filtREVs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- names(sfns$filtFWs)
names(derepRs) <- names(sfns$filtREVs)

# Core DADA2 Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Create Sequence Table
seqtab <- makeSequenceTable(mergers)

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
this_log(sprintf("Chimeras: %f %%", 100*(sum(seqtab)-sum(seqtab.nochim))/sum(seqtab)))
saveRDS(seqtab.nochim, file.path(opts$src_dir, "V1_seqtab.nochim.rds"))
seqtab.nochim <- readRDS(file.path(opts$src_dir, "V1_seqtab.nochim.rds"))
this_log("Saved V1_seqtab.nochim.rds")


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- names(dadaFs)
this_log("Read counts in the pipeline:", track)


# Assign Taxonomy
dna <- DNAStringSet(getSequences(seqtab.nochim))
load(file.path(opts$shared_dir, "SILVA_SSU_r132_March2018.RData")) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE) # use all processors
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
saveRDS(taxa, file.path(opts$src_dir, "V2_taxa.rds"))
this_log("Saved V2_taxa.rds")

# Construct Phylogenetic Tree
library(phangorn)
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) # Long running
fittedPhylTree <- fitGTR$tree
saveRDS(fittedPhylTree, file.path(opts$src_dir, "V3_fittedPhylTree.rds"))
this_log("Saved V3_fittedPhylTree.rds")
detach("package:phangorn", unload=TRUE)


# **********  phyloseq **********
library(Biostrings); packageVersion("Biostrings")
#Load sample names from seqtab.nochim as metadata:
samdf <- data.frame("name" = rownames(seqtab.nochim))
rownames(samdf) <- samdf$name
# Create and save phyloseq object
library(phyloseq); packageVersion("phyloseq")
library(ape); packageVersion("ape")
ps <- phyloseq(otu_table(t(seqtab.nochim), taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(taxa),
               phy_tree(fittedPhylTree))
saveRDS(ps, file.path(opts$src_dir, "V4_ps.rds"))
this_log("Saved V4_ps.rds")
