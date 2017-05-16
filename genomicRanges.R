###################
# Extract exons for a gene
###################
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exonsByGene <- exonsBy(txdb, by="gene", use.names=FALSE) # get exon list grouped by gene
exonsByGene_seqs <- getSeq(Hsapiens,exonsByGene) # get DNA sequence for each exon for each gene
entrezID <- names(exonsByGene) # get the entrez gene ID of all the genes in the list

## Bimap interface:
x <- org.Hs.egENSEMBL
mapped_genes <- mappedkeys(x) # Get the entrez gene IDs that are mapped to an Ensembl ID
xx <- as.list(x[mapped_genes]) # Convert to a list

# Given an index of the list of genes, from 1 to length(entrezID), one can obtain the entrez ID, 
# the ensemble ID, the range and the seqs of the gene transcripts
for(ind in seq(1,length(entrezID))){
print(ind)
entrezID[ind] # the entrez ID
# xx[entrezID[ind]] # get the ensemble ID
# exonsByGene[entrezID[ind]] # the range on the genome
# exonsByGene_seqs[[ind]] # the sequence
# Write to fasta format
out.fasta <- paste("~/Work/dataset/gene2vec/genedoc/", entrezID[ind],".fasta", sep="")
writeXStringSet(exonsByGene_seqs[[ind]], out.fasta, format="fasta")
}
