# Accomac_Sediment_Testing

## Load Packages
```{r}
#All the packages you will need for entire code, some of these could be unnecessary pending what you want to graph
library(dada2)
library(rmarkdown)
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(devtools)
library(vegan)
library(dbplyr)
library(microbiome)
library(tidyverse)
library(DECIPHER)
```

```{r}
formatPvalues <- function(pvalue) {
  ra<- ""
  if(pvalue <= 0.1) ra<- "."
  if(pvalue <= 0.05) ra<- "*"
  if(pvalue <= 0.01) ra<- "**"
  if(pvalue <= 0.001) ra<- "***"
  return(ra)
}
```

## DADA2 Pipeline
```{r}
# Path needs to have fastq files unzipped!
path <- "/Users/maggieshostak/Desktop/Accomac_Sediment/ShostakV4V5"
#list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)

#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
trimLeft = c(FWD, REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280,200), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE, trimLeft = c(19,20))
#head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)
```

```{r}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#dadaFs[[1]]

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#dadaRs[[1]]
```

```{r}
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

```{r}
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "/Users/maggieshostak/Desktop/Accomac_Sediment/data/track_sequences.csv", sep=",", quote=F, col.names=NA)
```

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "/Users/maggieshostak/Desktop/Accomac_Sediment/ShostakV4V5/silva_nr99_v138.1_train_set.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
```

```{r}
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
asv_headers[i] <- paste(">ASV", i, sep="_")
}

asv_fasta <- c(rbind(asv_headers, asv_seqs))

asv_otu <- t(seqtab.nochim)
row.names(asv_otu) <- sub(">", "", asv_headers)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)

otu_tax_table <- merge(asv_otu, asv_tax, by=0)

write(asv_fasta, "asv_fasta_AC_Sed.fa")
write.table(asv_otu, "asv_otu_AC_Sed.csv", sep=",", quote=F, col.names=NA)
write.table(asv_tax, "asv_tax_AC_Sed.csv", sep=",", quote=F, col.names=NA)
write.table(otu_tax_table, "otu_tax_table_AC_Sed.csv", sep=",", quote=F, col.names=NA)
```

```{r}
otu_counts <- read.csv("/Users/maggieshostak/Desktop/Accomac_Sediment/data/asv_otu_AC_Sed.csv") %>%
  pivot_longer(-ASV, names_to="sample_id", values_to = "count")
otu_counts

taxonomy <- read.csv("/Users/maggieshostak/Desktop/Accomac_Sediment/data/asv_tax_AC_Sed.csv")
taxonomy

metadata <- read.csv("/Users/maggieshostak/Desktop/Accomac_Sediment/data/accomac_sed_metadata.csv")
metadata
```

# Relative Abundance
```{r}
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="ASV") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "ASV"),
         names_to="level",
         values_to="taxon")
otu_rel_abund

write.table(otu_rel_abund, "/Users/maggieshostak/Desktop/Accomac_Sediment/data/otu_rel_abund.csv", sep=",", quote=F, col.names=NA)
```

# Barchart
```{r}
## Phylum
otu_rel_abund %>%
  filter(level=="Phylum") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/phylum_stacked_barchart_Acc_sed.tiff", width=20, height=10)

## Class
otu_rel_abund %>%
  filter(level=="Class") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/class_stacked_barchart_Acc_sed.tiff", width=25, height=10)

## Order
otu_rel_abund %>%
  filter(level=="Order") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/order_stacked_barchart_Acc_sed.tiff", width=55, height=10, limitsize = FALSE)

## Family
otu_rel_abund %>%
  filter(level=="Family") %>%
  group_by(sample_id, location, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  group_by(location, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund)) %>%
  ggplot(aes(x=location, y=mean_rel_abund, fill=taxon)) +
  geom_col(aes(x=location, y=mean_rel_abund), colour="black", stroke=10) +
    labs(x=NULL, 
         y="Mean Relative Abundance (%)") +
    theme_classic()
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/family_stacked_barchart_Acc_sed.tiff", width=50, height=10, limitsize = FALSE)
```

# NMDS
```{r}
# All Samples: Biofilm Separated
pc1 = read.csv("/Users/maggieshostak/Desktop/Accomac_Sediment/data/nmds_asv_otu_AC_Sed.csv")
pc1

com1 = pc1[,4:ncol(pc1)]
com1

m_com1 <- as.matrix(com1)

set.seed(1000)
nmds1 = metaMDS(m_com1, distance = "bray")
```

```{r}
data.scores1 <- as.data.frame(scores(nmds1)$sites)
data.scores1$location = pc1$location
data.scores1$sample_id = pc1$sample_id
head(data.scores1)

xx1 = ggplot(data.scores1, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = location))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "location", y = "NMDS2")
xx1
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/data/NMDS_AC_sediment_location.tiff", width = 10, height = 10)
```

```{r}
data.scores2 <- as.data.frame(scores(nmds1)$sites)
data.scores2$depth = pc1$depth
data.scores2$sample_id = pc1$sample_id
head(data.scores2)

xx2 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = depth))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "depth", y = "NMDS2")
xx2
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/data/NMDS_AC_sediment_depth.tiff", width = 10, height = 10)
```

```{r}
data.scores3 <- as.data.frame(scores(nmds1)$sites)
data.scores3$direction = pc1$direction
data.scores3$sample_id = pc1$sample_id
head(data.scores3)

xx3 = ggplot(data.scores3, aes(x = NMDS1, y = NMDS2)) +
 geom_point(size = 3, aes(colour = direction))+
  scale_fill_discrete()+
  ggtitle("NMDS Ordination - Samples Across Site")+
 theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
       axis.text.x = element_text(colour = "black", face = "bold", size = 12),
       legend.text = element_text(size = 12, face ="bold", colour ="black"),
       legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
       axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
       legend.title = element_text(size = 14, colour = "black", face = "bold"),
       panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
       legend.key=element_blank()) +
 labs(x = "NMDS1", colour = "direction", y = "NMDS2")
xx3
ggsave("/Users/maggieshostak/Desktop/Accomac_Sediment/data/NMDS_AC_sediment_direction.tiff", width = 10, height = 10)
```
