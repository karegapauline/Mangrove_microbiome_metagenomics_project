getwd()
#### add biome table, tree and metadata
biom_data <- import_biom(BIOMfilename = "table-with-taxonomyv2.biom", 
                         treefilename = "tree.nwk")
mapping_file <- import_qiime_sample_data(mapfilename = "sample-metadata.tsv")

#### Merge the OTU and mapping data into a phyloseq object
phylo <- merge_phyloseq(biom_data, mapping_file)

#### Add names to biom table and check phyloseq objects
colnames(tax_table(phylo))= c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
rank_names(phylo)

## Phyloseq components 
head(otu_table(phylo))
head(tax_table(phylo))
head(tax_table(phylo))
head(sample_data(phylo))
phy_tree(phylo)

###### Accessors
ntaxa(phylo)
nsamples(phylo)
sample_names(phylo)
taxa_names(phylo)
sample_sums(phylo)
taxa_sums(phylo)
rank_names(phylo)
rank_names(phylo)
sample_variables(phylo)
get_variable(phylo)
get_variable(phylo, varName = "X..single.end.PHRED.33.fastq.manifest.file")

#### number of samples
print ('Number of Samples in our Biom Table')
nsamples(phylo)

#### number of sequence variants
print ('Number of Sequence variants we have.')
ntaxa(phylo)

#### summary statistics of sampling depth
print ('Sequencing depth.')
depths <- sample_sums(phylo)
summary(depths)

#### data analysis and visualization
## plot_bar (v1)
p <- plot_bar(phylo)
plot(p)

## plot_bar (v2)
p <- plot_bar(phylo, fill = "Phylum") + facet_wrap(~X..single.end.PHRED.33.fastq.manifest.file, scales = "free_x", nrow = 1)
plot(p)

##### load the libraries needed
library(ggplot2)
library(ape)
library(phyloseq)
library(dplyr)
library(grid)
library(gridExtra)
library(readr)
library(scales)
library(tidyr)
library(reshape2)
library(purrr)
library(here)
library(plyr)
##### plot_composition manipulation
plot_composition <- function(phylo, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                             taxaRank2 = "Family", numberOfTaxa = 9, fill = NULL,
                             x = "Sample", y = "Abundance", facet_grid = NULL) 
#### plot_composition (Phylum within Bacteria)
p <- plot_composition(phylo, "Kingdom", "Bacteria", "Phylum", 5, fill = "Phylum") 
p <- p + facet_wrap(~X..single.end.PHRED.33.fastq.manifest.file, scales = "free_x", nrow = 1)
plot(p)

####### install the devtools package
install.packages("devtools")
library(devtools)

## plot_composition (Family within Proteobacteria)
p <- plot_composition(phylo, "Phylum", "Proteobacteria", "Family", 9, fill = "Family") 
p <- p + facet_wrap(~X..single.end.PHRED.33.fastq.manifest.file, scales = "free_x", nrow = 1)
plot(p)
###################### alpha diversity richness
### plot_richness
plot_richness(phylo, x = "X..single.end.PHRED.33.fastq.manifest.file", color = "X..single.end.PHRED.33.fastq.manifest.file", 
              measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(fill = EnvType), alpha = 0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())