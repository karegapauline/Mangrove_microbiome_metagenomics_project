setwd("/home/babydoll/Downloads/phyloseqfinal")
getwd()
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(biomformat)
#### add biome table, tree and metadata
biom_data <- import_biom(BIOMfilename = "table-with-taxonomyv2.biom", 
                         treefilename = "tree.nwk")
mapping_file <- import_qiime_sample_data(mapfilename = "metadata.tsv")

#### Merge the OTU and mapping data into a phyloseq object
phylo <- merge_phyloseq(biom_data, mapping_file)

#### Add names to biom table and check phyloseq objects
colnames(tax_table(phylo))= c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
rank_names(phylo)

## Phyloseq components 
head(otu_table(biom_data))
head(tax_table(biom_data))
head(sample_data(phylo))
phy_tree(biom_data)
rank_names(biom_data)
sample_variables(phylo)
View(phylo)
plot_bar(phylo, fill = "Phylum") + theme(legend.position="bottom")
### filtering our data set
# first generate prevalence table
prevelancedf = apply(X = otu_table(phylo),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})
#### Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(phylo),
                          tax_table(phylo))
prevelancedf[1:10,]
phylo1 <- subset_taxa(phylo, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
View(phylo)
library(dplyr)
library(plyr)
plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})
phyla2Filter = c("p__OP3", "p__OP8", "p__SAR406",
                 "p__GN04","p__")
# Filter entries with unidentified Phylum.
phylo1 = subset_taxa(phylo1, !Phylum %in% phyla2Filter)
qplot(colSums(otu_table(phylo1)),bins=30) +
  xlab("Logged counts-per-sample")
plot_richness(phylo1, measures=c("Observed","Chao1", "shannon"))
hell.tip.labels <- as(get_variable(phylo1, "X.2"), "character")
d <- distance(phylo1, method="bray", type="samples")
hell.hclust     <- hclust(d, method="average")
plot(hell.hclust)
#Lets write out a plot
pdf("My_dendro.pdf", width=7, height=7, pointsize=8)
plot(hell.hclust)
dev.off()
png("My_dendro.png", width = 7, height = 7, res=300, units = "in")
plot(hell.hclust)
dev.off()
v4.hell.ord <- ordinate(phylo1, "NMDS", "bray")
p1 = plot_ordination(phylo1, v4.hell.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 5)
p2 = plot_ordination(phylo1, v4.hell.ord, type="X.2", color="X.2", shape = "X.2") 
p2 + geom_polygon(aes(fill=X.2)) + geom_point(size=5) + ggtitle("samples")
