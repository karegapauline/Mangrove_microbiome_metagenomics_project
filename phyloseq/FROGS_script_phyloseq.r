######### installation of packages needed for analysis
#install.packages("devtools") # now installed
#library("devtools")
#library(phyloseq)
#library(ggplot2)
#library(plyr)
#library(reshape2)
#library(grid)
#library(gridExtra)
#library(ape)
#library(scales)

########## phyloseq analysis:importing data
import data
biomfile <- "/home/babydoll/Downloads/EANBIT_mangrove/data/feature-table.biom"
library(phyloseq)
## create a variable that contain the table  with the data
mangrove1 <-import_biom(biomfile)
mangrove1 <-importIntoEnv(biomfile)
food1
library(tidyverse)
sampledata <-read.csv("data/chaillou/sample_data.tsv", sep="\t", row.names = 1)
sample_data(food1) <- sampledata
food1

treefile <- read.tree("data/chaillou/tree.nwk")
phy_tree(food1) <- treefile
food1

## import data 
biomfile2 <- "data/chaillou/chaillou_with_sam_data.biom"
food <- import_biom(biomfile2, treefile, parseFunction = parse_taxonomy_greengenes)
food

#####################
##                 ##
##  DATA STRUCTURE ##
##                 ##
#####################

## Phyloseq components 
head(otu_table(food))
head(tax_table(food))
head(tax_table(food1))
head(sample_data(food))
phy_tree(food)

## Accessors
ntaxa(food)
nsamples(food)
sample_names(food)
taxa_names(food)
sample_sums(food)
taxa_sums(food)
rank_names(food)
rank_names(food1)
sample_variables(food)
get_variable(food)
get_variable(food, varName = "EnvType")
levels(get_variable(food, varName = "EnvType"))
get_sample(food, i = "otu_00520")
get_taxa(food, i = "MVT0.LOT07")

## sample library sizes
sample_sums(food)[c("MVT0.LOT01", "MVT0.LOT07", "MVT0.LOT09")]
## Otu overall abundances
taxa_sums(food)[c("otu_00520", "otu_00569", "otu_00527")]

#########################
##                     ##
##  DATA MANIPULATION  ##
## (MODIFICATION)      ##
##                     ##
#########################
## change taxonomic rank names
new_rank <- c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")
colnames(tax_table(food)) <- new_rank
rank_names(food1) ## or head(tax_table(food))

## change level order for EnvType
correct.order <- c("BoeufHache", "VeauHache", "DesLardons", 
                   "MerguezVolaille", "SaumonFume", "FiletSaumon", 
                   "FiletCabillaud", "Crevette")
sample_data(food)$EnvType <- factor(sample_data(food)$EnvType, 
                                    levels = correct.order)
levels(get_variable(food, "EnvType")) 

## Change otu count, or taxonomic affiliation
otu_table(food)["otu_00520", "DLT0.LOT08"] <- 0
tax_table(food)["otu_00520", "Species"] <- "Ornithinolytica"

#########################
##                     ##
##  DATA MANIPULATION  ##
## (FILTER, SMOOTHING) ##
##                     ##
#########################

## Filters, subset, prune 
samplesToKeep <- sample_names(food)[1:10]
prune_samples(samplesToKeep, food)

subset_samples(food, EnvType %in% c("DesLardons", "MerguezVolaille"))


small.food <- subset_taxa(food, Phylum == "Firmicutes" & Class == "Bacilli")
head(tax_table(small.food)[ , c("Phylum", "Class", "Order")])
## Unique combinations (Phylum, Class)
unique(tax_table(small.food)[ , c("Phylum", "Class")])

## aggregation
mergedData <- tax_glom(food, "Phylum")
ntaxa(mergedData) ## number of different phyla
tax_table(mergedData)[1:2, c("Phylum", "Order", "Class")]


## rarefaction and count transformation
foodRare <- rarefy_even_depth(food, rngseed = 1121983)
sample_sums(foodRare)[1:5]
count_to_prop <- function(x) { return( x / sum(x) )}
foodTrans <- transform_sample_counts(food, count_to_prop)
sample_sums(foodTrans)[1:5] ## should be 1

#########################
#########################
##                     ##
##  DATA VISUALIZATION ##
##                     ##
#########################
#########################

###################
##               ##
##  COMPOSITION  ##
##               ##
###################

## plot_bar (v1)
p <- plot_bar(food)
plot(p)

## plot_bar (v2)
p <- plot_bar(food, fill = "Phylum") + facet_wrap(~EnvType, scales = "free_x", nrow = 1)
plot(p)

## plot_composition (Phylum within Bacteria)
p <- plot_composition(food, "Kingdom", "Bacteria", "Phylum", 5, fill = "Phylum") 
p <- p + facet_wrap(~EnvType, scales = "free_x", nrow = 1)
plot(p)

## plot_composition (Family within Proteobacteria)
p <- plot_composition(food, "Phylum", "Proteobacteria", "Family", 9, fill = "Family") 
p <- p + facet_wrap(~EnvType, scales = "free_x", nrow = 1)
plot(p)

#####################
##                 ##
## ALPHA DIVERSITY ##
##    RICHNESS     ##
##                 ##
#####################

## plot_richness
plot_richness(food, x = "EnvType", color = "EnvType", 
              measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(fill = EnvType), alpha = 0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())

plot_richness(food, x = "FoodType", color = "FoodType", 
              measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(fill = FoodType), alpha = 0.2) + theme_bw() + geom_point() + theme(axis.text.x = element_blank())


## plot_richness (after trimming)
plot_richness(prune_taxa(taxa_sums(food) >= 500, food), color = "EnvType", 
              measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"), 
              x = "EnvType") + geom_boxplot(aes(fill = EnvType), alpha = 0.2) + theme_bw()  + geom_point() +  theme(axis.text.x = element_blank())

## Numeric values of alpha diversity indices
alpha.diversity <- estimate_richness(food, 
                                     measures = c("Observed", "Chao1", "Shannon"))
head(alpha.diversity)

## export diversity to text file
## write.table(alpha.diversity, "myfile.txt")

## Anova on observed richness
data <- cbind(sample_data(food), alpha.diversity)
food.anova <- aov(Observed ~ EnvType, data)
summary(food.anova) ## significant effect of environment type on richness

food.anova <- aov(Shannon ~ EnvType, data)
summary(food.anova)

####################
##                ##
## BETA DIVERSITY ##
##   DISTANCES    ##
##                ##
####################


## Generic method for distance
dist.bc <- distance(food, method = "bray") ## Bray-Curtis

## all available distances
distanceMethodList

## some code to visualize distances
sampleOrder <- levels(reorder(sample_names(food), as.numeric(get_variable(food, "EnvType")))) 
plot_dist_as_heatmap <- function(dist, order = sampleOrder, title = NULL) {
    data <- melt(as(dist, "matrix"))
    colnames(data) <- c("x", "y", "distance")
    if (!is.null(order)) {
        data$x <- factor(data$x, levels = order)
        data$y <- factor(data$y, levels = order)
    }
    p <- ggplot(data, aes(x = x, y = y, fill = distance)) + geom_tile() 
    p <- p + theme(axis.title.x = element_blank(), 
                   axis.title.y = element_blank(), 
                   axis.text.x = element_blank(), 
                   axis.text.y = element_blank()
                   )
    p <- p + scale_fill_continuous(limits = c(0, 1))
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    return(p)
}

## Compare Bray-Curtis (bray) and Jaccard (cc)
dist.bc <- distance(food, "bray")
p1 <- plot_dist_as_heatmap(dist.bc, title = "Bray-Curtis")
plot(p1)
dist.jac <- distance(food, "cc")
p2 <- plot_dist_as_heatmap(dist.jac, title = "Jaccard")
grid.arrange(p1, p2, ncol = 2, nrow = 1)

## Compare unifrac and wunifrac
dist.uf <- distance(food, method = "unifrac") ## Unifrac
dist.wuf <- distance(food, method = "wunifrac") ## Weighted Unifrac
p3 <- plot_dist_as_heatmap(dist.uf, title = "Unifrac")
p4 <- plot_dist_as_heatmap(dist.wuf, title = "wUnifrac")
grid.arrange(p3, p4, ncol = 2, nrow = 1)

## compare Bray-Curtis, Jaccard, Unifrac and wUnifrac
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

################
##            ##
## ORDINATION ##
##            ##
################


## Ordination (MDS/PCoA is the same thing)
ord <- ordinate(food, method = "MDS", distance = "bray")

## Alternative way to compute ordinations
## dist.bc <- distance(food, method = "bray")
## ord <- ordinate(food, method = "MDS", distance = dist.bc)

## plot_ordination
p <- plot_ordination(food, ord, color = "EnvType")
p <- p + theme_bw() + ggtitle("MDS + BC")
plot(p)

## plot_samples (add group means to the plot)
p <- plot_samples(food, ord, color = "EnvType", replicate = "EnvType")
p <- p + theme_bw() + ggtitle("MDS + BC")  
plot(p)

## Compare ordinations obtained from Bray-Curtis, Jaccard, Unifrac and wUnifrac distances
p1 <- plot_samples(food, ordinate(food, "MDS", "bray"), color = "EnvType") + ggtitle("MDS + BC") + theme_bw()
p2 <- plot_samples(food, ordinate(food, "MDS", "cc"), color = "EnvType") + ggtitle("MDS + Jaccard") + theme_bw()
p3 <- plot_samples(food, ordinate(food, "MDS", "unifrac"), color = "EnvType") + ggtitle("MDS + UF") + theme_bw()
p4 <- plot_samples(food, ordinate(food, "MDS", "wunifrac"), color = "EnvType") + ggtitle("MDS + wUF") + theme_bw()
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

################
##            ##
## CLUSTERING ##
##            ##
################

## General method for clustering samples 
## clustering <- hclust(distance.matrix, method = "linkage.function")
## plot(clustering)

## Clustering of samples (from UniFrac distances) using different linkage method
par(mfcol = c(1,3)) ## To plot the three clustering trees side-by-side
plot(hclust(dist.uf, method = "complete"))
plot(hclust(dist.uf, method = "ward.D2"))
plot(hclust(dist.uf, method = "single"))

## Better dendrograms, color leaves (samples) according to EnvType
## Env types
envtype <- get_variable(food, "EnvType")
## automatic color palette: one color per different sample type
palette <- hue_pal()(length(levels(envtype)))
## Map sample type to color
tipColor <- col_factor(palette, levels = levels(envtype))(envtype)
## Change hclust object to phylo object and plot
clust.uf <- as.phylo(hclust(dist.uf, method = "ward.D2"))
plot(clust.uf, tip.color = tipColor, direction = "downwards", main = "Ward linkage")

## To plot each linkage method side-by-side
par(mfcol= c(1,3))
clust.uf <- as.phylo(hclust(dist.uf, method = "single"))
plot(clust.uf, tip.color = tipColor, direction = "downwards", main = "Single linkage")
clust.uf <- as.phylo(hclust(dist.uf, method = "ward.D2"))
plot(clust.uf, tip.color = tipColor, direction = "downwards", main = "Ward linkage")
clust.uf <- as.phylo(hclust(dist.uf, method = "complete"))
plot(clust.uf, tip.color = tipColor, direction = "downwards", main = "Complete linkage")


################
##            ##
##  HEATMAP   ##
##            ##
################

## heatmap (v1)
p <- plot_heatmap(food)
plot(p)

## heatmap (v2, with more usual colors)
p <- plot_heatmap(food, low = "yellow", high = "red", na.value = "white") 
plot(p)

## heatmap (v3, with more usual colors and samples organized by EnvType)
p <- plot_heatmap(food, low = "yellow", high = "red", na.value = "white")
p <- p + facet_grid(~EnvType, scales = "free_x")
plot(p)

## heatmap (v4, alternative color scheme)
p <- plot_heatmap(food)
p <- p + scale_fill_gradient2(low = "#1a9850", mid = "#ffffbf", high = "#d73027", 
                              na.value = "white", trans = log_trans(4), 
                              midpoint = log(100, base = 4))
p <- p + facet_grid(~EnvType, scales = "free_x")
plot(p)

##################
##              ##
## MULTIVARIATE ##
##   ANOVA      ##
##              ##
##################

## MANOVA with CAP (Canonical Analysis of Principal Coordinates) 
## CAP through capscale
metadata <- as(sample_data(food), "data.frame") ## convert sample_data to data.frame
cap <- capscale(dist.uf ~ EnvType, data = metadata)
cap
anova <- anova(cap, permutations = 999)
print(anova)

## MANOVA with CAP (Canonical Analysis of Principal Coordinates) 
metadata <- as(sample_data(food), "data.frame")
adonis(dist.uf ~ EnvType, data = metadata, perm = 9999)

#####################
#####################
##                 ##
##  ADVANCED USE   ##
##                 ##
#####################
#####################

##################
##              ##
##    MANUAL    ##
##    IMPORT    ##
##              ##
##################


## Import raw data in R
sampledata <- read.csv("data/manual/sampledata.tsv", sep = "\t", row.names = 1)
taxtable <- read.csv("data/manual/taxtable.tsv", sep = "\t", row.names = 1)
otutable <- read.csv("data/manual/otutable.tsv", sep = "\t", row.names = 1)
tree <- read.tree("data/manual/tree.phy")

## Convert to R base data types
taxtable <- as.matrix(taxtable) 
otutable <- as.matrix(otutable)

## Convert to phyloseq base data types
sampledata <- sample_data(sampledata)
taxtable <- tax_table(taxtable) 
otutable <- otu_table(otutable, taxa_are_rows = TRUE)  

## Check name consistency (samples)
all(colnames(otutable) %in% rownames(sampledata)) ## OTU table and metadata

## Check name consistency (taxa)
all(rownames(otutable) %in% rownames(taxtable)) ## OTU table and Taxonomy table
all(rownames(otutable) %in% tree$tip.label) ## OTU table and tree leaves

## Construct phyloseq object 
manualData <- phyloseq(sampledata, otutable, taxtable, tree)
manualData

##################
##              ##
##   ADVANCED   ##
##   FILTERS    ##
##              ##
##################

## Sample wise condition (example I): filter taxa present in at least 5 samples
condition <- function(x) { x > 0 }
taxaToKeep <- genefilter_sample(food, condition, 5)
prune_taxa(taxaToKeep, food)

## Sample wise condition (example II): filter taxa in the top 250 taxa present in at least 3 samples
condition <- function(x) { order(x, decreasing = TRUE) <= 250 }
taxaToKeep <- genefilter_sample(food, condition, 3)
prune_taxa(taxaToKeep, food)

## sample wide condition (example I): filter taxa present in at least 5 samples
condition <- function(x) { sum(x > 0) >= 5 }
taxaToKeep <- filter_taxa(food, condition)
prune_taxa(taxaToKeep, food)

## sample wide condition (example II): filter taxa present with global abundance at least 100
condition <- function(x) { sum(x) >= 100 }
taxaToKeep <- filter_taxa(food, condition)
prune_taxa(taxaToKeep, food)

##################
##              ##
##   ADVANCED   ##
##  SMOOTHING   ##
##              ##
##################


## Merge all samples according to some metadata
mergedData <- merge_samples(food, "EnvType")
mergedData
sample_names(mergedData)
## Beware, some information is lost in the process
sample_data(mergedData)[1:2, ]
sample_data(food)[1:2, ]

## Merge taxa at at given height in the tree
mergedData <- tip_glom(food, h = 0.3)

##################
##              ##
## RAREFACTION  ##
##    CURVES    ##
##              ##
##################

## Rarefy samples to same depth and compute rarefaction curves
food <- rarefy_even_depth(food, rngseed = 1121983)
p <- ggrare(food, step = 1000, color = "EnvType", se = FALSE)
p <- p + facet_wrap(~EnvType, ncol = 4) + theme_bw()
plot(p)

###################
##               ##
## CUSTOM COLOR  ##
##  PALETTE      ##
##               ##
###################
# check EnvType order
levels(sample_data(food)$EnvType)

# if you want to change the order
## correct.order <- c("BoeufHache", "VeauHache", "DesLardons", 
##                   "MerguezVolaille", "SaumonFume", "FiletSaumon", 
##                   "FiletCabillaud", "Crevette")
##sample_data(food)$EnvType <- factor(sample_data(food)$EnvType, 
##                                    levels = correct.order)
# desired color palette
foodPalette <- c('#61001f','#b2182b','#d6604d','#f4a582','#92c5de','#4393c3','#2166ac','#053061')
# associated EnvType level
names(foodPalette) <- levels(sample_data(food)$EnvType)
foodPalette

## Try this palette on clustering
## Map sample type to color
tipColor <- col_factor(foodPalette, levels = levels(envtype))(envtype)
## Change hclust object to phylo object and plot
clust.uf <- as.phylo(hclust(dist.uf, method = "ward.D2"))
plot(clust.uf, tip.color = tipColor, direction = "downwards", main = "Ward linkage")

## Try this on any ggplot
food <- rarefy_even_depth(food, rngseed = 1121983)
p <- ggrare(food, step = 1000, color = "EnvType", se = FALSE)
p <- p + facet_wrap(~EnvType, ncol = 4) + theme_bw()
p <- p + scale_color_manual(values = foodPalette)
plot(p)

###################
##               ##
## CUSTOM ORDER  ##
## OTU IN CLUST  ##
##               ##
###################
# custom otu order in function of prevalence
prevalence <- estimate_prevalence(food,"EnvType")

# get in a different format to estimate correlations
correlationData <- estimate_prevalence(food, "EnvType", format = "wide")
correlationData <- t(correlationData)
correlation <- cor(correlationData, method="pearson")

##clustering sample and order otu according to the tree
otu.clust <- hclust(as.dist(1-c), method = "complete")
otuOrder <- otu.clust$labels[otu.clust$order]

## plot
p <- plot_heatmap(food, taxa.order = otuOrder)
p <- p + facet_grid(~EnvType, scales = "free", space = "free")
p <- p + scale_fill_gradient2 (low="#1a9850", mid = "#ffffbf", high="#d73027", na.value = "white", trans = log_trans(4), midpoint = log(100,base= 4 ))
plot(p)
