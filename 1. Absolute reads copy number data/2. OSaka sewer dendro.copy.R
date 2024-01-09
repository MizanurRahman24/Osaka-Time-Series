####Osaka sewer Timeseries data####
#####Kazuaki Matsui @ Rahman Mizanur####
############Change the directory according to the file location"##############
setwd("R:/Osaka Timeseries Master/V 1.5 copy num and reads to otu/3. Absolute reads copy number data/OTUClusteringOut")


######Load library
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(OTUtable)
library(reshape2)
library(phyloseq)
library(ape)
library(phyloseq)
library(ggplot2)

######Rdata preparation
# abundance table ????????????
Abundance_table.1 <- read.csv("otu_table.csv",header=T)
rownames(Abundance_table.1) <- Abundance_table.1$X
Abundance_table.1 <- Abundance_table.1[,-1]
Abundance_table <- Abundance_table.1[,-c(1,29,51,689,1376,2593,2765,2826)]
Abundance_table <- Abundance_table[-c(53:56),]
Abundance_table <- as.matrix(Abundance_table)

# taxa table
Taxa_table.1 <- read.csv("taxa_table.csv")
rownames(Taxa_table.1) <- Taxa_table.1$X
Taxa_table.1 <- Taxa_table.1[,-1]
Taxa_table<- Taxa_table.1[-c(1,29,51,689,1376,2593,2765,2826),]
Taxa_table <- as.matrix(Taxa_table)

# load sample ????????????
Sample_info <- read.csv("sample_data.csv")
rownames(Sample_info) <- Sample_info$X
Sample_info <- Sample_info[,-1]
Sample_info<- Sample_info[,-c(1,2)]
Sample_info<- Sample_info[-c(1:80,133:136),]
rownames(Sample_info) <- Sample_info$Sample_Name2
Sample_info <- Sample_info[,-1]

#####################??????#############################
OTU=otu_table(Abundance_table, taxa_are_rows = FALSE)
TAX=tax_table(Taxa_table)
SAMPLE=sample_data(Sample_info)

#?????????????????????????????????
physeq=phyloseq(OTU,TAX,SAMPLE)
physeq
saveRDS(physeq, "physeq.RData")

sewage_object <- readRDS("physeq.RData")
sewage_object
####Sewage relative abundance
sewage_object_relabun <- transform_sample_counts(sewage_object, function(x) x/sum(x))
# extract relative abundance info
sewage_relabun <- data.frame(sewage_object_relabun@otu_table@.Data)
sewage_relabun
# remove empty OTUs
sewage_relabun <- sewage_relabun[, colSums(sewage_relabun) > 0]

# load sample info
Sample_info <- read.csv("sample_data.csv")

Sample_info <- Sample_info[,-2]
Sample_info<- Sample_info[,-c(2,3)]
Sample_info<- Sample_info[-c(1:80,133:136),]
rownames(Sample_info) <- Sample_info$X
Sample_info

# remove OTUs with less than 1% as their maximum relative abundance
maxtax <- apply(sewage_relabun, 2, max)
mintax <- names(which(maxtax < 0.01))
sewage_filt <- sewage_relabun[, -which(colnames(sewage_relabun) %in% mintax)]
dim(sewage_filt)
# [1]  52 89

#############################
### sewage dendro and heat map ###
#############################

# convert to z score  = (x-μ)/σ
sewage.z <- t(zscore(t(sewage_filt)))
sewage.z



# change to chronological order
Sample_info <- Sample_info[order(Sample_info$order), ]
Sample_info

sewage.z <- sewage.z[match(rownames(Sample_info), rownames(sewage.z)),]
identical(rownames(sewage.z), rownames(Sample_info))
######Taxonomy al extract
Taxa_table.1 <- read.csv("taxa_table.csv")
rownames(Taxa_table.1) <- Taxa_table.1$X
Taxa_table<- Taxa_table.1[-c(1,29,51,689,1376,2593,2765,2826),]

# add all tax
temp_sewage <- data.frame(t(sewage.z))
temp_sewage$X <- rownames(temp_sewage)
temp_sewage <- merge(Taxa_table, temp_sewage, by = "X")
rownames(temp_sewage) <- paste0(temp_sewage$Phylum, "__", temp_sewage$Genus, "__", temp_sewage$X)
genera <- as.character(temp_sewage$Genus)
temp_sewage <- as.matrix(temp_sewage[,9:ncol(temp_sewage)])
sewage.z <- t(temp_sewage)
sewage.z

# change names
rownames(sewage.z) <- paste(Sample_info$month, Sample_info$Description)


# euclidian distances on z score matrices
Sewer_OTU_dist <- vegdist(t(sewage.z), method = "euclidian")


# cluster OTUs based on euclidian distances
Sewer_OTU_clus <- hclust(Sewer_OTU_dist, method = "average")


# add taxa info
OTUs.df <- data.frame(dendro_order = 1:length(colnames(sewage.z)), OTU = colnames(sewage.z))
OTUs.df$Phylum <- sapply(strsplit(as.character(OTUs.df$OTU), "__"), `[`, 1)
OTUs.df$Genus <- sapply(strsplit(as.character(OTUs.df$OTU), "__"), `[`, 2)


# define colors: blue is low, red is high
heat_colors <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(299))


# extract order of dendro tree leaves
Sewer_OTU_clus_order <- Sewer_OTU_clus$labels[c(Sewer_OTU_clus$order)]
Sewer_OTU_clus_order

Sewer_OTU_clus_order.df <- data.frame(t(data.frame(strsplit(Sewer_OTU_clus_order, "__"))))
colnames(Sewer_OTU_clus_order.df) <- c("Phylum", "Genus", "OTU_ID")
Sewer_OTU_clus_order.df$OTU <- Sewer_OTU_clus_order
rownames(Sewer_OTU_clus_order.df) <- Sewer_OTU_clus_order.df$OTU_ID


# heatmap with seasonal OTU colors
png("heat.def.png", units = "in", width = 8, height = 11, res = 600)
heatmap(as.matrix(t(sewage.z)),
        Rowv = as.dendrogram(Sewer_OTU_clus),
        Colv = NA, margins = c(2, 4),
        col = heat_colors)
dev.off()

# make same heatmap in ggplot2, so i can actually modify it
Osaka_sewer.z.df <- data.frame(sewage.z)
Osaka_sewer.z.df$Sample_name <- rownames(Osaka_sewer.z.df)
Osaka_sewer.z.df <- Osaka_sewer.z.df[match(paste(Sample_info$month, Sample_info$Description), Osaka_sewer.z.df$Sample_name), ]
Osaka_sewer.z.df$month_order <- 1:nrow(Osaka_sewer.z.df)
Osaka_sewer.z.df.m <- melt(Osaka_sewer.z.df, id.vars = c("month_order", "Sample_name"), variable.name = "OTU", value.name = "Z")
Osaka_sewer.z.df.m$OTU_ID <- sapply(strsplit(as.character(Osaka_sewer.z.df.m$OTU), "__"), '[', 3)


# need to keep order of OTUs in dendrogram
Sewer_OTU_clus_order.df$dendro_order <- 1:nrow(Sewer_OTU_clus_order.df)
Osaka_sewer.z.df.m <- merge(Osaka_sewer.z.df.m, Sewer_OTU_clus_order.df[c(3,5)], by = "OTU_ID")
Osaka_sewer.z.df.m <- Osaka_sewer.z.df.m[order(Osaka_sewer.z.df.m$dendro_order), ]


# make OTUs an ordered factor
Osaka_sewer.z.df.m$OTU_ID <- factor(Osaka_sewer.z.df.m$OTU_ID, levels = unique(Osaka_sewer.z.df.m$OTU_ID))


# specify color breaks (roughly using quantiles)
heat_colors <- c(rev(brewer.pal(9, "Blues")), "white", brewer.pal(9, "YlOrRd")[1:4])
heat_values <- c(seq(0, 0.2, length.out = 9), seq(0.21, 0.4, length.out = 2), seq(0.41, 1, length.out = 2))


# custom y axis labels
sewer_labels <- as.character(Sewer_OTU_clus_order.df$Genus)
sewer_labels[sewer_labels == ""] <- "sp."
sewer_labels <- paste0(Sewer_OTU_clus_order.df$Phylum, " (", sewer_labels, ")")


### sewer dendro###
png("heat.png", units = "in", width = 10, height = 11, res = 600)

sewer_dendro <- 
  ggplot(Osaka_sewer.z.df.m, aes(x = month_order, y = OTU_ID, fill = Z)) +
  geom_tile() +
  scale_fill_gradientn(colors = heat_colors, values = heat_values, breaks = c(-2, 0, 2, 4, 6)) +
  theme_classic() +
  scale_y_discrete(labels = sewer_labels, position = "right") +
  scale_x_continuous(breaks = c(1, 12, 24, 36, 48, 60)) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 7, color = "black", face = "italic"),
        axis.title.y = element_text(size = 9, color = "black", face = "bold"),
        axis.title.x = element_text(size = 9, color = "black", face = "bold"),
        legend.text = element_text(size = 9, color = "black"),
        legend.title = element_text(size = 9, color = "black", face = "bold"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.y = element_blank(),
        panel.border = element_blank()) +
  guides(fill = guide_colorbar(barheight = 2, barwidth = 0.3, units = "in")) +
  labs(x = "Weeks", y = "Sewer-associated OTUs\nPhylum (genus)", fill = "Z score")
sewer_dendro
dev.off()
ggsave("sewer heatmap.pdf", plot = sewer_dendro, device = "pdf", width = 4.5, height = 3.5, units = "in")

#ggsave("./Plots/sewer heatmap.pdf", plot = sewer_dendro, device = "pdf", width = 4.5, height = 3.5, units = "in")

