
###Osaka sewer timeseris data####
####Kazuaki Matsui @ Rahman Mizanur####
######change the working directory according to the location of the file########
setwd("R:/Osaka Timeseries Master/V 1.5 copy num and reads to otu/3. Absolute reads copy number data/OTUClusteringOut")
######Rdata preparation
#######Load Library
library(ape)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(OTUtable)
library(reshape2)
library(ape)
library(ggrepel)
library(indicspecies)

#######Load Data
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

# load sample info????????????
Sample_info <- read.csv("sample_data.csv")
rownames(Sample_info) <- Sample_info$X
Sample_info <- Sample_info[,-1]
Sample_info<- Sample_info[,-c(1,2)]
Sample_info<- Sample_info[-c(1:80,133:136),]
rownames(Sample_info) <- Sample_info$Sample_Name2
Sample_info <- Sample_info[,-1]

#####################??????#############################

#????????????phyloseq???????????????????????????
OTU=otu_table(Abundance_table, taxa_are_rows = FALSE)
TAX=tax_table(Taxa_table)
SAMPLE=sample_data(Sample_info)

#?????????????????????????????????
physeq=phyloseq(OTU,TAX,SAMPLE)
physeq
saveRDS(physeq, "physeq.RData")



################################################
### Analysis for sewage 2021-2022, Osaka Hirano
################################################


# load phyloseq object 
sewage_object <- readRDS("physeq.RData")
sewage_object

# convert to relative abundance
sewage_object_relabun <- transform_sample_counts(sewage_object, function(x) x/sum(x))

sewage_object_relabun
# extract relative abundance info
sewage_relabun <- data.frame(sewage_object_relabun@otu_table@.Data)

sewage_relabun
# remove empty ASVs
sewage_relabun <- sewage_relabun[, colSums(sewage_relabun) > 0]

# load sample info

Sample_info

Sample_info <- cbind(rownames(Sample_info), Sample_info)
rownames(Sample_info) <- NULL
colnames(Sample_info) <- c("Sample_name", "details", "month", "location", "Description", "I7_Index_ID", "index", "I5_Index_ID", "index2", "I7_Index_ID_.2", "index3", "I5_Index_ID_I_2", "index4")



# extract taxonomy
Taxonomy_all <- data.frame(sewage_object@tax_table@.Data)
Taxonomy_all$X <- rownames(Taxonomy_all)

Taxonomy_all
# add sample IDs to sample info
info_order <- c("S081", "S082", "S083", "S084", "S085", "S086", "S087", "S088", 
                "S089", "S090", "S091", "S092", "S093", "S094", "S095", "S096", "S097", "S098", "S099", "S100", 
                "S101", "S102", "S103", "S104", "S105", "S106", "S107", "S108", "S109", "S110", "S111", "S112", 
                "S113", "S114", "S115", "S116", "S117", "S118", "S119", "S120", "S121", "S122", "S123", 
                "S124", "S125", "S126", "S127", "S128", "S129", "S130", "S131", "S132")
Sample_info <- Sample_info[match(info_order, Sample_info$Sample_name),]
Sample_info$Sample_num <-  1:nrow(Sample_info)

#######################
### alpha diversity ###
#######################

# measure alpha diversity
shannon <- data.frame(Sample_name = rownames(sewage_relabun), shannon = diversity(sewage_relabun, "shannon"))


# add sample info
alpha.div <- merge(shannon, Sample_info, by = "Sample_name")
alpha.div


### figure ###
png("Alpha12.png", units = "in", width = 8, height = 5, res = 600)

plot_richness(sewage_object, x="month", measures="Shannon", color = "month")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

dev.off()

png("Alpha.location.png", units = "in", width = 8, height = 6.5, res = 600)
plot_richness(sewage_object, x="location", measures="Shannon", color = "location")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))
dev.off()

######################
### beta diversity ###
######################
#####Calculate axis weight
dist_bc <- as.matrix(vegdist(sewage_relabun, method = "bray"))
dist_bc
sewage.pcoa<- cmdscale(dist_bc, k=2, eig = TRUE, add = TRUE)
sewage.pcoa
positions<- sewage.pcoa$points
percent_explained <- 100*sewage.pcoa$eig / sum(sewage.pcoa$eig)
percent_explained[1:2]
#######[1] 23.50617 11.93962
# stat
sewage.pcoa <- pcoa(vegdist(sewage_relabun, method = "bray"))
sewage.pcoa
# extract values
sewage.pcoa.df <- data.frame(sewage.pcoa$vectors[,1:2])
sewage.pcoa.df

# add info
sewage.pcoa.df$Sample_name <- rownames(sewage.pcoa.df)
sewage.pcoa.df <- merge(sewage.pcoa.df, Sample_info, by = "Sample_name")

sewage.pcoa.df
### visualization  ###
png("PCoA.png", units = "in", width = 6, height = 4, res = 600)

months <-
  ggplot(sewage.pcoa.df, aes(x = Axis.2, y = Axis.1, color = month, shape = Description)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", size = 0.2) +
  geom_point(size = 4, alpha = 1.8) +
  theme_classic() +
  scale_color_manual(values = c("#5E4FA2", "#3288BD", "#66C2A5","#ABDDA4", "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142", "#4D004B"),
                     breaks = c("Jan","Feb", "Mar", "Apr",
                                "May", "Jun", "July", "Aug", "Sep",
                                "Oct", "Nov", "Dec")) +
  scale_shape_discrete(labels = c("direct", "pump")) +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.position = "left",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 2),
        panel.border = element_rect(color = "grey80", fill = NA, size = 1)) +
  guides(color = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1),
         shape = guide_legend(keyheight = 0.2, keywidth = 0.2, units = "in", ncol = 1)) +
  labs(y = "Axis 1\n(23.55%)", x = "Axis 2\n(11.95%)", color = "Sample\nmonth",
       shape = "sample\nlocation")
months
dev.off()
#ggsave("./Plots/pcoa.pdf", plot = pcoa, device = "pdf", width = 4.86, height = 2.6, units = "in")

#######THis is for further analysis to explain the figure######
# remove ASVs with less than 1% as their maximum relative abundance
maxtax <- apply(sewage_relabun, 2, max)
mintax <- names(which(maxtax < 0.01))
sewage_filt <- sewage_relabun[, -which(colnames(sewage_relabun) %in% mintax)]
dim(sewage_filt)
####[1] 52 89

#####cbind(Pos = data.frame(table(sewage.pcoa.df[sewage.pcoa.df$Axis.1 > 0, "month"])), 
   ######   Neg = data.frame(table(sewage.pcoa.df[sewage.pcoa.df$Axis.1 <= 0, "month"])))[-3]




################################
### determine indicator OTUs ###
################################
sewage_filt.df.t <- data.frame(t(sewage_filt))
sewage_filt.df.t$X <- rownames(sewage_filt.df.t)

# merge with tax info

sewage_filt.df.t <- merge(Taxonomy_all[-14], sewage_filt.df.t, by = "X")

# combine into one variable
sewage_filt.df.t$Names <- paste0(sewage_filt.df.t$Order, "__", sewage_filt.df.t$Family, "__", 
                                 sewage_filt.df.t$Genus, "__", sewage_filt.df.t$X)

# make that row names
rownames(sewage_filt.df.t) <- sewage_filt.df.t$Names

# transpose again, remove tax data
sewer_indic.1 <- data.frame(t(sewage_filt.df.t[,-c(1:8, ncol(sewage_filt.df.t))]))
sewer_indic <- sewer_indic.1[-c(53),]

# indic species analysis 3 month group
sewer_3mo <- multipatt(sewer_indic, as.vector(Sample_info$month), min.order = 3, max.order = 3, control = how(nperm = 999))
summary(sewer_3mo, indvalcomp = TRUE)

#########
