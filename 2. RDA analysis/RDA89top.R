
setwd("R:/Osaka Timeseries Master/CCA analysis/Data")

###########Load packages
library(tidyverse)
library(psych)
library(adespatial)
library(vegan)

# Import data
otu.table = read.csv("otu_table.top89.csv", row.names = 1)
sample.data = read.csv("sample_data.csv", row.names = 1)
data.frame = read.csv("data.frame.csv", row.names = 1)

# Set seed
set.seed(123)

# Perform RDA with all variables
rda1 = rda(otu.table ~ ., data = sample.data, scale = TRUE)
rda1

# Model summaries
RsquareAdj(rda1) # adjusted Rsquared 
vif.cca(rda1) # variance inflation factor (<10 OK)
anova.cca(rda1, permutations = 1000) # full model
anova.cca(rda1, permutations = 1000, by="margin") # per variable 

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "constrained"))
screeplot(rda1)


# Visualise results of RDA
png("rda89.png", width = 8, height = 7, units = "in", res = 600)
plot(rda1, type="n", scaling = 3)
title("Redundancy analysis")
# SITES
points(rda1, display="sp", pch=21, scaling=3, cex=1.5, col="blue") # sites
# text(rda1, display="sites", scaling = 3, col="black", font=2, pos=4)

text(rda1, display="bp", scaling=3, col="red1", cex=1, lwd=2)

###################
adj.R2 = round(RsquareAdj(rda1)$adj.r.squared, 3)
mtext(bquote(italic("R")^"2"~"= "~.(adj.R2)), side = 3, adj = 0.5)
dev.off()


