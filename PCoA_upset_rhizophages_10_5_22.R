setwd("C:/Users/Brandon Matsumoto/Documents/GitHub/ChickpeaGang/")

library(ape)
library(tidyverse)
library(vegan)
library(reshape2)
library(data.table)
library(ComplexUpset)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(tidyr)
library(janitor)
library(dplyr)
library(ggplot2)



#Load vOTU coverage table info and reformat

otu <- read.csv("combinedcoverm_w_paki_10_5_22.csv", header=T, stringsAsFactors = F, ) # load the table
otu <- otu %>%
  column_to_rownames("Contig") 
otu <- otu[rowSums(otu)> 0, ] # throw out rows (viruses) that don't appear in any sample

# Count vOTUs per sample
temp <- otu %>%
  mutate(across(.cols = everything(), ~ case_when(.x > 0 ~ 1, TRUE ~ 0)))
max(colSums(temp)) #26
which(colSums(temp) == 1)

otu_t <- otu %>% # Throw out samples (columns) that have no viruses
  t() %>% #transpose the table
  data.frame # turn it back from matrix to data frame

otu_t <- otu_t[rowSums(otu_t)> 0, ]

# Load SRA table and Greenlon sup. table 1
sra_run_table <- read.csv("SraRunTable.txt", header=T, stringsAsFactors = F) %>%
  select(c("Run", "BioSample"))
iso_to_loc <- read.csv("pnas.1900056116.sd01.csv", stringsAsFactors = F)[,1:38] %>% 
  dplyr::rename(BioSample = SRA_BioSample)
iso_to_loc <- iso_to_loc %>%
  left_join(sra_run_table, by="BioSample") %>%
  unique() %>%
  remove_rownames() %>%
  filter(!is.na(Run))

# Filter the vOTU table to sample that have a country assigned
otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run") %>%
  filter(Country.of.origin != "") %>% # Throwing out samples without strain specified
  filter(!is.na(Country.of.origin)) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")

#bar chart updated
#srr_filtered <- rownames(otu_t_valid) #pulling row names from out_t_valid
#srrdf <- as.data.frame(srr_filtered) #outputting rownames as dataframe to manipulate
#colnames(srrdf) <- c("Run")
#srr_smn <- merge(srrdf, sra_run_table,by="Run") #merged srrs that we have with their respective samn names
#srr_smn_geo <- merge(srr_smn, iso_to_loc, by="BioSample")
#bar_chart_data <- select(srr_smn_geo, c("Run.x","BioSample","isolation_source"))



# Create the distance matrix using the Jaccard algorithm (just presence/absence, not caring how high the coverage is)
pco_dist <- vegdist(otu_t_valid[, which(colnames(otu_t_valid) %like% "NZ")], method="jaccard")

pcoa_bray <- pcoa(pco_dist) # Creating the PCoA object
jac_variances <- data.frame(pcoa_bray$values$Relative_eig) %>% # Calculating the coordinates of each sample on a multidimensional PCoA plot
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

pcoa_bray_df <- data.frame(pcoa_bray$vectors) %>% # creating a data frame of the coordinates
  rownames_to_column(var = "Run") %>% # Creating a column from the rownames. The rownames are the SRR number
  mutate(Run = gsub("\\.Mean", "", Run)) %>% # Removing the .Mean from the end of the SR number
  data.frame %>%
  left_join(iso_to_loc, by="Run") %>% # Merging with the Greenlon sup. table 1 by Run (SRR...)
  filter(Country.of.origin != "") %>% # Throwing out samples without country of origin specified
  filter(!is.na(Country.of.origin)) %>%
  mutate(Country.of.origin = fct_relevel(Country.of.origin, c("AU", "CA", "US", "TU", "MR", "ET", "IN")))

eigenvalues<-round(jac_variances[,2], digits = 4)*100 # Rounding the variance explained by each PCoA axis to 2 decimals

pcoa_A <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.1, y = Axis.2, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
  #print(ggplot(pcoa_A))######
#plot(pcoa_A)
ggsave("Rhizophages_PCoA_Co1_Co2.pdf")

# Plotting coordinates 2 and 3
pcoa_B <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.2, y = Axis.3, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 3 ',eigenvalues[3],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
#ggsave(pcoa_B, "Rhizophages_PCoA_Co2_Co3.pdf")
ggsave("Rhizophages_PCoA_Co2_Co3.pdf")

# CAP analysis - by country of origin
iso_to_loc_2 <- iso_to_loc  %>%
  filter(Run %in% rownames(otu_t_valid)) %>%
  column_to_rownames("Run") %>%
  mutate(Country.of.origin = fct_relevel(Country.of.origin, c("AU", "CA", "US", "TU", "MR", "ET", "IN")))
cps <- capscale(pco_dist ~ Country.of.origin, data=iso_to_loc_2)
cps_gg <- scores(cps, tidy=T)
iso_to_loc_2 <- iso_to_loc_2 %>%
  rownames_to_column("Run")
cps_gg_sites <- cps_gg %>%
  filter(score == "sites") %>%
  rename(Run = label) %>%
  left_join(iso_to_loc_2[,c("Run", "Country.of.origin")], by="Run") %>%
  rename(label = Run)
cps_gg_cent <- cps_gg %>%
  filter(score == "centroids") %>%
  mutate(label = gsub(".*origin", "", label)) %>%
  mutate(label = fct_relevel(label, c("AU", "CA", "US", "TU", "MR", "ET", "IN")))
ggplot(data=cps_gg_sites, aes(x=CAP1, y=CAP2)) + 
  geom_point(aes(color=Country.of.origin)) + 
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)

# Coloring the PCoA by Mesorhizobium strain

# Filter the vOTU table to sample that have a strain assigned
otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run") %>%
  filter(ANI95.OTU != "") %>% # Throwing out samples without strain specified
  filter(!is.na(ANI95.OTU)) %>%
  filter(!is.na(Run)) %>%
  unique() %>%
  remove_rownames() %>%
  column_to_rownames("Run")

# Create the distance matrix using the Jaccard algorithm (just presence/absence, not caring how high the coverage is)
pco_dist <- vegdist(otu_t_valid[, which(colnames(otu_t_valid) %like% "NZ")], method="jaccard")

pcoa_bray <- pcoa(pco_dist) # Creating the PCoA object
jac_variances <- data.frame(pcoa_bray$values$Relative_eig) %>% # Calculating the coordinates of each sample on a multidimensional PCoA plot
  rownames_to_column(var = "PCaxis") %>% 
  data.frame

pcoa_bray_df <- data.frame(pcoa_bray$vectors) %>% # creating a data frame of the coordinates
  rownames_to_column(var = "Run") %>% # Creating a column from the rownames. The rownames are the SRR number
  mutate(Run = gsub("\\.Mean", "", Run)) %>% # Removing the .Mean from the end of the SR number
  data.frame %>%
  left_join(sra_run_table, by="Run") %>% # Merging with the SRR to BioSample table by SRR
  left_join(iso_to_loc, by="BioSample") %>% # Merging with the Greenlon sup. table 1 by BioSample
  filter(ANI95.OTU != "") %>% # Throwing out samples without country of origin specified
  filter(!is.na(ANI95.OTU)) %>%
  mutate(ANI95.OTU = fct_relevel(ANI95.OTU, c("1A", "1B", "1C", "1D", "1E", "1F", "1G", 
                                   "2A", "2B", "2C", "3A", "4A", "4B", 
                                   "5A", "5B", "5C", "5D", "6A", 
                                   "7A", "7B", "7D", 
                                   "8A", "8B", "9A", "10A")))

eigenvalues<-round(jac_variances[,2], digits = 4)*100 # Rounding the variance explained by each PCoA axis to 2 decimals
pcoa_C <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.1, y = Axis.2, color=ANI95.OTU), size = 1) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by strain') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_manual(values=c(brewer.pal(n=8, "Reds")[-1], brewer.pal(n=7, "YlOrRd")[-1], brewer.pal(n=6, "Greens")[-1], rev(brewer.pal(n=4, "Purples")[-1]), rev(brewer.pal(n=4, "Blues")[-1]), "black"))
ggsave("Rhizophages_PCoA_by_strain_Co1_Co2.pdf")

pcoa_D <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.2, y = Axis.3, color=ANI95.OTU), size = 1) +
  ylab(paste0('Co 3 ',eigenvalues[3],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by strain') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_manual(values=c(brewer.pal(n=8, "Reds")[-1], brewer.pal(n=7, "YlOrRd")[-1], brewer.pal(n=6, "Greens")[-1], rev(brewer.pal(n=4, "Purples")[-1]), rev(brewer.pal(n=4, "Blues")[-1]), "black"))
ggsave("Rhizophages_PCoA_by_strain_Co2_Co3.pdf")

ggarrange(pcoa_A, pcoa_B, pcoa_C, pcoa_D, ncol = 2, nrow = 2, labels=c("A", "B", "C", "D"), widths=c(1,1), heights=c(1,1))

# CAP analysis - by Mesorhizobium strain
iso_to_loc_2 <- iso_to_loc  %>%
  filter(Run %in% rownames(otu_t_valid)) %>%
  column_to_rownames("Run")
cps <- capscale(pco_dist ~ ANI95.OTU, data=iso_to_loc_2)
cps_gg <- scores(cps, tidy=T)
iso_to_loc_2 <- iso_to_loc_2 %>%
  rownames_to_column("Run") %>%
  mutate(ANI95.OTU = fct_relevel(ANI95.OTU, c("1A", "1B", "1C", "1D", "1E", "1F", "1G", 
                                              "2A", "2B", "2C", "3A", "4A", "4B", 
                                              "5A", "5B", "5C", "5D", "6A", 
                                              "7A", "7B", "7D", 
                                              "8A", "8B", "9A", "10A")))
cps_gg_sites <- cps_gg %>%
  filter(score == "sites") %>%
  rename(Run = label) %>%
  left_join(iso_to_loc_2[,c("Run", "ANI95.OTU")], by="Run") %>%
  rename(label = Run)
cps_gg_cent <- cps_gg %>%
  filter(score == "centroids")
ggplot(data=cps_gg_sites, aes(x=CAP1, y=CAP2)) + 
  geom_point(aes(color=ANI95.OTU)) + 
  theme_bw() + 
  scale_color_manual(values=c(brewer.pal(n=8, "Reds")[-1], brewer.pal(n=7, "YlOrRd")[-1], brewer.pal(n=6, "Greens")[-1], rev(brewer.pal(n=4, "Purples")[-1]), rev(brewer.pal(n=4, "Blues")[-1]), "black"))

# Creating the upset plot requires summing up coverage per country per virus, then converting the table to a presence/absence
# table, just like a Jaccard matrix. The way to do that is to change every value >0 to 1.

max.loc <- as.data.frame(otu_t) %>%
  rownames_to_column(var="Run") %>%
  mutate(Run = gsub("\\.Mean", "", Run))

# Reformat coverage and sum by country
max.loc <- merge(max.loc, sra_run_table, all.x=T) # Merging with the SRR to BioSample table by SRR
max.loc <- merge(max.loc, iso_to_loc[, c("BioSample", "Country.of.origin")], by="BioSample", all.x=T) %>% # Merging with the Greenlon sup. table 1 by BioSample
  reshape2::melt() %>% # converting from wide to long format
  mutate(Country.of.origin = replace(Country.of.origin, Run %like% "^P6", "PK")) %>% # Adding the Pakistan samples
  mutate(Country.of.origin = replace(Country.of.origin, is.na(Country.of.origin), "Unk")) %>% # Setting no country to unknown
  mutate(Country.of.origin = replace(Country.of.origin, Country.of.origin=="", "Unk")) %>% # Setting no country to unknown
  dplyr::select(c("Country.of.origin", "variable", "value")) %>% # Keeping only these three columns: country, vOTU and coverage
  group_by(variable, Country.of.origin) %>% # grouping by vOTU and country
  summarize_at("value", sum) %>% # Summing up coverage per vOTU per country
  reshape2::dcast(variable~Country.of.origin) # Converting back to wide format

# Reformat max.loc to presence/absence rather than coverage and remove samples with unknown country
max.loc[max.loc > 0] <- 1
max.loc <- max.loc[,-which(colnames(max.loc)=="Unk")]

# Add the lysogen/not lysogen info to a complex upset plot
lys <- scan("lysogenic_putative.tsv", what=character()) # Read the list of lysogenic viruses
lys <- gsub("_NZ_.*", "", lys) # DRAM duplicates the header. Fix that
lys <- gsub("\\.", "_", lys)

max.loc$lys <- FALSE # Create a column indicating whether a vOTU is lysogenic and setting it to FALSE as a default
max.loc$variable <- as.character(max.loc$variable) 
max.loc$variable <- gsub("\\.", "_", max.loc$variable)

max.loc$lys[max.loc$variable %in% lys] <- TRUE

sets=c("US", "CA", "AU", "IN", "PK", "ET", "MR", "TU")
upset(
  max.loc, min_degree=1,
  sets, sort_sets="descending", sort_intersections="descending",
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE,
      mapping=aes(fill=lys)
    ) + scale_fill_brewer(palette="Paired") + 
      ylab("Group size")
  ),
  width_ratio=0.2
)
ggsave("rhizophages_upset_complex.pdf")

# Figure out in how many Cicer hosts each vOTU appears
max.loc <- as.data.frame(otu_t) %>%
  rownames_to_column(var="Run") %>%
  mutate(Run = gsub("\\.Mean", "", Run))

# Reformat coverage and sum by country
max.loc <- merge(max.loc, sra_run_table, all.x=T) # Merging with the SRR to BioSample table by SRR
max.loc <- merge(max.loc, iso_to_loc[, c("BioSample", "Cicer.host")], by="BioSample", all.x=T) %>% # Merging with the Greenlon sup. table 1 by BioSample
  reshape2::melt() %>% # converting from wide to long format
  mutate(Cicer.host = replace(Cicer.host, is.na(Cicer.host), "Unk")) %>% # Setting no host plant to unknown
  mutate(Cicer.host = replace(Cicer.host, Cicer.host=="", "Unk")) %>% # Setting no host plant to unknown
  dplyr::select(c("Cicer.host", "variable", "value")) %>% # Keeping only these three columns: host plant, vOTU and coverage
  group_by(variable, Cicer.host) %>% # grouping by vOTU and host plant
  summarize_at("value", sum) %>% # Summing up coverage per vOTU per host plant
  reshape2::dcast(variable~Cicer.host) # Converting back to wide format

# Reformat max.loc to presence/absence rather than coverage and remove samples with unknown host plant
max.loc[max.loc > 0] <- 1
max.loc <- max.loc[,-which(colnames(max.loc)=="Unk")]
max.loc$numhosts <- rowSums(max.loc[,-1])
temp <- data.frame(c(nrow(max.loc[max.loc$numhosts==1,]), nrow(max.loc[max.loc$numhosts==2,]), nrow(max.loc[max.loc$numhosts==3,])))
colnames(temp)[1] <- "cicer_num"
temp$hosts <- c("one species", "two species", "three species")
ggbarplot(data=temp, x="hosts", y="cicer_num", fill="darkgreen") + xlab("Number of Cicer species") + 
  ylab("Number of vOTUs present")

sets=c("Ca", "Cr", "Ce")
temp <- max.loc[, -c(1, 5)]
upset(
  temp, min_degree=1,
  sets, sort_sets="descending", sort_intersections="descending",
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=TRUE) + 
      ylab("Number of vOTUs")
  ),
  width_ratio=0.2
)

# PERMANOVA on metadata
otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(sra_run_table, by="Run") %>%
  left_join(iso_to_loc, by="Run") %>%
  data.frame() %>%
  unique() %>%
  filter(!is.na("Run")) %>%
  column_to_rownames("Run")

# Adonis2 analysis of independent variables
ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  select(c(3, 9, 10:12, 16:21, 26, 30, 35, 38, 102))

# PERMANOVA requires no NAs, so choosing only the 275 complete cases
ad_comp <- ad[complete.cases(ad),]

pcoa_bray_df_num <- data.frame(pcoa_bray$vectors) 
pcoa_bray_df_num <- pcoa_bray_df_num[rownames(pcoa_bray_df_num) %in% ad_comp$Run, ]

ad_comp <- ad_comp %>%
  filter(Run %in% rownames(pcoa_bray_df_num))

# Leaves us with 94 samples, nothing is significant and total variance explained is only 32%
adonis2(pcoa_bray_df_num ~ Country.of.origin+Mean_annual_temp+Mean_annual_precipitation+soil.ph+soil.genus+soil.taxnwrb.group+altitude+ANI95.OTU, data=ad_comp, method="euclidian")

# Trying again with the parameters that were significant in Alex's paper
ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  select(c("Run", "soil.genus", "lat_lon", "Mean_annual_precipitation", "Mean_annual_temp", "soil.ph")) %>%
  mutate(lat_lon = gsub(" .*", "", lat_lon)) %>%
  mutate(lat_lon = as.numeric(lat_lon))

# PERMANOVA requires no NAs, so choosing only the 567 complete cases
ad_comp <- ad[complete.cases(ad),]

pcoa_bray_df_num <- data.frame(pcoa_bray$vectors) 
pcoa_bray_df_num <- pcoa_bray_df_num[rownames(pcoa_bray_df_num) %in% ad_comp$Run, ]

ad_comp <- ad_comp %>%
  filter(Run %in% rownames(pcoa_bray_df_num)) %>%
  unique()

# Leaves us with 144 samples, nothing is significant but latitude has a R^2=2% and p=0.055, soil pH has an R^2 of 1% and p=0.092
adonis2(pcoa_bray_df_num ~ soil.genus+lat_lon+Mean_annual_precipitation+Mean_annual_temp+soil.ph, data=ad_comp, method="euclidean")

# Trying again with only ANI
ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  select(c(3, 9, 102))

# PERMANOVA requires no NAs, so choosing only the 1089 complete cases
ad_comp <- ad[complete.cases(ad),]

pcoa_bray_df_num <- data.frame(pcoa_bray$vectors) 
pcoa_bray_df_num <- pcoa_bray_df_num[rownames(pcoa_bray_df_num) %in% ad_comp$Run, ]

ad_comp <- ad_comp %>%
  filter(Run %in% rownames(pcoa_bray_df_num)) %>%
  unique()

# Leaves us with 407 samples, total variance explained is only 7% p=0.003 999 permutations
adonis2(pcoa_bray_df_num ~ factor(ANI95.OTU), data=ad_comp, method="euclidean")

# Trying again with only latitude
ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  select(c("Run", "lat_lon")) %>%
  mutate(lat_lon = gsub(" .*", "", lat_lon)) %>%
  mutate(lat_lon = as.numeric(lat_lon))
ad$lat_lon[ad$lat_lon<10] = 10
ad$lat_lon[ad$lat_lon>=10 & ad$lat_lon<12] = 12
ad$lat_lon[ad$lat_lon>=12 & ad$lat_lon<14] = 14
ad$lat_lon[ad$lat_lon>=14 & ad$lat_lon<16] = 16
ad$lat_lon[ad$lat_lon>=16 & ad$lat_lon<18] = 18
ad$lat_lon[ad$lat_lon>=22 & ad$lat_lon<24] = 24
ad$lat_lon[ad$lat_lon>=24 & ad$lat_lon<26] = 26
ad$lat_lon[ad$lat_lon>=26 & ad$lat_lon<28] = 28
ad$lat_lon[ad$lat_lon>=30 & ad$lat_lon<32] = 32
ad$lat_lon[ad$lat_lon>=32 & ad$lat_lon<34] = 34
ad$lat_lon[ad$lat_lon>=34] = 36

# PERMANOVA requires no NAs, so choosing only the 845 complete cases
ad_comp <- ad[complete.cases(ad),]

pcoa_bray_df_num <- data.frame(pcoa_bray$vectors) 
pcoa_bray_df_num <- pcoa_bray_df_num[rownames(pcoa_bray_df_num) %in% ad_comp$Run, ]

ad_comp <- ad_comp %>%
  filter(Run %in% rownames(pcoa_bray_df_num)) %>%
  unique()

# Leaves us with 209 samples, total variance explained is only 2% p=0.001 999 permutations

adonis2(pcoa_bray_df_num ~ factor(lat_lon), data=ad_comp, method="euclidean")

# Trying directly with the Jaccard matrix
ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  select(c(3, 9, 102))

# PERMANOVA requires no NAs, so choosing only the 1089 complete cases
ad_comp <- ad[complete.cases(ad),]

pcoa_bray_df_num <- data.frame(pcoa_bray$vectors) 
pcoa_bray_df_num <- pcoa_bray_df_num[rownames(pcoa_bray_df_num) %in% ad_comp$Run, ]

otu_t_valid <- otu_t_valid %>%
  rownames_to_column("Run") %>%
  filter(Run %in% ad_comp$Run) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")

ad_comp <- ad_comp %>%
  filter(Run %in% rownames(otu_t_valid)) %>%
  unique()

pco_dist <- vegdist(otu_t_valid[, which(colnames(otu_t_valid) %like% "NZ")], method="jaccard")

# Leaves us with 209 samples p=0.3 999 permutations
adonis2(pco_dist ~ factor(ANI95.OTU), data=ad_comp, method="euclidean")

# Trying a Mantel test

ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  dplyr::select(Run, soil.genus, lat_lon, Mean_annual_precipitation, Mean_annual_temp, soil.ph) %>%
  mutate(lat_lon = gsub(" .*", "", lat_lon)) %>%
  mutate(lat_lon = as.numeric(lat_lon)) %>%
  filter(!is.na(Run)) %>%
  unique() %>%
  filter(Run %in% rownames(otu_t)) %>%
  remove_rownames() %>%
  column_to_rownames("Run") %>%
  filter_all(any_vars(!is.na(.))) %>%
  mutate(soil.genus.dummy = case_when(soil.genus == "Leptosols" ~ 1, 
                                      soil.genus == "Vertisols" ~ 2, 
                                      soil.genus == "Cambisols" ~3, 
                                      soil.genus == "Luvisols" ~ 4, 
                                      soil.genus == "Phaeozems" ~ 5, 
                                      soil.genus == "Regosols" ~ 6, 
                                      is.na(soil.genus) ~ 7)) %>%
  dplyr::select(-soil.genus)


otu_t_valid <- otu_t[rownames(otu_t) %in% rownames(ad),]

otu_dist <- vegdist(otu_t_valid, method="jaccard")
ad_dist <- vegdist(ad, method="bray", na.rm = T)

mantel(otu_dist, ad_dist)
# r=0.0375, p=0.012, num_perm=999
mant <- mantel.correlog(otu_dist, ad_dist)
mant_val <- mant$mantel.res[1:8,] %>% data.frame() %>% mutate(Pr.bin = ifelse(Pr.corrected. >= 0.1, "sig", "not_sig"))
ggplot(data=mant_val, aes(x=class.index, y=Mantel.cor)) + geom_point(aes(fill=Pr.bin), color="black", pch=21, size=2) + theme_classic() + geom_line() + scale_fill_manual(values=c("darkblue", "lightblue"))

# Add distance decay

require("geosphere")
require("geodist")
ad <- read.csv("pnas.1900056116.sd01.csv") %>%
  rename(BioSample = SRA_BioSample) %>%
  left_join(sra_run_table, by="BioSample") %>%
  dplyr::select(Run, lat_lon) %>%
  filter(!is.na(Run)) %>%
  filter(lat_lon != "not collected") %>%
  filter(lat_lon != "") %>%
  filter(lat_lon != "missing") %>%
  filter(Run %in% rownames(otu_t_valid)) %>%
  separate(lat_lon, c("lat", "lat_dir", "lon", "lon_dir"), sep=" ") %>%
  mutate(lat = as.numeric(lat)) %>%
  mutate(lon = as.numeric(lon))
ad$lon_cor <- ad$lon
ad$lon_cor[ad$lon_dir == "W"] <- ad$lon_cor[ad$lon_dir == "W"]*-1
ad_for_dist <- ad %>% dplyr::select(Run, lat, lon_cor) %>% unique() %>% remove_rownames %>% column_to_rownames("Run")

geo_dist <- geodist(ad_for_dist, measure="vincenty")

mantel.correlog(otu_dist, geo_dist)

# Trying to figure out which phages are in the symbiotic island based on the pacbio genomes
blst <- read.delim("vOTUs_confirmed_10k_vs_pacbio_blast_minqcov90.tsv", sep="\t", header=F, stringsAsFactors = F)
colnames(blst) <- c("qseqid","Refseq_accession","pident","qcovs","qstart","qend","sstart","send","length","evalue","bitscore","mismatch")
pacbio <- read.csv("Greenlon_PNAS_supp/pnas.1900056116.sd07.csv", header=T) %>%
  mutate(Refseq_accession = gsub("NZ_", "", Refseq_accession))
mrg <- pacbio %>%
  left_join(blst, by="Refseq_accession") %>%
  filter(!is.na(qseqid)) %>%
  filter(sstart >= sym.island.start) %>%
  filter(send <= sym.island.end)
unique(mrg$qseqid)

for (p in unique(mrg$qseqid)) {
   print(max.loc[max.loc$variable %like% p,])
}

# Identify AMGs in valid OTUs
orgN <- read.delim("DRAM/Distill/metabolism_summary_orgN.csv", header=T, stringsAsFactors = F, sep=";")
utilC <- read.delim("DRAM/Distill/metabolism_summary_C.csv", header=T, stringsAsFactors = F, sep=";")
misc <- read.delim("DRAM/Distill/metabolism_summary_MISC.csv", header=T, stringsAsFactors = F, sep=";")
transp <- read.delim("DRAM/Distill/metabolism_summary_transp.csv", header=T, stringsAsFactors = F, sep=";")
ener <- read.delim("DRAM/Distill/metabolism_summary_energy.csv", header=T, stringsAsFactors = F, sep=";")

otus <- sub("\\.", "_", sub("\\.", "_", colnames(otu_t)))
orgN_valid <- orgN %>%
  select(otus) %>%
  add_column(orgN[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
apply(orgN_valid, 1, function (x) {length(which(x>0))}) # Row 12 has 28 vOTUs with this gene: Catlytic type: Serine; processes the coat protein at the lysyl bond; involved in the maturation of the phage prohead
# Oh well, makes sense

utilC_valid <- utilC %>%
  select(otus) %>%
  add_column(utilC[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# Most common are GH24 lysozume and GH108 N-acetylmuramidase
# However, there are 4 chitinases:
colnames(utilC_valid)[which(utilC_valid[1,]>0)]

misc_valid <- misc %>%
  select(otus) %>%
  add_column(misc[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# Not interesting. All DNA related genes.

ener_valid <- ener %>%
  select(otus) %>%
  add_column(ener[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# 2 genes: reductive citratce cycle and methanogenesis

transp_valid <- transp %>%
  select(otus) %>%
  add_column(transp[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# Branched amino acid transport and peptide/nickel transport

shared_votu <- read_excel("3 shared votus.xlsx")
shared_votu <- t(shared_votu)

shared_votu_new <- shared_votu
colnames(shared_votu_new) <- shared_votu[1,]
shared_votu_new

#bar chart updated
srr <- out_t_valid[]

#bar chart for ubiquotus votu
srr <- read.csv("samnsrrgeowehave.csv") #need to change
all_data <- read.csv("pnas.1900056116.sd01.csv")
srr_samn_merge <- merge(shared_votu_new,srr, by="SRR")
all_merge <- merge(srr_samn_merge, all_data, by="SAMN")
stuff_I_want <- subset(all_merge, select = c("SAMN","SRR","NZ_SRUU01000009.1_Mesorhizobium_sp._M1C.F.Ca.ET.210.01.1.1_NODE_9_length_221169_cov_20.712__whole_genome_shotgun_sequence_fragment_2","NZ_SADM02000001.1_Mesorhizobium_sp._M2E.F.Ca.ET.166.01.1.1_NODE_1_length_2175095_cov_17.3373__whole_genome_shotgun_sequence_fragment_4","NZ_RZPU01000009.1_Mesorhizobium_sp._M7A.F.Ca.CA.004.02.1.1_NODE_9_length_43211_cov_3.77883__whole_genome_shotgun_sequence","isolation_source"))
stuff_I_want[3:5]<-apply(stuff_I_want[3:5],2,function(x) ifelse(x > 0,1,x))
#grabbing the first votu and all its data
votu1 <- subset(stuff_I_want, select = c("NZ_SRUU01000009.1_Mesorhizobium_sp._M1C.F.Ca.ET.210.01.1.1_NODE_9_length_221169_cov_20.712__whole_genome_shotgun_sequence_fragment_2","isolation_source"))
votu1[votu1==0] <- NA
votu1 <- votu1[complete.cases(votu1),]
votua <- votu1 %>% group_by(isolation_source) %>% tally()
#grabbing the second votu and all its data
votu2 <- subset(stuff_I_want, select = c("NZ_SADM02000001.1_Mesorhizobium_sp._M2E.F.Ca.ET.166.01.1.1_NODE_1_length_2175095_cov_17.3373__whole_genome_shotgun_sequence_fragment_4","isolation_source"))
votu2[votu2==0] <- NA
votu2 <- votu1[complete.cases(votu2),]
votub <- votu2 %>% group_by(isolation_source) %>% tally()
#grabbing the third votu and all its data
votu3 <- subset(stuff_I_want, select = c("NZ_RZPU01000009.1_Mesorhizobium_sp._M7A.F.Ca.CA.004.02.1.1_NODE_9_length_43211_cov_3.77883__whole_genome_shotgun_sequence","isolation_source"))
votu3[votu3==0] <- NA
votu3 <- votu1[complete.cases(votu3),]
votu3 %>% group_by(isolation_source) %>% tally()
votuc <- data.frame((votu3 %>% group_by(isolation_source) %>% tally()))
#coming all thre votu columns and data
combined_ab <- merge(votua,votub,by=c("isolation_source"))
combined_abc <- merge(combined_ab,votuc,by=c("isolation_source"))
names(combined_abc)[2] <- "a"
names(combined_abc)[3] <- "b"
names(combined_abc)[4] <- "c"
#removing unwanted rows
remove_row <- combined_abc[-c(4),]
melt_data <- melt(remove_row, id.var = c("isolation_source"), variable.name="vOTU")
#plotting
p <- ggplot(melt_data, aes(fill=isolation_source,x=vOTU,y=value)) + scale_fill_brewer(palette="Green") +  geom_bar(position="stack", stat="identity",colour="black") + geom_text(aes(label=value),size = 3, hjust = 0.5, vjust = 1, position ="stack") + ggtitle('Number of samples present in the 3 globally shared vOTUs in regards to Chickpea strain') + xlab('vOTU') + ylab('Number of Samples') + theme(plot.title = element_text(face='bold',hjust = 0.5))
p



