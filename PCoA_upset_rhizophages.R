setwd("~/Dropbox/UCDavis/Brandon/")

library(ape)
library(tidyverse)
library(vegan)
library(tidyverse)
library(reshape2)
library(data.table)
library(ComplexUpset)
library(RColorBrewer)
library(ggpubr)

#Load vOTU coverage table info and reformat

otu <- read.delim("coverm_output_mincov0.5.tsv", header=T, stringsAsFactors = F, sep="\t") # load the table
otu <- otu %>%
  column_to_rownames("Contig")
colnames(otu) <- gsub("\\.RPKM", "", colnames(otu))
otu <- otu[rowSums(otu)> 0, ] # throw out rows (viruses) that don't appear in any sample

# Count vOTUs per sample
temp <- otu %>%
  mutate(across(.cols = everything(), ~ case_when(.x > 0 ~ 1, TRUE ~ 0)))
max(colSums(temp)) #16
which(colSums(temp) == 1)

otu_t <- otu %>% # Throw out samples (columns) that have no viruses
  t() %>% #transpose the table
  data.frame # turn it back from matrix to data frame

otu_t[otu_t > 0] <- 1 # Convert to presence/absence
otu_t <- otu_t[rowSums(otu_t)> 2, ] # Remove communities with 2 or less viruses, this get rid of 81% of samples

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
# Adding Pakistan samples
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6B11_S135_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6B12_S136_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6C11_S147_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6D11_S159_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6E10_S170_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6E11_S171_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6E12_S172_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6F10_S182_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6F11_S183_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6F12_S184_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6G12_S196_L008"
iso_to_loc[nrow(iso_to_loc)+1, "sample_title"] <- "P6H10_S206_L008"
iso_to_loc[iso_to_loc$sample_title %like% "P6.*", "Country.of.origin"] <- "PK"
iso_to_loc[iso_to_loc$sample_title %like% "P6.*", "organism"] <- "Mesorhizobium"
iso_to_loc[(nrow(iso_to_loc)-11):nrow(iso_to_loc), "Run"] <- iso_to_loc[(nrow(iso_to_loc)-11):nrow(iso_to_loc), "sample_title"]

# Filter the vOTU table to sample that have a country assigned
otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run") %>%
  filter(Country.of.origin != "") %>% # Throwing out samples without strain specified
  filter(!is.na(Country.of.origin)) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")

# Create the distance matrix using the Jaccard algorithm (just presence/absence, not caring how high the coverage is)
# the reason we're using a Jaccard matrix is that the matrix is mostly 0's which can heavily bias other metrics
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
  mutate(Country.of.origin = fct_relevel(Country.of.origin, c("AU", "CA", "US", "TU", "MR", "ET", "IN", "PK")))

eigenvalues<-round(jac_variances[,2], digits = 4)*100 # Rounding the variance explained by each PCoA axis to 2 decimals

pcoa_A <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.1, y = Axis.2, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
ggsave("new_figs/Rhizophages_PCoA_Co1_Co2_min2.pdf", plot=pcoa_A)

# Plotting coordinates 2 and 3
pcoa_B <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.2, y = Axis.3, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 3 ',eigenvalues[3],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
ggsave("new_figs/Rhizophages_PCoA_Co2_Co3_min3.pdf.pdf", pcoa_B)

# CAP analysis - by country of origin
iso_to_loc_2 <- iso_to_loc  %>%
  filter(Run %in% rownames(otu_t_valid)) %>%
  column_to_rownames("Run") %>%
  mutate(Country.of.origin = fct_relevel(Country.of.origin, c("AU", "CA", "US", "TU", "MR", "ET", "IN", "PK")))
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
  mutate(label = fct_relevel(label, c("AU", "CA", "US", "TU", "MR", "ET", "IN", "PK")))
ggplot(data=cps_gg_sites, aes(x=CAP1, y=CAP2)) + 
  geom_point(aes(color=Country.of.origin)) + 
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
ggsave("new_figs/CAP_analysis_by_country_min2.pdf")

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
  mutate(ANI95.OTU = fct_relevel(ANI95.OTU, c("1A", "1B", "1C", "1D", "1E", "1G", 
                                   "2A", "2B", "2C", "2E", "3A", "4A", "4B", 
                                   "5A", "5B", "5C", "5D", "5E", "6A", 
                                   "7A", "7B", "7D", 
                                   "8A", "8B", "9A", "10A")))

eigenvalues<-round(jac_variances[,2], digits = 4)*100 # Rounding the variance explained by each PCoA axis to 2 decimals
pcoa_C <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.1, y = Axis.2, color=ANI95.OTU), size = 1) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by strain') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(palette="Spectral")
  #theme_bw() + scale_color_manual(values=c(brewer.pal(n=9, "Reds")[-1], brewer.pal(n=7, "YlOrRd")[-1], brewer.pal(n=6, "Greens")[-1], rev(brewer.pal(n=4, "Purples")[-1]), rev(brewer.pal(n=4, "Blues")[-1]), "black"))
ggsave("new_figs/Rhizophages_PCoA_by_strain_Co1_Co2_min2.pdf", pcoa_C)

pcoa_D <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.2, y = Axis.3, color=ANI95.OTU), size = 1) +
  ylab(paste0('Co 3 ',eigenvalues[3],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by strain') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(palette="Spectral")
  #theme_bw() + scale_color_manual(values=c(brewer.pal(n=9, "Reds")[-1], brewer.pal(n=7, "YlOrRd")[-1], brewer.pal(n=6, "Greens")[-1], rev(brewer.pal(n=4, "Purples")[-1]), rev(brewer.pal(n=4, "Blues")[-1]), "black"))
ggsave("new_figs/Rhizophages_PCoA_by_strain_Co2_Co3_min3.pdf", pcoa_D)

ggarrange(pcoa_A, pcoa_B, pcoa_C, pcoa_D, ncol = 2, nrow = 2, labels=c("A", "B", "C", "D"), widths=c(1,1), heights=c(1,1))
ggsave("new_figs/fig_2_PCoA_min2.pdf")

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
  theme_bw() + scale_color_brewer(palette="Spectral")
  #theme_bw() + 
  #scale_color_manual(values=c(brewer.pal(n=9, "Reds")[-1], brewer.pal(n=7, "YlOrRd")[-1], brewer.pal(n=6, "Greens")[-1], rev(brewer.pal(n=4, "Purples")[-1]), rev(brewer.pal(n=4, "Blues")[-1]), "black"))
ggsave("new_figs/CAP_analysis_by_strain_min2.pdf")

# Creating the upset plot requires summing up coverage per country per virus, then converting the table to a presence/absence
# table, just like a Jaccard matrix. The way to do that is to change every value >0 to 1.

otu_t_2 <- otu %>% # Throw out samples (columns) that have no viruses
  t() %>% #transpose the table
  data.frame # turn it back from matrix to data frame

otu_t_2[otu_t_2 > 0] <- 1 # Convert to presence/absence

max.loc <- as.data.frame(otu_t_2) %>%
  rownames_to_column(var="Run")

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
      ylab("Number of vOTUs")
  ),
  width_ratio=0.2
)
ggsave("new_figs/rhizophages_upset_complex_min2.pdf") # Note that Australia got thrown out, probably due to simple communities

# Which vOTUs are very common (in  or more countries)
max.loc[rowSums(max.loc[,-c(1, ncol(max.loc))]) >= 5,]
ubiq_vOTU <- max.loc$variable[rowSums(max.loc[,-c(1, ncol(max.loc))]) >= 5] # None when limiting to over 2 viruses per sample

# Figure out in how many Cicer hosts each vOTU appears
max.loc <- as.data.frame(otu_t) %>%
  rownames_to_column(var="Run")

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
ggsave("new_figs/Num_vOTUS_by_Cicer_min2.pdf")

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
ggsave("new_figs/upset_by_Cicer_min2.pdf")

# Which phages are in the symbiotic island based on the pacbio genomes - this is a relic of a previous iteration of the vOTU list.
# These vOTUs are no longer in the current list, hence they don't yield any results.
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
apply(orgN_valid, 1, function (x) {length(which(x>0))}) # Row 11 has 27 vOTUs with this gene: Catlytic type: Serine; processes the coat protein at the lysyl bond; involved in the maturation of the phage prohead
# Oh well, makes sense
# 16 vOTUs have a Catalytic type: Serine; cleaves.at an Ala-Gly or a Cys-Gly bond
# 15 vOTUs have a Catalytic type: Serine; hydrolyzes the N-blocked p-nitrophenyl esters of Gly, Ala, Phe, Val, Leu and Trp; prefers hydrophobic amino acids on either side of the scissile bond and will not cleave a peptide containing fewer than six amino acids; destruction of cleaved signal pepides in the periplasmic space

utilC_valid <- utilC %>%
  select(otus) %>%
  add_column(utilC[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# Most common are GH24 lysozume and GH108 N-acetylmuramidase
# However, there are 7 vOTUs with chitinases:
colnames(utilC_valid)[which(utilC_valid[1,]>0)]
colnames(utilC_valid)[which(utilC_valid[2,]>0)]

misc_valid <- misc %>%
  select(otus) %>%
  add_column(misc[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# Not interesting. All DNA related genes.

ener_valid <- ener %>%
  select(otus) %>%
  add_column(ener[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# 2 genes: reductive citrate cycle and methanogenesis

transp_valid <- transp %>%
  select(otus) %>%
  add_column(transp[,1:5], .before=1) %>%
  filter(rowSums(across(starts_with("NZ"))) > 0)
# Branched amino acid transport and peptide/nickel transport

# Accumulation curves for sampling effort effect by country

# Separate by country
ET.map <- filter(iso_to_loc, Country.of.origin == "ET")
ET.otu <- otu[,colnames(otu) %in% ET.map$Run]
ET.otu <- ET.otu[rowSums(ET.otu)>0,]

IN.map <- filter(iso_to_loc, Country.of.origin == "IN")
IN.otu <- otu[,colnames(otu) %in% IN.map$Run]
IN.otu <- IN.otu[rowSums(IN.otu)>0,]

PK.otu <- otu[,colnames(otu) %like% "P.*"]
PK.otu <- PK.otu[rowSums(PK.otu)>0,]

MR.map <- filter(iso_to_loc, Country.of.origin == "MR")
MR.otu <- otu[,colnames(otu) %in% MR.map$Run]
MR.otu <- MR.otu[rowSums(MR.otu)>0,]

TU.map <- filter(iso_to_loc, Country.of.origin == "TU")
TU.otu <- otu[,colnames(otu) %in% TU.map$Run]
TU.otu <- TU.otu[rowSums(TU.otu)>0,]

AU.map <- filter(iso_to_loc, Country.of.origin == "AU")
AU.otu <- otu[,colnames(otu) %in% AU.map$Run]
AU.otu <- AU.otu[rowSums(AU.otu)>0,]

CA.map <- filter(iso_to_loc, Country.of.origin == "CA")
CA.otu <- otu[,colnames(otu) %in% CA.map$Run]
CA.otu <- CA.otu[rowSums(CA.otu)>0,]

US.map <- filter(iso_to_loc, Country.of.origin == "US")
US.otu <- otu[,colnames(otu) %in% US.map$Run]
US.otu <- US.otu[rowSums(US.otu)>0,]

# Get the accumulation curves for all datasets and reformat
ET.sp <- specaccum(t(ET.otu), method = "random", permutations = 100)
ET.perm <- ET.sp$perm
ET.perm.tidy <- as_tibble(ET.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
ET.richness <- data.frame(Sites = ET.sp$sites, Species = ET.sp$richness)

IN.sp <- specaccum(t(IN.otu), method = "random", permutations = 100)
IN.perm <- IN.sp$perm
IN.perm.tidy <- as_tibble(IN.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
IN.richness <- data.frame(Sites = IN.sp$sites, Species = IN.sp$richness)

PK.sp <- specaccum(t(PK.otu), method = "random", permutations = 100)
PK.perm <- PK.sp$perm
PK.perm.tidy <- as_tibble(PK.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
PK.richness <- data.frame(Sites = PK.sp$sites, Species = PK.sp$richness)

MR.sp <- specaccum(t(MR.otu), method = "random", permutations = 100)
MR.perm <- MR.sp$perm
MR.perm.tidy <- as_tibble(MR.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
MR.richness <- data.frame(Sites = MR.sp$sites, Species = MR.sp$richness)

TU.sp <- specaccum(t(TU.otu), method = "random", permutations = 100)
TU.perm <- TU.sp$perm
TU.perm.tidy <- as_tibble(TU.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
TU.richness <- data.frame(Sites = TU.sp$sites, Species = TU.sp$richness)

AU.sp <- specaccum(t(AU.otu), method = "random", permutations = 100)
AU.perm <- AU.sp$perm
AU.perm.tidy <- as_tibble(AU.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
AU.richness <- data.frame(Sites = AU.sp$sites, Species = AU.sp$richness)

CA.sp <- specaccum(t(CA.otu), method = "random", permutations = 100)
CA.perm <- CA.sp$perm
CA.perm.tidy <- as_tibble(CA.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
CA.richness <- data.frame(Sites = CA.sp$sites, Species = CA.sp$richness)

US.sp <- specaccum(t(US.otu), method = "random", permutations = 100)
US.perm <- US.sp$perm
US.perm.tidy <- as_tibble(US.perm) %>% 
  mutate(Sites = 1:nrow(.)) %>% 
  gather(key = "Permutation", value = "Species", -Sites) 
US.richness <- data.frame(Sites = US.sp$sites, Species = US.sp$richness)

# Plot accumulation curves
perm.tidy <- rbind(mutate(ET.perm.tidy, Habitat = "Ethiopia"),
                   mutate(IN.perm.tidy, Habitat = "India"),
                   mutate(PK.perm.tidy, Habitat = "Pakistan"),
                   mutate(TU.perm.tidy, Habitat = "Turkey"),
                   mutate(MR.perm.tidy, Habitat = "Morocco"),
                   mutate(AU.perm.tidy, Habitat = "Australia"),
                   mutate(CA.perm.tidy, Habitat = "Canada"),
                   mutate(US.perm.tidy, Habitat = "USA")) %>% 
  mutate(Habitat = fct_relevel(Habitat, "Ethiopia", "India", "Pakistan", "Turkey", 
                               "Morocco", "Australia", "Canada", "USA")) 
richness <- rbind(mutate(ET.richness, Habitat = "Ethiopia"),
                  mutate(IN.richness, Habitat = "India"),
                  mutate(PK.richness, Habitat = "Pakistan"),
                  mutate(TU.richness, Habitat = "Turkey"),
                  mutate(MR.richness, Habitat = "Morocco"),
                  mutate(AU.richness, Habitat = "Australia"),
                  mutate(CA.richness, Habitat = "Canada"),
                  mutate(US.richness, Habitat = "USA")) %>% 
  mutate(Habitat = fct_relevel(Habitat, "Australia", "Canada", "USA", "Turkey", 
                               "Morocco", "Ethiopia", "India", "Pakistan")) 

acc.p <- ggplot(perm.tidy, aes(Sites, Species, color = Habitat)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_line(data = richness, size = 1) +
  scale_color_brewer(12, palette = "Paired", direction=-1) +
  xlim(0,100) +
  ylab("Cumulative richness\n(# vOTUs)") +
  xlab("Sampling effort\n(# samples)") +
  theme_light() +
  theme(text = element_text(size = 12),
        legend.position = "top",
        panel.border = element_blank()) 
acc.p
ggsave("new_figs/rhizophage_accumulation_curves_min2.pdf", acc.p, width=8, height=6)

# Viruses in plant metagenomes vs. MAGs

# Filter the vOTU table to sample that have a country assigned
otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run")

otu_t_valid <- otu_t_valid %>%
  filter(Country.of.origin != "") %>% # Throwing out samples without strain specified
  filter(!is.na(Country.of.origin)) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")

nrow(otu_t_valid[otu_t_valid$organism == "plant metagenome",]) # 114
nrow(otu_t_valid[otu_t_valid$organism != "plant metagenome",]) # 42

otu_t_valid_MG <- otu_t_valid %>%
  filter(organism == "plant metagenome")

otu_t_valid_Cul <- otu_t_valid %>%
  filter(organism != "plant metagenome")

# Create the distance matrix using the Jaccard algorithm (just presence/absence, not caring how high the coverage is)
pco_dist <- vegdist(otu_t_valid_MG[, which(colnames(otu_t_valid_MG) %like% "NZ")], method="jaccard")

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

pcoa_MG <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.1, y = Axis.2, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country in plant metagenomes') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
ggsave("new_figs/rhizophages_PCoA_MG_min2.pdf", pcoa_MG)

pco_dist <- vegdist(otu_t_valid_Cul[, which(colnames(otu_t_valid_Cul) %like% "NZ")], method="jaccard")

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

pcoa_Cul <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.1, y = Axis.2, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 1 ',eigenvalues[1],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
ggsave("new_figs/rhizophages_PCoA_Cul_Co1_Co2_min2.pdf", pcoa_Cul)

pcoa_Cul_2 <- ggplot(pcoa_bray_df) + # Plot the PCoA
  geom_point(aes(x = Axis.2, y = Axis.3, color=Country.of.origin), size = 1) +
  ylab(paste0('Co 3 ',eigenvalues[3],'%')) + #Extract y axis value from variance
  xlab(paste0('Co 2 ',eigenvalues[2],'%')) + #Extract x axis value from variance
  ggtitle('vOTU clustering by country') +
  coord_fixed(ratio = 1) +
  theme_bw() + scale_color_brewer(12, palette="Paired", direction=-1)
ggsave("new_figs/rhizophages_PCoA_Cul_Co2_Co3_min2.pdf", pcoa_Cul_2)

# ANOSIM
# By strain
temp2 <- otu_t_valid[!is.na(otu_t_valid$ANI95.OTU) & otu_t_valid$ANI95.OTU != "",] # Throw out samples without a dominant Mesorhizobium strain
com = temp2[,which(colnames(temp2) %like% "NZ_")] # Select only vOTU columns
m_com = as.matrix(com) # Convert to matrix
ano = anosim(m_com, temp2$ANI95.OTU, distance = "euclidean", permutations = 999) 
summary(ano) # min 4 vOTUs: R=0.32, p=0.001, when using only nodule MGs R=0.3, p=0.001
             # min 3 vOTUs: R=0.25, p=0.001
# By country
temp2 <- otu_t_valid[!is.na(otu_t_valid$Country.of.origin) & otu_t_valid$Country.of.origin != "",]
com = temp2[,which(colnames(temp2) %like% "NZ_")]
m_com = as.matrix(com)
ano = anosim(m_com, temp2$Country.of.origin, distance = "euclidean", permutations = 999) 
summary(ano) # min 4 vOTUs: R=0.06, p=0.02, when using only nodule MGs R=0.02, p=0.27
             # min 3 vOTUs: R=0.08, p=0.037
# By plant host
temp2 <- otu_t_valid[!is.na(otu_t_valid$Cicer.host),]
com = temp2[,which(colnames(temp2) %like% "NZ_")]
m_com = as.matrix(com)
ano = anosim(m_com, temp2$Cicer.host, distance = "euclidean", permutations = 999) # min 4 vOTUs: R=0.01, p=0.37
                                                                                  # min 3 vOTUs: R=0.08, p=0.16
# By soil genus
temp2 <- otu_t_valid[!is.na(otu_t_valid$soil.genus),]
com = temp2[,which(colnames(temp2) %like% "NZ_")]
m_com = as.matrix(com)
ano = anosim(m_com, temp2$soil.genus, distance = "euclidean", permutations = 999) # min 4 vOTUs: R=0.09, p=0.017
                                                                                  # min 3 vOTUs: R=0.04, p=0.21
p.adjust(c(0.001, 0.037, 0.16, 0.21), method="BH") # 0.004 0.074 0.21 0.21

# Pairwise dissimilarities within and between Mesorhizobium strains
temp2 <- otu_t_valid[!is.na(otu_t_valid$ANI95.OTU) & otu_t_valid$ANI95.OTU != "",]
com = temp2[,which(colnames(temp2) %like% "NZ_")]
jac_dist <- vegdist(com, method="jaccard")
df <- reshape2::melt(as.matrix(jac_dist), varnames = c("row", "col"))
df <- df[df$row != df$col,] # Remove self comparisons
df2 <- merge(df, iso_to_loc[, c("Run", "ANI95.OTU")], by.x="row", by.y="Run")
colnames(df2)[length(colnames(df2))] <- "row.strain"
df2 <- merge(df2, iso_to_loc[, c("Run", "ANI95.OTU")], by.x="col", by.y="Run")
colnames(df2)[length(colnames(df2))] <- "col.strain"
df2$grp[df2$row.strain == df2$col.strain] <- "within"
df2$grp[df2$row.strain != df2$col.strain] <- "between"
ggboxplot(data=df2, x="grp", y="value", fill="grp") + xlab("Comparison") + ylab("% dissimilarity")
ggsave("new_figs/dissimilarity_boxplot.pdf")
wilcox.test(df2$value[df2$grp=="within"], df2$value[df2$grp=="between"], alternative="less")

# Trying a Mantel test

ad <- iso_to_loc %>%
  dplyr::select(Run, lat_lon, Mean_annual_precipitation, Mean_annual_temp, altitude, soil.ph) %>%
  mutate(lat = as.numeric(gsub(" N.*", "", lat_lon))) %>%
  mutate(lon = as.numeric(gsub(" .*", "", gsub(".* N ", "", lat_lon)))) %>%
  mutate(lat_lon = na_if(lat_lon, "not collected")) %>%
  mutate(soil.ph = soil.ph/10) %>%
  filter(!is.na(Run)) %>%
  unique() %>%
  filter(Run %in% rownames(otu_t)) %>%
  remove_rownames() %>%
  column_to_rownames("Run") %>%
  filter_all(all_vars(!is.na(.))) %>%
  select(-lat_lon)

otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run")

otu_t_valid <- otu_t_valid %>%
  filter(Country.of.origin != "") %>% # Throwing out samples without strain specified
  filter(!is.na(Country.of.origin)) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")

otu_t_valid <- otu_t_valid %>%
  rownames_to_column("Samp") %>%
  filter(Samp %in% rownames(ad)) %>%
  column_to_rownames("Samp") %>%
  select(where(is.numeric))

otu_dist <- vegdist(otu_t_valid, method="jaccard", na.rm = T)
ad_dist <- vegdist(ad, method="euclidean", na.rm = T)

mantel(otu_dist, ad_dist)
# p=0.72, R=-0.1, num_perm=999
mant <- mantel.correlog(otu_dist, ad_dist)
mant_val <- mant$mantel.res[1:8,] %>% data.frame() %>% mutate(Pr.bin = ifelse(Pr.corrected. >= 0.1, "sig", "not_sig"))
ggplot(data=mant_val, aes(x=class.index, y=Mantel.cor)) + geom_point(aes(fill=Pr.bin), color="black", pch=21, size=2) + theme_classic() + geom_line() + scale_fill_manual(values=c("darkblue", "lightblue"))

# Add distance decay

require("geosphere")
require("geodist")
ad$lon_cor <- ad$lon
#ad$lon_cor[ad$lon_dir == "W"] <- ad$lon_cor[ad$lon_dir == "W"]*-1
ad_for_dist <- ad %>% dplyr::select(lat, lon_cor) %>% unique()

otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run")

otu_t_valid <- otu_t_valid %>%
  filter(Country.of.origin != "") %>% # Throwing out samples without strain specified
  filter(!is.na(Country.of.origin)) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")

otu_t_valid <- otu_t_valid[rownames(otu_t) %in% rownames(ad),] %>%
  select(where(is.numeric))
otu_t_valid <- otu_t_valid[rownames(otu_t_valid) %in% rownames(ad_for_dist),]

otu_dist <- vegdist(otu_t_valid, method="jaccard", na.rm = T)

geo_dist <- geodist(ad_for_dist, measure="vincenty")

mantel.correlog(otu_dist, geo_dist)  # No significance

# Common vOTUs (5 or more countries)
otu_t_valid <- otu_t %>%
  rownames_to_column(var="Run") %>%
  left_join(iso_to_loc, by="Run")

otu_t_valid <- otu_t_valid %>%
  filter(Country.of.origin != "") %>% # Throwing out samples without strain specified
  filter(!is.na(Country.of.origin)) %>%
  filter(!is.na(Run)) %>%
  column_to_rownames("Run")
colnames(otu_t_valid) <- gsub("\\.", "_", colnames(otu_t_valid))

ubiq_vOTU <- data.frame(ubiq_vOTU)
colnames(ubiq_vOTU)[1] <- "variable"
ubiq_vOTU$rand_name <- paste0("vOTU_", 1:2)

com_votus <- otu_t_valid %>%
  select(c(ubiq_vOTU$variable, "isolation_source")) %>%
  filter(!is.na(isolation_source)) %>%
  filter(isolation_source %like% "Cicer")
com_votus_sum <- com_votus %>%
  reshape2::melt() %>%
  group_by(variable, isolation_source) %>% summarize_at("value", sum) %>%
  left_join(ubiq_vOTU, by="variable")

ggbarplot(data=com_votus_sum, x = "rand_name", y = "value", group = "isolation_source", fill = "isolation_source") +
  scale_fill_brewer(palette = "Greens") + xlab("") + ylab("Number of samples")
ggsave("new_figs/common_vOTUs_Cicer_min2.pdf")

# What was the dominant Mesorhizobium strain in samples these vOTUs identified in?
ubiq_vOTU$variable
unique(otu_t_valid$ANI95_OTU[otu_t_valid$NZ_RZTI01000059_1_Mesorhizobium_sp__M6A_T_Ce_TU_002_03_1_1_NODE_59_length_35292_cov_4_27456__whole_genome_shotgun_sequence_fragment_1 == 1])
# 5A, 5B, 5D, 6A, 7A, 7B, 8A, 9A
unique(otu_t_valid$ANI95_OTU[otu_t_valid$NZ_RZSU01000038_1_Mesorhizobium_sp__M5C_F_Cr_IN_023_01_1_1_NODE_38_length_58474_cov_4_93071__whole_genome_shotgun_sequence_fragment_1 == 1])
# 1A, 2A, 5A, 5B, 5C, 5D, 5E, 7A, 8A, 9A