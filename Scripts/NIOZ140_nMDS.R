##%######################################################%##
#                                                          #
####              Exploratory Analysis of               ####
####            16S Amplicon Sequencing Data            ####
#                                                          #
##%######################################################%##

##%######################################################%##
#                                                          #
####           Project: foils statia NIOZ 140           ####
####            Adapted by Maaike April 2022            ####
####              (this is the right one)               ####
####               Corrected Primers                    ####
#                                                          #
##%######################################################%##
############################################################
##  version.string R version 4.2.1 (2022-06-23 ucrt)
##  nickname       Funny-Looking Kid   


##%######################################################%##
#                                                          #
####                      Beta Diversity              ####
#                                                          #
##%######################################################%##

####_______________________________________________________________________________________#%###
####                 Working Directory                 
####_______________________________________________________________________________________#%###
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")


####_______________________________________________________________________________________#%###
####                   Load libraries                                                      ####
####_______________________________________________________________________________________#%###
library(devtools)
library(phyloseq)
library(grid)
library(vegan)
library(compositions)
library("remotes")
library("ggh4x")
library(ggpubr)
library(ggrepel)
library("plotrix")
library("FactoMineR")
library("factoextra")
library(usedist)
library("ggthemes")
library(cowplot)
library(tidyverse)


####_______________________________________________________________________________________#%###
####                  Import Data &                    ####
####_______________________________________________________________________________________#%###
set.seed(42)

#Import earlier created tidy data tibble with pre-processed absolute abundance values
t2 <- read.csv("../Analysis/tidy_dataset_NIOZ140_nov22.csv", na.strings = c("")) 
head(t2)

#Make sure al T3 and negative controls are goooone
t2$revComp <- NULL
t2 <- t2 %>% filter( timepoint %in% c("T1", "T6")) %>% filter(surface != "nc") 
t2 <- t2 %>% mutate( timepoint = case_when(
  timepoint == 'T1' ~ "Day 1",
  timepoint == 'T6' ~ "Day 6" ,
)
)

#### Extra columns
time_UV <- str_c(t2$timepoint, "_", t2$treatment)
t2 <- t2 %>% 
  add_column(time_UV, .before ="Kingdom")

treat_time <- str_c(t2$treatment, "_", t2$timepoint)
t2 <- t2 %>% 
  add_column(treat_time, .before ="Kingdom")

time_mat <- str_c(t2$timepoint, "_", t2$Material)
t2 <- t2 %>% 
  add_column(time_mat, .before ="Kingdom")

t2 <- t2 %>% filter(timepoint != "T3")
head(t2)
###############################################################################################
###_______________________________________________________________________________________#%###
####  Colors!!                                                                  ####
####_______________________________________________________________________________________#%###
pal_isme <- c("#006d77", "#ffddd2", "#00C49A", "#e29578", "#83c5be")
Pal.plast <- c("#DDCC77","#117733", "#AA4499", "#88CCEE", "#332288" )
pal.time <- c("#44AA99", "#882255")
pal.uv <- c("#999933", "#CC6677")

##############################################################################################
###%#______________________________________________________________________________________#%###
####                NMDS based on previously calculated relative abundance data                 
###%#______________________________________________________________________________________#%###
ASV.RA.metadata <- t2 %>% dplyr::select(Sample, Description, timepoint, Material, treatment, surface, time_UV) %>% 
  distinct()  %>% arrange(Sample) %>% column_to_rownames(var = "Sample")
head(ASV.RA.metadata)

# Extract RA data from tidy tibble 
# Filter out all Abundance = 0 values, since this fucks up the standardization
ASV.RA.count <-  t2 %>%  dplyr::select(Sample, OTU, Sample_rel_abund)%>% filter (Sample_rel_abund > 0) %>%
  distinct() 
head(ASV.RA.count)

# Transform count table to get samples per row and ASV per column
# Not per se necessary, but this means we can use standard MARGIN in transfromation and ordination
ASV.RA.count.t <- ASV.RA.count %>% dplyr::select(Sample, OTU, Sample_rel_abund)  %>% mutate(across(c(OTU),factor)) %>% 
  pivot_wider(names_from = OTU, values_from = Sample_rel_abund) %>% replace(is.na(.), 0)  %>% arrange(Sample) %>% 
  column_to_rownames(var = "Sample") %>% as.data.frame() %>% sqrt()

# No additional data transformation needed
# Perform ordinations on transformed data
# We use Bray-Curtis distance in 2 dimensions (k=2)
# Autotransfrom = False to have control over transformations
# See above. When sample=rows and ASV=columns, we can use standard MARGINS
# Otherwise check the different ordination methods to 
ASV.RA.Bray <- metaMDS(ASV.RA.count.t , k = 3, autotransform = F, trymax = 999)

ASV.RA.Bray
 
# ASV_dist_Bray <-vegdist(ASV.RA.count, method = "bray")

nMDS.stress <- ASV.RA.Bray$stress

###_______________________________________________________________________________________#%###
####                      Plotting ordinations        ####
####______________________________________________________________________________________#%###

###____NMDS Bray-Curtis RA data___________________________________________________________####
###---------------------------------------------------------------------------------------#%###
# Prepare data frame to plot ordination
# extract NMDS scores (x and y coordinates) from the ordination of choice
data.scores = as.data.frame(scores(ASV.RA.Bray)$sites)

# Add relevant metadata to dataframe for plotting
data.scores$Samplename = ASV.RA.metadata$Description
data.scores$Timepoint= ASV.RA.metadata$timepoint
data.scores$Material = ASV.RA.metadata$Material
data.scores$treatment= ASV.RA.metadata$treatment
data.scores$surface= ASV.RA.metadata$surface
data.scores$Ti.UV = ASV.RA.metadata$time_UV
head(data.scores)

## Create NMDS ordination plot with stat ellipse
#POL:ISHED ONE
nmds_plot_ASV = 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 3.5, aes(shape = treatment, fill = Timepoint, stroke = 0.75))+ 
  geom_text_repel(aes(label= Material), size = 3, max.overlaps = 25) +
  annotate("text", x = -1.6, y = 1.55, label = paste0("k=2; stress= ", format(nMDS.stress, digits = 4)), hjust = 0, size = 3.5) +
  theme_pubr() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position="right", legend.box = "horizontal", 
        axis.title.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.line = element_blank(), 
        legend.title = element_text(size = 12),
        panel.background = element_blank(), panel.border = element_rect(colour = "darkgrey", fill = NA, linewidth = 1.2),
        legend.key=element_blank()
        ) +
  # Set axislabels and title 
  labs( title = "")  + 
  stat_ellipse(type = "t" , level = 0.5,
               geom = "path", linewidth = 1,
               aes(color = Timepoint,
                   linetype = Ti.UV),
               show.legend = T) +
  scale_linetype_manual( name = "Ellipses groups", values =  c(1,3,1,3)) +
  scale_shape_manual(values = c(21,24)) +
  scale_fill_manual(values = pal.time) +
  scale_color_manual(values = pal.time) +
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL, size = 4.5)),
         shape = guide_legend(title = "Treatment", override.aes = list(fill = "black", linetype = NULL, size = 4.5)),
         linetype = guide_legend(override.aes = list(color = c("#44AA99","#44AA99", "#882255", "#882255"))))


nmds_plot_ASV


legend.b <- get_legend(nmds_plot_ASV +
                       theme(legend.direction = "vertical",
                             legend.position = "bottom",
                             legend.justification = "center",
                             legend.spacing.x = unit(0.5, 'cm'),
                             legend.spacing.y = unit(0.1, 'cm')))

plot_grid(Chao1 + theme(legend.position ="none"),
          Simpson + theme(legend.position ="none"),
          Shannon + theme(legend.position ="none"),
          nmds_plot_ASV+ theme(legend.position ="none"),
          legend.a,
          legend.b,
          ncol = 2,
          align = 'v',
          axis = "l",
          labels = c("A", "B", "C", "D"),
          label_size = 12,
          rel_heights = c(1,1,0.35),
          rel_widths = c(1,1))

ggsave("Alpha_Beta_div_grid_3.tiff", width = 26, height  = 22, unit = "cm", dpi = 500)

# See groupings in the plot?
# Maybe add elipses


nmds_plot_ASV = 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 6, aes( shape = treatment, color = Timepoint))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  # Set axislabels and title 
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS EUX Foils BC RA-data")  + 
  scale_colour_manual(values = pal.time)

nmds_plot_ASV 

####_________Create Ellipses_____________________________________________________________####
# Generic function for ellipses
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 500) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Set up a data frame before running the function for ellipses
df_ell <- data.frame()

# Select ata for ellipse as a vector from the dataframe 
## Choose grouping factor for elipses in 4 places 
ell_data <-as.factor(data.scores$Ti.UV)       # <-------           Chose here grouping factor for elipse


# Run loop over data vector with ellipse function to get ellipse values
for(g in levels(ell_data)){
  
  df_ell <- rbind(df_ell, 
                  cbind(as.data.frame(with(data.scores[ell_data==g,],
                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))), 
                        
                        Ti.UV=g))        # <------------                   Chose here grouping factor for elipse
  
}

# data for labelling the ellipse
NMDS.mean.ASV=aggregate(data.scores[ ,c("NMDS1", "NMDS2")], 
                        list(group = ell_data), "mean")

## Now add elipses to plot
## Choose groupin factor for elipses in 4 places
## Select the correct stressvalue based on the NMDS ordination you plotted
## Check location of stressvalue based on NMDS axes
nmds_plot_ASV+
  
  geom_text_repel(aes(label= Material), fontface = "bold", size = 4, max.overlaps = 25) +
  
  
  # This part of the function creates outline for ellipse
  geom_path(data = df_ell, 
            aes(x = NMDS1, y = NMDS2,
                linetype = factor(Ti.UV)),
            size =1, color = "darkgrey"# <-----     Chose here grouping factor for elipse 
  ) +
  
  # This part of function creates fill for ellipse
  # geom_polygon(
  #   data=df_ell,
  #   aes(x=NMDS1,y=NMDS2,
  #     fill = Material),       # <-----------                    Chose here grouping factor for elipse
  #   alpha = 0.1)  +
  # 
  # Here we add stressvalue to the plot
  #Pick x and y coordinates based on what the plot looks like 
  annotate("text", x = -1.3, y = 1.5, label = paste0("k=3; stress= ", format(nMDS.stress, digits = 4)), hjust = 0, size = 5) +
  
  ##Choose colorvectors for elipses here
  scale_linetype_manual(values =  c(1,2,3,4)) +
  scale_colour_manual(values = pal.time) +
  guides(color = guide_legend(override.aes = list(shape = 15)))

##########################################################################################%###
####______________________________________________________________________________________#%###
####                NMDS Per timepoint              ####
####                             
####______________________________________________________________________________________#%###
ASV <- t2  %>%  select(Sample, timepoint, Material, treatment, surface, Sample, time_UV, treat_time, time_mat, OTU, Sample_rel_abund)%>% 
  distinct() 
ASV <- ASV %>% filter( timepoint %in% c("Day 1", "Day 6")) #%>% mutate(timepoint = ifelse(Material =="negative_c", "T1", timepoint))
unique(ASV$timepoint)

RA.T1 <- ASV %>% filter(timepoint == "Day 1") %>% select(Sample, OTU, Sample_rel_abund) %>% 
  filter (Sample_rel_abund > 0) %>%   distinct() 

RA.T1<- RA.T1 %>% select(Sample, OTU, Sample_rel_abund) %>% 
  pivot_wider(names_from = OTU, values_from = Sample_rel_abund) %>% replace(is.na(.), 0)  %>% 
  column_to_rownames(var = "Sample") %>% as.data.frame() %>% sqrt()


RA.T6 <-  ASV %>% filter(timepoint == "Day 6") %>% select(Sample, OTU, Sample_rel_abund) %>% 
  filter (Sample_rel_abund > 0) %>%   distinct() 


RA.T6<- RA.T6 %>% select(Sample, OTU, Sample_rel_abund) %>% 
  pivot_wider(names_from = OTU, values_from = Sample_rel_abund) %>% replace(is.na(.), 0)  %>% 
  column_to_rownames(var = "Sample") %>% as.data.frame() %>% sqrt()
head(RA.T6)

#Subset metadata per habitat
T1.metadata <- ASV.RA.metadata %>% filter(timepoint == "Day 1")
T6.metadata <- ASV.RA.metadata %>% filter(timepoint == "Day 6")

###%#_______________________________________________________________________________________#%###
####                      Plotting ordination T1              ####
####______________________________________________________________________________________#%###
T1.RA.Bray <- metaMDS(RA.T1 , k = 3, autotransform = F, trymax = 999)
T1.RA.Bray


###____NMDS Bray-Curtis RA data___________________________________________________________####
###%#---------------------------------------------------------------------------------------#%###
# Prepare data frame to plot ordination
# extract NMDS scores (x and y coordinates) from the ordination of choice
data.scores = as.data.frame(scores(T1.RA.Bray)$sites)

# Add relevant metadata to dataframe for plotting
data.scores$Samplename = T1.metadata$Description
data.scores$timepoint= T1.metadata$timepoint
data.scores$Material = T1.metadata$Material
data.scores$treatment= T1.metadata$treatment
data.scores$surface= T1.metadata$surface
head(data.scores)

nMDS.stress <- T1.RA.Bray$stress


## Create NMDS ordination plot
nmds_plot_ASV = 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = NMDS1, y = NMDS3)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 6, aes( shape = Material, color = treatment))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  # Set axislabels and title 
  labs( title = "nMDS T6",
       subtitle = "nMDS of relative abundance data of sample communities based on Bray-Curtis distances")  + 
  scale_colour_manual(values = rev(pal.uv)) +
  annotate("text", x = -1.3, y = 1, label = paste0("k=2; stress= ", format(nMDS.stress, digits = 4)), hjust = 0, size = 5) 

nmds_plot_ASV

## Create NMDS ordination plot
#POL:ISHED ONE
nmds_plot_T1_UV = 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 6, aes(shape = Material, fill = treatment, stroke = 1))+ # chose here what is shapes amd colors
  geom_text_repel(aes(label= treatment), size = 4, max.overlaps = 25) +
  annotate("text", x = -1.1, y = 0.65, label = paste0("k=3; stress= ", format(nMDS.stress, digits = 4)), hjust = 0, size = 4) +
  theme_pubr() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10, colour ="black"),
        legend.position="right", legend.box = "vertical", 
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) +
  labs( title = "nMDS of sqrt of RA data of sample on T1 based on Bray-Curtis distances")  + 
  stat_ellipse(type = "t" , level = 0.5,
               geom = "path", linewidth = 1.2,
               aes(color = treatment,  #Set color condition for elipses
                   alpha = 0.85)) +
  scale_linetype_manual( name = "Ellipse groups", values =  c(1,4,1,4)) +
  scale_shape_manual(values = c(22,21,23,24,25)) +
  scale_fill_manual(values =  pal.uv) + # select color palette, twice
  scale_color_manual(values = pal.uv) +
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),
         shape = guide_legend(override.aes = list(fill = "black")),
         alpha = "none")

nmds_plot_ASV

# See groupings in the plot?
# Maybe add elipses



###_________Create Ellipses_____________________________________________________________###
# Generic function for ellipses
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 500) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Set up a data frame before running the function for ellipses
df_ell <- data.frame()

# Select ata for ellipse as a vector from the dataframe 
## Choose grouping factor for elipses in 4 places 
ell_data <-as.factor(data.scores$treatment)       # <-------           Chose here grouping factor for elipse


# Run loop over data vector with ellipse function to get ellipse values
for(g in levels(ell_data)){
  
  df_ell <- rbind(df_ell, 
                  cbind(as.data.frame(with(data.scores[ell_data==g,],
                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))), 
                        
                        treatment=g))        # <------------                   Chose here grouping factor for elipse
  
}

# data for labelling the ellipse
NMDS.mean.ASV=aggregate(data.scores[ ,c("NMDS1", "NMDS2")], 
                        list(group = ell_data), "mean")

## Now add elipses to plot
## Choose groupin factor for elipses in 4 places
## Select the correct stressvalue based on the NMDS ordination you plotted
## Check location of stressvalue based on NMDS axes
nmds_plot_ASV+
  
  geom_text_repel(aes(label= treatment), fontface = "bold", size = 4, max.overlaps = 25) +
  
  
  # This part of the function creates outline for ellipse
  geom_path(data = df_ell, 
            aes(x = NMDS1, y = NMDS2, 
                color =treatment), # <-----     Chose here grouping factor for elipse 
            linetype = 1, size = 1, alpha = 1) +
  
  # This part of function creates fill for ellipse
  geom_polygon(
    data=df_ell,
    aes(x=NMDS1,y=NMDS2,
        fill = treatment),       # <-----------                    Chose here grouping factor for elipse
    alpha = 0.1)  +
  
  # Here we add stressvalue to the plot
  #Pick x and y coordinates based on what the plot looks like 
  
  
  ##Choose colorvectors for elipses here
  scale_colour_manual(values = rev(pal_darj2)) +
  scale_fill_manual(values = rev(pal_darj2)) 



###_______________________________________________________________________________________#%###
####                      Plotting ordination T6     ####
####______________________________________________________________________________________#%###
T6.RA.Bray <- metaMDS(RA.T6 , k = 3, autotransform = F, trymax = 150)
T6.RA.Bray

###____NMDS Bray-Curtis RA data___________________________________________________________####
###---------------------------------------------------------------------------------------#%###
# Prepare data frame to plot ordination
# extract NMDS scores (x and y coordinates) from the ordination of choice
data.scores = as.data.frame(scores(T6.RA.Bray)$sites)

# Add relevant metadata to dataframe for plotting
data.scores$Samplename = T6.metadata$Description
data.scores$timepoint= T6.metadata$timepoint
data.scores$Material = T6.metadata$Material
data.scores$treatment= T6.metadata$treatment
data.scores$surface= T6.metadata$surface
head(data.scores)

nMDS.stress <- T6.RA.Bray$stress



## Create NMDS ordination plot
nmds_plot_ASV = 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 6, aes( shape = Material, color = treatment))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5)) +
  # Set axislabels and title 
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS EUX Foils BC RA-data")  + 
  scale_colour_manual(values = Pal.plast) +
  
  annotate("text", x = -1, y = 0.6, label = paste0("k=3; stress= ", format(nMDS.stress, digits = 4)), hjust = 0, size = 5) 


nmds_plot_ASV
# See groupings in the plot?
# Maybe add elipses


nmds_plot_T6_pols = 
  # Create axis based on both NMDS values in 2 dimensions
  ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +  
  #Plot the points of the NMDS, select what you want as shape and color from metadata
  geom_point(size = 6, aes(shape = treatment, fill = Material, stroke = 1))+   #<-- choose variable for shape and fill
  geom_text_repel(aes(label= Material), size = 4, max.overlaps = 25) +
  annotate("text", x = -0.9, y = 1.2, label = paste0("k=3; stress= ", format(nMDS.stress, digits = 4)), hjust = 0, size = 4) +
  theme_pubr() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10, colour ="black"),
        legend.position="left", legend.box = "vertical", 
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11, colour = "black"),
        legend.title = element_text(size = 11, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1)) +
  labs( title = "nMDS sqrt of RA data of sample on T6 based on Bray-Curtis distances")  + 
  stat_ellipse(type = "t" , level = 0.5,
               geom = "path", size = 1.2,
               aes(color = Material,       #<-- what variable do you want to use for the ellipses
                   alpha = 0.85)) +
  scale_linetype_manual( name = "Ellipse groups", values =  c(1,4,1,4)) +
  scale_shape_manual(values = c(22,21,23,24,25)) +
  scale_fill_manual(values = Pal.plast) +       #<-- Pick color palette
  scale_color_manual(values = Pal.plast) +    
  guides(fill = guide_legend(override.aes = list( shape = 22, linetype = NULL)),
         shape = guide_legend( override.aes = list(fill = "black")),
         alpha = F)

nmds_plot_ASV





###_________Create Ellipses_____________________________________________________________###
# Generic function for ellipses
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 500) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Set up a data frame before running the function for ellipses
df_ell <- data.frame()

# Select ata for ellipse as a vector from the dataframe 
## Choose grouping factor for elipses in 4 places 
ell_data <-as.factor(data.scores$treatment)       # <-------           Chose here grouping factor for elipse


# Run loop over data vector with ellipse function to get ellipse values
for(g in levels(ell_data)){
  
  df_ell <- rbind(df_ell, 
                  cbind(as.data.frame(with(data.scores[ell_data==g,],
                                           veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))), 
                        
                        treatment=g))        # <------------                   Chose here grouping factor for elipse
  
}

# data for labelling the ellipse
NMDS.mean.ASV=aggregate(data.scores[ ,c("NMDS1", "NMDS2")], 
                        list(group = ell_data), "mean")

## Now add elipses to plot
## Choose groupin factor for elipses in 4 places
## Select the correct stressvalue based on the NMDS ordination you plotted
## Check location of stressvalue based on NMDS axes
nmds_plot_ASV +
  
  geom_text_repel(aes(label= Material), fontface = "bold", size = 4, max.overlaps = 25) +
  
  
  # This part of the function creates outline for ellipse
  geom_path(data = df_ell, 
            aes(x = NMDS1, y = NMDS2, 
                color = treatment), # <-----     Chose here grouping factor for elipse 
            linetype = 1, size = 1, alpha = 1) +
  
  # This part of function creates fill for ellipse
  geom_polygon(
    data=df_ell,
    aes(x=NMDS1,y=NMDS2,
        fill = treatment),       # <-----------                    Chose here grouping factor for elipse
    alpha = 0.1)  +
  
  # Here we add stressvalue to the plot
  #Pick x and y coordinates based on what the plot looks like 
  
  
  ##Choose colorvectors for elipses here
  scale_fill_manual(values = pal_alh)  


#### Plot Grid ----------------------------------------------------------------------------------------
legend.a <- get_legend(nmds_plot_T1_UV+
                         theme(legend.direction = "vertical",
                               legend.title.align = 0.5))

legend.b <- get_legend(nmds_plot_T1_pols+
                         theme(legend.direction = "vertical",
                               legend.title.align = 0.5))



plot_grid(nmds_plot_T1_UV + theme(legend.position ="none", plot.title = element_blank()),
          legend.a,
          nmds_plot_T6_UV + theme(legend.position ="none", plot.title = element_blank()),
          nmds_plot_T1_pols + theme(legend.position ="none", plot.title = element_blank()),
          legend.b,
          nmds_plot_T6_pols + theme(legend.position ="none", plot.title = element_blank()),
          labels = c("A", "", "B", "C", "", "D" ), 
          ncol = 3,
          nrow = 2,
          align = 'v',
          axis = "l",
          rel_heights = c(1,1),
          rel_widths = c(1,0.2,1))


