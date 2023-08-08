# install.packages("maps")
# library("maps")
setwd("C:/Users/mgoudriaan/Documents/R-files/Projects/NIOZ140-foils/R-project-files/amplicon-analysis_plots_alfa/scripts")

library("ggplot2")
library("ggspatial") 
library("sf")
library("ggpp")
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
require(rgdal)
library(broom)

sf_use_s2(FALSE)
carrib <- ne_countries(scale = 10, returnclass = "sf")
world_points<- st_centroid(carrib)
world_points <- cbind(carrib, st_coordinates(st_centroid(carrib$geometry)))

st.eux <- data.frame(longitude = -62.9666628, latitude = 17.499998)

carrib.p <- ggplot(data=carrib) + 
  geom_sf(fill= "grey30")+
  coord_sf(xlim = c(-61, -71), ylim = c(9,19)) +
  annotation_scale(bar_cols = c("black", "white"), location = "bl", width_hint = 0.4, pad_x = unit(1.4, "cm"), pad_y = unit(0.45, "cm")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(0, "cm"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = st.eux, aes(x = longitude, y = latitude),  
             size = 6, shape = "diamond filled", fill = "firebrick4", stroke = 2)+
  annotate(geom = "text", x = -69.2, y = 17, label = "Carribean Sea", 
           fontface = "italic", color = "grey20", size = 10) +
  annotate(geom = "text", x = -69.2, y = 18.6, label = "Greater Antilles", 
           fontface = "bold.italic", color = "black", size = 7, angle = 350) +
  annotate(geom = "text", x = -61.8, y = 17.7, label = "Leeward Islands", 
           fontface = "bold.italic", color = "black", size = 5.5, angle = 310) +
  annotate(geom = "text", x = -62.2, y = 13, label = "Windward Islands", 
           fontface = "bold.italic", color = "black", size = 5.5, angle = 70) +
  annotate(geom = "text", x = -67.7, y = 11.7, label = "Leeward Antilles", 
           fontface = "bold.italic", color = "black", size = 5.5, angle = 350) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "grey95"),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.text.x = element_text(colour = "black", size = 25),
        axis.title.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 30, colour = "black"),
        plot.title = element_text(size = 16)) +
   xlab("Longitude") + ylab("Latitude")


#Make an SP map of the island
steux <- readOGR("../topographyStatia/Contours_20n.shp") %>% broom::tidy(region = "NAME")
roads <- readOGR("../topographyStatia/main_road_20n.shp") %>% broom::tidy()   
sites <- readOGR("../topographyStatia/EUX_DiveSites_Moorings_from_coords.shp")  %>% data.frame() %>% filter(Name %in% c("Crook's Castle", "The Charles L Brown"))
site <- data.frame(longitude = -62.99, latitude = 17.483)
map <- readOGR("../topographyStatia/referenceStatia_AerialPhotos.shp", stringsAsFactors = F) 

summary(map@data)

EUX <- ggplot() + 
  geom_polygon(data=map,aes(x = long, y = lat, group = group),fill= "grey30", color = "black") +
  geom_point(data = site, aes(x = longitude, y = latitude),  
             size = 6, shape = "diamond filled", fill = "firebrick4", stroke = 2) +
  theme_classic()+ 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "grey95", color = "grey10", linewidth = 2),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
                axis.title.x = element_text(face = "bold", size = 18, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 18, colour = "black"),
        plot.title = element_text(size = 16), 
        plot.background = element_rect(fill = NA, color = NA)) +
  xlab("") + ylab("")

eux.tb <- tibble(x = -64, y =16, plot = list(EUX))

ggplot(data=carrib) + 
  geom_sf(fill= "grey30")+
  coord_sf(xlim = c(-61, -71), ylim = c(9,19)) +
  annotation_scale(bar_cols = c("black", "white"), location = "bl", width_hint = 0.4, pad_x = unit(1.4, "cm"), pad_y = unit(0.45, "cm")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0, "cm"), pad_y = unit(0, "cm"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = st.eux, aes(x = longitude, y = latitude),  
             size = 5, shape = "diamond filled", fill = "firebrick4", stroke = 2)+
  annotate(geom = "text", x = -69, y = 17, label = "Carribean Sea", 
           fontface = "italic", color = "grey10", size = 8) +
  geom_plot(data = eux.tb, aes(x,y, label = plot)) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.5), 
        panel.background = element_rect(fill = "grey95"),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.text.x = element_text(colour = "black", size = 25),
        axis.title.x = element_text(face = "bold", size = 30, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 30, colour = "black"),
        plot.title = element_text(size = 16)) +
  xlab("Longitude") + ylab("Latitude") 


ggdraw(carrib.p) +
  draw_plot(EUX, 
            x = 0.42,
            y = 0.36,
            width = 0.4,
            height = 0.4)
         