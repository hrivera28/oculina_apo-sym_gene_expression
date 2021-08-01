library(rnaturalearth)
library(sp)
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(ggspatial)
library(sf)
library(mapedit)
library(mapview)

countries<-ne_download(scale=10, type="countries", category = "cultural", returnclass="sf")

countries%>%mutate(Pop_fill=sqrt(as.numeric(POP_EST)))->countries
countries%>%mutate(Pop_log=log(as.numeric(POP_EST),base=2))->countries
countries[countries$Pop_log==-Inf,]$Pop_log<-0

state_points<- st_centroid(states)
state_points <- cbind(states, st_coordinates(st_centroid(states$geometry)))

cities<-ne_download(scale=10, type="populated_places", category="cultural", returnclass = "sf")
states<-ne_download(scale=50, type="admin_1_states_provinces", category="cultural", returnclass = "sf")


ggplot(data = countries)+ 
  geom_sf(fill= "antiquewhite")+ 
  geom_sf(data=states, fill="antiquewhite")+
  coord_sf(xlim = c(-90, -70), ylim = c(20, 40), expand = FALSE)+
  annotate(geom = "text", x = -87.3, y = 27.5, label = "Gulf of Mexico", 
           fontface = "italic", color = "grey22", size = 3) +
  annotation_scale(location = "bl", width_hint = 0.2, pad_y = unit(0.1, "in"), pad_x = unit(0.2, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.3, "in"), 
                         pad_y = unit(0.4, "in"), style = north_arrow_fancy_orienteering, 
                         height = unit(0.4, "in"), width = unit(0.4, "in")) + 
  xlab("") + ylab("") + 
  theme(panel.grid.major = element_line(color = gray(.75), linetype = "dashed", size = 0.3), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.text=element_text(face="bold", colour = "black"))

ggsave(filename="Figure_1.png", width=5, height=6, units = "in", dpi=300)
