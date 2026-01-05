library(sf)
library(tidyverse)
library(terra)
library(basemaps)
library(ggspatial)

# read in data
# data source: https://data.wsdot.wa.gov/arcgis/rest/services/Shared/MajorShorelinesData/FeatureServer/0
coastline <- sf::st_as_sf(vect("data/map/WSDOT_-_Major_Shorelines/WSDOT_-_Major_Shorelines.shp"))
coastline <- st_transform(coastline, crs = 4326)

coastline_plot <- ggplot() +
  geom_sf(data = coastline$geometry, fill = "lightblue", 
          color = "black", linewidth = 0.5) +
  geom_point(aes(x = -121.940857, y = 45.644077), size = 4,
             color = "firebrick") +
  geom_text(label = "Washington",
            aes(x = -123, y = 46.4),
            size = 3) +
  geom_text(label = "Oregon",
            aes(x = -123.5, y = 45.8),
            size = 3) +
  geom_text(label = "Bonneville",
            aes(x = -122.2, y = 45.85),
            size = 3) +
  geom_text(label = "dam",
            aes(x = -122.2, y = 45.75),
            size = 3) +
  geom_text(label = "Columbia River",
            aes(x = -120.75, y = 45.83),
            size = 3) +
  coord_sf(xlim = c(-125.5, -120),
           ylim = c(45.4, 46.522542)) +
  scale_x_continuous(breaks = c(-125, -123, -121),
                     labels = c("-125", "-123", "-121")) +
  scale_y_continuous(breaks = c(45.5, 46, 46.5),
                     labels = c("45.5", "46.0", "46.5")) +
  labs(x = "Longitude (\u00B0W)", y = "Latitude (\u00B0N)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "grey98", color = NA),
        panel.grid = element_blank(),
        axis.title = element_text(size = 9))
ggsave("final_code/figures/map.svg", coastline_plot, dpi = 400, 
       width = 4.75, height = 4, bg = "white")
  


