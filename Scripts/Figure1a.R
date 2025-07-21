#############  Figure 1a ###################

# Load packages
library("tidyverse");packageVersion("tidyverse")
library("maps");packageVersion("maps")
library("mapdata");packageVersion("mapdata")
library("sf");packageVersion("sf")
library("ggrepel");packageVersion("ggrepel")

# Load and format data
## *NOTE:* The shapefile and associated data for Scottish peatlands will have to be manually obtained from [ArcGIS Hub](https://hub.arcgis.com/datasets/snh::carbon-and-peatland-2016-map/about) since the files are too large to be stored in this repository.

# Get the map data for the UK
uk_map <- map_data("worldHires") %>%
  filter(subregion == "Great Britain" | subregion == "Scotland")

# Load the shapefiles for England and Wales
england_peatlands <- st_read("../Data/Figure1/Peaty_Soils_Location_(England)___BGS_&_NSRI.shp") # Source: https://hub.arcgis.com/datasets/Defra::peaty-soils-location-england/about
wales_peatlands <- st_read("../Data/Figure1/Unified_peat_map_for_Wales.shp") # Source: https://hub.arcgis.com/datasets/theriverstrust::unified-peat-map-for-wales/about

# The shapefile and associated data for 
scotland_peatlands <- st_read("../Data/Figure1/CARBONPEATLANDMAP_SCOTLAND.shp") # Source: https://hub.arcgis.com/datasets/snh::carbon-and-peatland-2016-map/about
scotland_peatlands <- scotland_peatlands %>%
  filter(grepl("peat", COMPSOIL) | grepl("Peat", COMPSOIL) | grepl("peat", MSSG_NAME) | grepl("Peat", MSSG_NAME) | grepl("peat", SMU_NAME) | grepl("Peat", SMU_NAME)) # Only peat areas

## Ensure all data are in the same coordinate system
england_peatlands <- st_transform(england_peatlands, crs = st_crs(4326))
wales_peatlands <- st_transform(wales_peatlands, crs = st_crs(4326))
scotland_peatlands <- st_transform(scotland_peatlands, crs = st_crs(4326))

## Create a data frame for the labels
labels <- data.frame(
  name = c("Balmoral", "Bowness", "Crocach", "Langwell", "Migneint", "Moor House", "Stean"),
  lat = c(56.92341, 54.93297, 58.39304, 58.19629, 52.99565, 54.69166, 54.1308),
  lon = c(-3.15831, -3.23945, -4.00182, -3.61489, -3.81652, -2.38228, -1.95015)
)

## Plot the map
peat_soil_map <- ggplot() +
  geom_polygon(data = uk_map, aes(x = long, y = lat, group = group), fill = "grey60", color = "grey20") +
  geom_sf(data = england_peatlands, fill = "#7B3F00", color = NA, alpha = 1) +
  geom_sf(data = wales_peatlands, fill = "#7B3F00", color = NA, alpha = 1) +
  geom_sf(data = scotland_peatlands, fill = "#7B3F00", color = NA, alpha = 1) +
  coord_sf(xlim = c(-8, 2), ylim = c(50, 60.75)) +
  theme_void() +
  geom_point(data = labels,
             aes(x = lon, y = lat),
             color = "#E41A1C",
             size = 2) +
  geom_text_repel(data = labels,
                  aes(x = lon, y = lat, label = name),
                  fontface = "bold",
                  size = 4,
                  color = "white",
                  bg.color = "black",
                  bg.r = .15,
                  box.padding = 1,
                  segment.color = "#E41A1C",
                  force = 4
  )

## Save the plot
ggsave(filename = "../Plots/Figure1/Figure_1a_map.png", plot = peat_soil_map, device = "png", dpi = 600, height = 6, width = 4, units = "in")


