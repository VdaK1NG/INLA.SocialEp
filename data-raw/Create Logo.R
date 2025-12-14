# Load packages
library("cowplot")
library("showtext")
library("mapSpain")
library("ggplot2")
library("dplyr")

DIF <- c("#133BF2", "#7189F7", "#FFFFFF", "#FF867A", "#FF2F1B")
CScale_dif <- colorRampPalette(DIF)

data <- mapSpain::esp_get_ccaa()

# Add a great Google Font
font_add_google("Do Hyeon")

# Import a image to R
png <- ggdraw()

subplot  <-  ggplot(data = data %>% filter(ine.ccaa.name!="Canarias", ine.ccaa.name!="Balears, Illes", ine.ccaa.name!="Melilla", ine.ccaa.name!="Ceuta")) +
  geom_sf(aes(fill=ine.ccaa.name)) +
  theme_void() +
  theme_transparent(legend.position="none") +
  scale_fill_manual(values = CScale_dif(15))

ggsave(
  "back.png",
  plot = subplot,
  scale = 0.5
)

img <- png::readPNG("back.png")
g <- grid::rasterGrob(img, interpolate=TRUE)

logo <- ggplot() +
  geom_hexagon(size = h_size, fill = CScale_dif(17)[4],  color = NA) +
  geom_hexagon(size = h_size, fill = NA, color = h_color)  +
  annotation_custom(g, xmin=0.45, xmax=1.65, ymin=-0.05, ymax=1.6) +
  annotate("text", x=1, y=1.75, size = 12, family   = "serif", label = "bold(INLA)", parse = TRUE, col="white") +
  annotate("text", x=1, y=1.5, size = 12, family   = "serif", label = "bold(SocialEp)", parse = TRUE, col="white") +
  theme_void() + theme_transparent()

ggsave(
  "logo.png",
  plot = logo,
  width = 2,
  height = 2
)
