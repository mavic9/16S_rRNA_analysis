#by Viktor Mamontov
library(ggradar)
library(tidyverse)
library(tidyquant)
library(scales)
library(corrr)
library(dplyr)
library(ggplot2)

df <- read.csv2('../kit_rating/kit_rating.tsv', header=TRUE,
                sep='\t', check.names = F)

reverse_values <- function(x){
  return (min(x) - x + max(x))
}

df_scores <- df %>% mutate_at(vars(-Type), reverse_values)

df_scores %>% mutate_at(vars(-Type), funs(rescale(., to=c(0, 1)))) -> df_radar
df_radar$RowSum <- rowSums(df_radar[, -which(colnames(df_radar) == "Type")])
df_radar <- df_radar %>% arrange(-RowSum)
df_radar$RowSum <- NULL


#######

colnames(df_radar) = c("Type", "v(DNA)",
                       "DIN", "Inh", "18S/16S", "P", "Rep", "α")

color_list <- c("PowerSoil" = "#ff5a5f", "PowerFecal" = "#ffb400", "Stool" = "#007a87",  
                "B&T" = "#8ce071", "LSBio" = "#7b0051", "PureLink" = "#00d1c1", 
                "Microbiome" = "#ffaa91", "Monarch" = "#b4a76c")

df_radar$Type <- factor(df_radar$Type,
                        levels = df_radar$Type)
df_radar$Type

r <- ggradar(df_radar, base.size = 10, axis.label.size = 3.9, grid.label = F,
             legend.position = "right", group.line.width = 1, group.point.size = 2,
             background.circle.colour = 'white',
             group.colours = color_list, fill = T, fill.alpha = 0.5)
r + facet_wrap(~ Type, ncol = 4) + theme_void() + 
  theme(
    strip.text = element_text(
      size = 13,
      colour = "Black",
      margin = margin(t = 5, b = 5)
    ),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )



