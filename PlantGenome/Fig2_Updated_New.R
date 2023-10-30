rm(list=ls())
df1 <- read.csv("CV_Results_5fold.csv",header=T)

colnames(df1)[3] <- "Par"

df1$Trait <- factor(df1$Trait, levels=c("TCH","CCS","Fibre"))
df1$Fold <- factor(df1$Fold, levels=c("1","2","3","4","5"))
df1$Par <- factor(df1$Par,levels= c("dip_26K", "Con_26K", "Con_58K"))
df1$Model <- factor(df1$Model, levels=c("GBLUP","RKHS","AK1","AK4"))

str(df1)
head(df1)

library(ggplot2)

# Define gradient color scheme
my_colors <- c("#3366CC", "#66A5D3", "#99C5E8", "#CCE5F2")  # Gradient of blue colors
my_font_size <- 12  # Font size for journal submission

# Create a new variable for facet labels
df1$Facet_Label <- ifelse(df1$Par == "dip_26K", "Diploid\n 26K Markers",
                          ifelse(df1$Par == "Con_26K", "Continuous\n26K Markers",
                                 ifelse(df1$Par == "Con_58K", "Continuous\n58K Markers",
                                        "Unknown")))
df1$Facet_Label <- factor(df1$Facet_Label,levels= c("Diploid\n 26K Markers", "Continuous\n26K Markers", "Continuous\n58K Markers"))

library(ggplot2)


# Create the plot with modified facet labels and no legends
plot1 <- ggplot(df1, aes(x = Model, y = PA, fill = Model)) +
  geom_boxplot(width = 0.7, alpha = 0.7, outlier.shape = NA) +

  facet_grid(Trait ~ Facet_Label, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Mean Prediction Accuracy") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    strip.placement = "outside",  # Place facet labels outside the plot area
    legend.position = "none" ) # Remove the legends
     # Adjust plot margins for better spacing
  

plot1
# Save the plot as a PNG file
ggsave("CV_Results_5fold_New_updated.png", plot = plot1, width = 8, height = 6, dpi = 300)


  

