rm(list=ls())
df1 <- read.csv("CV_Results_5fold.csv",header=T)

colnames(df1)[3] <- "Par"

df1$Trait <- factor(df1$Trait, levels=c("TCH","CCS","Fibre"))
df1$Fold <- factor(df1$Fold, levels=c("1","2","3","4","5"))
df1$Par <- factor(df1$Par,levels= c("dip_26K", "Con_26K", "Con_58K"))
df1$Model <- factor(df1$Model, levels=c("GBLUP","RKHS","AK1","AK4"))

str(df1)
head(df1)



#############################################################
library(ggplot2)

# Define gradient color scheme
my_colors <- c("#3366CC", "#66A5D3", "#99C5E8", "#CCE5F2")  # Gradient of blue colors
my_font_size <- 12  # Font size for journal submission

# Create the plot
plot1 <- ggplot(df1, aes(x = Model, y = PA, fill = Model)) +
  geom_boxplot(width = 0.7, alpha = 0.7, outlier.shape = NA) +
  geom_text(
    data = aggregate(PA ~ Model + Trait + Par, data = df1, FUN = mean),
    aes(label = sprintf("%.2f", PA)), vjust = 0.47, size = 3
  ) +
  facet_grid(Trait ~ Par, scales = "free_x", space = "free_x") +
  labs(x = "Model", y = "Mean Prediction Accuracy") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    strip.placement = "outside"  # Place facet labels outside the plot area
  )+guides(fill = FALSE)  # Remove the fill legend

plot1
# Save the plot as a PNG file
ggsave("cv_results_5fold_new.png", plot = plot1, width = 8, height = 6, dpi = 300)


#######################################################################
################# MSE ########################

rm(list=ls())
df1 <- read.csv("CV_Results_5fold.csv",header=T)

colnames(df1)[3] <- "Par"

df1$Trait <- factor(df1$Trait, levels=c("TCH","CCS","Fibre"))
df1$Fold <- factor(df1$Fold, levels=c("1","2","3","4","5"))
df1$Par <- factor(df1$Par,levels= c("dip_26K", "Con_26K", "Con_58K"))
df1$Model <- factor(df1$Model, levels=c("GBLUP","RKHS","AK1","AK4"))

str(df1)
head(df1)



#############################################################
library(ggplot2)

# Define gradient color scheme
my_colors <- c("#3366CC", "#66A5D3", "#99C5E8", "#CCE5F2")  # Gradient of blue colors
my_font_size <- 12  # Font size for journal submission

# Create the plot
plot2 <- ggplot(df1, aes(x = Model, y = RMSE, fill = Model)) +
  geom_boxplot(width = 0.7, alpha = 0.7, outlier.shape = NA) +
  geom_text(
    data = aggregate(RMSE ~ Model + Trait + Par, data = df1, FUN = mean),
    aes(label = sprintf("%.2f", RMSE)), vjust = 0.47, size = 3
  ) +
  facet_grid(Trait ~ Par, scales = "free_x", space = "free_x") +
  labs(x = "Model", y = "Prediction Accuracy") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    strip.placement = "outside"  # Place facet labels outside the plot area
  )
plot2

# Save the plot as a PNG file
ggsave("cv_results_5fold_mse.png", plot = plot2, width = 8, height = 6, dpi = 300)


#######################################################################
