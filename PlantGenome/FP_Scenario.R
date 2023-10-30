rm(list=ls())
df1 <- read.csv("FP_Scenario1.csv",header=T)
str(df1) # 36 by 7
colnames(df1)[3] <- "Par"
df1$Scenario <- as.factor(df1$Scenario)
df1$Trait <- factor(df1$Trait, levels=c("TCH","CCS", "Fibre"))

df1$Par <- factor(df1$Par,levels=c("dip_26K", "con_26K", "con_58K"))
df1$Model <- factor(df1$Model, levels=c("GBLUP","RKHS","AK4","BayesR"))

library(ggplot2)
# Define the color palette with a gradient of purple
color_palette <- c("#3366CC", "#66A5D3", "#99C5E8", "#CCE5F2")
my_font_size <- 12

plot1 <- ggplot(df1, aes(x = Model, y = PA, fill = Par)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = PA - SE, ymax = PA + SE), width = 0.2, position = position_dodge(width = 0.9)) +
  #geom_text(aes(label = PA), position = position_dodge(width = 0.9), vjust =-0.02, size = 3,angle = 90)+
  facet_wrap(~Trait, scales = "free_y") +
  labs(x = "", y = "Prediction Accuracy") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) + ggtitle("A) Scenario 1")


##########################################################
df2 <- read.csv("FP_Scenario2.csv",header=T)
str(df2) # 36 by 7
colnames(df2)[3] <- "Par"
df2$Scenario <- as.factor(df2$Scenario)
df2$Trait <- factor(df2$Trait, levels=c("TCH","CCS", "Fibre"))

df2$Par <- factor(df2$Par,levels=c("dip_26K", "con_26K", "con_58K"))
df2$Model <- factor(df2$Model, levels=c("GBLUP","RKHS","AK4","BayesR"))

library(ggplot2)
# Define the color palette with a gradient of purple
color_palette <- c("#3366CC", "#66A5D3", "#99C5E8", "#CCE5F2")
my_font_size <- 12

plot2 <- ggplot(df2, aes(x = Model, y = PA, fill = Par)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = PA - SE, ymax = PA + SE), width = 0.2, position = position_dodge(width = 0.9)) +
  
  #geom_text(aes(label = PA), position = position_dodge(width = 0.9), vjust = -0.02, size = 3,angle = 90)+
  facet_wrap(~Trait, scales = "free_y") +
  labs(x = "", y = "Prediction Accuracy") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + ggtitle("B) Scenario 2")

########################################################
library(patchwork)


# Combine plots vertically
combined_plots <- plot1 / plot2

# Set the layout to two rows and one column
combined_plots <- combined_plots + plot_layout(ncol = 1, nrow = 2)

# Display the combined plots
combined_plots

# Save the plot as a PNG file
ggsave("FP_Scenario1.png", plot = plot1, width = 8, height = 6, dpi = 300)
# Save the plot as a PNG file
ggsave("FP_Scenario2.png", plot = plot2, width = 8, height = 6, dpi = 300)
# Save the plot as a PNG file
ggsave("FP_bothScenario.png", plot = combined_plots, width = 8, height = 6, dpi = 300)
