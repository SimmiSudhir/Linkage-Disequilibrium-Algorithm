rm(list=ls())
df1 <- read.csv("AK1toAK4_selection.csv",header=T)
df1$Model <- as.factor(df1$Model)
library(ggplot2)


# # Scatter plot with highlighted point for AK4
# plot_df1 <- ggplot(df1, aes(x = DIC, y = mean_PA)) +
#   geom_point() +
#   geom_point(data = df1[df1$Model == "AK4", ], color = "red", size = 3) +
#   labs(x = "DIC", y = "Prediction Accuracy") +
#   theme_minimal()
# 
# # Display the plot
# plot_df1
# Set the desired y-axis limits
#y_limits <- c(min(df1$mean_PA) - 0.1, max(df1$mean_PA) + 0.1)
y_limits <- seq(0.40, 0.50, by = 0.10)
my_font_size <- 12
# Calculate the breaks for DIC

dic_breaks <- seq(7100, 7400, by = 100)
# Define the desired cell size for grid lines
grid_size <- 0.2  # Adjust the value to decrease the cell size


# Scatter plot with highlighted point for AK4 and labels for all points
plot_df1 <- ggplot(df1, aes(x = DIC, y = mean_PA, label = Model)) +
  geom_point() +
  geom_point(data = df1[df1$Model == "AK4", ], color = "red", size = 3) +
  geom_line() +
  geom_text(position = position_nudge(x = 0.2, y = 0.01), size = 3) +
  labs(x = "DIC", y = "Prediction Accuracy") +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    panel.grid.major.x = element_line(size = 0.5),  # Adjust the size as needed
    panel.grid.minor.x = element_line(size = 0.2)
  ) +
  scale_x_continuous(breaks = dic_breaks)

# Display the plot
plot_df1



# Save the plot with adjusted width and height
ggsave("AK1toAK4.png", plot = plot_df1, width = 8, height = 6)
