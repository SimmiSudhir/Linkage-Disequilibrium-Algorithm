setwd("C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Research Chapter-Allele dosages/Data/Output/15Plates_data/total_batch_file/Analysis_Review/Fig_Review")
rm(list=ls())
df1 <- read.csv("CV_Results.csv",header=T)

colnames(df1)[2] <- "Par"

df1$Trait <- as.factor(df1$Trait)
df1$Par <- as.factor(df1$Par)
df1$Model <- as.factor(df1$Model)

str(df1)
head(df1)


library(ggplot2)


# Define gradient color scheme
my_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00")   # Gradient of blue colors
my_font_size <- 12  # Font size for journal submission

# Reorder levels of Trait variable
df1$Trait <- factor(df1$Trait, levels = c("TCH", "CCS", "Fibre"))

df1$Par <- factor(df1$Par, levels = c("Dip_26K", "Con_26K", "Con_58K"))
df1$Model <- factor(df1$Model, levels=c("GBLUP","RKHS","AK1","AK4"))
# Create the plot
# Create separate line plots for each trait
plot_tch <- ggplot(df1[df1$Trait == "TCH", ], aes(x = Model, y = RMSE, color = Par, group = Par)) +
  geom_line() +
  geom_point() +
  labs(x = "", y = "Root Mean Square Error", color = "Parameterization", title="TCH") +
  scale_color_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    #legend.text = element_text(size = my_font_size),
    #legend.title = element_text(size = my_font_size)
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "gray80", size = 0.5)
  )
# Save the plot as a PNG file
#ggsave("TCH_mse.png", plot = plot_tch, width = 8, height = 6, dpi = 300)


plot_ccs <- ggplot(df1[df1$Trait == "CCS", ], aes(x = Model, y = RMSE, color = Par, group = Par)) +
  geom_line() +
  geom_point( ) +
  labs(x = "", y = " ", color = "Parameterization", title="CCS") +
  scale_color_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    #legend.text = element_text(size = my_font_size),
    #legend.title = element_text(size = my_font_size)
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "gray80", size = 0.5)
  )
plot_fibre <- ggplot(df1[df1$Trait == "Fibre", ], aes(x = Model, y = RMSE, color = Par, group = Par)) +
  geom_line() +
  geom_point() +
  labs(x = "", y = " ", color = "Parameterization", title="Fibre") +
  scale_color_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size, angle = 45, hjust = 1),
    axis.text.y = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "gray80", size = 0.5)
    #legend.position = "none"
  )
########################################################

library(patchwork)

# Combine the plots in a single row
combined_plots <- plot_tch + plot_ccs + plot_fibre +
  plot_layout(ncol = 3)

# Display the combined plots
combined_plots
ggsave("MSE_traits_cv.png", plot = combined_plots, width = 8, height = 6, dpi = 300)


# # Function to extract legend
# get_legend <- function(p) {
#   g <- ggplotGrob(p)
#   leg <- g$grobs[[which(g$layout$name == "guide-box")]]
#   leg
# }
# 
# 
# # Create legends for each plot
# legend_tch <- get_legend(plot_tch)
# legend_ccs <- get_legend(plot_ccs)
# legend_fibre <- get_legend(plot_fibre)
# #install.packages("grid")
# library(grid)
# # Arrange the plots in a grid with shared legends
# grid_arrange_shared_legend <- function(...) {
#   plots <- list(...)
#   num_plots <- length(plots)
#   
#   
#   # Create grid layout
#   grid_layout <- grid.layout(num_plots, 2)
#   grid.newpage()
#   pushViewport(viewport(layout = grid_layout))
#   
#   # Add plots and legends to the grid
#   for (i in 1:num_plots) {
#     plot_vp <- viewport(layout.pos.row = i, layout.pos.col = 1)
#     legend_vp <- viewport(layout.pos.row = i, layout.pos.col = 2)
#     
#     pushViewport(plot_vp)
#     print(plots[[i]], newpage = FALSE)
#     upViewport()
#     
#     pushViewport(legend_vp)
#     if (i == 1) print(legend_tch, newpage = FALSE)
#     if (i == 2) print(legend_ccs, newpage = FALSE)
#     if (i == 3) print(legend_fibre, newpage = FALSE)
#     upViewport()
#   }
# }
# 
# # Arrange the plots in a grid with shared legends
# grid_arrange_shared_legend(plot_tch, plot_ccs, plot_fibre)
