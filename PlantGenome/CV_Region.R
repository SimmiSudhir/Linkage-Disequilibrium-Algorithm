rm(list=ls())
df <- read.csv("CV_Region.csv",header=T)
str(df)
df$Trait <- factor(df$Trait,levels=c("TCH","CCS","Fibre"))

# Add error bars to the chart (optional)

my_font_size <- 12

plot <- ggplot(df, aes(x = Region, y = PA, fill = Trait)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = PA - SD, ymax = PA + SD),
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(x = "Region", y = " Mean Prediction Accuracy") +
  
  scale_fill_brewer(palette = "Set3") + # Adjust the c
theme_minimal() +
  theme(
    axis.text.x = element_text(size = my_font_size),
    axis.text.y = element_text(size = my_font_size),
    strip.text = element_text(size = my_font_size),
    legend.text = element_text(size = my_font_size),
    legend.title = element_text(size = my_font_size),
    strip.placement = "outside" ,
      panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()# Place facet labels outside the plot area
  )

plot
# Save the plot as a PNG file
ggsave("cv_region_10Rep.png", plot = plot, width = 8, height = 6, dpi = 300)

####################################################
###### check the significance

# Load the necessary library (if not already loaded)
# # install.packages("dplyr")
# library(dplyr)
# 
# # Perform ANOVA for each trait
# results <- df %>%
#   group_by(Trait) %>%
#   summarise(p_value = anova(lm(PA ~ Region, data = .))
#             $`Pr(>F)`[1])
# 
# # Display the results
# print(results)
# #############################
# # Perform Tukey HSD post hoc test
# post_hoc_results <- df %>%
#   group_by(Trait) %>%
#   do(tukey = TukeyHSD(aov(PA ~ Region, data = .)))
# 
# # Display the post hoc test results
# print(post_hoc_results)



# Assuming you have your data in a data frame called 'df'
# Perform ANOVA for each trait
anova_tch <- aov(PA ~ Region, data = df[df$Trait == "TCH", ])
anova_ccs <- aov(PA ~ Region, data = df[df$Trait == "CCS", ])
anova_fibre <- aov(PA ~ Region, data = df[df$Trait == "Fibre", ])

# Check the results of the ANOVA tests
summary(anova_tch)
summary(anova_ccs)
summary(anova_fibre)

# Perform post hoc tests if ANOVA is significant
posthoc_tch <- TukeyHSD(anova_tch)
posthoc_ccs <- TukeyHSD(anova_ccs)
posthoc_fibre <- TukeyHSD(anova_fibre)

# Check post hoc test results
posthoc_tch
posthoc_ccs
posthoc_fibre

