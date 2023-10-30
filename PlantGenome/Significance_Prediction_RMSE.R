#################################################
################################################
################################################
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Create an empty data frame to store ANOVA results
df_with_all_scenarios <- df1
# Create an empty data frame to store ANOVA results
# Create an empty data frame to store Pr(>F) values
result_summary_df <- data.frame(
  Trait = character(),
  Par = character(),
  Pr_gt_F_Model = numeric(),
  stringsAsFactors = FALSE  # Prevent automatic conversion of strings to factors
)

# Perform ANOVA within each combination of Trait and Par and store Pr(>F) values
combinations <- unique(df_with_all_scenarios[c("Trait", "Par")])

for (i in 1:nrow(combinations)) {
  subset_data <- subset(df_with_all_scenarios, Trait == combinations[i, "Trait"] & Par == combinations[i, "Par"])
  
  # Check if there are data points for this combination
  if (nrow(subset_data) > 0) {
    aov_result <- aov(PA ~ Model, data = subset_data)
    result_summary <- summary(aov_result)
    
    # Extract Pr(>F) value
    pr_gt_f <- result_summary[[1]][1, "Pr(>F)"]
    
    # Append the result to the data frame
    result_summary_df <- rbind(result_summary_df, data.frame(
      Trait = combinations[i, "Trait"],
      Par = combinations[i, "Par"],
      Pr_gt_F_Model = pr_gt_f,
      stringsAsFactors = FALSE
    ))
  }
}

# View the updated result_summary_df
print(result_summary_df)
# ##result_summary[[1]] (one example)
# Df  Sum Sq   Mean Sq F value Pr(>F)
# Model        3 0.00202 0.0006733  0.1618 0.9205
# Residuals   16 0.06660 0.0041625   
#print(result_summary_df)
## these results are the within scenario checking model performance for each trait
# #Trait     Par Pr_gt_F_Model
# 1   TCH dip_26K     0.9205148
# 2   CCS dip_26K     0.7448120
# 3 Fibre dip_26K     0.5247013
# 4   TCH Con_26K     0.9160779
# 5   CCS Con_26K     0.9902073
# 6 Fibre Con_26K     0.7837308
# 7   TCH Con_58K     0.9281862
# 8   CCS Con_58K     0.9364144
# 9 Fibre Con_58K     0.5555734

######################################################
######
## now checking statistical tests that compare the same model's performance across different scenarios
## .	For example, you can compare how the GBLUP model performs across the three marker parameterizations (dip_26K, Con_26K, Con_58K) using appropriate statistical tests. This would help determine if the choice of marker parameterization significantly impacts prediction accuracy.
# Create an empty data frame to store p-values
p_value_df <- data.frame(
  Trait1 = character(),
  Trait2 = character(),
  Model = character(),
  Scenario1 = character(),
  Scenario2 = character(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Create a list of traits, models, and scenarios
traits <- c("TCH", "CCS", "Fibre")
models <- c("GBLUP", "RKHS", "AK1", "AK4")
scenarios <- c("dip_26K", "Con_26K", "Con_58K")

# Iterate through all combinations
for (trait in traits) {
  for (model in models) {
    for (i in 1:(length(scenarios) - 1)) {
      for (j in (i + 1):length(scenarios)) {
        scenario1 <- scenarios[i]
        scenario2 <- scenarios[j]
        
        # Subset the data for the current trait, model, and scenarios
        subset_data1 <- subset(df1, Trait == trait & Model == model & Par == scenario1)$PA
        subset_data2 <- subset(df1, Trait == trait & Model == model & Par == scenario2)$PA
        
        # Perform t-test
        t_test_result <- t.test(subset_data1, subset_data2)
        
        # Extract p-value
        p_value <- t_test_result$p.value
        
        # Append results to the data frame
        p_value_df <- rbind(p_value_df, data.frame(
          Trait1 = trait,
          Trait2 = trait,
          Model = model,
          Scenario1 = scenario1,
          Scenario2 = scenario2,
          P_Value = p_value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# View the p-value data frame
print(p_value_df)
# ###################################
# In the code provided earlier, I used an independent samples t-test, also known as a two-sample t-test. Here's an explanation of the t-test used in the code:
# 
# Independent Samples T-Test (Two-Sample T-Test):
# 
# In the context of the code, the independent samples t-test is used to compare the means of two different scenarios (parameterizations) within the same model and trait. It assesses whether the choice of marker parameterization significantly impacts prediction accuracy for a specific model and trait.
# 
# Here's how the code implements the t-test:
#   
#   It starts by iterating through all combinations of traits, models, and scenarios.
# For each combination, it selects the data corresponding to the two scenarios being compared (e.g., dip_26K vs. Con_26K) within the same model and trait.
# It then performs a t-test on these two sets of data to determine whether there is a statistically significant difference between their means.
# The t-test results include a p-value, which indicates the probability that the observed difference in means occurred by chance. A small p-value (typically less than a significance level, e.g., 0.05) suggests that there is a significant difference between the means.
# The code collects the p-values for all pairwise comparisons of scenarios within each combination of trait and model, allowing you to assess which scenario pairs exhibit significant differences in prediction accuracy.
# 
# You can adapt this code to your specific data by replacing the data frame and variables with the appropriate names from your dataset. Make sure to customize it according to your data structure and research questions.
# 

print(p_value_df)
# ## Trait1 Trait2 Model Scenario1 Scenario2   P_Value
# 1     TCH    TCH GBLUP   dip_26K   Con_26K 0.9624358
# 2     TCH    TCH GBLUP   dip_26K   Con_58K 0.8089347
# 3     TCH    TCH GBLUP   Con_26K   Con_58K 0.8442677
# 4     TCH    TCH  RKHS   dip_26K   Con_26K 0.9545213
# 5     TCH    TCH  RKHS   dip_26K   Con_58K 0.9096927
# 6     TCH    TCH  RKHS   Con_26K   Con_58K 0.9574501
# 7     TCH    TCH   AK1   dip_26K   Con_26K 0.9291043
# 8     TCH    TCH   AK1   dip_26K   Con_58K 0.8938658
# 9     TCH    TCH   AK1   Con_26K   Con_58K 0.8218965
# 10    TCH    TCH   AK4   dip_26K   Con_26K 0.9286954
# 11    TCH    TCH   AK4   dip_26K   Con_58K 0.8206950
# 12    TCH    TCH   AK4   Con_26K   Con_58K 0.8887433
# 13    CCS    CCS GBLUP   dip_26K   Con_26K 0.4474715
# 14    CCS    CCS GBLUP   dip_26K   Con_58K 0.2225419
# 15    CCS    CCS GBLUP   Con_26K   Con_58K 0.6489416
# 16    CCS    CCS  RKHS   dip_26K   Con_26K 0.8314993
# 17    CCS    CCS  RKHS   dip_26K   Con_58K 0.5372045
# 18    CCS    CCS  RKHS   Con_26K   Con_58K 0.7001439
# 19    CCS    CCS   AK1   dip_26K   Con_26K 0.3727494
# 20    CCS    CCS   AK1   dip_26K   Con_58K 0.2089379
# 21    CCS    CCS   AK1   Con_26K   Con_58K 0.7783130
# 22    CCS    CCS   AK4   dip_26K   Con_26K 0.4288043
# 23    CCS    CCS   AK4   dip_26K   Con_58K 0.2964601
# 24    CCS    CCS   AK4   Con_26K   Con_58K 0.8571268
# 25  Fibre  Fibre GBLUP   dip_26K   Con_26K 1.0000000
# 26  Fibre  Fibre GBLUP   dip_26K   Con_58K 0.5893516
# 27  Fibre  Fibre GBLUP   Con_26K   Con_58K 0.5771587
# 28  Fibre  Fibre  RKHS   dip_26K   Con_26K 1.0000000
# 29  Fibre  Fibre  RKHS   dip_26K   Con_58K 0.3298351
# 30  Fibre  Fibre  RKHS   Con_26K   Con_58K 0.3594106
# 31  Fibre  Fibre   AK1   dip_26K   Con_26K 1.0000000
# 32  Fibre  Fibre   AK1   dip_26K   Con_58K 0.5770082
# 33  Fibre  Fibre   AK1   Con_26K   Con_58K 0.5306399
# 34  Fibre  Fibre   AK4   dip_26K   Con_26K 0.4937293
# 35  Fibre  Fibre   AK4   dip_26K   Con_58K 0.3001381
# 36  Fibre  Fibre   AK4   Con_26K   Con_58K 0.7273641
# > 
#######################################################
#Calculate NRMSE for each scenario (combination of Trait, Model, and Marker Parameterization):
# Calculate NRMSE ( no need)
df1$NRMSE <- df1$RMSE / max(df1$RMSE)
# Load necessary libraries
library(ggplot2)
#Create a plot to visualize NRMSE values for each scenario:
# Create a plot
ggplot(df1, aes(x = Model, y = NRMSE, fill = Par)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Trait ~ Par, scales = "free") +
  labs(x = "Model", y = "NRMSE") +
  theme_minimal()

# Assuming you have a data frame named 'df1' with your data
# Create an empty data frame to store the results
# Assuming you have your data loaded in a data frame called df1

library(dplyr)

############################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!! second approach
# Assuming you have a data frame df1 with columns Trait, Par, Model, mean_PA, RMSE, NRMSE

# 1. Calculate NRMSE for each combination
df1$NRMSE <- df1$RMSE / (max(df1$RMSE) - min(df1$RMSE))

# 2. Summarize NRMSE values separately for each trait
summary_by_trait <- df1 %>%
  group_by(Trait, Model, Par) %>%
  summarize(mean_NRMSE = mean(NRMSE), se_NRMSE = sd(NRMSE) / sqrt(n()))

# 3. Compare summarized NRMSE values for each trait across different models and parameterizations

# 4. Evaluate overall model performance based on the combined results from all traits
overall_summary <- df1 %>%
  group_by(Model, Par) %>%
  summarize(mean_NRMSE = mean(NRMSE), se_NRMSE = sd(NRMSE) / sqrt(n()))

# You can further analyze and visualize the results as needed.



# Example for Trait TCH
anova_result_model_tch <- aov(mean_NRMSE ~ Model, data = subset(summary_by_trait, Trait == "TCH"))
summary(anova_result_model_tch)

# Pairwise comparisons
pairwise_comp_model_tch <- TukeyHSD(anova_result_model_tch)
print(pairwise_comp_model_tch)
##############################################
## TCH
summary(anova_result_model_tch)

# #Df    Sum Sq   Mean Sq F value  Pr(>F)    
# Model        3 3.022e-04 0.0001007   29.66 0.00011 ***
#   Residuals    8 2.718e-05 0.0000034                    
# ---
#   Signif. codes:  
#   0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > # Pairwise comparisons
#   > pairwise_comp_model <- TukeyHSD(anova_result_model)
# > print(pairwise_comp_model)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = mean_NRMSE ~ Model, data = subset(summary_by_trait, Trait == "TCH"))
# 
# $Model
# diff          lwr           upr     p adj
# RKHS-GBLUP -0.008844296 -0.013663522 -4.025070e-03 0.0016657
# AK1-GBLUP  -0.010182088 -0.015001314 -5.362862e-03 0.0006493
# AK4-GBLUP  -0.013600892 -0.018420118 -8.781666e-03 0.0000834
# AK1-RKHS   -0.001337793 -0.006157019  3.481433e-03 0.8108422
# AK4-RKHS   -0.004756596 -0.009575822  6.262995e-05 0.0530171
# AK4-AK1    -0.003418803 -0.008238029  1.400423e-03 0.1840319
# 
# # # Example for Trait CCS
anova_result_model_ccs <- aov(mean_NRMSE ~ Model, data = subset(summary_by_trait, Trait == "CCS"))
 summary(anova_result_model_ccs)
# 
# # Pairwise comparisons
 pairwise_comp_model_ccs <- TukeyHSD(anova_result_model_ccs)
 print(pairwise_comp_model_ccs)
# 
 #summary(anova_result_model_ccs)
# Df    Sum Sq   Mean Sq F value Pr(>F)
# Model        3 1.781e-07 5.938e-08   0.261  0.852
# Residuals    8 1.823e-06 2.278e-07               
# > # Pairwise comparisons

# # 
#  print(pairwise_comp_model_ccs)
#  Tukey multiple comparisons of means
#  95% family-wise confidence level
#  
#  Fit: aov(formula = mean_NRMSE ~ Model, data = subset(summary_by_trait, Trait == "CCS"))
#  
#  $Model
#                 diff           lwr         upr     p adj
#  RKHS-GBLUP -7.432181e-05 -0.0013224285 0.001173785 0.9973343
#  AK1-GBLUP  -7.432181e-05 -0.0013224285 0.001173785 0.9973343
#  AK4-GBLUP   2.229654e-04 -0.0010251413 0.001471072 0.9376868
#  AK1-RKHS    0.000000e+00 -0.0012481067 0.001248107 1.0000000
#  AK4-RKHS    2.972873e-04 -0.0009508194 0.001545394 0.8688651
#  AK4-AK1     2.972873e-04 -0.0009508194 0.001545394 0.8688651
 
 ## For Fibre 
 # # # Example for Trait CCS
 anova_result_model_fibre <- aov(mean_NRMSE ~ Model, data = subset(summary_by_trait, Trait == "Fibre"))
 summary(anova_result_model_fibre)
 # 
 # # Pairwise comparisons
 pairwise_comp_model_fibre <- TukeyHSD(anova_result_model_fibre)
 print(pairwise_comp_model_fibre)
 
 # #summary(anova_result_model_fibre)
 # Df    Sum Sq   Mean Sq F value Pr(>F)  
 # Model        3 1.005e-05 3.352e-06   7.223 0.0115 *
 #   Residuals    8 3.712e-06 4.640e-07                 
 # ---
 #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

 ###############################################
 ## comparison among parametrisation
 # Load necessary libraries
 library(stats)
 library(dplyr)
 library(tidyr)
 library(multcomp)
 
 # Get unique combinations of Model and Trait
 unique_models <- unique(df1$Model)
 unique_traits <- unique(df1$Trait)
 
 # Initialize a list to store the results
 results_list <- list()
 
 # Loop through all combinations of Model and Trait
 for (model_name in unique_models) {
   for (trait_name in unique_traits) {
     # Filter the dataframe for the current combination
     filtered_data <- df1 %>%
       filter(Model == model_name, Trait == trait_name)
     
     # Perform ANOVA to compare NRMSE among parametrizations
     anova_result <- aov(NRMSE ~ Par, data = filtered_data)
     
     # Perform Tukey's HSD post-hoc test
     posthoc_result <- TukeyHSD(anova_result)
     
     # Store the results in a list
     result_key <- paste(model_name, trait_name, sep = "_")
     results_list[[result_key]] <- list(
       anova_summary = summary(anova_result),
       posthoc_results = posthoc_result
     )
   }
 }
 
 # Access the results for a specific combination (e.g., GBLUP_TCH)
 specific_result <- results_list[["GBLUP_TCH"]]
 
 # Access ANOVA summary for the specific combination
 anova_summary <- specific_result$anova_summary
 
 # Access post-hoc results for the specific combination
 posthoc_results <- specific_result$posthoc_results
 