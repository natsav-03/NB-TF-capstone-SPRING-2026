#==================================================================================================
#   Transcription factor control of ADRN/MES identity, viability and prognosis in neuroblastoma
#   41030 Engineering Capstone Spring 2026 - Natalija Savkovic
#
#   Run from top to bottom in a fresh R session. Section numbers match those cited in the Methodology.

#   Data files must be in a folder called data/ inside the project folder:
#     data/depmap_GD.txt, data/depmap_GE_scaled.txt, data/depmap_meta.txt, data/depmap_TF Activity.txt
#     data/KOCAK_meta.txt, data/KOCAK_survival.txt, data/KOCAK_TF Activity.txt
#==================================================================================================


                         ### 0 - Setup and Data Import ###

#one time set up

      # install.packages(c("tidyverse", "survival", "survminer", "circlize"))
      #
      # # Bioconductor packages need a different installer:
      # if (!requireNamespace("BiocManager", quietly = TRUE))
      #   install.packages("BiocManager")
      #
      # BiocManager::install("ComplexHeatmap")


#Setup
  library(tidyverse)        #data import, reshaping, summaries and ggplot2 figures
  library(survival)         #Kaplan-Meier, log-rank and Cox models (Aim 3)
  library(survminer)        #Kaplan-Meier plots (Aim 3)
  library(ComplexHeatmap)   #heatmaps (Aim 4)
  library(circlize)         #colour scale for heatmaps (Aim 4)


#Read the tab-separated files
#DepMap files (Aims 1 and 2)
  depmap_GD <- read_tsv("data/depmap_GD.txt")
  depmap_GE_scaled <- read_tsv("data/depmap_GE_scaled.txt")
  depmap_meta <- read_tsv("data/depmap_meta.txt")
  depmap_TF_activity <- read_tsv("data/depmap_TF Activity.txt")

#Kocak files (Aim 3)
  kocak_meta <- read_tsv("data/KOCAK_meta.txt")
  kocak_survival <- read_tsv("data/KOCAK_survival.txt")
  kocak_TF_activity <- read_tsv("data/KOCAK_TF Activity.txt")


#check the files have loaded with the following
      # dim(depmap_GD)
      # dim(depmap_GE_scaled)
      # dim(depmap_meta)
      # dim(depmap_TF_activity)
      # dim(kocak_meta)
      # dim(kocak_survival)
      # dim(kocak_TF_activity)



#==================================================================================================
#                                         AIM 1: IDENTITY
#==================================================================================================

                      ### Aim 1.1 TF activity - Differential Testing ###

#confirm ID column =
      # colnames(depmap_meta)[1:3]
      # [1] "ModelID"   "MYCN_cn"   "condition"
      # > colnames(depmap_TF_activity)[1:3]
      # [1] "...1"  "ABL1"  "AEBP1"

#rename depmap_TF_activity col
  depmap_TF_activity <- rename(depmap_TF_activity, condition = ...1)

#combining the AM.score of depmap_meta with depmap_TF_activity
  depmap_TF_long <- left_join(depmap_TF_activity, select(depmap_meta, condition, AM.score, AM_class), by = "condition")

#to check the 2 columns were added at the end
      # select(depmap_TF_long, condition, AM.score, AM_class)

#depmap_TF_long is too wide, to stack the TF downward its split into 2 columns (TF and Activity)
      # check with dim() and head(), the coloums should have increased significantly
  depmap_TF_long <- pivot_longer(
    depmap_TF_long,
    cols = -c(condition, AM.score, AM_class),
    names_to = "TF",
    values_to = "activity"
  )

#remove rows with no AM classification - as these cant be compared, check with dim and table
  depmap_TF_long <- filter(depmap_TF_long, !is.na(AM_class))

#master list of the 763 TFs of interest (the TFs with inferred activity), used by every aim
  tf_list <- unique(depmap_TF_long$TF)


#test placeholder TF = MYCN
  mycn_TF_data <- filter(depmap_TF_long, TF == "MYCN", AM_class %in% c("ADRN", "MES"))
    #wilcox.test
      wilcox.test(
        mycn_TF_data$activity[mycn_TF_data$AM_class == "ADRN"],
        mycn_TF_data$activity[mycn_TF_data$AM_class == "MES"]
        )

      #group medians
        #adrn_median <- median(mycn_TF_data$activity[mycn_TF_data$AM_class == "ADRN"])
        #adrn_median
        #mes_median  <- median(mycn_TF_data$activity[mycn_TF_data$AM_class == "MES"])
        #mes_median

      #effect size
        #effect_size <- adrn_median - mes_median
        #effect_size


#scaling up to all TFs, check with dim and head
  aim1.1_TF_results <- summarise(
    group_by(depmap_TF_long, TF),
    p_value = wilcox.test(activity[AM_class == "ADRN"], activity[AM_class == "MES"])$p.value,
    adrn_median = median(activity[AM_class == "ADRN"]),
    mes_median = median(activity[AM_class == "MES"]),
    effect_size = adrn_median - mes_median
  )

  #sort by p-value
    aim1.1_TF_results_sorted <- arrange(aim1.1_TF_results, p_value)
  #search for specific TFs
      # filter(aim1_TF_results, TF == "PHOX2B")



              ### Aim 1.1 GE Expression - Differential Testing ###


#confirm ID column
  # colnames(depmap_GE_scaled)[1:3]
  # [1] "...1" "A1BG" "A1CF"

#rename depmap_GE_scaled col
  depmap_GE_scaled <- rename(depmap_GE_scaled, condition = ...1)

#combining the AM.score of depmap_meta with depmap_GE_scaled
  depmap_GE_long <- left_join(depmap_GE_scaled, select(depmap_meta, condition, AM.score, AM_class), by = "condition")

#depmap_GE_long is too wide, stack the genes downward into 2 columns (TF and expression)
  depmap_GE_long <- pivot_longer(
      depmap_GE_long,
      cols = -c(condition, AM.score, AM_class),
      names_to = "TF",
      values_to = "expression"
    )

#filter to just the 763 TFs of interest (depmap_GE_scaled is genome-wide), tf_list made in Aim 1.1
  depmap_GE_long <- filter(depmap_GE_long, TF %in% tf_list)

#remove rows with no AM classification
  depmap_GE_long <- filter(depmap_GE_long, !is.na(AM_class))



#test placeholder TF = MYCN
  mycn_GE_data <- filter(depmap_GE_long, TF == "MYCN", AM_class %in% c("ADRN", "MES"))

  #wilcox.test
  wilcox.test(
    mycn_GE_data$expression[mycn_GE_data$AM_class == "ADRN"],
    mycn_GE_data$expression[mycn_GE_data$AM_class == "MES"]
  )

  #group medians
  # adrn_median <- median(mycn_GE_data$expression[mycn_GE_data$AM_class == "ADRN"])
  # adrn_median
  # mes_median <- median(mycn_GE_data$expression[mycn_GE_data$AM_class == "MES"])
  # mes_median

  #effect size
  # effect_size <- adrn_median - mes_median
  # effect_size


#scaling up to all TFs, check with dim and head
  aim1.1_GE_results <- summarise(
    group_by(depmap_GE_long, TF),
    p_value = wilcox.test(expression[AM_class == "ADRN"], expression[AM_class == "MES"])$p.value,
    adrn_median = median(expression[AM_class == "ADRN"]),
    mes_median = median(expression[AM_class == "MES"]),
    effect_size = adrn_median - mes_median
  )

#sort by p-value
  aim1.1_GE_results_sorted <- arrange(aim1.1_GE_results, p_value)
  #search for specific TFs
  # filter(aim1_TF_results, TF == "PHOX2B")



                    ### Aim 1.2 TF activity - Correlation Analysis ###

#test placeholder TF = MYCN
  mycn_TF_cor_data <- filter(depmap_TF_long, TF == "MYCN")

  #runs correlation between activity values and AM.score
  #spearmean method ranks the values first rather then using raw numbers
  cor.test(mycn_TF_cor_data$activity, mycn_TF_cor_data$AM.score, method = "spearman")


#scaling up to all TFs, check with dim and head
  aim1.2_TF_correlation_results <- summarise(
    group_by(depmap_TF_long, TF),
    cor_result = list(cor.test(activity, AM.score, method = "spearman")),
    rho = cor_result[[1]]$estimate,
    p_value = cor_result[[1]]$p.value
  )

  #removing the cor_result column
  aim1.2_TF_correlation_results <- select(aim1.2_TF_correlation_results, -cor_result)


#combine TF aim1.1 and aim1.2, join by the TF column and rename the p_value column as there are 2
  aim1_TF_combined <- rename(
    left_join(aim1.1_TF_results, aim1.2_TF_correlation_results, by = "TF"),
    wilcox_p_value = p_value.x,
    cor_p_value = p_value.y
    )


                ### Aim 1.2 GE expression - Correlation Analysis ###

#test placeholder TF = MYCN
  mycn_GE_cor_data <- filter(depmap_GE_long, TF == "MYCN")
  cor.test(mycn_GE_cor_data$expression, mycn_GE_cor_data$AM.score, method = "spearman")

 #scaling up to all TFs, check with dim and head
  aim1.2_GE_correlation_results <- summarise(
    group_by(depmap_GE_long, TF),
    cor_result = tryCatch(
      list(cor.test(expression, AM.score, method = "spearman")),
      error = function(e) list(NULL)
    ),
    rho = if (!is.null(cor_result[[1]])) cor_result[[1]]$estimate else NA,
    p_value = if (!is.null(cor_result[[1]])) cor_result[[1]]$p.value else NA
  )

  #removing the cor_result column
  aim1.2_GE_correlation_results <- select(aim1.2_GE_correlation_results, -cor_result)

#combine GE aim1.1 and aim1.2,
  aim1_GE_combined <- rename(
    left_join(aim1.1_GE_results, aim1.2_GE_correlation_results, by = "TF"),
    wilcox_p_value = p_value.x,
    cor_p_value = p_value.y
  )



                         ### Aim 1.3 TF activity - Visualisation ###

# Volcano plot: effect size vs statistical significance
  ggplot(aim1_TF_combined, aes(x = effect_size, y = -log10(wilcox_p_value))) +
    geom_point(alpha = 0.6, color = "grey40", size = 0.8) +
                  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
                  # geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = -1, color = "Effect size threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 1, color = "Effect size threshold"), linetype = "dashed") +
    scale_color_manual(
      name = "Reference lines",
      values = c("p = 0.05 threshold" = "red", "Effect size threshold" = "blue")) +
    scale_y_continuous(breaks = seq(0, 6, by = 0.5)) +
    scale_x_continuous(breaks = seq(-3, 3, by = 0.5)) +
    labs(
      title = "Aim 1: Differential TF Activity, ADRN vs MES",
      x = "Effect Size",
      y = "-log10(p-value)"
    ) +
    theme_minimal()

#Figure 3A: volcano plot of differential TF activity
#Volcano plot with coloured
  aim1_TF_volcano_coloured <- mutate(aim1_TF_combined,
                                  significance = case_when(
                                    wilcox_p_value < 0.05 & effect_size > 1 ~ "ADRN-associated",
                                    wilcox_p_value < 0.05 & effect_size < -1 ~ "MES-associated",
                                    TRUE ~ "Not significant"
                                  )
  )

  ggplot(aim1_TF_volcano_coloured , aes(x = effect_size, y = -log10(wilcox_p_value))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = -1, color = "Effect size threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 1, color = "Effect size threshold"), linetype = "dashed") +
    scale_color_manual(
      name = "Legend",
      values = c(
        "ADRN-associated" = "palegreen4",
        "MES-associated" = "mediumorchid3",
        "Not significant" = "grey70",
        "p = 0.05 threshold" = "red",
        "Effect size threshold" = "blue"
      ),
      breaks = c(
        "ADRN-associated",
        "MES-associated",
        "Not significant",
        "p = 0.05 threshold",
        "Effect size threshold"
      )
    )+
    scale_y_continuous(breaks = seq(0, 6, by = 0.5)) +
    scale_x_continuous(breaks = seq(-3, 3, by = 0.5)) +
    labs(
      title = "Aim 1: Differential TF Activity, ADRN vs MES",
      x = "Effect Size",
      y = "-log10(p-value)"
    ) +
    theme_minimal()



                      ### Aim 1.3 GE expression - Visualisation ###

# Volcano plot: effect size vs statistical significance
  ggplot(aim1_GE_combined, aes(x = effect_size, y = -log10(wilcox_p_value))) +
    geom_point(alpha = 0.6, color = "grey40", size = 0.8) +
    # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    # geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = -1, color = "Effect size threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 1, color = "Effect size threshold"), linetype = "dashed") +
    scale_color_manual(
      name = "Reference lines",
      values = c("p = 0.05 threshold" = "red", "Effect size threshold" = "blue")) +
    scale_y_continuous(breaks = seq(0, 6, by = 0.5)) +
    scale_x_continuous(breaks = seq(-3, 3, by = 0.5)) +
    labs(
      title = "Aim 1: Differential Gene Expression, ADRN vs MES",
      x = "Effect Size",
      y = "-log10(p-value)"
    ) +
    theme_minimal()

#Figure 3B: volcano plot of differential gene expression
  #Volcano plot with coloured
  aim1_GE_volcano_coloured <- mutate(aim1_GE_combined,
                                  significance = case_when(
                                    wilcox_p_value < 0.05 & effect_size > 1 ~ "ADRN-associated",
                                    wilcox_p_value < 0.05 & effect_size < -1 ~ "MES-associated",
                                    TRUE ~ "Not significant"
                                  )
  )

  ggplot(aim1_GE_volcano_coloured , aes(x = effect_size, y = -log10(wilcox_p_value))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = -1, color = "Effect size threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 1, color = "Effect size threshold"), linetype = "dashed") +
    scale_color_manual(
      name = "Legend",
      values = c(
        "ADRN-associated" = "palegreen4",
        "MES-associated" = "mediumorchid3",
        "Not significant" = "grey70",
        "p = 0.05 threshold" = "red",
        "Effect size threshold" = "blue"
      ),
      breaks = c(
        "ADRN-associated",
        "MES-associated",
        "Not significant",
        "p = 0.05 threshold",
        "Effect size threshold"
      )
    )+
    scale_y_continuous(breaks = seq(0, 6, by = 0.5)) +
    scale_x_continuous(breaks = seq(-3, 3, by = 0.5)) +
    labs(
      title = "Aim 1: Differential Gene Expression, ADRN vs MES",
      x = "Effect Size",
      y = "-log10(p-value)"
    ) +
    theme_minimal()



                ### Aim 1.4 - Consistent TFs: Significant in BOTH Differential and Correlation Tests ###

## --- TF Activity ---

#TFs significant in BOTH the differential test AND the correlation test
  consistent_tfs <- filter(aim1_TF_combined,
                           wilcox_p_value < 0.05,
                           cor_p_value < 0.05
  )
  arrange(consistent_tfs, desc(abs(rho)))
  dim(consistent_tfs)

#Summary counts / percentage
  sum(aim1_TF_combined$wilcox_p_value < 0.05, na.rm = TRUE)                                            # total significant differential
  sum(aim1_TF_combined$wilcox_p_value < 0.05 & aim1_TF_combined$cor_p_value < 0.05, na.rm = TRUE)      # consistent
  mean(aim1_TF_combined$cor_p_value[aim1_TF_combined$wilcox_p_value < 0.05] < 0.05, na.rm = TRUE) * 100 # % consistent


## --- Gene Expression ---

  consistent_ge_tfs <- filter(aim1_GE_combined,
                              wilcox_p_value < 0.05,
                              cor_p_value < 0.05
  )
  arrange(consistent_ge_tfs, desc(abs(rho)))
  dim(consistent_ge_tfs)

  sum(aim1_GE_combined$wilcox_p_value < 0.05, na.rm = TRUE)
  sum(aim1_GE_combined$wilcox_p_value < 0.05 & aim1_GE_combined$cor_p_value < 0.05, na.rm = TRUE)
  mean(aim1_GE_combined$cor_p_value[aim1_GE_combined$wilcox_p_value < 0.05] < 0.05, na.rm = TRUE) * 100


#Checking whether the key identity TFs (PHOX2B, ASCL1, FEV, ZEB2, MYCN) fall in this consistent group

  filter(consistent_tfs, TF %in% c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN"))
  filter(consistent_ge_tfs, TF %in% c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN"))


#Checking overlap: TFs consistent in BOTH activity and expression

  overlap_consistent <- intersect(consistent_tfs$TF, consistent_ge_tfs$TF)
  overlap_consistent
  length(overlap_consistent)



                ### Aim 1.5 - Inconsistent TFs: Significant Differential Test but Weak Correlation ###

## --- TF Activity ---

#TFs significant in the differential test but NOT in the correlation test
  inconsistent_tfs <- filter(aim1_TF_combined,
                             wilcox_p_value < 0.05,
                             cor_p_value >= 0.05
  )
  arrange(inconsistent_tfs, wilcox_p_value)
  dim(inconsistent_tfs)

#Reverse pattern: significant correlation but NOT significant differential test
  reverse_pattern_tfs <- filter(aim1_TF_combined,
                                wilcox_p_value >= 0.05,
                                cor_p_value < 0.05
  )
  dim(reverse_pattern_tfs)

#Summary counts
  sum(aim1_TF_combined$wilcox_p_value < 0.05, na.rm = TRUE)                                    # total significant differential
  sum(aim1_TF_combined$wilcox_p_value < 0.05 & aim1_TF_combined$cor_p_value >= 0.05, na.rm = TRUE)  # inconsistent
  sum(aim1_TF_combined$wilcox_p_value < 0.05 & aim1_TF_combined$cor_p_value <  0.05, na.rm = TRUE)  # consistent (both significant)


## --- Gene Expression ---

#TFs significant in the differential expression test but NOT in correlation
  inconsistent_ge_tfs <- filter(aim1_GE_combined,
                                wilcox_p_value < 0.05,
                                cor_p_value >= 0.05
  )
  arrange(inconsistent_ge_tfs, wilcox_p_value)
  dim(inconsistent_ge_tfs)

#Summary counts
  sum(aim1_GE_combined$wilcox_p_value < 0.05, na.rm = TRUE)                                       # total significant differential
  sum(aim1_GE_combined$wilcox_p_value < 0.05 & aim1_GE_combined$cor_p_value >= 0.05, na.rm = TRUE)   # inconsistent


#Checking specific standout TFs across both activity and expression

  filter(aim1_TF_combined, TF %in% c("MYCN", "TSHZ3", "HCFC1", "POU3F1"))
  filter(aim1_GE_combined, TF %in% c("MYCN", "HSF2", "ASCL1", "EBF1", "MYC"))


#Checking overlap: TFs showing this pattern in BOTH activity and expression

  overlap_inconsistent <- intersect(inconsistent_tfs$TF, inconsistent_ge_tfs$TF)
  overlap_inconsistent
  length(overlap_inconsistent)


#==================================================================================================
#                                         AIM 2: VIABILITY
#==================================================================================================

                  ### Aim 2.1 - Differential Testing ###

#confirm ID column =
    # colnames(depmap_GD)[1:3]
    # [1] "...1" "A1BG" "A1CF"

#rename depmap_GD
  depmap_GD <- rename(depmap_GD, condition = ...1)

#combining the AM.score of depmap_meta with depmap_GD
  depmap_GD_long <- left_join(depmap_GD, select(depmap_meta, condition, AM.score, AM_class), by = "condition")

#to check the 2 columns were added at the end
      #select(depmap_GD_long, condition, AM.score, AM_class)


#reshape = depmap_GD_long is too wide, to stack the TF downward its split into 2 columns (TF and Dependency)
      # check with dim() and head(), the coloums should have increased significantly

  depmap_GD_long <- pivot_longer(
      depmap_GD_long,
      cols = -c(condition, AM.score, AM_class),
      names_to = "TF",
      values_to = "dependency"
  )

#filter to the 763 TFs of interest (tf_list made in Aim 1.1)
  depmap_GD_long <- filter(depmap_GD_long, TF %in% tf_list)

#remove rows with no AM classification - as these cant be compared, check with dim and table
  depmap_GD_long <- filter(depmap_GD_long, !is.na(AM_class))


#test placeholder TF = MYCN
  mycn_GD_data <- filter(depmap_GD_long, TF == "MYCN", AM_class %in% c("ADRN", "MES"))

    #wilcox.test
      wilcox.test(
        mycn_GD_data$dependency[mycn_GD_data$AM_class == "ADRN"],
        mycn_GD_data$dependency[mycn_GD_data$AM_class == "MES"]
      )

    #group medians
      # adrn_median <- median(mycn_GD_data$dependency[mycn_GD_data$AM_class == "ADRN"])
      # adrn_median
      # mes_median <- median(mycn_GD_data$dependency[mycn_GD_data$AM_class == "MES"])
      # mes_median

      #effect size
      # effect_size <- adrn_median - mes_median
      # effect_size


#scaling up to all TFs, check with dim and head, failed initially as there were TFs that didn't have enough data.

  aim2.1_GD_results <- summarise(
    group_by(depmap_GD_long, TF),
    p_value = tryCatch(
      wilcox.test(dependency[AM_class == "ADRN"], dependency[AM_class == "MES"])$p.value,
      error = function(e) NA
    ),
    adrn_median = median(dependency[AM_class == "ADRN"],na.rm = TRUE),
    mes_median = median(dependency[AM_class == "MES"],na.rm = TRUE),
    effect_size = adrn_median - mes_median
    )
  sum(is.na(aim2.1_GD_results$effect_size) & !is.na(aim2.1_GD_results$p_value))
  #sum(is.na(aim2.1_GD_results$p_value)) to check how many TFs were skipped with tryCatch
    #to check which TFs were skipped (no NA dependency Values)
                # na_tfs <- filter(aim2.1_GD_results, is.na(p_value))
                # na_tfs
                # TF    p_value adrn_median mes_median effect_size
                # <chr>   <dbl>       <dbl>      <dbl>       <dbl>
                #   1 FOXD1      NA          NA         NA          NA
                # 2 HIPK2      NA          NA         NA          NA
                # 3 MLXIP      NA          NA         NA          NA
                # 4 NR2E3      NA          NA         NA          NA
                # 5 SMAD5      NA          NA         NA          NA
                # 6 ZBED1      NA          NA         NA          NA

    #TFs entirely abset as columns in depmap_GD
                # missing_tfs <- setdiff(tf_list, aim2.1_GD_results$TF)
                # missing_tfs
                # [1] "AP1"    "FOXO6"  "NFKB"   "NKRF"   "NKX3.2" "NKX6.3" "TBPL2"
                # length(missing_tfs)
                # [1] 7

#sort by p-value
  aim2.1_GD_results_sorted <- arrange(aim2.1_GD_results, p_value)
  #search for specific TFs
  #filter(aim2.1_GD_results, TF %in% c("MYCN", "PHOX2B", "HAND2"))



                ### Aim 2.2 - Correlation Analysis ###

#test placeholder TF = MYCN
  mycn_GD_cor_data <- filter(depmap_GD_long, TF == "MYCN")

  #runs correlation between activity values and AM.score
  #spearmean method ranks the values first rather then using raw numbers
  cor.test(mycn_GD_cor_data$dependency, mycn_GD_cor_data$AM.score, method = "spearman")


#scaling up to all TFs, check with dim and head
  aim2.2_correlation_results <- summarise(
    group_by(depmap_GD_long, TF),
    cor_result = tryCatch(
      list(cor.test(dependency, AM.score, method = "spearman")),
      error = function(e) list(NULL)
    ),
    rho = if (!is.null(cor_result[[1]])) cor_result[[1]]$estimate else NA,
    p_value = if (!is.null(cor_result[[1]])) cor_result[[1]]$p.value else NA
  )

  #removing the cor_result column
  aim2.2_correlation_results <- select(aim2.2_correlation_results, -cor_result)


  #combine aim2.1 and aim2.2, join by the TF column and rename the p_value column as there are 2
  aim2_combined <- rename(
    left_join(aim2.1_GD_results, aim2.2_correlation_results, by = "TF"),
    wilcox_p_value = p_value.x,
    cor_p_value = p_value.y
  )


#consistent and inconsistent dependency TFs (same definitions as Aim 1.4 and 1.5), used for Table 4
  #consistent: significant in BOTH the differential test and the correlation test
  consistent_gd_tfs <- filter(aim2_combined,
                              wilcox_p_value < 0.05,
                              cor_p_value < 0.05
  )
  arrange(consistent_gd_tfs, desc(abs(rho)))
  dim(consistent_gd_tfs)

  sum(aim2_combined$wilcox_p_value < 0.05, na.rm = TRUE)
  sum(aim2_combined$wilcox_p_value < 0.05 & aim2_combined$cor_p_value < 0.05, na.rm = TRUE)

  #inconsistent: significant in the differential test but NOT the correlation test
  inconsistent_gd_tfs <- filter(aim2_combined,
                                wilcox_p_value < 0.05,
                                cor_p_value >= 0.05
  )
  arrange(inconsistent_gd_tfs, wilcox_p_value)
  dim(inconsistent_gd_tfs)

  sum(aim2_combined$wilcox_p_value < 0.05, na.rm = TRUE)
  sum(aim2_combined$wilcox_p_value < 0.05 & aim2_combined$cor_p_value >= 0.05, na.rm = TRUE)


                       ### Aim 2.3 - Visualisation ###

#Volcano plot: effect size vs statistical significance
   ggplot(aim2_combined, aes(x = effect_size, y = -log10(wilcox_p_value))) +
    geom_point(alpha = 0.6, color = "grey40", size = 0.8) +
     # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
     # geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "blue") +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = -0.3, color = "Effect size threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 0.3, color = "Effect size threshold"), linetype = "dashed") +
    scale_color_manual(
      name = "Reference lines",
      values = c("p = 0.05 threshold" = "red", "Effect size threshold" = "blue")
    ) +
    scale_y_continuous(breaks = seq(0, 4, by = 0.5)) +
    scale_x_continuous(breaks = seq(-2, 2, by = 0.25)) +
    labs(
      title = "Aim 2: Differential TF Dependency, ADRN vs MES",
      x = "Effect Size",
      y = "-log10(p-value)"
    ) +
    theme_minimal()

#Figure 4: volcano plot of differential dependency
#Volcano plot with coloured
  aim2_volcano_coloured <- mutate(aim2_combined,
                                  significance = case_when(
                                    wilcox_p_value < 0.05 & effect_size > 0.3 ~ "MES-dependent",
                                    wilcox_p_value < 0.05 & effect_size < -0.3 ~ "ADRN-dependent",
                                    TRUE ~ "Below threshold"
                                  )
  )


  ggplot(aim2_volcano_coloured , aes(x = effect_size, y = -log10(wilcox_p_value))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = -0.3, color = "Effect size threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 0.3, color = "Effect size threshold"), linetype = "dashed") +
    scale_color_manual(
      name = "Legend",
      values = c(
        "ADRN-dependent" = "palegreen4",
        "MES-dependent" = "mediumorchid3",
        "Below threshold" = "grey70",
        "p = 0.05 threshold" = "red",
        "Effect size threshold" = "blue"
      ),
      breaks = c(
        "ADRN-dependent",
        "MES-dependent",
        "Below threshold",
        "p = 0.05 threshold",
        "Effect size threshold"
      )
    )+
    scale_y_continuous(breaks = seq(0, 4, by = 0.5)) +
    scale_x_continuous(breaks = seq(-2, 2, by = 0.25)) +
    labs(
      title = "Aim 2: Differential TF Dependency, ADRN vs MES",
      x = "Effect Size (ADRN − MES median dependency)",
      y = "-log10(p-value)"
    ) +
    theme_minimal()


           ### Aim 2.4 - Differential Dependencies by MYCN Status ###

#MYCN amplification is unevenly spread across cell states (most ADRN lines are amplified,
#most MES lines are not), so the differential dependencies are also summarised by MYCN status
#to check whether a result reflects cell state or MYCN amplification

#cell lines with dependency data and an ADRN/MES class
  gd_lines <- filter(depmap_meta, condition %in% depmap_GD$condition, !is.na(AM_class))

#number of lines in each cell state by MYCN status
  table(gd_lines$AM_class, gd_lines$MYCN_status)

#the five differentially dependent TFs from Aim 2.1, plus MYCN for comparison
  mycn_check_tfs <- c("MYC", "KAT5", "EGR2", "SOX11", "MLX", "MYCN")

#add MYCN status to the long dependency table
  candidate_GD <- filter(depmap_GD_long, TF %in% mycn_check_tfs)
  candidate_GD <- left_join(candidate_GD, select(depmap_meta, condition, MYCN_status), by = "condition")

#median dependency for each TF by cell state and MYCN status (Table 6)
#descriptive only - groups are too small (1 to 4 lines in three of four groups) to test
  aim2.4_mycn_summary <- summarise(
    group_by(candidate_GD, TF, AM_class, MYCN_status),
    n_lines = sum(!is.na(dependency)),
    median_dependency = median(dependency, na.rm = TRUE),
    .groups = "drop"
  )

  print(aim2.4_mycn_summary, n = 30)


#==================================================================================================
#                                         AIM 3: PROGNOSIS
#==================================================================================================

                     ### Aim 3.0 - Preprocessing ###

#confirm ID column in kocak_TF_activity
  # colnames(kocak_TF_activity)[1:3]
  # [1] "...1"  "ABL1"  "AEBP1"

#confirm ID column and survival columns in kocak_survival
  # colnames(kocak_survival)
    # [1] "condition" "EFS_d"     "EFS_bin"   "OS_d"      "OS_bin"
  # head(kocak_survival)

#confirm ID column in kocak_meta
  # colnames(kocak_meta)[1:3]
  # [1] "condition"   "MYCN_status" "MYCN_bin"

#rename kocak_TF_activity
  kocak_TF_activity <- rename(kocak_TF_activity, condition = ...1)

#combine EFS and OS survival data with kocak_TF_activity
  kocak_TF_long <- left_join(
    kocak_TF_activity,
    select(kocak_survival, condition, EFS_d, EFS_bin, OS_d, OS_bin),
    by = "condition"
  )


#kocak_TF_long is too wide, stack the TF downward into 2 columns (TF and activity)
  #same pivot_longer pattern used for depmap_TF_long in Aim 1
  kocak_TF_long <- pivot_longer(
    kocak_TF_long,
    cols = -c(condition, EFS_d, EFS_bin, OS_d, OS_bin),
    names_to = "TF",
    values_to = "activity"
  )

#filter to just the 763 TFs of interest (tf_list made in Aim 1.1)
  kocak_TF_long <- filter(kocak_TF_long, TF %in% tf_list)


#remove rows with no survival data - separate for EFS and OS since a sample could be missing one but not the other
  kocak_TF_long_EFS <- filter(kocak_TF_long, !is.na(EFS_d), !is.na(EFS_bin))
  kocak_TF_long_OS  <- filter(kocak_TF_long, !is.na(OS_d), !is.na(OS_bin))


              ### Aim 3.1 - Kaplan-Meier Survival Analysis ###

#Figure 5: Kaplan-Meier curves for MYCN activity
#stratify patients into high vs low activity groups using a median split
    #used to gather the log-rank p-value of all TFs and visualization for individual TFs

#--- EFS ---
  mycn_EFS_data <- filter(kocak_TF_long_EFS, TF == "MYCN")
  mycn_EFS_data <- mutate(
    mycn_EFS_data,
    activity_group = if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
  )
  mycn_EFS_km_fit <- survfit(Surv(EFS_d, EFS_bin) ~ activity_group, data = mycn_EFS_data)
  #survdiff(Surv(EFS_d, EFS_bin) ~ activity_group, data = mycn_EFS_data)  # log-rank p-value

  ggsurvplot(
    mycn_EFS_km_fit,
    data = mycn_EFS_data,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (days)",
    title = "Aim 3: Kaplan-Meier Curve - MYCN Activity (High vs Low) - Event-Free Survival"
  )

#--- OS ---
  mycn_OS_data <- filter(kocak_TF_long_OS, TF == "MYCN")
  mycn_OS_data <- mutate(
    mycn_OS_data,
    activity_group = if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
  )
  mycn_OS_km_fit <- survfit(Surv(OS_d, OS_bin) ~ activity_group, data = mycn_OS_data)
  # survdiff(Surv(OS_d, OS_bin) ~ activity_group, data = mycn_OS_data)  # log-rank p-value

  ggsurvplot(
    mycn_OS_km_fit,
    data = mycn_OS_data,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (days)",
    title = "Aim 3: Kaplan-Meier Curve - MYCN Activity (High vs Low) - Overall Survival"
  )


#for each TF: split into high/low activity groups by median, run log-rank test

#--- EFS ---
  aim3.1_EFS_logrank_results <- summarise(
    group_by(kocak_TF_long_EFS, TF),
    logrank_p = tryCatch(
      {
        grp <- if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
        survdiff(Surv(EFS_d, EFS_bin) ~ grp)$pvalue
      },
      error = function(e) NA
    )
  )

#--- OS ---
  aim3.1_OS_logrank_results <- summarise(
    group_by(kocak_TF_long_OS, TF),
    logrank_p = tryCatch(
      {
        grp <- if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
        survdiff(Surv(OS_d, OS_bin) ~ grp)$pvalue
      },
      error = function(e) NA
    )
  )


#Kaplan-Meier curves for STOX1 and ZNF143, the strongest prognostic TFs that were not identity TFs (Figures 7 and 8)

#STOX1 - EFS
  stox1_EFS_data <- filter(kocak_TF_long_EFS, TF == "STOX1")
  stox1_EFS_data <- mutate(
    stox1_EFS_data,
    activity_group = if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
  )
  stox1_EFS_km_fit <- survfit(Surv(EFS_d, EFS_bin) ~ activity_group, data = stox1_EFS_data)

  ggsurvplot(
    stox1_EFS_km_fit,
    data = stox1_EFS_data,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (days)",
    title = "Aim 3: Kaplan-Meier Curve - STOX1 Activity (High vs Low) - Event-Free Survival"
  )

#STOX1 - OS
  stox1_OS_data <- filter(kocak_TF_long_OS, TF == "STOX1")
  stox1_OS_data <- mutate(
    stox1_OS_data,
    activity_group = if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
  )
  stox1_OS_km_fit <- survfit(Surv(OS_d, OS_bin) ~ activity_group, data = stox1_OS_data)

  ggsurvplot(
    stox1_OS_km_fit,
    data = stox1_OS_data,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (days)",
    title = "Aim 3: Kaplan-Meier Curve - STOX1 Activity (High vs Low) - Overall Survival"
  )


#ZNF143 - EFS
  znf143_EFS_data <- filter(kocak_TF_long_EFS, TF == "ZNF143")
  znf143_EFS_data <- mutate(
    znf143_EFS_data,
    activity_group = if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
  )
  znf143_EFS_km_fit <- survfit(Surv(EFS_d, EFS_bin) ~ activity_group, data = znf143_EFS_data)

  ggsurvplot(
    znf143_EFS_km_fit,
    data = znf143_EFS_data,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (days)",
    title = "Aim 3: Kaplan-Meier Curve - ZNF143 Activity (High vs Low) - Event-Free Survival"
  )

#ZNF143 - OS
  znf143_OS_data <- filter(kocak_TF_long_OS, TF == "ZNF143")
  znf143_OS_data <- mutate(
    znf143_OS_data,
    activity_group = if_else(activity > median(activity, na.rm = TRUE), "High", "Low")
  )
  znf143_OS_km_fit <- survfit(Surv(OS_d, OS_bin) ~ activity_group, data = znf143_OS_data)

  ggsurvplot(
    znf143_OS_km_fit,
    data = znf143_OS_data,
    pval = TRUE,
    risk.table = TRUE,
    xlab = "Time (days)",
    title = "Aim 3: Kaplan-Meier Curve - ZNF143 Activity (High vs Low) - Overall Survival"
  )



           ### Aim 3.2 - Cox Proportional Hazards Modelling ###

#treat activity as a continuous predictor

#--- EFS ---
  mycn_EFS_cox_fit <- coxph(Surv(EFS_d, EFS_bin) ~ activity, data = mycn_EFS_data)
  summary(mycn_EFS_cox_fit)   # gives hazard ratio, CI, and p-value

#--- OS ---
  mycn_OS_cox_fit <- coxph(Surv(OS_d, OS_bin) ~ activity, data = mycn_OS_data)
  summary(mycn_OS_cox_fit)


#for each TF: fit Cox model, extract hazard ratio, CI and p-value
  #tryCatch used as in Aim 2, since some TFs may fail to converge or have insufficient data

#--- EFS ---
  aim3.2_EFS_cox_results <- summarise(
    group_by(kocak_TF_long_EFS, TF),
    cox_result = tryCatch(
      list(coxph(Surv(EFS_d, EFS_bin) ~ activity)),
      error = function(e) list(NULL)
    ),
    hazard_ratio = if (!is.null(cox_result[[1]])) exp(coef(cox_result[[1]])) else NA,
    p_value = if (!is.null(cox_result[[1]])) summary(cox_result[[1]])$coefficients[, "Pr(>|z|)"] else NA,
    ci_lower = if (!is.null(cox_result[[1]])) exp(confint(cox_result[[1]]))[1] else NA,
    ci_upper = if (!is.null(cox_result[[1]])) exp(confint(cox_result[[1]]))[2] else NA
  )
  aim3.2_EFS_cox_results <- select(aim3.2_EFS_cox_results, -cox_result)
  aim3.2_EFS_cox_results_sorted <- arrange(aim3.2_EFS_cox_results, p_value)

#--- OS ---
  aim3.2_OS_cox_results <- summarise(
    group_by(kocak_TF_long_OS, TF),
    cox_result = tryCatch(
      list(coxph(Surv(OS_d, OS_bin) ~ activity)),
      error = function(e) list(NULL)
    ),
    hazard_ratio = if (!is.null(cox_result[[1]])) exp(coef(cox_result[[1]])) else NA,
    p_value = if (!is.null(cox_result[[1]])) summary(cox_result[[1]])$coefficients[, "Pr(>|z|)"] else NA,
    ci_lower = if (!is.null(cox_result[[1]])) exp(confint(cox_result[[1]]))[1] else NA,
    ci_upper = if (!is.null(cox_result[[1]])) exp(confint(cox_result[[1]]))[2] else NA
  )
  aim3.2_OS_cox_results <- select(aim3.2_OS_cox_results, -cox_result)
  aim3.2_OS_cox_results_sorted <- arrange(aim3.2_OS_cox_results, p_value)

#check how many TFs failed / were skipped
  sum(is.na(aim3.2_EFS_cox_results$p_value))
  sum(is.na(aim3.2_OS_cox_results$p_value))



#combine Cox and log-rank results for each outcome, join by TF
  aim3_EFS_combined <- left_join(aim3.2_EFS_cox_results, aim3.1_EFS_logrank_results, by = "TF")
  aim3_OS_combined  <- left_join(aim3.2_OS_cox_results, aim3.1_OS_logrank_results, by = "TF")

  #sort by Cox p-value
  aim3_EFS_combined_sorted <- arrange(aim3_EFS_combined, p_value)
  aim3_OS_combined_sorted  <- arrange(aim3_OS_combined, p_value)



#search for specific TFs (the Aim 1 identity TFs, Table 8)
  # filter(aim3_EFS_combined, TF %in% c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN"))
          # # A tibble: 5 × 6
          # TF     hazard_ratio  p_value ci_lower ci_upper logrank_p
          # <chr>         <dbl>    <dbl>    <dbl>    <dbl>     <dbl>
          #   1 ASCL1         0.873 1.17e- 1    0.737    1.03   4.71e- 2
          # 2 FEV           0.739 3.95e- 4    0.625    0.873  1.23e- 2
          # 3 MYCN          2.09  1.24e-19    1.78     2.45   1.99e-13
          # 4 PHOX2B        0.933 4.05e- 1    0.793    1.10   5.67e- 1
          # 5 ZEB2          1.56  1.61e- 8    1.34     1.82   3.84e- 5
  # filter(aim3_OS_combined, TF %in% c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN"))
          # A tibble: 5 × 6
          # TF     hazard_ratio  p_value ci_lower ci_upper logrank_p
          # <chr>         <dbl>    <dbl>    <dbl>    <dbl>     <dbl>
          #   1 ASCL1         0.943 5.78e- 1    0.765    1.16   3.04e- 1
          # 2 FEV           0.645 3.60e- 5    0.524    0.794  1.16e- 3
          # 3 MYCN          2.74  5.77e-21    2.22     3.38   5.44e-11
          # 4 PHOX2B        1.00  9.84e- 1    0.819    1.23   9.84e- 1
          # 5 ZEB2          1.74  7.56e- 9    1.44     2.10   1.19e- 5


#How correlated are EFS and OS hazard ratios?
  cor.test(aim3_EFS_combined$hazard_ratio, aim3_OS_combined$hazard_ratio, method = "spearman")


                    ### Aim 3.3 - Visualisation ###
#Figure 6: volcano plots of Cox hazard ratios
  #Volcano plot: log(hazard ratio) vs significance, coloured - EFS
  aim3_EFS_volcano_coloured <- mutate(aim3_EFS_combined,
                                      significance = case_when(
                                        p_value < 0.05 & hazard_ratio > 1 ~ "Poor prognosis (higher risk)",
                                        p_value < 0.05 & hazard_ratio < 1 ~ "Favourable prognosis (lower risk)",
                                        TRUE ~ "Not significant"
                                      )
  )



  ggplot(aim3_EFS_volcano_coloured, aes(x = log(hazard_ratio), y = -log10(p_value))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 0, color = "HR = 1 (no effect)"), linetype = "dashed") +
    scale_color_manual(
      name = "Legend",
      values = c(
        "Poor prognosis (higher risk)" = "indianred1",
        "Favourable prognosis (lower risk)" = "turquoise",
        "Not significant" = "grey70",
        "p = 0.05 threshold" = "red",
        "HR = 1 (no effect)" = "blue"
      ),
      breaks = c(
        "Poor prognosis (higher risk)",
        "Favourable prognosis (lower risk)",
        "Not significant",
        "p = 0.05 threshold",
        "HR = 1 (no effect)"
      )
    ) +
    labs(
      title = "Aim 3: TF Activity & Event-Free Survival (Cox Model), Kocak Cohort",
      x = "log(Hazard Ratio)",
      y = "-log10(p-value)"
    ) +
    theme_minimal()


#Volcano plot: log(hazard ratio) vs significance, coloured - OS
  aim3_OS_volcano_coloured <- mutate(aim3_OS_combined,
                                     significance = case_when(
                                       p_value < 0.05 & hazard_ratio > 1 ~ "Poor prognosis (higher risk)",
                                       p_value < 0.05 & hazard_ratio < 1 ~ "Favourable prognosis (lower risk)",
                                       TRUE ~ "Not significant"
                                     )
  )


  ggplot(aim3_OS_volcano_coloured, aes(x = log(hazard_ratio), y = -log10(p_value))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    geom_hline(aes(yintercept = -log10(0.05), color = "p = 0.05 threshold"), linetype = "dashed") +
    geom_vline(aes(xintercept = 0, color = "HR = 1 (no effect)"), linetype = "dashed") +
    scale_color_manual(
      name = "Legend",
      values = c(
        "Poor prognosis (higher risk)" = "indianred1",
        "Favourable prognosis (lower risk)" = "turquoise",
        "Not significant" = "grey70",
        "p = 0.05 threshold" = "red",
        "HR = 1 (no effect)" = "blue"
      ),
      breaks = c(
        "Poor prognosis (higher risk)",
        "Favourable prognosis (lower risk)",
        "Not significant",
        "p = 0.05 threshold",
        "HR = 1 (no effect)"
      )
    ) +
    labs(
      title = "Aim 3: TF Activity & Overall Survival (Cox Model), Kocak Cohort",
      x = "log(Hazard Ratio)",
      y = "-log10(p-value)"
    ) +
    theme_minimal()



          ### Aim 3.4 - MYCN-Adjusted Cox Models for Key TFs ###

  #MYCN amplification is the single biggest known clinical confounder in NB survival, and is also
  #strongly linked to ADRN identity. A TF's univariate hazard ratio could partly reflect MYCN status
  #rather than an independent prognostic effect, so this section checks whether top candidate TFs
  #hold up once MYCN_bin is added as a covariate. Not run across all 756 TFs (too heavy for the
  #scope of this project) - only for a selected set of key candidates.

  #Data validation identified that the pre-existing MYCN_bin variable in the
  #Kocak metadata was coded with 1 representing non-amplified status;
  #a relabeled variable (MYCN_amplified) was created for clarity in adjusted models
    kocak_meta_mycn_amp <- mutate(kocak_meta, MYCN_amplified = if_else(MYCN_status == "AMP", 1, 0))

#join MYCN_amplified from kocak_meta onto the EFS/OS long tables
  kocak_TF_long_EFS_mycn <- left_join(kocak_TF_long_EFS, select(kocak_meta_mycn_amp, condition, MYCN_amplified), by = "condition")
  kocak_TF_long_OS_mycn  <- left_join(kocak_TF_long_OS, select(kocak_meta_mycn_amp, condition, MYCN_amplified), by = "condition")

  #remove rows with missing MYCN_amplified
    kocak_TF_long_EFS_mycn_adj <- filter(kocak_TF_long_EFS_mycn, !is.na(MYCN_amplified))
    kocak_TF_long_OS_mycn_adj  <- filter(kocak_TF_long_OS_mycn, !is.na(MYCN_amplified))

#define the 12 key TFs to check
  mycn_adjust_tfs <- c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN", "MYC", "HAND2", "STOX1", "ZNF395", "ZNF143", "E2F1", "RB1")

#test placeholder TF = PHOX2B

  #--- EFS, unadjusted vs adjusted ---
  phox2b_EFS_adj_data <- filter(kocak_TF_long_EFS_mycn_adj, TF == "PHOX2B")

  phox2b_EFS_unadjusted <- coxph(Surv(EFS_d, EFS_bin) ~ activity, data = phox2b_EFS_adj_data)
  phox2b_EFS_adjusted   <- coxph(Surv(EFS_d, EFS_bin) ~ activity + MYCN_amplified, data = phox2b_EFS_adj_data)

  # summary(phox2b_EFS_unadjusted)
              # Call:
              #   coxph(formula = Surv(EFS_d, EFS_bin) ~ activity, data = phox2b_EFS_adj_data)
              #
              # n= 475, number of events= 140
              #
              # coef exp(coef) se(coef)      z Pr(>|z|)
              # activity -0.06914   0.93320  0.08304 -0.833    0.405
              #
              # exp(coef) exp(-coef) lower .95 upper .95
              # activity    0.9332      1.072     0.793     1.098
              #
              # Concordance= 0.522  (se = 0.024 )
              # Likelihood ratio test= 0.69  on 1 df,   p=0.4
              # Wald test            = 0.69  on 1 df,   p=0.4
              # Score (logrank) test = 0.69  on 1 df,   p=0.4

  # summary(phox2b_EFS_adjusted)
              # Call:
              #   coxph(formula = Surv(EFS_d, EFS_bin) ~ activity + MYCN_amplified, data = phox2b_EFS_adj_data)
              #
              # n= 475, number of events= 140
              #
              # coef exp(coef) se(coef)     z Pr(>|z|)
              # activity       0.02213   1.02237  0.08596 0.257    0.797
              # MYCN_amplified 1.43193   4.18678  0.18861 7.592 3.15e-14 ***
              #   ---
              #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
              #
              # exp(coef) exp(-coef) lower .95 upper .95
              # activity           1.022     0.9781    0.8639     1.210
              # MYCN_amplified     4.187     0.2388    2.8929     6.059
              #
              # Concordance= 0.602  (se = 0.027 )
              # Likelihood ratio test= 48.13  on 2 df,   p=4e-11
              # Wald test            = 59.42  on 2 df,   p=1e-13
              # Score (logrank) test = 69.78  on 2 df,   p=7e-16



#scaling up to the 12 key TFs

#--- EFS: unadjusted
  aim3.4_EFS_unadjusted <- summarise(
    group_by(filter(kocak_TF_long_EFS_mycn_adj, TF %in% mycn_adjust_tfs), TF),
    cox_result = tryCatch(
      list(coxph(Surv(EFS_d, EFS_bin) ~ activity)),
      error = function(e) list(NULL)
    ),
    hazard_ratio_unadj = if (!is.null(cox_result[[1]])) exp(coef(cox_result[[1]])) else NA,
    p_value_unadj = if (!is.null(cox_result[[1]])) summary(cox_result[[1]])$coefficients[, "Pr(>|z|)"] else NA
  )
  aim3.4_EFS_unadjusted <- select(aim3.4_EFS_unadjusted, -cox_result)

#--- EFS: adjusted for MYCN_amplified ---
  aim3.4_EFS_adjusted <- summarise(
    group_by(filter(kocak_TF_long_EFS_mycn_adj, TF %in% mycn_adjust_tfs), TF),
    cox_result = tryCatch(
      list(coxph(Surv(EFS_d, EFS_bin) ~ activity + MYCN_amplified)),
      error = function(e) list(NULL)
    ),
    hazard_ratio_adj = if (!is.null(cox_result[[1]])) exp(coef(cox_result[[1]])["activity"]) else NA,
    p_value_adj = if (!is.null(cox_result[[1]])) summary(cox_result[[1]])$coefficients["activity", "Pr(>|z|)"] else NA
  )
  aim3.4_EFS_adjusted <- select(aim3.4_EFS_adjusted, -cox_result)

  #combine unadjusted and adjusted EFS results side by side
  aim3.4_EFS_comparison <- left_join(aim3.4_EFS_unadjusted, aim3.4_EFS_adjusted, by = "TF")
  aim3.4_EFS_comparison_sorted <- arrange(aim3.4_EFS_comparison, p_value_adj)

 #--- OS: unadjusted ---
  aim3.4_OS_unadjusted <- summarise(
    group_by(filter(kocak_TF_long_OS_mycn_adj, TF %in% mycn_adjust_tfs), TF),
    cox_result = tryCatch(
      list(coxph(Surv(OS_d, OS_bin) ~ activity)),
      error = function(e) list(NULL)
    ),
    hazard_ratio_unadj = if (!is.null(cox_result[[1]])) exp(coef(cox_result[[1]])) else NA,
    p_value_unadj = if (!is.null(cox_result[[1]])) summary(cox_result[[1]])$coefficients[, "Pr(>|z|)"] else NA
  )
  aim3.4_OS_unadjusted <- select(aim3.4_OS_unadjusted, -cox_result)

#--- OS: adjusted for MYCN_amplified ---
  aim3.4_OS_adjusted <- summarise(
    group_by(filter(kocak_TF_long_OS_mycn_adj, TF %in% mycn_adjust_tfs), TF),
    cox_result = tryCatch(
      list(coxph(Surv(OS_d, OS_bin) ~ activity + MYCN_amplified)),
      error = function(e) list(NULL)
    ),
    hazard_ratio_adj = if (!is.null(cox_result[[1]])) exp(coef(cox_result[[1]])["activity"]) else NA,
    p_value_adj = if (!is.null(cox_result[[1]])) summary(cox_result[[1]])$coefficients["activity", "Pr(>|z|)"] else NA
  )
  aim3.4_OS_adjusted <- select(aim3.4_OS_adjusted, -cox_result)

  #combine unadjusted and adjusted OS results side by side
  aim3.4_OS_comparison <- left_join(aim3.4_OS_unadjusted, aim3.4_OS_adjusted, by = "TF")
  aim3.4_OS_comparison_sorted  <- arrange(aim3.4_OS_comparison, p_value_adj)



          ### Aim 3.5 - FEV and ZEB2 Activity by MYCN Amplification Status ###

#add MYCN status to the Kocak activity table and summarise median activity (expect n = 68 AMP, 407 NAMP)
  fev_zeb2_mycn <- left_join(
    filter(kocak_TF_long, TF %in% c("FEV", "ZEB2")),
    select(kocak_meta, condition, MYCN_status),
    by = "condition"
  )

  summarise(
    group_by(fev_zeb2_mycn, TF, MYCN_status),
    n_patients = n(),
    median_activity = round(median(activity, na.rm = TRUE), 2),
    .groups = "drop"
  )


#==================================================================================================
#                                        AIM 4: INTEGRATION
#==================================================================================================


#objects used from earlier aims (all created above in this script)
  # aim1_TF_combined   (Aim 1.1 + 1.2, TF activity: effect_size, rho, p-values)
  # aim1_GE_combined   (Aim 1.1 + 1.2, GE expression: effect_size, rho, p-values)
  # aim2_combined      (Aim 2.1 + 2.2, dependency: effect_size, rho, p-values)
  # aim3_EFS_combined  (Aim 3.1 + 3.2, EFS: hazard_ratio, p_value, logrank_p)
  # aim3_OS_combined   (Aim 3.1 + 3.2, OS: hazard_ratio, p_value, logrank_p)

                 ### Aim 4.0 - Consolidate & Standardise ###

#master TF list - aim1_TF_combined already reflects the full 763 TFs of interest
  master_tf_list <- unique(aim1_TF_combined$TF)

#direction convention check:
  #effect_size = adrn_median - mes_median
  #  Aim 1 (activity/expression): positive = higher in ADRN
  #  Aim 2 (dependency): positive = MES more dependent (lower dependency score = more dependent)
  #rho is a correlation against AM.score (MES - ADRN), so it has the opposite polarity to effect_size
  #flip rho so it points the same way as effect_size within each aim

  aim1_TF_combined <- mutate(aim1_TF_combined, rho_flipped = -1 * rho)
  aim1_GE_combined <- mutate(aim1_GE_combined, rho_flipped = -1 * rho)
  aim2_combined    <- mutate(aim2_combined,    rho_flipped = -1 * rho)

#Aim 3 hazard ratio is multiplicative and has no ADRN/MES polarity - convert to log scale so it's
#symmetric around 0 and combinable with the other z-scored metrics, keeping its own natural sign
#(positive = poor prognosis, negative = favourable prognosis)
  aim3_EFS_combined <- mutate(aim3_EFS_combined, log_hr = log(hazard_ratio))
  aim3_OS_combined  <- mutate(aim3_OS_combined,  log_hr = log(hazard_ratio))

#z-score helper - standardises a metric to mean 0, sd 1 across all TFs so effect sizes,
#correlations and hazard ratios become comparable on the same scale
  z_score <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }

#standardise each metric within its own table
  aim1_TF_combined <- mutate(aim1_TF_combined,
                             z_effect_size = z_score(effect_size),
                             z_rho = z_score(rho_flipped)
  )

  aim1_GE_combined <- mutate(aim1_GE_combined,
                             z_effect_size = z_score(effect_size),
                             z_rho = z_score(rho_flipped)
  )

  aim2_combined <- mutate(aim2_combined,
                          z_effect_size = z_score(effect_size),
                          z_rho = z_score(rho_flipped)
  )

  aim3_EFS_combined <- mutate(aim3_EFS_combined, z_log_hr = z_score(log_hr))
  aim3_OS_combined  <- mutate(aim3_OS_combined,  z_log_hr = z_score(log_hr))

#build a single master table, joining all standardised metrics by TF onto the full TF list
  aim4_master <- data.frame(TF = master_tf_list)

  aim4_master <- left_join(
    aim4_master,
    select(aim1_TF_combined, TF, aim1_TF_z_effect = z_effect_size, aim1_TF_z_rho = z_rho),
    by = "TF"
  )
  aim4_master <- left_join(
    aim4_master,
    select(aim1_GE_combined, TF, aim1_GE_z_effect = z_effect_size, aim1_GE_z_rho = z_rho),
    by = "TF"
  )
  aim4_master <- left_join(
    aim4_master,
    select(aim2_combined, TF, aim2_z_effect = z_effect_size, aim2_z_rho = z_rho),
    by = "TF"
  )
  aim4_master <- left_join(
    aim4_master,
    select(aim3_EFS_combined, TF, aim3_EFS_z = z_log_hr),
    by = "TF"
  )
  aim4_master <- left_join(
    aim4_master,
    select(aim3_OS_combined, TF, aim3_OS_z = z_log_hr),
    by = "TF"
  )

#check how many TFs are missing Aim 3 data (expect 7: FIGLA, HOXD8, IRF2BPL, NACC1, NPAS4, SP6, VSX1)
  sum(is.na(aim4_master$aim3_EFS_z))

  #exclude TFs with incomplete data - the consensus score needs to compare TFs on equal footing
  #across all three aims, so a TF missing Aim 3 data isn't silently averaged over fewer aims and
  #treated as if it had the same evidence base as a TF with full 3-aim support
  incomplete_tfs <- filter(aim4_master, is.na(aim3_EFS_z) | is.na(aim3_OS_z))$TF
  incomplete_tfs
  length(incomplete_tfs)

  aim4_master <- filter(aim4_master, !is.na(aim3_EFS_z) & !is.na(aim3_OS_z))

  #confirm the exclusion worked - should now be 756 TFs, all with complete data across all three aims
  nrow(aim4_master)
  sum(is.na(aim4_master$aim3_EFS_z))


                      ### Aim 4.1 - Aim-level Scores ###

#Aim 1 score: average of TF activity effect size + rho and GE effect size + rho (4 metrics, equal weight)
  #signed: positive = ADRN-associated, negative = MES-associated
  aim4_master <- mutate(aim4_master,
                        aim1_score = rowMeans(
                          cbind(aim1_TF_z_effect, aim1_TF_z_rho, aim1_GE_z_effect, aim1_GE_z_rho),
                          na.rm = TRUE
                        )
  )

#Aim 2 score: average of dependency effect size + rho (2 metrics, equal weight)
  #signed: negative = ADRN-dependent, positive = MES-dependent
  aim4_master <- mutate(aim4_master,
                        aim2_score = rowMeans(
                          cbind(aim2_z_effect, aim2_z_rho),
                          na.rm = TRUE
                        )
  )

#Aim 3 score: average of EFS + OS log(HR) z-scores (2 metrics, equal weight)
  #kept in its own natural sign: positive = poor prognosis, negative = favourable prognosis
  aim4_master <- mutate(aim4_master,
                        aim3_score = rowMeans(
                          cbind(aim3_EFS_z, aim3_OS_z),
                          na.rm = TRUE
                        )
  )

#Cap Outlier Z-scores

  #winsorize helper - caps any z-score beyond +/- 3 SD to +/- 3, so a single extreme value
  #can't dominate an aim-level score on its own. TFs are kept in the table

  outlier_cap <- 3

  winsorize_z <- function(x, cap = outlier_cap) {
    x_capped <- pmin(x, cap)
    x_capped <- pmax(x_capped, -cap)
    x_capped
  }

  #check which TFs actually get capped, and on which metric, before applying it
  filter(aim4_master,
         abs(aim1_TF_z_effect) > outlier_cap | abs(aim1_TF_z_rho) > outlier_cap |
           abs(aim1_GE_z_effect) > outlier_cap | abs(aim1_GE_z_rho) > outlier_cap |
           abs(aim2_z_effect) > outlier_cap | abs(aim2_z_rho) > outlier_cap |
           abs(aim3_EFS_z) > outlier_cap | abs(aim3_OS_z) > outlier_cap
  )$TF

  #second version of the master table, with the 8 individual metrics capped and the aim-level
  #scores recomputed from those capped metrics
  aim4_master_capped <- mutate(aim4_master,
                               aim1_TF_z_effect = winsorize_z(aim1_TF_z_effect),
                               aim1_TF_z_rho    = winsorize_z(aim1_TF_z_rho),
                               aim1_GE_z_effect = winsorize_z(aim1_GE_z_effect),
                               aim1_GE_z_rho    = winsorize_z(aim1_GE_z_rho),
                               aim2_z_effect    = winsorize_z(aim2_z_effect),
                               aim2_z_rho       = winsorize_z(aim2_z_rho),
                               aim3_EFS_z       = winsorize_z(aim3_EFS_z),
                               aim3_OS_z        = winsorize_z(aim3_OS_z)
  )

  aim4_master_capped <- mutate(aim4_master_capped,
                               aim1_score = rowMeans(
                                 cbind(aim1_TF_z_effect, aim1_TF_z_rho, aim1_GE_z_effect, aim1_GE_z_rho),
                                 na.rm = TRUE
                               ),
                               aim2_score = rowMeans(
                                 cbind(aim2_z_effect, aim2_z_rho),
                                 na.rm = TRUE
                               ),
                               aim3_score = rowMeans(
                                 cbind(aim3_EFS_z, aim3_OS_z),
                                 na.rm = TRUE
                               )
  )



#Weighted Aim-level Scores

  #rho is weighted higher than effect_size within Aim 1 and Aim 2, because rho uses the full
  #continuous AM-score axis while effect_size collapses cell lines into two groups at an
  #arbitrary ADRN/MES threshold. Aims stay equally weighted against each other (identity, viability and
  #prognosis are co-equal objectives of the integration), and Aim 3's EFS/OS pair stays equal
  #weight since neither one relies on a threshold split the way effect_size does.

  effect_weight <- 0.4
  rho_weight    <- 0.6


  weighted_avg2 <- function(x1, x2, w1 = effect_weight, w2 = rho_weight) {
    numerator   <- ifelse(is.na(x1), 0, w1 * x1) + ifelse(is.na(x2), 0, w2 * x2)
    denominator <- ifelse(is.na(x1), 0, w1) + ifelse(is.na(x2), 0, w2)
    numerator / denominator
  }

  #third version of the master table, with rho weighted higher than effect_size within Aim 1 and Aim 2.
  aim4_master_weighted <- mutate(aim4_master,
                                 aim1_TF_component = weighted_avg2(aim1_TF_z_effect, aim1_TF_z_rho),
                                 aim1_GE_component = weighted_avg2(aim1_GE_z_effect, aim1_GE_z_rho)
  )

  aim4_master_weighted <- mutate(aim4_master_weighted,
                                 aim1_score = rowMeans(
                                   cbind(aim1_TF_component, aim1_GE_component),
                                   na.rm = TRUE
                                 ),
                                 aim2_score = weighted_avg2(aim2_z_effect, aim2_z_rho),
                                 aim3_score = rowMeans(
                                   cbind(aim3_EFS_z, aim3_OS_z),
                                   na.rm = TRUE
                                 )
  )

#Capped AND Weighted Aim-level Scores

  aim4_master_capped_weighted <- mutate(aim4_master_capped,
                                        aim1_TF_component = weighted_avg2(aim1_TF_z_effect, aim1_TF_z_rho),
                                        aim1_GE_component = weighted_avg2(aim1_GE_z_effect, aim1_GE_z_rho)
  )

  aim4_master_capped_weighted <- mutate(aim4_master_capped_weighted,
                                        aim1_score = rowMeans(
                                          cbind(aim1_TF_component, aim1_GE_component),
                                          na.rm = TRUE
                                        ),
                                        aim2_score = weighted_avg2(aim2_z_effect, aim2_z_rho),
                                        aim3_score = rowMeans(
                                          cbind(aim3_EFS_z, aim3_OS_z),
                                          na.rm = TRUE
                                        )
  )



                ### Aim 4.2 - Consensus Score & Sensitivity Check ###

#run on  ORIGINAL master (outliers included)

  #consensus score: unsigned, built from the magnitude of each aim-level score so a TF ranks
  #highly for showing strong signal across identity, viability and prognosis, regardless of
  #which direction (ADRN/MES, poor/favourable) it points
  aim4_master <- mutate(aim4_master,
                        consensus_score = rowMeans(
                          cbind(abs(aim1_score), abs(aim2_score), abs(aim3_score)),
                          na.rm = TRUE
                        )
  )

  #direction annotations for readability - kept separate from the ranking itself
  aim4_master <- mutate(aim4_master,
                        identity_direction = case_when(
                          aim1_score > 0 ~ "ADRN-associated",
                          aim1_score < 0 ~ "MES-associated",
                          TRUE ~ NA_character_
                        ),
                        prognosis_direction = case_when(
                          aim3_score > 0 ~ "Poor prognosis",
                          aim3_score < 0 ~ "Favourable prognosis",
                          TRUE ~ NA_character_
                        )
  )

  #sort by consensus score, highest = strongest cross-aim evidence
    aim4_consensus_sorted <- arrange(aim4_master, desc(consensus_score))

    aim_grouped_top20 <- head(arrange(aim4_master, desc(consensus_score)), 20)$TF



#run on CAPPED master (outliers capped)

  aim4_master_capped <- mutate(aim4_master_capped,
                               consensus_score = rowMeans(
                                 cbind(abs(aim1_score), abs(aim2_score), abs(aim3_score)),
                                 na.rm = TRUE
                               )
  )

  aim4_master_capped <- mutate(aim4_master_capped,
                               identity_direction = case_when(
                                 aim1_score > 0 ~ "ADRN-associated",
                                 aim1_score < 0 ~ "MES-associated",
                                 TRUE ~ NA_character_
                               ),
                               prognosis_direction = case_when(
                                 aim3_score > 0 ~ "Poor prognosis",
                                 aim3_score < 0 ~ "Favourable prognosis",
                                 TRUE ~ NA_character_
                               )
  )

  aim4_consensus_sorted_capped <- arrange(aim4_master_capped, desc(consensus_score))

  aim_grouped_top20_capped <- head(arrange(aim4_master_capped, desc(consensus_score)), 20)$TF



#Compare the ORIGINAL top 20 against the CAPPED top 20
  #this is the direct answer to "how much are the outliers changing the ranking"

  intersect(aim_grouped_top20, aim_grouped_top20_capped)
  length(intersect(aim_grouped_top20, aim_grouped_top20_capped))

  #TFs that were in the original top 20 but dropped out once outliers were capped
  setdiff(aim_grouped_top20, aim_grouped_top20_capped)

  #side-by-side rank comparison for every TF: original rank vs capped rank
  aim4_rank_comparison <- data.frame(
    TF = aim4_master$TF,
    consensus_score = aim4_master$consensus_score,
    consensus_rank = rank(desc(aim4_master$consensus_score), ties.method = "min"),
    consensus_score_capped = aim4_master_capped$consensus_score,
    consensus_rank_capped = rank(desc(aim4_master_capped$consensus_score), ties.method = "min")
  )
  aim4_rank_comparison <- mutate(aim4_rank_comparison,
                                 rank_change = consensus_rank - consensus_rank_capped
  )

  aim4_rank_comparison_sorted <- arrange(aim4_rank_comparison, consensus_rank)

  #TFs whose rank moved the most once outliers were capped - these are the ones an outlier metric was propping up
  aim4_rank_comparison_sorted_by_change <- arrange(aim4_rank_comparison, desc(abs(rank_change)))



#run on WEIGHTED master (outliers included)

  aim4_master_weighted <- mutate(aim4_master_weighted,
                                 consensus_score = rowMeans(
                                   cbind(abs(aim1_score), abs(aim2_score), abs(aim3_score)),
                                   na.rm = TRUE
                                 )
  )

  aim4_master_weighted <- mutate(aim4_master_weighted,
                                 identity_direction = case_when(
                                   aim1_score > 0 ~ "ADRN-associated",
                                   aim1_score < 0 ~ "MES-associated",
                                   TRUE ~ NA_character_
                                 ),
                                 prognosis_direction = case_when(
                                   aim3_score > 0 ~ "Poor prognosis",
                                   aim3_score < 0 ~ "Favourable prognosis",
                                   TRUE ~ NA_character_
                                 )
  )

  aim4_consensus_sorted_weighted <- arrange(aim4_master_weighted, desc(consensus_score))

  aim_grouped_top20_weighted <- head(arrange(aim4_master_weighted, desc(consensus_score)), 20)$TF

#Compare the ORIGINAL top 20 against the WEIGHTED top 20
  #this is the direct answer to "how much DOES weighting rho higher change the ranking"

  intersect(aim_grouped_top20, aim_grouped_top20_weighted)
  length(intersect(aim_grouped_top20, aim_grouped_top20_weighted))

  #TFs that were in the original top 20 but dropped out once rho was weighted higher
  setdiff(aim_grouped_top20, aim_grouped_top20_weighted)

  #side-by-side rank comparison for every TF: original rank vs weighted rank
  aim4_rank_comparison_weighting <- data.frame(
    TF = aim4_master$TF,
    consensus_score = aim4_master$consensus_score,
    consensus_rank = rank(desc(aim4_master$consensus_score), ties.method = "min"),
    consensus_score_weighted = aim4_master_weighted$consensus_score,
    consensus_rank_weighted = rank(desc(aim4_master_weighted$consensus_score), ties.method = "min")
  )
  aim4_rank_comparison_weighting <- mutate(aim4_rank_comparison_weighting,
                                           rank_change = consensus_rank - consensus_rank_weighted
  )

  aim4_rank_comparison_weighting_sorted <- arrange(aim4_rank_comparison_weighting, consensus_rank)

  #TFs whose rank moved the most once rho was weighted higher
  aim4_rank_comparison_weighting_sorted_by_change <- arrange(aim4_rank_comparison_weighting, desc(abs(rank_change)))



#run on CAPPED AND WEIGHTED master (outliers capped)


  aim4_master_capped_weighted <- mutate(aim4_master_capped_weighted,
                                        consensus_score = rowMeans(
                                          cbind(abs(aim1_score), abs(aim2_score), abs(aim3_score)),
                                          na.rm = TRUE
                                        )
  )

  aim4_master_capped_weighted <- mutate(aim4_master_capped_weighted,
                                        identity_direction = case_when(
                                          aim1_score > 0 ~ "ADRN-associated",
                                          aim1_score < 0 ~ "MES-associated",
                                          TRUE ~ NA_character_
                                        ),
                                        prognosis_direction = case_when(
                                          aim3_score > 0 ~ "Poor prognosis",
                                          aim3_score < 0 ~ "Favourable prognosis",
                                          TRUE ~ NA_character_
                                        )
  )

  aim4_consensus_sorted_capped_weighted <- arrange(aim4_master_capped_weighted, desc(consensus_score))

  aim_grouped_top20_capped_weighted <- head(arrange(aim4_master_capped_weighted, desc(consensus_score)), 20)$TF


#Compare the ORIGINAL top 20 against the CAPPED AND WEIGHTED top 20
  #this is the direct answer to "how much do the two adjustments change the ranking together"

  intersect(aim_grouped_top20, aim_grouped_top20_capped_weighted)
  length(intersect(aim_grouped_top20, aim_grouped_top20_capped_weighted))

  #TFs that were in the original top 20 but dropped out once both adjustments were applied
  setdiff(aim_grouped_top20, aim_grouped_top20_capped_weighted)

  #side-by-side rank comparison for every TF: original rank vs capped+weighted rank
  aim4_rank_comparison_capped_weighted <- data.frame(
    TF = aim4_master$TF,
    consensus_score = aim4_master$consensus_score,
    consensus_rank = rank(desc(aim4_master$consensus_score), ties.method = "min"),
    consensus_score_capped_weighted = aim4_master_capped_weighted$consensus_score,
    consensus_rank_capped_weighted = rank(desc(aim4_master_capped_weighted$consensus_score), ties.method = "min")
  )
  aim4_rank_comparison_capped_weighted <- mutate(aim4_rank_comparison_capped_weighted,
                                                 rank_change = consensus_rank - consensus_rank_capped_weighted
  )

  aim4_rank_comparison_capped_weighted_sorted <- arrange(aim4_rank_comparison_capped_weighted, consensus_rank)

  #TFs whose rank moved the most once both adjustments were applied together
  aim4_rank_comparison_capped_weighted_sorted_by_change <- arrange(aim4_rank_comparison_capped_weighted, desc(abs(rank_change)))


#Compare the CAPPED-ONLY ranking against the CAPPED AND WEIGHTED ranking
  #isolates what weighting adds on top of capping, since both already share the same capped base

  intersect(aim_grouped_top20_capped, aim_grouped_top20_capped_weighted)
  length(intersect(aim_grouped_top20_capped, aim_grouped_top20_capped_weighted))



#Compare all FOUR approaches together

  #rank of every TF under each of the four approaches, side by side
  aim4_all_ranks <- data.frame(
    TF = aim4_master$TF,
    rank_baseline = rank(desc(aim4_master$consensus_score), ties.method = "min"),
    rank_capped = rank(desc(aim4_master_capped$consensus_score), ties.method = "min"),
    rank_weighted = rank(desc(aim4_master_weighted$consensus_score), ties.method = "min"),
    rank_capped_weighted = rank(desc(aim4_master_capped_weighted$consensus_score), ties.method = "min")
  )

  aim4_all_ranks_sorted <- arrange(aim4_all_ranks, rank_baseline)

  #TFs that make the top 20 under every single approach - the most robust candidates regardless of methodology choice
  robust_top20_all_methods <- Reduce(
    intersect,
    list(aim_grouped_top20, aim_grouped_top20_capped, aim_grouped_top20_weighted, aim_grouped_top20_capped_weighted)
  )
  robust_top20_all_methods
  length(robust_top20_all_methods)

  #union of every TF that makes the top 20 under at least one approach, flagged by which approach(es) included it
  all_top20_union <- unique(c(
    aim_grouped_top20, aim_grouped_top20_capped, aim_grouped_top20_weighted, aim_grouped_top20_capped_weighted
  ))

  aim4_top20_membership <- data.frame(
    TF = all_top20_union,
    in_baseline = all_top20_union %in% aim_grouped_top20,
    in_capped = all_top20_union %in% aim_grouped_top20_capped,
    in_weighted = all_top20_union %in% aim_grouped_top20_weighted,
    in_capped_weighted = all_top20_union %in% aim_grouped_top20_capped_weighted
  )

  aim4_top20_membership <- mutate(aim4_top20_membership,
                                  n_methods = rowSums(cbind(in_baseline, in_capped, in_weighted, in_capped_weighted))
  )

  #sorted so TFs robust across the most methods appear first, least-agreed-upon TFs last
  aim4_top20_membership_sorted <- arrange(aim4_top20_membership, desc(n_methods))

  #pairwise overlap count between every pair of the four approaches, as a matrix
  top20_by_method <- list(
    baseline = aim_grouped_top20,
    capped = aim_grouped_top20_capped,
    weighted = aim_grouped_top20_weighted,
    capped_weighted = aim_grouped_top20_capped_weighted
  )

  pairwise_overlap_matrix <- matrix(
    NA,
    nrow = 4, ncol = 4,
    dimnames = list(names(top20_by_method), names(top20_by_method))
  )

  for (row_method in names(top20_by_method)) {
    for (col_method in names(top20_by_method)) {
      pairwise_overlap_matrix[row_method, col_method] <- length(
        intersect(top20_by_method[[row_method]], top20_by_method[[col_method]])
      )
    }
  }

  pairwise_overlap_matrix


#Capping Threshold Robustness Test

  #tests whether the chosen approach (capping + weighting) is robust to the exact cap threshold,
  #or whether the specific choice of 3 SD is doing unjustified work. Re-runs the capped+weighted
  #pipeline at 2.5, 3, and 4 SD and compares the resulting rankings, using aim4_master (uncapped,
  #unweighted) as the common starting point each time.

  build_capped_weighted_scores <- function(base_data, cap) {

    winsorize_at_cap <- function(x) {
      pmin(pmax(x, -cap), cap)
    }

    result <- mutate(base_data,
                     aim1_TF_z_effect = winsorize_at_cap(aim1_TF_z_effect),
                     aim1_TF_z_rho    = winsorize_at_cap(aim1_TF_z_rho),
                     aim1_GE_z_effect = winsorize_at_cap(aim1_GE_z_effect),
                     aim1_GE_z_rho    = winsorize_at_cap(aim1_GE_z_rho),
                     aim2_z_effect    = winsorize_at_cap(aim2_z_effect),
                     aim2_z_rho       = winsorize_at_cap(aim2_z_rho),
                     aim3_EFS_z       = winsorize_at_cap(aim3_EFS_z),
                     aim3_OS_z        = winsorize_at_cap(aim3_OS_z)
    )

    result <- mutate(result,
                     aim1_TF_component = weighted_avg2(aim1_TF_z_effect, aim1_TF_z_rho),
                     aim1_GE_component = weighted_avg2(aim1_GE_z_effect, aim1_GE_z_rho)
    )

    result <- mutate(result,
                     aim1_score = rowMeans(cbind(aim1_TF_component, aim1_GE_component), na.rm = TRUE),
                     aim2_score = weighted_avg2(aim2_z_effect, aim2_z_rho),
                     aim3_score = rowMeans(cbind(aim3_EFS_z, aim3_OS_z), na.rm = TRUE)
    )

    result <- mutate(result,
                     consensus_score = rowMeans(cbind(abs(aim1_score), abs(aim2_score), abs(aim3_score)), na.rm = TRUE)
    )

    result
  }

  #run the same pipeline at three different thresholds
  aim4_cap2_5 <- build_capped_weighted_scores(aim4_master, 2.5)
  aim4_cap3   <- build_capped_weighted_scores(aim4_master, 3)
  aim4_cap4   <- build_capped_weighted_scores(aim4_master, 4)

  #top 20 under each threshold
  top20_cap2_5 <- head(arrange(aim4_cap2_5, desc(consensus_score)), 20)$TF
  top20_cap3   <- head(arrange(aim4_cap3, desc(consensus_score)), 20)$TF
  top20_cap4   <- head(arrange(aim4_cap4, desc(consensus_score)), 20)$TF

  #TFs in the top 20 no matter which threshold is used
  robust_top20_all_thresholds <- Reduce(intersect, list(top20_cap2_5, top20_cap3, top20_cap4))
  robust_top20_all_thresholds
  length(robust_top20_all_thresholds)

  #full rank comparison across the three thresholds, for every TF
  aim4_threshold_ranks <- data.frame(
    TF = aim4_master$TF,
    rank_cap2_5 = rank(desc(aim4_cap2_5$consensus_score), ties.method = "min"),
    rank_cap3   = rank(desc(aim4_cap3$consensus_score), ties.method = "min"),
    rank_cap4   = rank(desc(aim4_cap4$consensus_score), ties.method = "min")
  )

  aim4_threshold_ranks_sorted <- arrange(aim4_threshold_ranks, rank_cap3)

  #biological plausibility check: track known/candidate TFs across thresholds specifically,
  #rather than relying on overlap counts alone - does the ranking still make sense
  key_candidate_tfs <- c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN", "MYC", "HAND2", "SOX11", "EGR2")
  filter(aim4_threshold_ranks_sorted, TF %in% key_candidate_tfs)

  #TFs whose rank is most sensitive to the exact threshold chosen (largest spread across the three)
  aim4_threshold_ranks <- mutate(aim4_threshold_ranks,
                                 rank_range = pmax(rank_cap2_5, rank_cap3, rank_cap4) - pmin(rank_cap2_5, rank_cap3, rank_cap4)
  )

  aim4_threshold_ranks_sorted_by_range <- arrange(aim4_threshold_ranks, desc(rank_range))



#SUPPLEMENTARY: Naive vs Aim-Grouped Comparison

  #secondary check only: compares the aim-grouped consensus score with a naive score
  #(unweighted average of all eight metrics, ignoring aim grouping). Reported briefly in the
  #Methodology; it was not used to change the method, which remains the aim-grouped approach.
  
  aim4_master <- mutate(aim4_master,
                        naive_score = rowMeans(
                          cbind(
                            abs(aim1_TF_z_effect), abs(aim1_TF_z_rho),
                            abs(aim1_GE_z_effect), abs(aim1_GE_z_rho),
                            abs(aim2_z_effect), abs(aim2_z_rho),
                            abs(aim3_EFS_z), abs(aim3_OS_z)
                          ),
                          na.rm = TRUE
                        )
  )

  aim4_naive_sorted <- arrange(aim4_master, desc(naive_score))
  naive_top20        <- head(aim4_naive_sorted, 20)$TF

  cor.test(aim4_master$consensus_score, aim4_master$naive_score, method = "spearman")
  intersect(aim_grouped_top20, naive_top20)
  length(intersect(aim_grouped_top20, naive_top20))

  aim4_master_capped <- mutate(aim4_master_capped,
                               naive_score = rowMeans(
                                 cbind(
                                   abs(aim1_TF_z_effect), abs(aim1_TF_z_rho),
                                   abs(aim1_GE_z_effect), abs(aim1_GE_z_rho),
                                   abs(aim2_z_effect), abs(aim2_z_rho),
                                   abs(aim3_EFS_z), abs(aim3_OS_z)
                                 ),
                                 na.rm = TRUE
                               )
  )

  aim4_naive_sorted_capped <- arrange(aim4_master_capped, desc(naive_score))
  naive_top20_capped        <- head(aim4_naive_sorted_capped, 20)$TF

  cor.test(aim4_master_capped$consensus_score, aim4_master_capped$naive_score, method = "spearman")
  intersect(aim_grouped_top20_capped, naive_top20_capped)
  length(intersect(aim_grouped_top20_capped, naive_top20_capped))



  #final consensus ranking table, built from the FINAL chosen method (capped + weighted, 3 SD)
  #rather than the original equal-weighted/uncapped version - this is what all visuals below use
  aim4_final_table <- mutate(aim4_consensus_sorted_capped_weighted,
                             consensus_rank = rank(desc(consensus_score), ties.method = "min")
  )

  aim4_final_table <- select(
    aim4_final_table,
    TF, consensus_rank, consensus_score, identity_direction, prognosis_direction,
    aim1_score, aim2_score, aim3_score
  )

  #note for reporting: 756 TFs total (7 excluded for missing Aim 3 data - see Aim 4.0/4.1)

  #add a dependency direction label from the sign of the Aim 2 score
  #positive = MES lines more dependent, negative = ADRN lines more dependent
  aim4_final_table <- mutate(
    aim4_final_table,
    dependency_direction = ifelse(aim2_score > 0, "MES-dependent", "ADRN-dependent")
  )

  #top 30 with all three directions (Table 11)
  top30 <- arrange(aim4_final_table, consensus_rank)
  top30 <- select(top30, consensus_rank, TF, consensus_score,
                  identity_direction, dependency_direction, prognosis_direction)
  head(top30, 30)



                  ### Aim 4.3 - Visualisation ###
#common colour scale for every heatmap: white = 0, range -3 to 3 (matches the capping threshold)
#circlize is loaded in section 0
  heatmap_colours <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

#column labels state the direction of each colour, so the legend is not needed to read them
  heatmap_column_labels <- c(
    "Identity (Aim 1)\nblue = MES, red = ADRN",
    "Viability (Aim 2)\nblue = MES, red = ADRN",
    "Prognosis (Aim 3)\nblue = favourable, red = poor"
  )

#column labels drawn as a text row under each heatmap, so the two-line labels sit below the last TF
#height is set explicitly so the second line of each label is not cut off
  column_label_annotation <- HeatmapAnnotation(
    labels = anno_text(
      heatmap_column_labels,
      rot = 0,
      just = "center",
      location = 0.5,
      height = unit(2, "cm"),
      gp = gpar(fontsize = 12)
    )
  )

  #NOTE: in every heatmap below, viability (aim2_score) is multiplied by -1 for DISPLAY ONLY.
  #In the scores, positive Aim 2 = MES-dependent; flipping it makes red = ADRN-dependent,
  #so the identity and viability columns use the same colour direction. The ranking is unchanged.

#Figure 10: heatmap of signed aim-level scores for the top 30 TFs by consensus score

  top_tfs <- head(aim4_final_table, 30)$TF

  heatmap_data <- filter(aim4_master_capped_weighted, TF %in% top_tfs)
  heatmap_data <- arrange(heatmap_data, desc(consensus_score))
  heatmap_matrix <- as.matrix(select(heatmap_data, aim1_score, aim2_score, aim3_score))
  heatmap_matrix[, "aim2_score"] <- -1 * heatmap_matrix[, "aim2_score"]
  rownames(heatmap_matrix) <- heatmap_data$TF
  colnames(heatmap_matrix) <- heatmap_column_labels

  Heatmap(
    heatmap_matrix,
    col = heatmap_colours,
    heatmap_legend_param = list(
      title = "Direction\nand\nStrength",
      at = c(-3, -2, -1, 0, 1, 2, 3)
    ),
    na_col = "grey80",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_title = "Transcription Factor",
    column_title = "Aim 4: Top 30 TFs by Consensus Score Across Identity, Viability and Prognosis",
    show_column_names = FALSE,
    bottom_annotation = column_label_annotation
  )


#Figure 11: top 30 strongly ADRN and top 30 strongly MES TFs, ranked by identity score - as two separate heatmaps
  top30_adrn_tfs <- head(arrange(aim4_final_table, desc(aim1_score)), 30)$TF
  top30_mes_tfs  <- head(arrange(aim4_final_table, aim1_score), 30)$TF


  #heatmap 1 of 2: top 30 strongly ADRN TFs
  adrn30_heatmap_data <- filter(aim4_master_capped_weighted, TF %in% top30_adrn_tfs)
  adrn30_heatmap_data <- arrange(adrn30_heatmap_data, desc(aim1_score))
  adrn30_heatmap_matrix <- as.matrix(select(adrn30_heatmap_data, aim1_score, aim2_score, aim3_score))
  adrn30_heatmap_matrix[, "aim2_score"] <- -1 * adrn30_heatmap_matrix[, "aim2_score"]
  rownames(adrn30_heatmap_matrix) <- adrn30_heatmap_data$TF
  colnames(adrn30_heatmap_matrix) <- heatmap_column_labels

  Heatmap(
    adrn30_heatmap_matrix,
    col = heatmap_colours,
    heatmap_legend_param = list(
      title = "Direction\nand\nStrength",
      at = c(-3, -2, -1, 0, 1, 2, 3)
    ),
    na_col = "grey80",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_title = "Transcription Factor",
    column_title = "Aim 4: Top 30 Strongly ADRN TFs",
    show_column_names = FALSE,
    bottom_annotation = column_label_annotation
  )



  #heatmap 2 of 2: top 30 strongly MES TFs
  mes30_heatmap_data <- filter(aim4_master_capped_weighted, TF %in% top30_mes_tfs)
  mes30_heatmap_data <- arrange(mes30_heatmap_data, aim1_score)
  mes30_heatmap_matrix <- as.matrix(select(mes30_heatmap_data, aim1_score, aim2_score, aim3_score))
  mes30_heatmap_matrix[, "aim2_score"] <- -1 * mes30_heatmap_matrix[, "aim2_score"]
  rownames(mes30_heatmap_matrix) <- mes30_heatmap_data$TF
  colnames(mes30_heatmap_matrix) <- heatmap_column_labels

  Heatmap(
    mes30_heatmap_matrix,
    col = heatmap_colours,
    heatmap_legend_param = list(
      title = "Direction\nand\nStrength",
      at = c(-3, -2, -1, 0, 1, 2, 3)
    ),
    na_col = "grey80",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_title = "Transcription Factor",
    column_title = "Aim 4: Top 30 Strongly MES TFs",
    show_column_names = FALSE,
    bottom_annotation = column_label_annotation
  )


#Figure 12: heatmap of the nine key TFs from Aims 1 to 3, regardless of consensus rank
  key_candidate_tfs <- c("PHOX2B", "ASCL1", "FEV", "ZEB2", "MYCN", "MYC", "HAND2", "SOX11", "EGR2")

  candidate_data <- filter(aim4_master_capped_weighted, TF %in% key_candidate_tfs)
  candidate_matrix <- as.matrix(select(candidate_data, aim1_score, aim2_score, aim3_score))
  candidate_matrix[, "aim2_score"] <- -1 * candidate_matrix[, "aim2_score"]
  rownames(candidate_matrix) <- candidate_data$TF
  colnames(candidate_matrix) <- heatmap_column_labels

  Heatmap(
    candidate_matrix,
    col = heatmap_colours,
    heatmap_legend_param = list(
      title = "Direction\nand\nStrength",
      at = c(-3, -2, -1, 0, 1, 2, 3)
    ),
    na_col = "grey80",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_title = "Key Candidate TF",
    column_title = "Aim 4: Key Candidate TFs Across Identity, Viability and Prognosis",
    show_column_names = FALSE,
    bottom_annotation = column_label_annotation
  )


                  ### Aim 4.4 - Additional Visualisation ###

#Figure 9: rank-change slope chart, baseline rank vs final (capped + weighted) rank
  #shows every TF that made the top 20 under either method, connecting its two ranks with a
  #line - visualises the MYC/MYCN swap and SOX11's decline directly, rather than as a table

  slope_tfs <- unique(c(aim_grouped_top20, aim_grouped_top20_capped_weighted))

  slope_data <- filter(aim4_all_ranks, TF %in% slope_tfs)
  slope_data <- select(slope_data, TF, rank_baseline, rank_capped_weighted)

  slope_data_long <- pivot_longer(
    slope_data,
    cols = c(rank_baseline, rank_capped_weighted),
    names_to = "method",
    values_to = "rank"
  )
  slope_data_long <- mutate(slope_data_long,
                            method = case_when(
                              method == "rank_baseline" ~ "Baseline",
                              method == "rank_capped_weighted" ~ "Final (Capped + Weighted)"
                            ),
                            method = factor(method, levels = c("Baseline", "Final (Capped + Weighted)"))
  )

  ggplot(slope_data_long, aes(x = method, y = rank, group = TF)) +
    geom_line(color = "grey60") +
    geom_point(size = 2, color = "steelblue") +
    geom_text(
      data = filter(slope_data_long, method == "Final (Capped + Weighted)"),
      aes(label = TF), hjust = -0.1, size = 3
    ) +
    scale_y_reverse(breaks = seq(0, max(slope_data_long$rank), by = 5)) +
    labs(
      title = "Aim 4: Rank Change from Baseline to Final Method",
      x = NULL, y = "Consensus Rank"
    ) +
    theme_minimal()


#Figure 13: bubble plot of identity vs prognosis for all 756 TFs, sized and coloured by viability strength
  #shows all 756 TFs at once, rather than a pre-filtered top-N list, so the key TFs can
  #be seen in the context of the full landscape

  bubble_data <- mutate(aim4_master_capped_weighted,
                        is_key_candidate = TF %in% key_candidate_tfs,
                        viability_strength = abs(aim2_score)
  )

  ggplot(bubble_data, aes(x = aim1_score, y = aim3_score)) +
    geom_point(
      aes(size = viability_strength, fill = viability_strength),
      shape = 21, color = "white", stroke = 0.2, alpha = 0.75
    ) +
    geom_point(
      data = filter(bubble_data, is_key_candidate),
      aes(size = viability_strength, fill = viability_strength),
      shape = 21, color = "black", stroke = 0.8
    ) +
    geom_text(
      data = filter(bubble_data, is_key_candidate),
      aes(label = TF), size = 3, vjust = -1.3, fontface = "bold"
    ) +
    scale_fill_gradientn(
      name = "Viability\n|Z-score|",
      colors = c("#f7fbff", "#6baed6", "#08306b"),
      breaks = c(0, 1, 2, 3),
      limits = c(0, 3),
      guide = guide_legend(
        override.aes = list(size = c(1, 3.7, 6.3, 9))
      )
    ) +
    scale_size_continuous(
      name = "Viability\n|Z-score|",
      range = c(1, 9),
      breaks = c(0, 1, 2, 3),
      limits = c(0, 3)
    ) +
    labs(
      title = "Aim 4: Identity vs Prognosis Across All 756 TFs",
      subtitle = "*Black-outlined points are key candidate TFs",
      x = "MES <--- Identity Score (Aim 1) ---> ADRN",
      y = "Favourable <--- Prognosis Score (Aim 3) ---> Poor"
    ) +
    theme_minimal()



#==================================================================================================
#   Record of R and package versions used for this run (for reproducibility)
#==================================================================================================
  sessionInfo()

