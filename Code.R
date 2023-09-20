####Text analysis (Figure 2) ####

# Load required libraries
library(tidyverse)
library(tm)
library(wordcloud)
library(dplyr)
library(tidytext)

# Cleanig and processig text
process_text<-function(text){
  text<-gsub("https\\S*", "", text) 
  text<-gsub("@\\S*", "", text) 
  text<-gsub("[\r\n]", "", text)
  text<-gsub("[[:punct:]]", "", text)
  return(text)
}


# Define the GitHub URL of the CSV file
github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/text_antibiotic_IP.csv" # (use _OP postfix for outpatients)

# Read the CSV file from the GitHub URL
text_antibiotic <- read_csv(github_url)

# Tokenize the text data and count drug-word occurrences
drug_words <- text_antibiotic %>%
  unnest_tokens(word, text) %>%
  count(drug, word, sort = TRUE)

# Calculate the total word count per drug
total_words <- drug_words %>%
  group_by(drug) %>%
  summarize(total = sum(n))

# Join the total word count back to the drug_words dataset
drug_words <- left_join(drug_words, total_words)

# Calculate term frequency-inverse document frequency (TF-IDF)
drug_words <- drug_words %>%
  bind_tf_idf(word, drug, n)

# Display the drug_words dataset
drug_words

# Filter data for a specific drug (e.g., ETHAMBUTOL)
df_trimmed <- drug_words %>%
  filter(drug == "ETHAMBUTOL") %>%
  select(-total) %>%
  arrange(desc(tf_idf))

# Generate a word cloud for the selected drug
wordcloud(words = df_trimmed$word, freq = df_trimmed$tf_idf)

#####Causal Impact Analysis (Figure 4)####

# Load required libraries
library(CausalImpact)

#trend_total_comp_tot_IP.csv for inpatient 
#trend_total_comp_tot_OP.csv for outpatients
github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/trend_total_comp_tot_IP.csv"

# Read the CSV file containing the trend data
trend_total <- read_csv(github_url)

# Display unique values in the "antimicrobial" column
unique_antimicrobials <- unique(trend_total$antimicrobial)
cat("Unique Antimicrobials:", unique_antimicrobials, "\n")

# Define a list of example antimicrobials
antibiotics <- c("AZITHROMYCIN", "LINEZOLID", "LEVOFLOXACIN", "CHLOROQUINE")

# Initialize an empty dataframe to store results
output <- data.frame()

# Loop through each antimicrobial
for (k in 1:length(antibiotics)) {
  
  # Subset data for the current antimicrobial
  trend_total_short <- subset(trend_total, antimicrobial == antibiotics[k])
  
  # Loop through different time windows
  for (j in seq(3, 100, 3)) {
    cat("Time Window (weeks):", j, "\n")
    
    # Define pre-period and post-period for Causal Impact analysis (115 is Covid Start week)
    pre.period <- c(115 - j, 115)
    post.period <- c(116, 116 + j)
    
    # Create a 'no_effect' series to simulate no impact
    no_effect <- trend_total_short$trend
    no_effect[116:length(no_effect)] <- no_effect[1:length(no_effect[116:length(no_effect)])]
    
    # Fit ARIMA model and simulate 'no_effect'
    arima_sim <- arima(trend_total_short$trend[(115 - j):115], order = c(0, 0, 0))
    no_effect <- arima.sim(model = list(order = c(0, 0, 0)), n = length(trend_total_short$trend),
                           mean = arima_sim[[1]], sd = sqrt(arima_sim[[2]]))
    
    # Perform Causal Impact analysis
    impact <- CausalImpact(zoo(data.frame(cbind(trend_total_short$trend, no_effect))),
                           pre.period, post.period)
    
    # Plot the Causal Impact analysis
    plot(impact)
    
    # Store the results in the 'output' dataframe
    output <- rbind(output, c(j, impact$summary["Average", 10], impact$summary["Average", 11],
                              impact$summary["Average", 12], antibiotics[k]))
  }
}

# Rename columns in the 'output' dataframe
colnames(output) <- c("week", "value", "lower", "upper", "antibiotic")

# Convert percentages to numeric values
output$value <- as.numeric(as.character(output$value)) * 100
output$lower <- as.numeric(as.character(output$lower)) * 100
output$upper <- as.numeric(as.character(output$upper)) * 100

# Convert 'week' to a factor and reorder levels
output$week <- as.factor(output$week)
output$week <- factor(output$week, levels = as.character(sort(as.numeric(levels(output$week)))))

# Create a bar plot with error bars
library(ggplot2)
ggplot(output, aes(x = week, y = value, fill = antibiotic)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(0.9)) +
  ylab("Effect percentage") +
  xlab("Pre- & Post-Covid Time Window in Weeks") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1),
        axis.text.y = element_text(size = 13, hjust = 1),
        axis.title.x = element_text(color = "black", size = 15, face = "bold"),
        axis.title.y = element_text(color = "black", size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold")) +
  scale_fill_manual(values = c("AZITHROMYCIN" = "#80B1D3", "LINEZOLID" = "#D95F02",
                               "MICONAZOLE" = "#FFED6F", "LEVOFLOXACIN" = "#1B9E77",
                               "CHLOROQUINE" = "#9B4348"))



#### Death Odds Ratio Visulization (Figure 6) ####

# Load required libraries
library(RColorBrewer)
library(tidyverse)

# Read the CSV files
github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/DeathOddsRatio_noConfounding_tot.csv"
barplot_df <- read_csv(github_url)

github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/LogisticRegression_DeathOddsRatio_Strains.csv"
barplot_df_lg <- read_csv(github_url)

github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/SurvivalAnalysis_DeathOddsRatio_Strains.csv"
barplot_df_survival <- read_csv(github_url)

# Preprocess the dataframes
barplot_df_shortened <- barplot_df[barplot_df$pvalue < 0.05, ]
barplot_df_shortened$Antimicrobial <- toupper(barplot_df_shortened$Antimicrobial)

barplot_df_lg_shortened <- barplot_df_lg[which(barplot_df$pvalue < 0.05), ]
barplot_df_lg_shortened$Antimicrobial <- toupper(barplot_df_lg_shortened$Antimicrobial)
barplot_df_lg_shortened <- barplot_df_lg_shortened[complete.cases(barplot_df_lg_shortened), ]

barplot_df_lg_survival <- barplot_df_survival[which(barplot_df_survival$pvalue < 0.05), ]
barplot_df_lg_survival$Antimicrobial <- toupper(barplot_df_lg_survival$Antimicrobial)
barplot_df_lg_survival <- barplot_df_lg_survival[complete.cases(barplot_df_lg_survival), ]

# Define lists and dataframes
ants <- unique(barplot_df_shortened$Antimicrobial)
orgs <- unique(barplot_df_shortened$Organism)
output_df_processed_oddsratio <- c()
output_df_processed_estimate <- c()

# Loop through organisms
for (i in orgs) {
  # Prepare data for odds ratio
  holder_dimension <- dim(barplot_df_shortened[barplot_df_shortened$Organism == i, ])[1]
  ants_tmp <- ants[which(!ants %in% barplot_df_shortened$Antimicrobial[barplot_df_shortened$Organism == i])]
  holder_df <- list(
    OddsRatio = rep(0, (10 - holder_dimension)),
    upper = rep(0, (10 - holder_dimension)),
    lower = rep(0, (10 - holder_dimension)),
    pvalue = rep(0, (10 - holder_dimension)),
    Antimicrobial = sample(ants_tmp, (10 - holder_dimension)),
    Organism = rep(i, (10 - holder_dimension))
  )
  output_df_processed_oddsratio <- rbind(output_df_processed_oddsratio, barplot_df_shortened[barplot_df_shortened$Organism == i, ], data.frame(holder_df))
  
  # Prepare data for estimate
  tmp_df <- barplot_df_shortened[barplot_df_shortened$Organism == i, ]
  estimate_tmp <- barplot_df_lg[barplot_df_lg$Organism == i, ]
  estimate_tmp <- estimate_tmp[match(tmp_df$Antimicrobial, estimate_tmp$Antimicrobial), ]
  estimate_tmp$Antimicrobial <- tmp_df$Antimicrobial
  estimate_tmp$Organism <- i
  holder_df_estimate <- list(
    Estimate = rep(0, (10 - dim(estimate_tmp)[1])),
    upper = rep(0, (10 - dim(estimate_tmp)[1])),
    lower = rep(0, (10 - dim(estimate_tmp)[1])),
    Antimicrobial = holder_df$Antimicrobial,
    Organism = rep(i, (10 - dim(estimate_tmp)[1]))
  )
  estimate_tmp <- rbind(estimate_tmp, holder_df_estimate)
  estimate_tmp[is.na(estimate_tmp)] <- 0
  output_df_processed_estimate <- rbind(output_df_processed_estimate, estimate_tmp)
}

# Define a color palette
color_palette <- c(
  "#8DD3C7", "#80B1D3", "#FFFFB3", "#BEBADA", "#FB8072", "#FDB462", "#D95F02", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#1B9E77", "#FFED6F", "#BC80BD", "#CCEBC5", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"
)

# Define terms for antimicrobials
terms_antimicrobials <- c(
  "AMIKACIN", "AMOXICILLIN/CLAVUL.", "AMPICILLIN", "CEFEPIME", "CEFTAZIDIME",
  "CEFUROXIME", "CIPROFLOXACIN", "CLINDAMYCIN", "COTRIMOXAZOLE", "ERTAPENEM",
  "ERYTHROMYCIN", "FUSIDIC ACID", "GENTAMICIN", "IMIPENEM", "MEROPENEM",
  "NITROFURANTOIN", "NORFLOXACIN", "PIPERACILLIN/TAZOBACTAM", "TETRACYCLINE"
)

# Create shortened organism names
shortened_orgs <- c(
  "E. COLI", "ESBL POSITIVE E. COLI", "S. AUREUS", "K. PNEUMONIAE", "ENT. EXCEPT E. COLI & K. PNEUMONIAE",
  "P. AERUGINOSA", "MRSA", "GBS", "ESBL POSITIVE K. PNEUMONIAE", "E. FAECALIS"
)

output_df_processed_oddsratio$Organism <- shortened_orgs[match(output_df_processed_oddsratio$Organism, unique(output_df_processed_oddsratio$Organism))]

# Create bar plots for odds ratio
ggplot(data = output_df_processed_oddsratio, aes(x = Organism, y = OddsRatio, fill = Antimicrobial)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = upper, ymax = lower), position = position_dodge(width = 0.9), width = 0.2, color = "black") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title.x = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    strip.text.x = element_text(size = 18, color = "black", face = "bold"),
    legend.position = "bottom",  # Set the legend position to the bottom
    plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, unit = "cm")  # Increase the top margin to provide more space
  ) +
  geom_vline(xintercept = seq(0.5, length(unique(output_df_processed$Organism)) - 0.5), color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")

# Create bar plots for estimate
ggplot(data = output_df_processed_estimate, aes(x = Organism, y = Estimate, fill = Antimicrobial)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.9), width = 0.2, color = "black") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, hjust = 1),
    axis.title.x = element_text(color = "black", size = 18, face = "bold"),
    axis.title.y = element_text(color = "black", size = 18, face = "bold"),
    strip.text.x = element_text(size = 18, color = "black", face = "bold"),
    legend.position = "bottom",  # Set the legend position to the bottom
    plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, unit = "cm")  # Increase the top margin to provide more space
  ) +
  geom_vline(xintercept = seq(0.5, length(unique(output_df_processed$Organism)) - 0.5), color = "black", size = 0.5)


#####Survival and logistic regression analysis (Figure 6)####

# Load required libraries
library(epitools)
library(tidyverse)
library(readxl)
library(survey)
library(survival)
library(car)
library(tidycmprsk)

# Partial header for GeneralData, AMR data (MicroCultureData) and the Terms files may be found in the File folder
github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/"
gen_info_file <- file.path(data_dir, "GeneralData.csv")
terms_file <- file.path(data_dir, "Terms_total.csv")
organism_files <- list.files(path = github_url, pattern = "MicroCultureData_\\d{4}\\.xlsx", full.names = TRUE)
organisms_name_file <- file.path(data_dir, "Organisms.csv")

# Load data
gen_info <- read_csv( paste0(github_url,gen_info_file), col_names = c("PATIENT_ID", "DEATH_INDICATOR", "PROBLEM_LIST",
                                                  "ADMISSION_DIAGNOSIS", "CLINICAL_DX", "BMI",
                                                  "BMI_DATE", "ADMISSION_DATE", "DISCHARGE_DATE",
                                                  "TREATMENT_OUTCOME"))
gen_info$PATIENT_ID <- as.character(gen_info$PATIENT_ID)

terms <- read_csv(paste0(github_url,terms_file))

# Function to read and preprocess AMR data
read_and_preprocess_amr_data <- function(file_path) {
  amr_data <- read_xlsx(file_path) %>%
    select(PATIENT_ID, SPECIMEN_DATE_FINAL, Drug, ORGANISM_NAME, AST_Result_CAT, 
           PATIENT_GENDER, PATIENT_NATIONALITY, PATIENT_AGE)
  amr_data$ORGANISM_NAME <- organisms_name$UpdatedOrganismName[match(amr_data$ORGANISM_NAME, organisms_name$Oragnism)]
  amr_data$ORGANISM_NAME <- gsub("[[:punct:]]", "", amr_data$ORGANISM_NAME)
  amr_data <- amr_data[amr_data$ORGANISM_NAME != "NOI", ]
  return(amr_data)
}

# Read and preprocess AMR data for multiple years
amr_data_list <- lapply(organism_files, read_and_preprocess_amr_data)
amr_tot <- do.call(rbind, amr_data_list)

# Analysis for top organisms
org_df <- names(rev(tail(sort(table(amr_tot$ORGANISM_NAME)), 10)))
barplot_df <- c()
barplot_df_lg <- c()
barplot_df_survival <- c()

for (j in 1:length(org_df)) {
  amr_tot_filetered_organism <- amr_tot %>%
    filter(grepl(paste0("^", org_df[j]), ORGANISM_NAME))
  
  names_ab <- names(tail(sort(table(amr_tot_filetered_organism$Drug)), 10))
  print(j)
  
  for (k in 1:length(names_ab)) {
    print(k)
    amr_tot_filetered <- amr_tot %>%
      filter(grepl(paste0("^", org_df[j]), ORGANISM_NAME) & grepl(names_ab[k], Drug)) %>%
      mutate(PHENOTYPE = ifelse(grepl("R", AST_Result_CAT), "R", "S"), 
             WEEK = year_week(SPECIMEN_DATE_FINAL, "2017-01-01")) %>%
      arrange(WEEK)
    
    amr_tot_filetered$COMBINED <- paste0(amr_tot_filetered$PATIENT_ID, "_", amr_tot_filetered$WEEK)
    
    gen_info$PATIENT_ID <- as.character(gen_info$PATIENT_ID)
    amr_tot_filetered <- inner_join(amr_tot_filetered, gen_info, by = c("PATIENT_ID" = "PATIENT_ID"))
    
    if (length(unique(amr_tot_filetered$PHENOTYPE)) > 1) {
      table_freq <- table(amr_tot_filetered$PHENOTYPE, amr_tot_filetered$DEATH_INDICATOR)[c(2, 1),]
      if (length(which(table_freq == 0)) == 0) {
        # Perform logistic regression
        predictor_columns <- colnames(amr_tot_filetered)[-1]
        predictor_columns_short <- predictor_columns[c(8, 5, 6, 7, 15, 20, 21, 22:length(predictor_columns))]
        
        formula_str <- paste("DEATH_INDICATOR ~", paste(predictor_columns_short, collapse = " + "))
        model_formula <- as.formula(formula_str)
        
        model <- glm(model_formula, family = binomial(link = 'logit'), data = amr_tot_filetered)
        
        # Print summary and other statistics
        print(summary(model)$coefficients[c(1, 2, 3),])
        
        # Perform survival analysis
        formula_str <- paste("Surv(PATIENT_AGE, DEATH_INDICATOR) ~", 
                             paste(predictor_columns_short[-4], collapse = " + "))
        model_formula <- as.formula(formula_str)
        
        tb <- ""
        tryCatch(
          {
            tb <- crr(model_formula, data = amr_tot_filetered)
            print(tb$tidy)
            barplot_df_survival <- rbind(barplot_df_survival, c(tb$tidy$estimate[1], 
                                                                tb$tidy$estimate[1] - tb$tidy$std.error[1], 
                                                                tb$tidy$estimate[1] + tb$tidy$std.error[1], 
                                                                names_ab[k], org_df[j]))
          },
          error = function(e) {
            # Handle any potential errors during survival analysis
            print("Error in survival analysis")
          }
        )
      }
    }
  }
}

# The rest of your code for analysis goes here, if there are any additional steps or calculations you need to perform.

# Finally, save the results to a CSV file
barplot_df_survival <- data.frame(barplot_df_survival)
barplot_df_survival$Estimate <- as.numeric(as.character(barplot_df_survival$Estimate))
barplot_df_survival$upper <- as.numeric(as.character(barplot_df_survival$upper))
barplot_df_survival$lower <- as.numeric(as.character(barplot_df_survival$lower))
colnames(barplot_df_survival) <- c("Estimate", "upper", "lower", "Antimicrobial", "Organism")
result_file <- file.path(data_dir, "SurvivalAnalysis_DeathOddsRatio_Strains_21Aug2023.csv")
write_csv(barplot_df_survival, result_file)


####Cross Covariance Analysis (Figure 5)####
library(dplyr)
library(dtw)
library(ggplot2)

# Filter AMR data
amr_tot_filtered_organism <- amr_tot %>%
  filter(grepl(paste0("^", org_df[3]), ORGANISM_NAME))

# Get the top 10 drug names based on frequency
names_ab_amr <- names(tail(sort(table(amr_tot_filtered_organism$Drug)), 10))

# Initialize a matrix for correlation results
correlation_res_drug <- matrix(0, nrow = length(names_ab), ncol = length(names_ab_amr))

# Loop through drug names and weeks
values <- data.frame()
for (j in 1:length(names_ab_amr)) {
  for (k in 1:15) {
    trend_total_red <- trend_total[trend_total$antimicrobial == names_ab[k], ]
    amr_tot_filtered <- amr_tot %>%
      filter(grepl(paste0("^", org_df[3]), ORGANISM_NAME) & grepl(names_ab_amr[j], Drug)) %>%
      mutate(PHENOTYPE = ifelse(grepl("R", AST_Result_CAT), "R", "S"), 
             WEEK = year_week(SPECIMEN_DATE_FINAL, "2017-01-01")) %>%
      arrange(WEEK)
    
    amr_tot_filtered <- amr_tot_filtered[which(amr_tot_filtered$WEEK %in% trend_total_red$weeks), ]
    
    amr_tot_filtered_sum <- table(amr_tot_filtered$WEEK, amr_tot_filtered$PHENOTYPE)
    tmp <- apply(amr_tot_filtered_sum, 1, sum)
    amr_tot_filtered_weekly_sum <- as.data.frame.matrix(table(amr_tot_filtered$WEEK, amr_tot_filtered$PHENOTYPE) / tmp)
    amr_tot_filtered_weekly_sum$total_tests <- tmp
    
    ccf_values <- ccf(amr_tot_filtered_weekly_sum$R, trend_total_red$trend, type = "correlation")
    ccf_values <- ccf_values$acf[, , 1][20:40]
    
    if (length(which(ccf_values > 0.14 | ccf_values < -0.14)) > 1) {
      values <- rbind(values, cbind(ccf_values, rep(names_ab_amr[j], length(ccf_values)), rep(names_ab[k], length(ccf_values))))
    }
  }
}

values <- data.frame(values)
values$ccf_values <- as.numeric(as.character(values$ccf_values))

# Create a boxplot
ggplot(values, aes(x = V2, y = ccf_values, fill = V3)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0.14, linetype = "dotted", color = "blue", size = 1.5) +
  geom_hline(yintercept = -0.14, linetype = "dotted", color = "blue", size = 1.5) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 13, hjust = 1),
        axis.title.x = element_text(color = "black", size = 15, face = "bold"),
        axis.title.y = element_text(color = "black", size = 15, face = "bold"),
        strip.text.x = element_text(size = 15, color = "black", face = "bold")
  ) +
  ylab("ccf") +
  xlab("Antimicrobial")

