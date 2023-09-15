####Text analysis (Figure 2) ####

# Load required libraries
library(tidyverse)
library(tm)
library(wordcloud)

# Define the GitHub URL of the CSV file
github_url <- "https://raw.githubusercontent.com/DaneshMoradigaravand/DubaiAMRProject/main/File/text_antibiotic.csv"

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





