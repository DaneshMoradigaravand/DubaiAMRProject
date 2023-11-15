# Antimicrobial Utilization and Resistance Analysis in a Hospital Network

This repository contains code and data related to the research article titled "Unveiling the dynamics of antimicrobial utilization and resistance in a large hospital network over five years: Insights from health record data analysis" [DOI: 10.21203/rs.3.rs-3014899/v1](https://www.researchsquare.com/article/rs-3014899/v1). The study explores the impact of the COVID-19 pandemic on antimicrobial resistance (AMR) using population-level data from clinical, laboratory, and prescription records. We are adding the data resulting from the analysis to this repository (File folder). Original EHR data may be directly accessed with the permission from Dubai Health Authority. 

## Abstract
Antimicrobial Resistance (AMR) poses a significant global public health challenge, which has been further compounded by the COVID-19 pandemic. However, there is a lack of comprehensive population-level data integrating clinical, laboratory, and prescription data to understand the impact of the pandemic on AMR evolution. In this study, we present an analysis of data extracted from a centralized electronic platform that captures health records of 60,551 patients across a network of public healthcare facilities in Dubai, United Arab Emirates. Our analysis employs various analytical methods, including time-series analysis, natural language processing (NLP), and unsupervised clustering algorithms, to investigate trends in antimicrobial usage and resistance over time, assess the impact of prescription practices on resistance rates, and explore the effects of COVID-19 on antimicrobial usage and resistance.

## Methodology
Data was extracted from the centralized electronic platform, encompassing inpatient and outpatient records of patients diagnosed with bacterial infections between 01/01/2017 and 31/05/2022. The dataset includes structured and unstructured Electronic Health Record data, microbiological laboratory data (antibiogram, molecular typing, and COVID-19 testing information), as well as antibiotic prescribing data. We utilized various analytical techniques, such as time-series analysis, NLP, and unsupervised clustering algorithms, to analyze the data and derive insights into antimicrobial utilization and resistance patterns.

## Results
Our findings identified a significant impact of COVID-19 on antimicrobial prescription practices, with short-term and long-lasting over-prescription of these drugs. Resistance to antimicrobials increased the odds ratio of all mortality to an average of 2.18 (95% CI: 1.87-2.49) for the most commonly prescribed antimicrobials. Moreover, the effects of antimicrobial prescription practices on resistance were observed within one week of initiation. Significant trends in antimicrobial resistance, exhibiting fluctuations for various drugs and organisms, with an overall increasing trend in resistance levels, particularly post-COVID-19 were identified.


Please refer to the research article for detailed findings and discussion.

## File Content

The Code file contains R codes for predocing the figures for causal impact analysis and natural language processing. The repository contains the following files:

| File Name                                           | Description                           |
| --------------------------------------------------- | ------------------------------------- |
| Correlation_Resistance_Prescription_CCF_Inpatient.csv | Correlation between resistance and prescription for inpatient data |
| Correlation_Resistance_Prescription_CCF_Outpatient.csv | Correlation between resistance and prescription for outpatient data |
| Count_Antimicrobial_Prescription_Inpatient.csv       | Count of antimicrobial prescriptions in inpatient data |
| Count_Antimicrobial_Prescription_Outpatient.csv      | Count of antimicrobial prescriptions in outpatient data |
| DeathOddsRatio.csv                                  | Odds ratio of death data from resistance |
| Drugs_Abbreviation_Class.csv                        | Abbreviation and classification of drugs |
| inpatient_prescription_resistance_odds_ratio.csv     | Odds ratio of prescription on resistance in inpatient data |
| outpatient_prescription_resistance_odds_ratio.csv    | Odds ratio of prescription on resistance in outpatient data |
| IP_season_frequency.csv                         | Seasonal fluctuation in inpatients prescriptions |
| OP_season_frequency.csv                         | Seasonal fluctuation in outpatients prescriptions |
| LogisticRegression_DeathOddsRatio_Strains.csv   | Logistic regression results for the impact of resistance on death |
| SurvivalAnalysis_DeathOddsRatio_Strains.csv     | Survival analysis results                      |
| UpdatedOrganisms.csv                            | List of organisms and standardized names       |
| drugs.csv                                       | List of drugs and frequency                    |
| inpatient_prescription_resistance_odds_ratio.csv | Odds ratio for the effect of prescription on resistance for inpatients |
| outpatient_prescription_resistance_odds_ratio.csv | Odds ratio for the effect of prescription on resistance for outpatients |
| text_antibiotic_IP.csv                          | The extracted text for the inpatient prescriptions |
| text_antibiotic_OP.csv                          | The extracted text for the outpatient prescriptions |
| trend_total_comp_tot_IP.csv                     | The trend of inpatient prescriptions            |
| trend_total_comp_tot_OP.csv                     | The trend of outpatient prescriptions           |


## Citation
If you use the code or findings from this study, please cite the following research article:
```



| Type of Cost Items                       | Cost per Unit (SAR) | Total (SAR)   |
|------------------------------------------|---------------------|---------------|
| 1. **Requested Funds**                   |                     | 1,600,000     |
|   - University Overhead (10%)            |                     | 160,000       |
|   - PI’s Salary                          | 60,000 (per year)  | 240,000       |
|   - Co-PI’s Salary                       | 50,000 (per year)  | 200,000       |
| 2. **Human Resources**                   |                     |               |
|   - Co-Researcher Consultant             | 24,000 (per year)     | 96,000        |
| 3. **Sample Collection, Sequencing, and Data Storage** |             |               |
|   - Short-read Sequencing (Illumina)     | 212,600 (per sample) X3000 | 637,800       |
|   - ONT Sequencing                       | 375 (per sample) X100   | 37,500        |
|   - Single-particle Analysis (Cryo-EM)   | 200 (per session)  | 20,000        |
| 4. **Lab Consumables**                   |                     |               |
|   - Microbiology Lab                      | 24,046.75 (per year)     | 96,187        |
|   - Structural Biology (3rd and 4th year)| 15,000 (per year)            | 60,000        |
|   - Grids for Data Acquisition (Aim 2)   |               | 7,500         |
| 5. **Publication Costs**                 | 11,253.12 (per publication) | 45,013 |



This paper is still under review
```

For any questions or inquiries, please contact the authors.




| Type of Cost Items                       | Cost per Unit/per Year (SAR) | Total (SAR)   |
|------------------------------------------|---------------------|---------------|
| 1. **PIs and Overhead**                   |                     |  |
|   - University Overhead (10%)            |                     | 160,000       |
|   - PI’s Salary                          | 60,000 (per year)  | 240,000       |
|   - Co-PI’s Salary                       | 50,000 (per year)  | 200,000       |
| 2. **Human Resources**                   |                     |               |
|   - Co-Researcher Consultant             | 24,000 (per year)     | 96,000        |
| 3. **Sample Collection, Sequencing, and Data Storage** |             |               |
|   - Short-read Sequencing (Illumina)     | 212,600 (per sample) X3000 | 637,800       |
|   - ONT Sequencing                       | 375 (per sample) X100   | 37,500        |
|   - Single-particle Analysis (Cryo-EM)   | 200 (per session)  | 20,000        |
| 4. **Lab Consumables**                   |                     |               |
|   - Microbiology Lab                      | 24,046.75 (per year)     | 96,187        |
|   - Structural Biology (3rd and 4th year)| 15,000 (per year)            | 60,000        |
|   - Grids for Data Acquisition (Aim 2)   |               | 7,500         |
| 5. **Publication Costs**                 | 11,253.12 (per publication) | 45,013 |



This paper is still under review
