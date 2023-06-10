# Antimicrobial Utilization and Resistance Analysis in a Hospital Network

This repository contains code and data related to the research article titled "Unveiling the dynamics of antimicrobial utilization and resistance in a large hospital network over five years: Insights from health record data analysis" [DOI: 10.21203/rs.3.rs-3014899/v1](https://www.researchsquare.com/article/rs-3014899/v1). The study explores the impact of the COVID-19 pandemic on antimicrobial resistance (AMR) using population-level data from clinical, laboratory, and prescription records. We are adding the data resulting from the analysis to this repository (File folder). Original EHR data may be directly accessed with the permission from Dubai Health Authority. 

## Abstract
Antimicrobial Resistance (AMR) poses a significant global public health challenge, which has been further compounded by the COVID-19 pandemic. However, there is a lack of comprehensive population-level data integrating clinical, laboratory, and prescription data to understand the impact of the pandemic on AMR evolution. In this study, we present an analysis of data extracted from a centralized electronic platform that captures health records of 60,551 patients across a network of public healthcare facilities in Dubai, United Arab Emirates. Our analysis employs various analytical methods, including time-series analysis, natural language processing (NLP), and unsupervised clustering algorithms, to investigate trends in antimicrobial usage and resistance over time, assess the impact of prescription practices on resistance rates, and explore the effects of COVID-19 on antimicrobial usage and resistance.

## Methodology
Data was extracted from the centralized electronic platform, encompassing inpatient and outpatient records of patients diagnosed with bacterial infections between 01/01/2017 and 31/05/2022. The dataset includes structured and unstructured Electronic Health Record data, microbiological laboratory data (antibiogram, molecular typing, and COVID-19 testing information), as well as antibiotic prescribing data. We utilized various analytical techniques, such as time-series analysis, NLP, and unsupervised clustering algorithms, to analyze the data and derive insights into antimicrobial utilization and resistance patterns.

## Results
Our analysis revealed a significant impact of the COVID-19 pandemic on antimicrobial prescription practices, with evidence of short-term and long-lasting over-prescription of these drugs. Furthermore, resistance to antimicrobials increased the odds ratio of mortality to an average of 2.5, and the effects of prescription practices on resistance were observed within one week of initiation. We also identified significant trends in antimicrobial resistance, with fluctuations observed for various drugs and organisms. Overall, there was an increasing trend in resistance levels, particularly in the post-COVID-19 period.

Please refer to the research article for detailed findings and discussion.

## File Content

The repository contains the following files:

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


## Citation
If you use the code or findings from this study, please cite the following research article:
```
This paper is still under review
```

For any questions or inquiries, please contact the authors.
