# An-R-package-for-of-High-Dimensional-Error-Prone-Data-Analysis
This repository consists of the leukemia dataset and the analysis result by using the SIMEXBoost package in R CRAN

Two excel files are contained:
1. leukemia_big.csv: The full leukemia dataset.
2. SIMEXBoost_analysis.csv: Results of sensitivity analyses with selected covariates and corresponding estimates. In this data analysis, we examine different magnitude of measurement error effects and explore linear and quadratic extrapolation functions when implementing the simulation and extrapolation (SIMEX) method. Here "Correct-L-0.1", "Correct-L-0.3", and "Correct-L-0.5" represents the SIMEXBoost method with the linear extrapolation function and measurement error variance 0.1, 0.3, and 0.5, respectively; "Correct-Q-0.1", "Correct-Q-0.3", and "Correct-Q-0.5" represents the SIMEXBoost method with the quadratic extrapolation function and measurement error variance 0.1, 0.3, and 0.5, respectively.
3. SIMEXBoost_data_analysis.r: Programming code and the demonstration of data analysis
