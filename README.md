# ralid
Jamovi module providing CCC, ICC, Bland-Altman, TOST, error metrics,   within-session repeatability, and related plots.

## Statistics

**Lin's Concordance Correlation Coefficient:** A single index of agreement that combines correlation with a bias correction, providing an alternative to ICC when you want a single concordance measure. Lin's CCC is moment-based, not model-based, using means, variances, and covariances of paired measurements, and does not rely on explicit variance components underpinning ANOVA-based ICC models. [Paper](https://github.com/ebonfowl/ralid/blob/main/Papers/Lin_1989_Concordance_Correlation_Coefficient.pdf)

**Liao's Concordance Correlation Coefficient:** An enhanced concordance measure that incorporates the joint geometry of paired observations, yielding a more accurate assessment of agreement when methods differ in mean, variance, or scaling. Unlike Lin’s moment-based CCC, Liao’s CCC explicitly adjusts for differences in precision and slope, reducing the upward bias in agreement estimates when variability is unequal. [Paper](https://github.com/ebonfowl/ralid/blob/main/Papers/Liao_2003_An_Improved_CCC.pdf)

**Intraclass Correlation Coefficients:** A flexible assessment of either agreement or consistency between scores. The ralid package includes one-way random effects, two-way random effects, and two-way mixed effects models that assess both single and average of K raters. Validity and reliability assessment are both provided via ICC in this module. [Paper](https://github.com/ebonfowl/ralid/blob/main/Papers/Koo_2015_Selecting_and_Reporting_ICC.pdf)

**Bland-Altman Limits of Agreement:** An intuitive assessment of agreement between measures which is very useful for visualizing systematic bias, proportionate bias, and clinical significance. The ralid module includes an option to plot regression lines to evaluate proportionate bias visually, and outputs the slope of this line and tests the hypothesis that β<sub>1</sub> ≠ 0. It also includes options to plot clinical thresholds defined by the user and the TOST equivalence interval discussed below. [Paper](https://github.com/ebonfowl/ralid/blob/main/Papers/Altman_1983_Limits_of_Agreement.pdf)

**TOST Equivalence Assessment:** An inferential procedure that tests whether the mean bias between two measures falls within user-specified equivalence bounds (i.e., practically negligible differences). It reverses the usual null hypothesis by considering differences outside the bounds as the null, and rejects non-equivalence only if both one-sided tests (lower and upper) are significant, which is equivalent to a confidence interval (90% CI when α = .05) lying entirely within the equivalence region. [Paper](https://github.com/ebonfowl/ralid/blob/main/Papers/Schuirmann_1987_Equivalence_Testing.pdf)

**Error Metrics:** Several absolute and relative error metrics are provided by the ralid module.
- Mean bias
- Standard deviation of difference scores
- **MAE (Mean Absolute Error):** The average absolute difference between two measures, representing typical unsigned error in the original units.
- **RMSE (Root Mean Square Error):** The square root of the average squared differences between measures, weighting larger errors more heavily.
- **NRMSE (Normalized RMSE):** RMSE expressed relative to a scale factor (e.g., the mean or range) to allow comparison across variables with different units or magnitudes.
- **MAPE (Mean Absolute Percentage Error):** The average absolute error expressed as a percentage of the reference value, indicating relative error magnitude.
- **sMAPE (Symmetric MAPE):** A percentage error metric that uses the average magnitude of the two measures in the denominator to reduce asymmetry present in standard MAPE.
- **TE (Typical Error):** The within-subject standard error of measurement calculated as the standard deviation of paired differences divided by √2.
- **TE% (Percent Typical Error):** Typical Error expressed as a percentage of the reference measure’s mean, providing a unitless index of relative measurement error.

## Plots

### Bland-Altman Plots

**With Regression Lines**

<img width="720" height="540" alt="image" src="https://github.com/user-attachments/assets/4142e708-df89-4bb1-a835-43cc58ed1f8d" />

**With Clinical Thresholds**

<img width="720" height="540" alt="image" src="https://github.com/user-attachments/assets/9b182629-4340-4ebe-9908-82cab3330b18" />

### Mountain Plots

<img width="720" height="540" alt="image" src="https://github.com/user-attachments/assets/d86143b7-2795-46f4-af92-cd33a6cb66bf" />

### TOST Equivalance Region Plots

<img width="720" height="450" alt="image" src="https://github.com/user-attachments/assets/89af8c03-2158-4846-be08-f0d44f230393" />
