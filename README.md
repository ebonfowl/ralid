# ralid
Jamovi module providing CCC, ICC, Bland-Altman, TOST, error metrics,   within-session repeatability, and related plots.

## Statistics

**Lin's Concordance Correlation Coefficient:** A single measure of agreement between scores without the assumptions inherent to other, more commonly used metrics, like intraclass correlation. [Paper](https://www.jstor.org/stable/2532051?seq=1)

**Intraclass Correlation Coefficients:** A flexible assessment of either agreement or consistency between scores. The ralid package includes one-way random effects, two-way random effects, and two-way mixed effects models that assess both single and average of K raters. Validity and reliability assessment are both provided via ICC in this module. [Paper](https://www.sciencedirect.com/science/article/pii/S1556370716000158)

**Bland-Altman Limits of Agreement:** An intuitive assessment of agreement between measures which is very useful for visualizing systematic bias, proportionate bias, and clinical significance. The ralid module includes an option to plot regression lines to evaluate systematic bias visually, and outputs the slope of this line and tests the hypothesis that β<sub>1</sub> ≠ 0. It also includes options to plot clinical thresholds defined by the user and the TOST equivalance interval discussed below. [Paper](https://www.jstor.org/stable/2987937?seq=1)

**Error Metrics:** Several absolute and relative error metrics are provided by the ralid module.
- Mean bias
- Standard deviation of difference scores
- Mean absolute error (MAE)
- Root mean square error (RMSE)
- Normalized RMSE
- Mean absolute percent error (MAPE)
- Symmetric MAPE (sMAPE)
- Typical error (TE)
- Typical error percent (TE%)

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
