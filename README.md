# survSens

This package performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding in observational studies with either survival or competing risks outcomes.

## Release Note

**12/27/2022 Version 1.1.0**

Motivated by Cinelli, C., & Hazlett, C. (2020), we introduce partial $R^2_{T \sim U | X, Z}$ and partial $R^2_{Z \sim U | X}$ to measure the dependency between the unmeasured confounder $U$ and the outcome $T$, and between $U$ and the treatment $Z$, respectively.

In the outcome model, suppose the full model
$$T \sim X + Z + U$$
has a log likelihood $l_{1}$, and the reduced model
$$T \sim X + Z$$
has a log likelihood $l_{0}$, with $k$ events, we calculate partial $R^2_{T \sim U | X, Z}$ similar to $R^2$ in O'Quigley, J., Xu, R., & Stare, J. (2005)
$$R^2_{T \sim U | X, Z} = 1 - e^{-2(l_1 - l_0)/k}.$$

In the treatment model, according to Cox, D. R., & Snell, E. J. (1989), the full model
$$Z \sim X + U$$
has deviance $DEV_1$,
$$DEV_1 = -2 \sum[z_i \log(\hat p_i) + (1-z_i)\log(1-\hat p_i)]$$
where
$$\hat p_i = P(z_i = 1 | x_i, u_i)$$
and the reduced model
$$Z \sim X$$
has deviance $DEV_0$,
$$DEV_0 = -2 \sum[z_i \log(\hat p_i) + (1-z_i)\log(1-\hat p_i)]$$
where
$$\hat p_i = P(z_i = 1 | x_i)$$
then
$$R^2_{Z \sim U | X} = 1 - e^{(DEV_1-DEV_0)/n}.$$

The returned partial $R^2$ has a sign, which indicates whether the association is positive or negative. Its absolute value is calculated as above.

**04/29/2020 Version 0.1.0**

Performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding in observational studies with either survival or competing risks outcomes as described in Huang, R., Xu, R., & Dulai, P. S. (2020).

## Install package from GitHub

```r
if(!require(devtools))install.packages("devtools")
devtools::install_github("Rong0707/survSens")
```

## Usage examples

For survival outcomes, 
```r
# Load the dataset included in the package.
data(survdata)
# Stochastic EM with regression
tau.res = survSensitivity(survdata$t, survdata$d, survdata$Z, survdata$X, "stoEM_reg", B = 5)
# EM with regression
tau.res = survSensitivity(survdata$t, survdata$d, survdata$Z, survdata$X, "EM_reg", Bem = 50)

# Contour plot with coefficients as axes.
plotsens(tau.res, coeff0 = 1.131)
# Contour plot with partial R-squared as axes.
plotsens(tau.res, coeff0 = 1.131, TRUE)
```

For competing risks outcomes,
```r
# Load the dataset included in the package
data(comprdata)
# Stochastic EM with regression
tau.res = comprSensitivity(comprdata$t, comprdata$d, comprdata$Z, comprdata$X, "stoEM_reg", B = 5)
# EM with regression
tau.res = comprSensitivity(comprdata$t, comprdata$d, comprdata$Z, comprdata$X, "EM_reg", Bem = 50)

# Contour plot with coefficients as axes.
plotsens(tau.res$tau1, coeff0 = 1.244)
# Contour plot with partial R-squared as axes.
plotsens(tau.res$tau1, coeff0 = 1.244, TRUE)
```

Output consists of dataframe(s) for estimated treatment effect(s), as well as a contour plot for visualization.

## Reference
* Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending omitted variable bias. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(1), 39-67.
* Huang, R., Xu, R., & Dulai, P. S. (2020). Sensitivity analysis of treatment effect to unmeasured confounding in observational studies with survival and competing risks outcomes. Statistics in Medicine, 39(24), 3397-3411.
* O'Quigley, J., Xu, R., & Stare, J. (2005). Explained randomness in proportional hazards models. Statistics in medicine, 24(3), 479-489.
* Cox, D. R., & Snell, E. J. (1989). The analysis of binary data, 2nd ed. London: Chapman and Hall.
