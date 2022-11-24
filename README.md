# survSens

This package performs a dual-parameter sensitivity analysis of treatment effect to unmeasured confounding in observational studies with either survival or competing risks outcomes.

## Release Note

## Install package from GitHub

```r
if(!require(devtools))install.packages("devtools")
devtools::install_github("Rong0707/survSens")
```

## Usage examples

For survival outcomes, 
```r
#load the dataset included in the package.
data(survdata)
#stochastic EM with regression
tau.res = survSensitivity(survdata$t, survdata$d, survdata$Z, survdata$X, "stoEM_reg", B = 5)
#EM with regression
tau.res = survSensitivity(survdata$t, survdata$d, survdata$Z, survdata$X, "EM_reg", Bem = 50)

plotsens(tau.res, coeff0 = 1.131)
```

For competing risks outcomes,
```r
#load the dataset included in the package
data(comprdata)
#stochastic EM with regression
tau.res = comprSensitivity(comprdata$t, comprdata$d, comprdata$Z, comprdata$X, "stoEM_reg", B = 5)
#EM with regression
tau.res = comprSensitivity(comprdata$t, comprdata$d, comprdata$Z, comprdata$X, "EM_reg", Bem = 50)

plotsens(tau.res$tau1, coeff0 = 1.244)
```

Output consists of dataframe(s) for estimated treatment effect(s), as well as a contour plot for visualization.
