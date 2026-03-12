# MSFR

author: Roberta De Vito, Alejandra Avalos-Pacheco

Fit the Multi-Study Factor Regression model via the [ECM algorithm](#1-fitting-a-msfa-model-via-the-ecm-algorithm).

## 1 Fitting the MSFR model via the ECM Algorithm

The following example illustrates how to fit the MSFR model via the ECM Algorithm, 
using a data set generating from a simulation scenario. 

## Getting the data

```{r help2, echo = TRUE, results = TRUE, tidy = TRUE}
library(MSFR)
data(Scenario1_MSFR.rda)
```

## Obtaining suitable starting values for model parameters
Then we get inizialization for model parameters, with q=3 common
factors and q_s=1 study-specific factors for two studies (S=2).

```{r, starting values, messages = FALSE}
start_value <- start_msfa(X_s, B_s, p_b, k, j_s, constraint = "block_lower2", method = "adhoc")
```

## Fitting the model via ECM
Now we can proceed for estimating the model parameters via the ECM algorithm

```{r get estimate, results = FALSE}
ECM_MSFR <- ecm_msfa(X_s, B_s, start=start_value,  nIt = 10000, trace = FALSE)
```

The estimated matrix of common loadings can be visualized
```{r get common, results = FALSE}
Phi <- ECM_MSFR$Phi

The estimated matrix of study-specific loadings can be visualized
```{r get spec, results = FALSE}
Lambda_1 <- ECM_MSFR$Lambda_s[[1]]
Lambda_2 <- ECM_MSFR$Lambda_s[[2]]

The estimated matrix of regression coefficients for covariates effect can be visualized
```{r get cov, results = FALSE}
beta <- ECM_MSFR$$beta
