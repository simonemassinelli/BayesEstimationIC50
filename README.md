# Bayesian-Dose-Response: Trametinib $IC_{50}$ Estimation

## About the Project

This repository contains the source code, datasets, and Bayesian statistical architecture developed for modeling the dose-response dynamics of Trametinib, a targeted cancer drug, across various cancer cell lines.
The system acts as an advanced inferential framework: rather than relying on fragile, single-point estimates, it leverages Markov Chain Monte Carlo (MCMC) simulations to extract full posterior distributions for the parameter of interest. This provides a mathematically robust evaluation of drug sensitivity that naturally incorporates uncertainty, biological variability, and experimental noise.

## Objectives & Biological Context

To model the relationship between drug concentration and biological response, pharmacology relies on the Hill equation, generalized as the four-parameter log-logistic model (LL.4). The primary objective of this study is to estimate the half-maximal inhibitory concentration ($IC_{50}$)—the exact drug concentration required to reduce cancer cell viability by 50%.
Because Trametinib is evaluated based on its toxicity to cancer cells, the model operates "in reverse" compared to standard efficacy models: as the drug concentration increases, cell survival decreases, dropping from maximum baseline viability ($E_{max}$) toward complete cell death ($E_{min}$).

## Dataset & Genomic Filtering

The project relies on high-throughput viability measurements and biological annotations gathered from the Cancer Dependency Map (DepMap) Portal, maintained by the Broad Institute.
A critical component of the data engineering workflow was genomic filtering. To ensure that the observed biological response was truly driven by the drug's primary mechanism of action and not by overlapping genetic variables, the pharmacological measurements were cross-referenced with somatic mutation data. Cell lines presenting known off-target mutations were explicitly excluded to isolate the true pharmacological signal.

## Methodology & Modeling

Because estimating all four parameters of the LL.4 equation simultaneously would drastically increase posterior dimensionality and computational complexity, the project employs a hybrid parameter-isolation strategy implemented in JAGS.

* **Frequentist Pre-computation:** Preliminary frequentist fits were performed using the `drc` package in R to estimate the upper asymptote ($E_{max}$), the lower asymptote ($E_{min}$), and the Hill coefficient ($h$).


* **Bayesian Core:** These estimates are supplied to the Bayesian model as fixed inputs. This strategic dimensionality reduction allows the MCMC sampler to concentrate all of its resolving power entirely on exploring the posterior distribution of the $IC_{50}$.



## Uncertainty Quantification & Shrinkage Engine

The system rejects a rigid "Generic Pooled" approach in favor of a multi-level Hierarchical Patient-Specific Model. It estimates a distinct $IC_{50}$ for each individual cell line while assuming these values are drawn from a shared, population-level Gamma distribution.

* **Hierarchical Shrinkage:** Well-measured cell lines collectively establish a population norm. A regularization effect—shrinkage—pulls unstable or noisy individual estimates toward that population mean without forcing artificial homogeneity across genuinely different cancer lines.


* **Robustness to Outliers:** Unlike empirical procedures that rely heavily on a small number of observations, Bayesian inference combines information across the entire dose-response profile. The estimated curve is only moderately influenced by extreme isolated outliers.



## Estimation Strategy & Diagnostics

The final layer of the project involves rigorous model comparison and MCMC validation to ensure statistical soundness:

* **Model Selection:** The hierarchical architecture drastically outperformed the baseline pooled model, producing a massive reduction in the Deviance Information Criterion (DIC) (-1884.8 vs 162.0).


* **Posterior Predictive Checks (PPCs):** Replicated datasets generated from the posterior predictive distribution confirm the hierarchical model's ability to perfectly capture individual variability and overall population structure.


* **MCMC Diagnostics:** The reliability of the generated samples was confirmed via visual traceplots (demonstrating rapid mixing and stationarity) and formal statistics, achieving a perfect Gelman-Rubin ($\hat{R}$) score of $\approx$ 1.0.



## Author

Simone Massinelli
