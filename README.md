# Bayesian Data Analysis

Working through Gelman's [*Bayesian Data Analysis*](http://www.stat.columbia.edu/~gelman/book/) and practicing in Julia.

A lot of the derivations checks did not need complex code so only giving a short snippet here.

## `MathsChecks.ipynb`

Contains few checks of integrals using [SageMath](https://wiki.sagemath.org/)

## `chapter_5_hierarchical_models`

Contains some coded explaration of the conjugate problem from Chap 5.1. Plotting in Julia was taking too long, so all the data processing was done with Julia (`sec_5_3_cmd.jl`), whilst plotting was coded up in a Python notebook (`Orchestrate_5_3.ipynb`). INterface was via the CMD calls
