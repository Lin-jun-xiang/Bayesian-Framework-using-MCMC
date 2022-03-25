# MCMC
Uncertainty Analysis by Markov Chain Monte Carlo

---
## Introduction

The contaminated site model exist some uncertainty (ex: model, boundary condition, contaminated parameters and **hydraulic conductivity** uncertainty...), therefore, 
reduce the uncertainty and show the range of uncertainty will be our target.

Uncertainty target:
>1. Hydraulic conductivity will be our uncertainty target.

>2. In the demo case, we suppose our model is homogeneous conductivity field which only have 2 values of conductivuty.

>3. In the Mingshuan site case, we suppose our model is heterogeneous conductivity field (so called random field), and the uncertainty target will be randomFieldGenerator's parameters (eg:mean, variance, len scale)

---
## Markov chain Monte Carlo

Using observation data (ex: contamination, conductivity...) to updated our uncertainty parameters (Baye's theorem).

Prior vs Posterior:
>Prior : the realizations before MCMC which will have more uncertainty.

>Posterior : the realizations after MCMC which will have lower uncertainty.

Markov chain explained:

As following image, after M-H algorithm iteration, we will get posterior (stable) parameters (𝜽).

<p align="center">
<img src="https://user-images.githubusercontent.com/63782903/160072144-1a8d6d96-5fb3-418a-b9f4-966232f4bff5.png" alt="Cover" width="50%"/>
</p>
