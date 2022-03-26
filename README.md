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

* Notice : the uncertainty parameter target only one - hydraulic conductivity in "case-2."  
           However, the uncertainty parameter target have 5 - gs_mean, gs_variance, ls_x, ls_y and ls_z in "case-3.".
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

---
## Heterogeneous demo case

Create a random field and run the transport simulation as "true" site, then
use the true site contaminated data as our observation.

Geostatistical model parameter (ref : Troldborg et al., 2010) which are the randomFieldGenerator's parameters (eg:mean, variance, len scale):

<p align="center">
  <img src="https://user-images.githubusercontent.com/63782903/160170618-3d152419-1a8e-40cd-93e4-88e98d53cc63.png" width=50%>
</p>

First, the available observation data are used to update the
covariance parameters via Bayesian inference. We use these
covariance parameters to generate random conductivity
fields (unconditioned).

The "true" site model:
>1. gs_mean : -10.981
>2. gs_var : 3.2305442
>3. gs_ls : [11.1, 40.73, 1.53]

The hydraulic conductivity field (Kxx=Kyy=Kzz) of true site in FEFLOW as following image:

<p align="center">
  <image src="https://user-images.githubusercontent.com/63782903/160131301-ffb38b51-fda7-4ce6-a7f3-cbe390b4dc62.png" width=50%/>
</p>

Then we setting the Mass concentration Dirichlet BC, then start 100 days simulation.

As following image, the left image show the concentration distribution in slice-1, and
we will take the control plane data as our concentration of observation.

Using Voronoi method to calculate the control plane mass flux is "0.00286233".

<p align="center">
  <image src="https://user-images.githubusercontent.com/63782903/160133801-5e558cd4-e2ba-4b93-a4b3-64ed8135b8dc.png" width=50%>
</p>

  ..
