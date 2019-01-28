Shiny Applications for Epidemiology
============

epiapps.com is a public hosted repository of useful statistical tools in the area of epidemiology.
These applications are hosted and maintained by [Fellows Statistics Inc.](http://www.fellstat.com/).
To report an issue or bug, please create a GitHub issue on [the repository](https://github.com/fellstat/epiapps).

---

Consensus Estimation
------------

This tool assists in synthesizing multiple independent estimates of a 
quantity (e.g. population size or prevalence). Stakeholders may add additional information
regarding the methodological quality of the studies and prior knowledge of the metric.


<h3 style="text-align: center;">
  <a href="https://epiapps.com/apps/combine_estimates/"> Launch Application</a> 
</h3>


### Details: 

|   |   |
|--:|---|
| ***Authors:*** | Ian E. Fellows |
| ***Github Repository:*** | <https://github.com/fellstat/combine_estimates> |

---


Multiple Source Capture Recapture
------------

Implements user interfaces for log-linear models, Bayesian model 
averaging and Bayesian Dirichlet process mixture models.


<h3 style="text-align: center;">
  <a href="https://epiapps.com/apps/rcapture/"> Analysis Application</a> | <a href="https://epiapps.com/apps/capture_power/"> Power Calculator</a>
</h3>


### Details: 

|   |   |
|--:|---|
| ***Authors:*** | Ian E. Fellows |
| ***Video Tutorial:*** | <https://www.youtube.com/watch?v=PgmyUnFlo5Y&feature=youtu.be> |
|    ***Manual:***  |  <https://fellstat.github.io/shinyrecap/> |
| ***Github Repository:*** | <https://github.com/fellstat/shinyrecap> |

---

Population Size Estimation Using Multiple Data Sources
------------

Implements a user interface for an algorithm for presenting a Bayesian hierarchical model 
for estimating the sizes of 
local and national populations. The model incorporates 
multiple commonly used data sources including mapping data, surveys, interventions, 
capture-recapture data, estimates or guesstimates from organizations, and expert opinion.

<h3 style="text-align: center;">
  <a href="https://epiapps.com/apps/size_estimation/"> Launch Application</a> 
</h3>


### Details: 

|   |   |
|--:|---|
| ***Authors:*** | Jacob Parsons using the algorithm developed by Le Bao, Adrian E. Raftery and Kyongwon Kim|
| ***CRAN Repository:*** | <https://CRAN.R-project.org/package=SizeEstimation> |
| ***Reference:*** | Bao, L., Raftery, A. E., & Reddy, A. (2015). Estimating the sizes of populations at risk of HIV infection from multiple data sources using a Bayesian hierarchical model. Statistics and its Interface, 8(2), 125-136. | 

---



Incidence Estimation In Cross-sectional Surveys Using Testing History
------------

Utilizes crosssectional survey data containing information on participants' testing history and diagnosis to estimate incidence.


<h3 style="text-align: center;">
  <a href="https://epiapps.com/apps/testing_incidence/"> Launch Application</a> 
</h3>


### Details: 

|   |   |
|--:|---|
| ***Authors:*** | Ian E. Fellows |
| ***Video Tutorial:*** | <https://www.youtube.com/watch?v=YVPcLLs9zxc&t=08s> |
|    ***Manual:***  |  <https://github.com/fellstat/TestingHistoryIncidence/wiki/Shiny-App-Documentation> |
| ***Github Repository:*** | <https://github.com/fellstat/TestingHistoryIncidence> |
| ***Example Data:*** | <https://raw.githubusercontent.com/fellstat/TestingHistoryIncidence/master/inst/shiny_ui/tstdat.csv> |

---