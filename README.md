**Nonparametric Bayesian Instrumental Variable Analysis: Evaluating Heterogeneous Effects of Arterial Access Sites for Opening Blocked Blood Vessels**

**Samrachana Adhikari, Sherri Rose, Sharon-Lise Normand**

Percutaneous coronary interventions (PCIs) are nonsurgical procedures to open blocked blood vessels to the heart, frequently using a catheter to place a stent. The catheter can be inserted into the blood vessels using an artery in the groin or an artery in the wrist. Because clinical trials have indicated that access via the wrist may result in fewer post procedure complications, shortening the length of stay, and ultimately cost less than groin access, adoption of access via the wrist has been encouraged. However, patients treated in usual care are likely to differ from those participating in clinical trials, and there is reason to believe that the effectiveness of wrist access may differ between males and females. Moreover, the choice of artery access strategy is likely to be influenced by patient or physician unmeasured factors. To study the effectiveness of the two artery access site strategies on hospitalization charges, we use data from a state-mandated clinical registry including 7,963 patients undergoing PCI. A hierarchical Bayesian likelihood-based instrumental variable analysis under a latent index modeling framework is introduced to jointly model outcomes and treatment status. Our approach accounts for unobserved heterogeneity via a latent factor structure, and permits nonparametric error distributions with Dirichlet process mixture models. Our results demonstrate that artery access in the wrist reduces hospitalization charges compared to access in the groin, with higher mean reduction for male patients.

https://arxiv.org/abs/1804.08055


**Instructions to install the R package:**

1. Directly from github repository

- install.packages('devtools')
- library(devtools)
- install_github('SamAdhikari/BayesIV_0.1')

2.  By downloading the folder from GitHub

- download (fork) the repository in a local directory
- build .tar.zip file of  the R package using "R CMD build BayesIV"  
- install in R from a local  repository


