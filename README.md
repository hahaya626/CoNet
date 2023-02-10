# CoNet
We develop a cox proportional hazard model for network regression in TWAS, CoNet to detect the association between a given network and survival phenotype. CoNet is developed under a two-stage TWAS framework. In the first stage, the nonparametric Dirichilet process regression(DPR) model is used to estimat the SNP effect size for each specific gene within onde network in the eQTL study. In the second stage, CoNet adopts PMI to quantify the edge of the network to describe the general relationship among the nodes, and conducts the association analysis with all the nodes with edges in the model. CoNet can simultaneously identify the protential nodes and edges that are associated with survival time. 
