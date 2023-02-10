# CoNet
We develop a cox proportional hazard model for network regression in TWAS, CoNet to detect the association between a given network and survival phenotype. CoNet is developed under a two-stage TWAS framework. In the first stage, the nonparametric Dirichilet process regression(DPR) model is used to estimat the SNP effect size for each specific gene within onde network in the eQTL study. In the second stage, CoNet adopts PMI to quantify the edge of the network to describe the general relationship among the nodes, and conducts the association analysis with all the nodes with edges in the model. CoNet can simultaneously identify the protential nodes and edges that are associated with survival time. \<br>
This approach is described in, \<br>
>Jiayi Han 1,2, Liye Zhang 1,2, Ran Yan1,2, Tao Ju 1,2, Xiuyuan Jin 1,2, Shukang Wang1,2, Zhongshang Yuan 1,2 and Jiadong Ji 3,*
>CoNet: Efficient network regression for survival analysis in transcriptome-wide association studies—with applications to studies of breast cancer

#Installation
It is easy to intall the development version of CoNet package using the 'devtools' package.\<br>
#install.packages("devtools")
library(devtools)
install_github("hanjiayi/CoNet")
#Usage



