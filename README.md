## Modeling the complete dose-time drug sensitivity surface

**Reference:** [Recursive model for dose-time responses in pharmacological studies](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2831-4)

Majority of the drug sensitivity prediction models attempt to predict a single representative metric of the complete dose-response curve such as IC<sub>50</sub> or AUC at steady-state. This can potentially fail to provide some crucial information such as the trend of change in sensitivity as dose increases or the difference in sensitivity trends between two dose-response curves with similar AUC and/or IC<sub>50</sub> values. These information can often be used for a particular patient receiving precision therapy to select the most effective drug dosage or to avoid potential drug toxicity at the current time point. 

![Dose-time-sensitivity](https://github.com/dhruba018/Dose_time_Response_Recursive_Model/blob/master/3D_dose_time_resp_curve_example.png)

**Description.**  
We develop a recursive methodology to model the complete dose-time drug response profile post drug administration that follows [Gompertz law](https://en.wikipedia.org/wiki/Gompertz%E2%80%93Makeham_law_of_mortality) in time and [sigmoidal model](https://en.wikipedia.org/wiki/Sigmoid_function) in dose, and can be used to predict the sensitivity at a certain dose and/or time point. The details of the model is described in the 2019 paper: [Recursive model for dose-time responses in pharmacological studies](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2831-4). 

This repository contains the necessary code to reproduce the results described in the paper and the corresponding data for the synthetic experiment.  

    RecursiveHybridModel.m: Model class with all necessary methods  
    synthetic_data_analysis_SRD_v2.m: Synthetic experiment code
    synthetic_recursive.mat: Data file for synthetic experiment
