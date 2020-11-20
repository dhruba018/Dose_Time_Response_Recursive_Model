## Modeling the complete dose-time drug sensitivity surface

**Reference:** [Recursive model for dose-time responses in pharmacological studies](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2831-4)

Majority of the drug sensitivity predictive modeling approaches attempt to predict a single representative metric of the complete dose-response curve such as the *[half maximal inhibitory concentration (IC50)](https://en.wikipedia.org/wiki/IC50)* or the *area under the curve (AUC)* at the steady-state of 72 hours. This can potentially fail to provide some crucial information such as the trend of change in sensitivity as dose increases or the difference in sensitivity trends between two dose-response curves with similar *AUC* and/or *IC<sub>50</sub>* values. These information can often be used for a particular patient receiving precision therapy to select the most effective drug dosage or to avoid potential drug toxicity at the current time point. 

![Dose-time-sensitivity](https://github.com/dhruba018/Dose_time_Response_Recursive_Model/blob/master/3D_dose_time_resp_curve_example.png)

### Description

#### Data
We use *in vitro* dose-time proteomic and cellular viability information post drug administration for 10 BRAF Melanoma V<sup>600E/D</sup> cell lines from the HMS-LINCS database which is part of the NIH Library of Integrated Network-based Cellular Signatures (LINCS) Program at Harvard Medical School. We were forced to limit our analysis to a single dataset since, to our knowledge, HMS-LINCS is the only publicly available source offering functional responses as well as predictors. To know more about the data, please look into our investigation in modeling drug sensitivity using proteomic features in our previous work from BHI 2018: (An investigation of proteomic data for application in precision medicine)[https://ieeexplore.ieee.org/abstract/document/8333447]. 


#### Description
We develop a recursive methodology to model the complete dose-time drug response profile post drug administration that follows [Gompertz law](https://en.wikipedia.org/wiki/Gompertz%E2%80%93Makeham_law_of_mortality) in time and the [4-parameter sigmoidal model](https://en.wikipedia.org/wiki/Sigmoid_function) in dose

![recursive_eqn](https://latex.codecogs.com/gif.latex?y_%7Bt%2C%20d%2C%20i%7D%20%3D%20%5Cunderbrace%7B%5Cleft%5B%20a_%7B0%2C%20i%7D%20&plus;%20%5Cfrac%7Bb_%7B0%2C%20i%7D%20-%20a_%7B0%2C%20i%7D%7D%7B1%20&plus;%20%5Cleft%28%20%5Cdfrac%7Bc_%7B0%2C%20i%7D%7D%7Bd%7D%20%5Cright%29%5E%7B%5Ctheta_%7B0%2C%20i%7D%7D%7D%20%5Cright%5D%7D_%7B%5Cbf%20%5Ctext%7BSigmoidal%20Model%7D%7D%20%5Cunderbrace%7Be%5E%7B%5Cgamma_%7Bd%2C%20i%7D%20%5Cleft%28%201%20%5C%2C%20-%20%5C%2C%20e%5E%7B%5Calpha_%7Bd%2C%20i%7Dt%7D%20%5Cright%29%7D%7D_%7B%5Cbf%20%5Ctext%7BGompertz%20Model%7D%7D)

, and can be used to predict the sensitivity at a certain dose and/or time point. The details of the model is described in the 2019 paper: [Recursive model for dose-time responses in pharmacological studies](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2831-4). 

This repository contains the necessary code to reproduce the results described in the paper and the corresponding data for the synthetic experiment.  

    RecursiveHybridModel.m: Model class with all necessary methods  
    synthetic_data_analysis_SRD_v2.m: Synthetic experiment code
    synthetic_recursive.mat: Data file for synthetic experiment
