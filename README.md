# CPS
Simulated Data which was used in the "Cross potential selection: a proposal for optimizing crossing combinations in recurrent selection using the usefulness criterion of future inbred lines" was uploaded here. This paper has been published in G3 Genes|Genomes|Genetics (https://doi.org/10.1093/g3journal/jkae224)
You can do the same breeding simulation to my manuscript when you douwnload these files.

---

* R
    * 1.functionCode.R : registering the functions which were used in the scripts
    * 2.0.simulationSetting198.R : simulating the initial population
    * 2.1.simulationSetting150.R : simulating the population in generation 0 by four-way cross
    * 3.0.GS.R : breeding strategy using normal genomic selection
    * 3.1.OCS.R : breeding strategy using optimal cross selection
    * 3.2.CPS.R : breeding strategy using cross potential selection
    * 4.0.Comparison.R : comparing the results among 3 breeding strategies
* results
    * seedInd.csv : setting random numbers
    * snpInfoDf.rds : genome data of initial population filtered by MAF and LD
    * genomeMat.rds : genome marker data of initial population