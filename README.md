# CPS
Simulated Data which was used in the "Cross Potential Selection: A Proposal for Optimizing Crossing Combinations in Recurrent Selection Based on the Ability of Futurely Derived Inbred Lines" was uploaded here.
You can do the same breeding simulation to my manuscript when you douwnload these files.

---

* R
    * 1.functionCode.R : registering the functions which were used in the scripts
    * 2.0.simulationSetting198.R : simulating the initial population
    * 2.1.simulationSetting150.R : simulating the population in generation 0 by four-way cross
    * 3.0.GS.R : breeding strategy using normal genomic selection
    * 3.1.OCS.R : breeding strategy using optimal cross selection
    * 3.2.CPS.R : breeding strategy using cross potential selection
    * 4.0.ComparisonOCS.R : comparing the results among OCS
    * 4.1.Comparison.R : comparing the results among 4 breeding strategies
    * 4.2.Comparison50.R : comparing the results among 2 breeding strategies in 25-year breeding program
* results
    * seedInd.csv : setting random numbers
    * snpInfoDf.rds : genome data of initial population filtered by MAF and LD