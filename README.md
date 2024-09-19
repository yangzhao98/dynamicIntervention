# dynamicIntervention

R codes for implementing data analyses in the manuscript entitled "Dynamic lipid-lowering interventions reduce cardiovascular disease and all-cause mortality risk in the community-dwelling population: a target trial emulation using the Chinese Multi-provincial Cohort Study" (In principal acception). 

1. _dynamicIntervention()_ defines the dynamic lipid-lowering interventions, including LDL-C, non-LDL-C, and lipid-lowering medications.
2. _Utilities()_ generates the Tables and Figures provided in the manuscript.
3. _lmtp_shiftDat_toGitHub()_ validates the parametric g-formula-based estimators using the subsequent doubly robust estimators, in which all info related to the Chinese Multi-provincial Cohort Study was removed. Please replace the info with your own data to run this function successfully. 
