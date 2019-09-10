# Jlab JLCEE96 Experiment Proton and Deuteron(per-nucleon) Reduced Cross Section.

* Data files: 
  * sig_r    proton   : [xlsx](../data/JAM/10055.xlsx), [csv](../data/JAM/csv/10055.csv)   
  * sig_r    deuteron : [xlsx](../data/JAM/10056.xlsx), [csv](../data/JAM/csv/10056.csv)   

## Reference: 

Maria Ioana Niculescu's thesis: https://misportal.jlab.org/ul/publications/view_pub.cfm?pub_id=5397

## Source file: 

Deuteron: https://hallcweb.jlab.org/resdata/database/jlabd2.txt (dueteron cross section, not per-nucleon)

Proton: https://hallcweb.jlab.org/resdata/database/jlabh2.txt


## Uncertainty:
The source file listed only the total uncertainty. 

*The total systematic uncertainty is 3.1% from table 3.5, ref p107, among which 
	*norm_c: the beam energy, target density, and charge uncertainties are pulled out to get the normalization uncertainty = sqrt(0.3^2+0.7^2+1^2) = 1.26%. 
	*syst_u: the remained is denoted as uncorrelated systematic uncertainty = sqrt(3.1^2 - 1.26^2)=2.83%

*stat: statistical = sqrt(total_uncertainty^2-3.1%^2). The results match the statistical uncertainty listed in A.1/A.2 of reference 


