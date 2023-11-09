# SLAC Proton and Deuteron F2 (reanalyzed)
## Data files: 
  * F2  proton     : [xlsx](../dataframe/10010.xlsx), [csv](../dataframe/csv/10010.csv)   
  * F2  deuteron   : [xlsx](../dataframe/10011.xlsx), [csv](../dataframe/csv/10011.csv)   

## Reference:  
L.W.Whitlow, SLAC-Report-357, Ph.D. Thesis, Stanford University, March 1990.     

## Source:
http://www.slac.stanford.edu/exp/e140/F2.H2_357         

http://www.slac.stanford.edu/exp/e140/F2.D2_357

* E0, E', theta are not available in the source file. They were borrowed from the corresponding cross section data files slac_sp / slac_sd

## Uncertainties
Source: http://www.slac.stanford.edu/exp/e140/HELP.DOC_357

**norm**: 
overall normalization uncertainty. 1.7% for Deuteron and 2.1% for Hydrogen. See p76 in ref.

**dFST**: 
fractional statistical COUNTING uncertainty in F_2. This is the same as the dST uncertainty in the measured cross sections.       
          

**dFSR**: 
fractional random systematic uncertainty in F_2. This is the same as the dSR uncertainty in the measured cross sections. [This uncertainty is due to FLUCTUATIONS in beam charge, in scattering kinematics, in target density, and in detector efficiencies -- and so, is random in nature.]       

**dFSY**: 
fractional systematic uncertainty in F_2 due to background contamination. This is the same as the dSY uncertainty in the measured cross sectionsWhen propagated within an experiment, assume perfectly correlated. When propagated between experiments, assume uncorrelated. [This uncertainty   is due to uncertainties in calibrations of scattering kinematics and uncertainties due to background contamination.]   

**dFSE**: 
fractional systematic uncertainty in F_2 due to the E' dependence of the spectrometer acceptance. This is the same as the dSE uncertainty in the measured cross sections. Because these errors are perfectly correlated with E', they carry a (+/-) sign to keep track of (+/anti) correlations. When propagated within an experiment, assume correlated and respect the sign.  When propagated between experiments, assume uncor-  related and ignore the sign.  [This uncertainty is due to the  possibility that the spectrometer acceptance is dependent upon the E' setting of the central momentum.]        

**dFRC**: 
fractional systematic uncertainty in F_2 due to the radiative corrections. This is **NOT** the same as dRC because the epsilon-correlation of dRC effects the extraction of F_2. The general formula for dFRC is given by Equation 5.39 of Reference: 

     dFRC = .023(epsilon-.85+(1+.5R)(Rfac-1)/(Rfac(1+R)),

where

     Rfac = 1+(1-epsilon)/(epsilon(1+R)).

For x<.1 or Q¢2<1GeV, this uncertainty is increased by a factor of 1.5.  When propagated within or between experiments, assume perfectly correlated.           

**dFN1**: 
fractional statistical uncertainty in F_2 due to the relative normalization uncertainties of the experiments. This is the same as the dN1 uncertainty in the measured cross sectionsWhen propagated within an experiment, assume perfectly correlated. When propagated between experiments, assume uncorrelated although they are distinctly, positively correlated as indicated in Table 5.4 of Reference. [This uncertainty reflects possible statistical errors in our global normalization procedure.]      

**dFN2**: 
fractional systematic uncertainty in F_2 due to the relative normalizations of the experiments. This is the same as the dN2 uncertainty in the measured cross sections of Files E.2 and E.3, and so, the columns are labeled with both. When propagated within an experiment, assume perfectly correlated When propagated between experiments, assume uncorrelated. [This uncertainty reflects possible systematic errors in our global normalization procedure.]      

**dFSZ**: 
fractional systematic uncertainty in F_2 due to the experimental uncertainty in R=R¢{1990}, and so, does not correlate to an uncertainty in the measured cross sections. This error does not include the uncertainty in F_2 due to [the uncertainty in R¢{1990} due to (the uncertainty in radiative corrections)], which is explicitly included in dFRC above. When propagated within or between experiments assume perfec correlated.          

**stat**: 
total fractional random error, uncorrelated point to point stat=sqrt(dFST¢2+dFSR¢2). Note that this is the same quantity defined as "stat" in the measured cross section.

**syst**: 
total fractional experimental systematic error, correlated between neighboring points.  
     
     syst=sqrt(dFSY¢2+dFSE¢2+dFRC¢2+dFN1¢2+dFN2¢2+dFSZ¢2). 

Note that this is **NOT** the same quantity defined as "syst" in the measured cross section.
               

## Relative Normalization

**scale** is the relative shift between SLAC experiments introduced when combining the dataset. See p76 in ref for details.

|EXP    |J      |N_H2  |N_D2 | N_d/p|
|:--:   |:--:|:--:  |:--: |  :--:|
|e49a     |1    |1.012   |1.001|   0.922|
|e49b     |2    |0.981   |0.981|   1    |
|e61 |3    |1.011   |1.033|   1.92 |
|e87 |4    |0.982   |0.986|   1.013|
|e89a     |5    |0.989   |0.985|   0.995|
|e89b     |6    |0.953   |0.949|   0.995|
|e139     |7    |       |1.088|         |
|e140     |8    |       |1      |      |
              
