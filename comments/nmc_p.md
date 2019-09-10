# NMC Proton and Deuteron F2

## Reference: 
  Nucl.Phys. B483 (1997) 3-43（hep-ph/9610231）, table 3 and 4

## Data files: 
  * F2  proton     : [xlsx](../data/JAM/10020.xlsx), [csv](../data/JAM/csv/10020.csv)   
  * F2  deuteron   : [xlsx](../data/JAM/10039.xlsx), [csv](../data/JAM/csv/10039.csv)   

## Note on R: 

  For x < 0.12, the target-averaged R in table 3 and 4 were extracted from NMC data as a function of x at an averaged Q2. -> runs in 1989
  
  For x >= 0.12, R1990 [ref: PRB 250 (1990) 193] was used -> runs from 1986-1989
  
  Error on R are NOT propagated to sig_r. Since the NMC didn't consider the R errors when getting F2 ??.


## Uncertainties:
__dE__:    systematic uncertainty (%) from incident muon energy. Correlated between p and d targets, but independent between different energies.

__dE'__:   systematic uncertainty (%) from scattered muon energy.Correlated between p and d targets, but independent between different energies.

__dAC__:   systematic uncertainty (%) from acceptance. Fully correlated for all energies and targets.

__dRC__:   systematic uncertainty (%) from radiative correction. Correlated between all energies but independent for p and d targets.

__dRE__:   systematic uncertainty from the reconstruction efficiency. Fully correlated for all energies and targets.

__dST__:   statistical uncertainty (%). dST = dFstat/F2p * 100

__norm__:  2.5% normalization uncertainty including relative normalization
