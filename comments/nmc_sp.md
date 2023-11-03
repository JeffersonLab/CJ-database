# NMC proton and Deuteron Per-nucleon reduced cross section 

## Reference: 
Nucl.Phys. B483 (1997) 3-43（hep-ph/9610231）, table 3 and 4

* The cross section showed in Table 3 and 4 in ref are **NOT** corrected for higher order electroweak processes. The reduced cross section provided here are reconstructed from R and F2 in the same table.

## Data files: 
  * sigr  proton     : [xlsx](../data/dataframe/10022.xlsx), [csv](../data/dataframe/csv/10022.csv)   
  * sigr  deuteron   : [xlsx](../data/dataframe/10040.xlsx), [csv](../data/dataframe/csv/10040.csv)   

## Notes:
*R=L/T: 

  For x < 0.12, the target-averaged R in table 3 and 4 were extracted from NMC data as a function of x at an averaged Q2. -> runs in 1989
  
  For x >= 0.12, R1990 [ref: PRB 250 (1990) 193] was used -> runs from 1986-1989
  
  Error on R are NOT propagated to sig_r. Since the NMC didn't consider the R errors when getting F2.

*Muon Triggers:
  * T1: scattering angle > 10 mrad, 0.006 < xbj < 0.6, 0.5<Q2<75 GeV2
  * T2: scattering angle between 5 - 25 mrad, 0.002< xbj <0.14, 0.5<Q2<25 GeV2


## Relative normalization:
Data are simultanously fitted and shifted wrt BCDMS and SLAC: see table 2 in reference:

| data set | proton%  |  deuteron % |
| :--:     | :--:     |    :--:     |
| SLAC     | -0.4     |  0.9        |  
| BCDMS| -1.8|-0.7|
|NMC 90 GeV| -2.7|-2.7|
|NMC 120 GeV| 1.1|1.1|
|NMC 200 GeV T1| 1.1|1.1|
|NMC 280 GeV T1| -1.7|-1.7|
|NMC 200 GeV T2| -2.9|-2.9|
|NMC 280 GeV T2|2|2|




## Uncertainties:
__dE__:    systematic uncertainty (%) from incident muon energy. Correlated between p and d targets, but independent between different energies.

__dE'__:   systematic uncertainty (%) from scattered muon energy.Correlated between p and d targets, but independent between different energies.

__dAC__:   systematic uncertainty (%) from acceptance. Fully correlated for all energies and targets.

__dRC__:   systematic uncertainty (%) from radiative correction. Correlated between all energies but independent for p and d targets.

__dRE__:   systematic uncertainty from the reconstruction efficiency. Fully correlated for all energies and targets.

__dFsyst__: overall absolute systematic uncertainty, which is the above errors added in quadrature.

__dST__:   statistical uncertainty (%). dST = dFstat/F2p * 100.

__norm__:  2.5% normalization uncertainty including relative normalization.
