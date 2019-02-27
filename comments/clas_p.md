# F2 proton and Deuteron from CLAS6

## Reference:
Phys.Rev.C73:045205,2006,
CLAS-NOTE 2003-001, 2005-013.

## Source: 
Emails from Mikhail Osipenko.

## Words on Uncertainties (from Osipenko):

F2p has even systematic uncertainties split into different sources. Unfortunately there is no "overall normalizzation" part, however because CLAS is not very precise detector the dominant uncertainty is point-to-point correlated (due to wide acceptance and complex detector structure). The overall part is perhaps 1-3%, but we never separated it.

1) eff - efficiency+acceptance, the data are corrected by means of GEANT3 Monte Carlo simulation of CLAS detector, these are rough estimate of possible mismatch in detector description and reality,

2) eep - e+e- contamination correction, i.e. mostly e- coming from pi0 decay,

3) phe - photo-electron correction, means efficiency of Cherenkov counter of CLAS used to identify electrons,

4) rdc - radiative correction (i.e. model used in calculations),

5) mom - momentum correction, small adjustments of electron momentum due to non-ideal knowledge of CLAS magnetic field,

6) rlt - for F2 only, model dependence due to assumed R=sigma_L/sigma_T ratio used in the extraction.

All these corrections were applied to the data in tables. Uncertainties, being an estimate of unknown error, are just rough guess of how much these corrections could be wrong.

All these uncertainties are point-to-point correlated, but in unknown way, what means that if eg. the model of radiative correction or R_LT was wrong by 10% all data points will be affected together, the magnitude of effect will be similar to the given uncertainty, but the specific kinematic function could be different.

These uncertainties are independent of each other, you can sum them in quadrature."
