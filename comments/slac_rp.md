# SLAC R proton and deuteron (reanalyzed)

## Reference: 
L.W.Whitlow, SLAC-Report-357, Ph.D. Thesis, Stanford University,March 1990. 

## Source:
http://www.slac.stanford.edu/exp/e140/R.H2_357

http://www.slac.stanford.edu/exp/e140/R.D2_357

88 data points each.

* E140 R data is not included in this dataset (see e140_rp)

## Uncertainties:
  As suggested by the help document below, 3 types of uncertainties are considered:
  * stat is the absolute point-to-point statistical uncertainty,
  * syst is the absolute total systematic uncertainty which is treated as uncorrelated.
  * dRRC is the absolute uncertainty due to radiative correlations which is fully correlated.



------------------HELP.DOCUMENT.-------------------------- 

E.15, E.16 from http://www.slac.stanford.edu/exp/e140/HELP.DOC_357


Files contain the values of R extracted from the global reanalysis of SLAC deep inelastic data. For most purposes, the user is suggested to utilize the average of the hydrogen and deuterium values as presented in File E.17. This exploits our conclusion that R¢p=R¢d (=R¢n) and reduces the statistical scatter of the data. Also, some parts of the systematic error in R¢p are expected to be correlated to some parts of the syste- 
matic error in R¢d. The average results, given in File E.17, display the correct propagation of all such errors. 
 
**Note that there is one misalignment between these files.** While each file contains 88 datalines, there was no hydrogen extraction at (x,Q¢2)=(.7,5) and no deuterium extraction at (x,Q¢2)=(.625,4). 
 
 
x = Bjorken scaling variable. 
Q¢2 = 4-momentum transfer squared (GeV¢2). 
R = the extracted value of R¢p or R¢d. 
stat = total ABSOLUTE random error in R, uncorrelated point to point. 
syst = total ABSOLUTE experimental systematic error, complexly correlated between neighboring points. syst=sqrt(dRSY¢2+dRSE¢2+dRN1¢2+dRN2¢2), where dRSY, dRSE, dRN1, and dRN2 are defined below, and in more detail in Chapters 4 and 5 of Reference. 
dRSY = absolute systematic uncertainty in R due to background contamination and kinematic calibration uncertainties. 
dRSE = absolute systematic uncertainty in R due to the E' dependence of the spectrometer acceptances. 
dRN1 = absolute statistical uncertainty in R due to the statistical uncertainties in the relative normalizations of the experiments. 
dRN2 = absolute systematic uncertainty in R due to the systematic uncertainties in the relative normalizations of the experiments. 
 
General remarks: 
^^^^^^^^^^^^^^^^ 
The errors dRSY, dRSE, dRN1, and dRN2 are VERY complexly correlated between data points. The correct propagation of these errors is discussed in Section 5.3.1 of Reference, and is not reconstructible from 
the information provided here.**When making fits to the data, simply assume that syst is uncorrelated between points, and let the weights for the fitting procedure be determined by 1/(stat¢2+syst¢2).** When averaging these values together, on the other hand, make the conservative assumption that each of these errors is perfectly correlated across all measurements. 
 
**NOTE: When manipulating these data, beware that the uncertainty distributions are NOT gaussian. See Appendix B.2 of Reference for the best way to deal with this fact.** 
 
The uncertainty dRRC due to radiative corrections uncertainties is not included in this table and is given by Equation 4.5 of Reference: 
 
 dRRC = .023(1+.5R) , 
 
also defined absolutely and perfectly correlated across all measurements on both targets. 
 
