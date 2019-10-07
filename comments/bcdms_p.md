# BCDMS proton and deuteron F2 structure functions and reduced cross section

## Data files: 
  * F2       proton   : [xlsx](../data/JAM/10016.xlsx), [csv](../data/JAM/csv/10016.csv)   
  * F2       deuteron : [xlsx](../data/JAM/10017.xlsx), [csv](../data/JAM/csv/10017.csv)   
  * sig_r    proton   : [xlsx](../data/JAM/10018.xlsx), [csv](../data/JAM/csv/10018.csv)   
  * sig_r    deuteron : [xlsx](../data/JAM/10019.xlsx), [csv](../data/JAM/csv/10019.csv)  

## Sources:
  * data converted from CJ database
  * This dataset is a collection of data from different beam energies(100,120,200,280 GeV). Apart from correlated errors, the data (at each energy, and also combined energy) are also presented in the [HepData on-line data review"](http://hepdata.cedar.ac.uk/review/f2/BCDMS.shtml). Details of correlated systematic uncertainties can be found in corresponding CERN reports.
  * F2 files contain the experiment extracted F2 assuming R=R_QCD.
  * The experiment also provided F2 calculated with R=0. This is used as the reduced cross section:  sig_r = F2(R=0)

## References:
  * Proton F2 and R: Phys.Lett. B223 (1989) 485-489， [Paper](https://inspirehep.net/record/276661?ln=en)
  and [CERN report](http://cds.cern.ch/record/185732/files/cer-000097167.pdf) 
  * Deuteron F2 and R:Phys.Lett. B237 (1990) 592-598， [Paper](https://inspirehep.net/record/285497?ln=en)
  and [CERN report](http://cds.cern.ch/record/203765/files/199001439.pdf)
  

## Uncertainties:
* Correlated systematic errors:
  The published paper gave statistical and total systematic uncertainties. 5 sources of correlated factors are obtained from the CERN report: 
  * __cor1_c, cor2_c, cor3_c:__ uncertainties due to beam momentum calibration, spectrometer magnetic field calibration, and spectrometer resolution, respectively;
  * __cor4_c:__  systematic error due to detector and trigger inefficiencies;
  * __cor5_c:__  uncertainty in the relative normalization of data from external and internal targets.

* Normalization error: 0.030

