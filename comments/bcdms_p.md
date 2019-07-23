#BCDMS proton and deuteron F2 structure functions and reduced cross section

## Data files: 
  * F2       proton   : [xlsx](../data/JAM/10016.xlsx), [csv](../data/JAM/csv/10016.csv)   
  * F2       deuteron : [xlsx](../data/JAM/10017.xlsx), [csv](../data/JAM/csv/10017.csv)   
  * sig_r    proton   : [xlsx](../data/JAM/10018.xlsx), [csv](../data/JAM/csv/10018.csv)   
  * sig_r    deuteron : [xlsx](../data/JAM/10019.xlsx), [csv](../data/JAM/csv/10019.csv)  
## Source: CJ database
  * F2 is a combination of cross sections at 4 different energies, assuming R=Q_QCD, with the letter explained in the published paper.
<<<<<<< HEAD
  * The experiment also provided F2 calculated with R=0. This is used as the reduced cross section:  sig_r = F2(R=0)
=======
  * The experiment also provided data directly measured for R, and F2 calculated at R=0. This is used as the reduced cross section:  sig_r = F2(R=0)
>>>>>>> b032bd2299e9fb2a5af185fd993967c4752f3f58


## References:
  * Proton F2 and R: [Paper](https://inspirehep.net/record/276661?ln=en)
  and [CERN report](http://cds.cern.ch/record/185732/files/cer-000097167.pdf) 
  * Deuteron F2 and R: [Paper](https://inspirehep.net/record/285497?ln=en)
  and [CERN report](http://cds.cern.ch/record/203765/files/199001439.pdf)
  * Apart from correlated errors, data correspond to the "combined structre functions" presented in the "HepData on-line data review":
[F2(proton)](http://hepdata.cedar.ac.uk/review-cgi/struct3/BCDMS/PL223B485/f2protcomb) 
and [F2(deuteron)](http://hepdata.cedar.ac.uk/review-cgi/struct3/BCDMS/PL237B592/f2nucldeutcomb)
  * Details of correlated systematic uncertainties can be found in corresponding CERN reports.

## Uncertainties:
* Correlated systematic errors:
  The published paper gave statistical and total systematic uncertainties. 5 sources of correlated factors are obtained from the CERN report: 
  * __cor1_c, cor2_c, cor3_c:__ uncertainties due to beam momentum calibration, spectrometer magnetic field calibration, and spectrometer resolution, respectively;
  * __cor4_c:__  systematic error due to detector and trigger inefficiencies;
  * __cor5_c:__  uncertainty in the relative normalization of data from external and internal targets.

* Normalization error 0.030

