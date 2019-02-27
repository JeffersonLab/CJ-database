__BCDMS proton and deuteron F2 structure functions and reduced cross section__

* Data files: 
  * xlsx: [proton](../format/10016.xlsx), [deuteron](../format/10017.xlsx)
  * csv : [proton](../csv/10016.csv), [deuteron](../csv/10017.csv)

* Source: CJ database

* References:
  * Proton F2 and R: [Paper](https://inspirehep.net/record/276661?ln=en)
  and [HEPData](https://hepdata.net/record/ins285519) 
  * Deuteron F2 and R: [Paper](https://inspirehep.net/record/285497?ln=en)
  and [HEPData](https://hepdata.net/record/ins285497)
  * Apart from correlated errors, data correspond to the "combined structre functions" presented in the "HepData on-line data review":
[F2(proton)](http://hepdata.cedar.ac.uk/review-cgi/struct3/BCDMS/PL223B485/f2protcomb) 
and [F2(deuteron)](http://hepdata.cedar.ac.uk/review-cgi/struct3/BCDMS/PL237B592/f2nucldeutcomb)

* Reduced cross section sig_r = F2(R=0)

* Correlated systematic errors:
  * [AA Oct 2016] BUT: not sure where the correlated columns were taken from; it seems not from the above references. See also comment in the [old HepData page](http://hepdata.cedar.ac.uk/view/ins285497).
  * However, the sum in quadrature of the 5 sources of correlated errors corresponds to the syst error found in the combined table from teh previous point.

* Normalization error 0.030
   * [AA Oct 2016] __Already included in the correlated vectors???? Seems not, therefore I added a calculated `*norm_c` column__

* __NOTE:__ this is F2 from a combination of cross sections at 4 different energies, assuming R=Q_QCD, with the letter explained in the published paper.
Also available are data directly measured for R, and F2 calculated at R=0.
