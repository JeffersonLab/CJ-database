# CJ Unpolarized DIS Database Homepage

__Reference__: [arXiv:2309.16851](https://arxiv.org/abs/2309.16851). 

See also
* CTEQ-JLab collaboration [website](https://www.jlab.org/theory/cj/).
* [note](src/cj-notes.pdf) for reduced cross section and F2 calculation.


## World DIS data tables
World __proton__ and __deuteron__ data of unploarized DIS cross sections, F2 structure functions, and the longitudinal to transverse cross section ratio R are collected or extracted from various experiments. Data were collected for the CJ global fit and related analysis. Now open for general use. See details under the [data](./data) directory.


## Neutron F2 extraction
Based on the collected F2 data, we performed a data-driven extraction of __neutron F2__ and __neutron-to-proton F2n/F2p ratio__ within the CJ15 framework (see eq. 7-9 in reference for details). Data from all experiemnts are cross-normalized and combined into a single Excel file, both in the original kineamtics, as well as rebinned in Q^2. Check the [f2n](./f2n) directory.

## Structure function grids
Within CJ framework, we calculated various structure functions (F2, F3, FL, etc) at given x, Q^2 grids. Results are provided under folder [SFN_grids](./SFN_grids) in the [LHAPDF](https://lhapdf.hepforge.org) format. An example plotting script is available at ```src/plot_sfn.py```


** This database is maintained by Alberto Accardi (accardi at jlab org) and Shujie Li (shujieli at lbl gov). Please help us improve (e.g. request a bug fix or adding new data) by sumitting a github issue.
