# Structure function grids

We provide pre-calculated structure function values at given x, Q^2 grids in the [LHAPDF](https://lhapdf.hepforge.org) format. 

## Usage:
* 1. install LHAPDF library from the [official website](https://lhapdf.hepforge.org) 
* 2. Set up $LHAPDF_PATH, and put downloaded grids under the path.
* 3. load LHAPDF module in python to read grids, and calculate structure functions at desired x, Q^2.  See ```src/plot_sfn.py```


## SFn Index (see also [txgrids](https://jeffersonlab.github.io/txgrids/_build/html/grids.html)):

| sfn_index | observable  |
|-----------|-------------|
| 900       | F2 g        |
| 901       | FL g        |
| 902       | F2 gZ       |
| 903       | FZ gZ       |
| 904       | F3 gZ       |
| 905       | F2 Z        |
| 906       | FL Z        |
| 907       | F3 Z        |
| 908       | F2 NC       |
| 909       | FL NC       |
| 910       | F3 NC       |
| 930       | F2 W-       |
| 931       | FL W-       |
| 932       | F3 W-       |
| 940       | F2 W+       |
| 941       | FL W+       |
| 942       | F3 W+       |



## Available grids:

__CJ15nlo_mod_SFn__: within the [CJ15 framework](https://inspirehep.net/literature/1420566), structure functions are calculated with modified uncertainties (eigendirections increased from 19 to 24).  
* Reference: appendix B in [arXiv:2309.16851](https://arxiv.org/abs/2309.16851)
* 49 sets in total. Set 0000 is the central value, 1-48 are error sets.
* Only F2 are available with higher twist. Others are all leading twist (LT)

