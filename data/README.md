
# World DIS data tables

__Reference__: [arXiv:2309.16851](https://arxiv.org/abs/2309.16851). 

World DIS data for proton and deuteron F2 structure functions are provided in both xlsx/csv format ( in folder [dataframe](./dataframe)) as well as a plain [text](./text) format (deprecated).
* A five digit database ID (```100xx```) is assigned to each experiment dataset as shown in the table below. 
* Detailed explanation of data source and variables are stored under the [comments](./comments) folder. Please click the database ID in the table below for quick navigation.
* The data format is explained on the [wikipage](https://github.com/JeffersonLab/CJ-database/wiki).

**Notes**
* All data are __per-nucleon__.
* __(\*)__ = The actual measurement was d/p cross section ratios. F2 d/p was extracted by the experiment under the assumption Rp=Rd.
* __(A)__ = The data is available but not yet collected
* Extraction of reduced cross section and F2 is discussed in the [note][note] 

|Experiment                 | &sigma;<sub>r</sub>        |  F2                       |   R                    |
| :--:                      | :--:                       | :--:                      | :--:                   | 
|SLAC-Whitlow               | p:  [10014][slac_sp]       | p: [10010][slac_p]        | p: [10064][slac_rp]    |
|                           | d:  [10015][slac_sd]       | d: [10011][slac_d]        | d: [10065][slac_rp]    |
|                           | d/p: [10034][slac_dp]      | d/p __(\*)__: [10034][slac_dp] |                   |
|SLAC-Whitlow(rebinned)     |                            | rebinned p: [10012][slac_p_rebin] |                |
|                           |                            | rebinned d: [10013][slac_d_rebin] |                |
|SLAC-E140                  |                            |                           | d: [10066][e140_r]     |    
|SLAC-E140x                 | p: [10037][e140x_sp]       | p: [10035][e140x_p]       | p: [10067][e140x_rp]   |
|                           | d: [10038][e140x_sd]       | d: [10036][e140x_d]       | d: [10068][e140x_rd]   |
|NMC                        | p: [10022][nmc_sp]         | p: [10020][nmc_p]         |                        |
|                           | d: [10040][nmc_sd]         | d: [10039][nmc_d]         |                        |
|                           | d/p:[10021][nmc_dp]        | d/p __(\*)__:[10021][nmc_dp] |                     |
|BCDMS                      | p: [10018][bcdms_sp]       | p: [10016][bcdms_p]       | p: [10069][bcdms_rp]   |
|                           | d: [10019][bcdms_sd]       | d: [10017][bcdms_d]       | d: [10070][bcdms_rd]   |
|JLab E06-009               | d: [10042][e06009_sd]      | d: [10041][e06009_d]      | d: [10071][e06009_d]   |
|(includes E04-001, E02-109)|                            |                           |                        |
|JLab E94-110               | p: [10044][e94110_sp]      | p: [10043][e94110_p]      | p: [10074][e94110_rp]  |
|[JLab E03-103][e03103]     | p:10047                    | p:10045                   |                        |
|                           | d:10048                    | d:10046                   |                        |
|JLab E99-118               | p:[10052][e99118_sp]       | p:[10049][e99118_p]       | p:   __(A)__           |
|                           | d:[10053][e99118_sd]       | d:[10050][e99118_d]       | p-d: __(A)__           |
|                           | d/p:[10054][e99118_sdp]    | d/p:[10051][e99118_dp]    |                        |
|JLab JLCEE96               | p: [10055][ioana_sp]       | p: [10072][ioana_p]       |                        |
|                           | d: [10056][ioana_sd]       | d: [10073][ioana_d]       |                        |
|[JLab E00-116][e00116]     | p: 10003                   | p:  10001                 |                        |
|                           | d: 10004                   | p:  10002                 |                        |
|CLAS6                      | p: [10059][clas_sp]        | p: [10057][clas_p]        |                        |
|                           | d: [10060][clas_sd]        | d: [10058][clas_d]        |                        |
|BONUS                      |                            | n: [10061][bonus_n]       |                        |
|                           |                            | n/d: [10033][bonus_nd]    |                        |
|HERA I+II                  | p: [10026 - 10032][hera]   |                           |                        |
|[HERMES][hermes]           | p: 10007                   | p: 10005                  |                        |
|                           | d: 10008                   | d: 10006                  |                        |
|                           | d/p: 10009                 |                           |                        |
|E665                       |                            | p: [10062][e665_p]        |                        |
|                           |                            | d: [10063][e665_d]        |                        |

[note]:../src/cj-notes.pdf
[slac_sp]:comments/slac_sp.md
[slac_sd]:comments/slac_sp.md
[slac_p]:comments/slac_p.md
[slac_d]:comments/slac_p.md
[slac_p_rebin]:comments/slac_rebinned.md
[slac_d_rebin]:comments/slac_rebinned.md
[slac_dp]:comments/slac_dp.md
[slac_rp]:comments/slac_rp.md
[slac_rd]:comments/slac_rp.md
[e140_r]:comments/e140_r.md
[e140x_sp]:comments/e140x_sp.md
[e140x_sd]:comments/e140x_sp.md
[e140x_p]:comments/e140x_p.md
[e140x_d]:comments/e140x_p.md
[e140x_rp]:comments/e140x_r.md
[e140x_rd]:comments/e140x_r.md
[nmc_sp]:comments/nmc_sp.md
[nmc_sd]:comments/nmc_sp.md
[nmc_p]:comments/nmc_p.md
[nmc_d]:comments/nmc_p.md
[nmc_rp]:comments/nmc_rp.md
[nmc_rd]:comments/nmc_rp.md
[nmc_dp]:comments/nmc_dp.md
[e06009]:comments/e06009_sd.md
[e06009_d]:comments/e06009_d.md
[e06009_sd]:comments/e06009_sd.md
[e03103]:comments/e03103.md
[e02109]:comments/e02109.md
[e94110_sp]:comments/e94110_sp.md
[e94110_p]:comments/e94110_p.md
[e94110_rp]:comments/e94110_rp.md
[ioana_sp]:comments/ioana_sd.md
[ioana_sd]:comments/ioana_sd.md
[ioana_p]:comments/ioana_d.md
[ioana_d]:comments/ioana_d.md
[e99118_p]:comments/e99118_p.md
[e99118_d]:comments/e99118_p.md
[e99118_dp]:comments/e99118_p.md
[e99118_sp]:comments/e99118_sp.md
[e99118_sd]:comments/e99118_sp.md
[e99118_sdp]:comments/e99118_sdp.md
[bonus_n]:comments/bonus_n.md
[slac101_d]:comments/slac101_d.md
[e00116]:comments/e00116.md
[hermes]:comments/HERMES_DIS.md
[hera]:comments/HERA2.md
[bcdms_p]:comments/bcdms_p.md
[bcdms_d]:comments/bcdms_p.md
[bcdms_sp]:comments/bcdms_p.md
[bcdms_sd]:comments/bcdms_p.md
[bcdms_rp]:comments/bcdms_r.md
[bcdms_rd]:comments/bcdms_r.md
[clas_p]:comments/clas_p.md
[clas_d]:comments/clas_p.md
[clas_sp]:comments/clas_p.md
[clas_sd]:comments/clas_d.md
[e665_p]: comments/e665_p.md
[e665_d]: comments/e665_p.md
[bonus_nd]: comments/bonus_nd.md
