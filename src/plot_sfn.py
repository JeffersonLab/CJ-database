# Example python code to plot CJ15nlo_mod structure functions from LHAPDF format grids
# Shujie Li, Oct 2023 
# 
# Prerequisites: 
# 1. Install LHAPDF, follow instructions at https://lhapdf.hepforge.org
# 2. Download CJ15nlo_mod SFN, visit https://github.com/JeffersonLab/CJ-database
#

import numpy as np
np.seterr(divide='ignore', invalid='ignore') ## supress division by zero warning

import lhapdf
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize']=(9,6)
plt.rcParams['font.size'] = 22

# dictionary for SFN index
sfn_index = {
900:"F2 g",  901:"FL g",
902:"F2 gZ", 903:"FZ gZ", 904:"F3 gZ",
905:"F2 Z",  906:"FL Z",  907:"F3 Z",
908:"F2 NC", 909:"FL NC", 910:"F3 NC",
930:"F2 W-", 931:"FL W-", 932:"F3 W-",
940:"F2 W+", 941:"FL W+", 942:"F3 W+"
}
# swap key and value to get a dict for SFN name
sfn_name = {v: k for k, v in sfn_index.items()}


## calculate structure functions with the LHAPDF grid ##
#  Inputs:
#      x, Q2 can can be arrays
#      iset: SFN index. see https://jeffersonlab.github.io/txgrids/_build/html/grids.html
#      tabname: check available table name with lhapdf.availablePDFSets()
#      kerr = 1: calculate uncertainties
#           = 0: error set to zero
# Outputs: 
#       central values, uncertainties
def calc_sfn(x,Q2,tabname="CJ15nlo_mod_SFp",iset=908,kerr=0):
    err = np.zeros(len(x))
    # central set
    set0     = lhapdf.mkPDF(tabname,0)
    sfn      = np.array(set0.xfxQ2(iset,x,Q2))
    # uncertainties
    if kerr==1:
        sets = lhapdf.mkPDFs(tabname)
        nerr = int((len(sets)-1)/2) #number of error eigen-directions = (total sets - central sets)/2
        if nerr==0:
            return sfn, err
        ## CJ style error calculation with all error sets. see See Eq.(11) in [arXiv:2309.16851]
        if "CJ" in tabname.upper():
            for ii in np.arange(0,nerr):
                sfn1    = np.array(sets[2*ii+1].xfxQ2(iset,x,Q2))
                sfn2    = np.array(sets[2*ii+2].xfxQ2(iset,x,Q2))
                err    += (sfn1-sfn2)**2
            err = np.sqrt(err)/2.0
            
        else:
        ## for all other sets, take mean and sigma
            n   = len(sets)
            sfn = []
            for ii in np.arange(0,n):
                sfn.append(np.array(sets[ii].xfxQ2(iset,x,Q2)))
            sfn = np.array(sfn)
            err = sfn.std(axis=0)
            sfn = sfn.mean(axis=0)
    return sfn, err

## calculate p-n, or n/p structure functions with the LHAPDF grid
#  Inputs:
#      x, Q2 can can be arrays
#      iset: SFN index. see https://jeffersonlab.github.io/txgrids/_build/html/grids.html
#      tabp and tabn: proton and neutron set name
#      kerr = 1: calculate uncertainties;
#           = 0: error set to zero
#      p_n  = 1: return proton - neutron; 
#           = 0: return neutron/proton (default)
#  Outputs:
#      central values, uncertainties
def calc_sfn_np(x,Q2,tabp = "CJ15nlo_mod_SFp", tabn = "CJ15nlo_mod_SFn", iset=908,kerr=0,p_n=0):
    err = np.zeros(len(x))
    # central values
    if kerr==0:
        setp0      = lhapdf.mkPDF(tabp,0)
        setn0      = lhapdf.mkPDF(tabn,0)
        if p_n:
            sfn0       = np.array(setp0.xfxQ2(iset,x,Q2))-np.array(setn0.xfxQ2(iset,x,Q2))
        else:
            sfn0       = np.array(setn0.xfxQ2(iset,x,Q2))/np.array(setp0.xfxQ2(iset,x,Q2))
        
    else:
        setp=lhapdf.mkPDFs(tabp)
        setn=lhapdf.mkPDFs(tabn)
        if not len(setp) == len(setn):
            print(f"WARNING(calc_sfn_np): tabp and tabn number of sets not match: {len(setp)} vs {len(setn)}, can not calculate error, set to zero" )
        elif len(setp)>1: ## CJ15 style error calculation. See Eq.(11) in [arXiv:2309.16851] 
            nerr = int((len(setp)-1)/2)
            print("nerr=",nerr)
            if p_n:
                sfn0  = np.array(setp[0].xfxQ2(iset,x,Q2))-np.array(setn[0].xfxQ2(iset,x,Q2))
            else:
                sfn0  = np.array(setn[0].xfxQ2(iset,x,Q2))/np.array(setp[0].xfxQ2(iset,x,Q2))
            for ii in np.arange(0,nerr):
                setp1    = setp[2*ii+1]
                setp2    = setp[2*ii+2]
                setn1    = setn[2*ii+1]
                setn2    = setn[2*ii+2]
                if p_n:
                    sfn1    = np.array(setp1.xfxQ2(iset,x,Q2))-np.array(setn1.xfxQ2(iset,x,Q2))
                    sfn2    = np.array(setp2.xfxQ2(iset,x,Q2))-np.array(setn2.xfxQ2(iset,x,Q2))   
                else:
                    sfn1    = np.array(setn1.xfxQ2(iset,x,Q2))/np.array(setp1.xfxQ2(iset,x,Q2))
                    sfn2    = np.array(setn2.xfxQ2(iset,x,Q2))/np.array(setp2.xfxQ2(iset,x,Q2)) 
                err += (sfn1-sfn2)**2
        err = np.sqrt(err)/2.0
    return sfn0, err


## -------------------------------------------------
##    Example: plot F2p and F2n/p at Q2=10 GeV2
## -------------------------------------------------

if __name__ == "__main__":

    path_to_grids ="../SFN_grids/CJ15nlo_mod_SFn/"
    if len(path_to_grids)<1:
        print("!!!\nWARNING: please make sure the dataset path are included in $LHAPDF_PATH, or add the path in the code\n!!!")
    lhapdf.pathsPrepend(path_to_grids)
    # print(lhapdf.paths())
    
    x         = np.linspace(0,0.99,1000)
    qq        = np.ones(len(x))*10 # 10 GeV2
    kerr      = 1 # calculate error or not
    
    ## plot F2 p 
    fig       = plt.figure() 
    ax        = fig.add_subplot(1, 1, 1)

    iset      = sfn_name["F2 NC"] #908
    targ      = "p"
    setname   = "CJ15nlo_mod_SFp"

    sfn1,err1 = calc_sfn(x, qq,setname,iset=iset,kerr=kerr)
    plt.fill_between(x,sfn1-err1,sfn1+err1,color='b', alpha=0.5,label = setname)#setname.replace('_', '\_'))
    plt.legend()
    plt.xlabel("x")
    plt.ylabel(sfn_index[iset]+" "+targ)
    plt.show()
    
    
    ## plot F2 n/p 
    fig       = plt.figure() 
    ax        = fig.add_subplot(1, 1, 1)
    sfn2,err2 = calc_sfn_np(x, qq,"CJ15nlo_mod_SFp","CJ15nlo_mod_SFn",iset=iset,kerr=kerr,p_n=0)
    plt.fill_between(x,sfn2-err2,sfn2+err2,color='r',label="CJ15nlo_mod",alpha=0.3)
    plt.legend()
    plt.ylabel("$F_2^n/F_2^p$")
    plt.ylim(0,1)
    plt.xlabel("x")
    plt.show()


