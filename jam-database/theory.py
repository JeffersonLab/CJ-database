#!/usr/bin/env python
import sys,os
import numpy as np
import pandas as pd
import time
from obslib.idis.aux import AUX
from tools.config import conf
from tools.tools import checkdir
from scipy.special import gamma

class HT:

    def __init__(self):

        #--target mass corrections
        if   'idis tmc' in conf: self.tmc = conf['idis tmc']
        elif 'idis_tmc' in conf: self.tmc = conf['idis_tmc']
        elif 'tmc'      in conf: self.tmc = conf['tmc'] #--for backward compatibility
        else:                    self.tmc = False

        self.mell = conf['mellin']
        self.params={}
        self.params['F2p']=np.array([0.0,0.5,3,0,0])
        self.params['FLp']=np.array([0.0,0.5,3,0,0])
        self.params['F3p']=np.array([0.0,0.5,3,0,0])
        self.params['W2p']=np.array([0.0,0.5,3,0,0])
        self.params['WLp']=np.array([0.0,0.5,3,0,0])
        self.params['W3p']=np.array([0.0,0.5,3,0,0])

        self.params['F2n']=np.array([0.0,0.5,3,0,0])
        self.params['FLn']=np.array([0.0,0.5,3,0,0])
        self.params['F3n']=np.array([0.0,0.5,3,0,0])
        self.params['W2n']=np.array([0.0,0.5,3,0,0])
        self.params['WLn']=np.array([0.0,0.5,3,0,0])
        self.params['W3n']=np.array([0.0,0.5,3,0,0])

        self.FLAV = ['F2p','FLp','F3p','F2n','FLn','F3n']
        self.FLAV.extend(['W2p','WLp','W3p','W2n','WLn','W3n'])
        self.PAR = ['N','a','b','c','d']
        
    def setup(self):
        pass

    def get_state(self):
        return (self.params)

    def set_state(self,state):
        self.params = state

    def get_ht(self,x,Q2,nucleon,stf):
        #--if TMCs are active, evaluate at xi instead of x
        if self.tmc!=False:
            mu2=conf['aux'].M**2/Q2
            rho=(1+4*mu2*x**2)**0.5
            xi=2*x/(1+rho)
            x = xi
        N,a,b,c,d = self.params['%s%s'%(stf,nucleon)]
        return N*x**a*(1-x)**b*(1+c*x**0.5+d*x)/Q2 

    def beta(self,a,b):
        return gamma(a)*gamma(b)/gamma(a+b)
  
    def get_moments(self,flav,N=None,b_shift=0):
        """
        if N==None: then parametrization is to be use to compute moments along mellin contour
        else the Nth moment is returned
        """
        if N==None: N=self.mellin.N
        M,a,b,c,d=self.params[flav]
        b += b_shift
        mom=self.beta(N+a,b+1)+c*self.beta(N+a+0.5,b+1)+d*self.beta(N+a+1.0,b+1)
        norm=self.beta(2+a,b+1)+c*self.beta(2+a+0.5,b+1)+d*self.beta(2+a+1.0,b+1)
  
        result = M*mom/norm
  
        #--add admixture
        if self.choice=='valence'  and flav=='dv1':  result += self.get_mix(N)
 
        return result 

    def get_ht_mell(self,Q2,nucleon,stf):
        N = self.mell.N
        M,a,b,c,d = self.params['%s%s'%(stf,nucleon)]
        mom=self.beta(N+a,b+1)+c*self.beta(N+a+0.5,b+1)+d*self.beta(N+a+1.0,b+1)

        result = M*mom
        return np.einsum('i,x->xi',result,1/Q2)

class OFFSHELL_MODEL:

    def __init__(self,pol=False):

        if 'offpdf model' in conf: self.model = conf['offshell model']
        else:                      self.model = 1

    def get_model(self,q,nucleon,nucleus):
  
        if self.model==1:

            if 'factor_of_two' in conf: two=conf['factor_of_two']
            else: two = True
            qup=np.copy(q[1])
            qdn=np.copy(q[2])
            if nucleon=='p':
                #--in helium the average is taken (extra factor of 2 for down since up is normalized to 2)
                if nucleus=='h':
                    q[1] = 0.5*(qup+2.0*qdn)
                    if two: q[2] = 0.5*(qup+2.0*qdn)/2.0
                    else:   q[2] = 0.5*(qup+2.0*qdn)
            if nucleon=='n':
                if nucleus=='d' or nucleus=='h':
                    q[1]=qdn
                    q[2]=qup
                if nucleus=='t':
                    if two: q[1] = 0.5*(qup+2.0*qdn)/2.0
                    else:   q[1] = 0.5*(qup+2.0*qdn)
                    q[2] = 0.5*(qup+2.0*qdn)

        return q




#--with Mellin precalculation
class THEORY(AUX):

    def __init__(self):

        self.mell=conf['mellin']
        self.mellnpts = self.mell.density
        self.N=self.mell.N
        self.M=self.mell.N
        self.pdf=conf['pdf']
        if 'ht4' in conf: self.ht4 = conf['ht4']


        if 'dsmf' in conf: self.dsmf=conf['dsmf']
        if 'hsmf' in conf: self.hsmf=conf['hsmf']

        if 'dsmf_type' in conf: self.dsmf_type = conf['dsmf_type']
        if 'hsmf_type' in conf: self.hsmf_type = conf['hsmf_type']
 
        self.setup()
  
        if conf['order']=='LO':  self.order=0
        if conf['order']=='NLO': self.order=1
  
        # nuclear smearing
        if   'idis nuc' in conf and conf['idis nuc']: self.nuc = True
        elif 'nuc'      in conf and conf['nuc']:      self.nuc = True #--for backward compatibility
        else:                                         self.nuc = False
 
        self.Nf     = lambda Q2: conf['alphaS'].get_Nf(Q2)
        self.alphaS = lambda Q2: conf['alphaS'].get_alphaS(Q2)

        #--einsum storage
        self.storage = {}

        #--target mass corrections
        if   'idis tmc' in conf: self.tmc = conf['idis tmc']
        elif 'idis_tmc' in conf: self.tmc = conf['idis_tmc']
        elif 'tmc'      in conf: self.tmc = conf['tmc'] #--for backward compatibility
        else:                    self.tmc = False

        #--higher twist
        if   'idis ht' in conf: self.ht_flag = conf['idis ht']
        elif 'ht'      in conf: self.ht_flag = conf['ht']
        else:                   self.ht_flag = False

        if self.nuc:
            self.setup_nucleus()
            self.grids = {}
            self.load_grids()

        #--offshell corrections
        if   'idis offpdf' in conf: self.offpdf = conf['idis offpdf']
        elif 'idis_offpdf' in conf: self.offpdf = conf['idis_offpdf']
        elif 'offpdf'      in conf: self.offpdf = conf['offpdf'] #--for backward compatibility
        else: self.offpdf = False

        #--cannot have offshell corrections without nuclear corrections
        if self.offpdf: self.nuc = True

        if self.offpdf: self.off_model = OFFSHELL_MODEL()

    ###################################
    #--nuclear smearing grid generation
    ###################################
    #--set up nuclear smearing grids
    def setup_nucleus(self):
        eps_d       = -0.00222
        eps_h       = -0.00772
        eps_t       = -0.00848
        mp          = 0.93827231
        mn          = 0.93956563
        mN          = (mp + mn)/2.0
        mD          =   mp +   mn + eps_d
        mH          = 2*mp +   mn + eps_h
        mT          =   mp + 2*mn + eps_t
        self.ymax = {}
        self.ymax['d']  = mD/mN
        self.ymax['h']  = mH/mN
        self.ymax['t']  = mT/mN
        self.ng = 100
        self.gX,self.gW=np.polynomial.legendre.leggauss(self.ng)

    def load_grids(self,process='idis'):

        if '%s tabs'%process in conf: tabs = conf['%s tabs'%process]
        else: return

        data_dir = '%s/database/%s/expdata' %(os.environ['FITPACK'],process)
        grid_dir = '%s/grids/idis'%(os.environ['FITPACK'])

        self.grids[process] = {}

        checkdir(grid_dir)
        done_grids = os.listdir(grid_dir)

        for idx in tabs:
            #--load full dataset prior to cuts
            tab = pd.read_excel('%s/%s.xlsx'%(data_dir,idx)).to_dict(orient='list')
            x   = np.array(tab['X'])
            Q2  = np.array(tab['Q2'])
            obs = tab['obs'][0]
            tar = tab['target'][0].strip()

            tars = []
            if '/' in tar:
                num,den = tar.split('/')
                tars = [num,den]
            else:
                tars = [tar]

            if tar=='p' or tar=='n': continue

            dsmf_type = self.dsmf_type
            hsmf_type = self.hsmf_type

            name = '%s_%s_ng=%s_mellnpts=%s_dsmftype=%s_hsmftype=%s_tmc=%s.npy'%(process,idx,self.ng,self.mellnpts,dsmf_type,hsmf_type,self.tmc)
            filename = '%s/%s'%(grid_dir,name)
            if name in done_grids:
                self.grids[process][idx] = np.load(filename,allow_pickle=True).item()
                print('Loading nuclear smearing grid %s'%filename)
                continue
            grid = {}
            print('generating nuclear smearing grid for %s dataset %s with mell npts = %s,dsmf type=%s, hsmf type=%s, tmc=%s'%(process,idx,self.mellnpts,dsmf_type,hsmf_type,self.tmc))

            for tar in tars:
 
                if tar in ['p', 'n']: continue
                grid[tar] = {}

                #--generate both onshell and offshell components
                for kind in ['onshell','offshell']:
                    grid[tar][kind] = {}
                    if tar in ['d']:
                        grid[tar][kind]['f22_F2']  = self.gen_smf(x,Q2,'d','f22', 'F2',kind)
                        grid[tar][kind]['fLL_FL']  = self.gen_smf(x,Q2,'d','fLL', 'FL',kind)
                        grid[tar][kind]['fL2_F2']  = self.gen_smf(x,Q2,'d','fL2', 'F2',kind)
                        grid[tar][kind]['f33_F3']  = self.gen_smf(x,Q2,'d','f33', 'F3',kind)
                        grid[tar][kind]['f22']     = self.gen_smf(x,Q2,'d','f22', None,kind)
                        grid[tar][kind]['fLL']     = self.gen_smf(x,Q2,'d','fLL', None,kind)
                        grid[tar][kind]['fL2']     = self.gen_smf(x,Q2,'d','fL2', None,kind)
                        grid[tar][kind]['f33']     = self.gen_smf(x,Q2,'d','f33', None,kind)
                    if tar in ['h']:
                        grid[tar][kind]['f22p_F2'] = self.gen_smf(x,Q2,'h','f22p','F2',kind)
                        grid[tar][kind]['f22n_F2'] = self.gen_smf(x,Q2,'h','f22n','F2',kind)
                        grid[tar][kind]['fLLp_FL'] = self.gen_smf(x,Q2,'h','fLLp','FL',kind)
                        grid[tar][kind]['fLLn_FL'] = self.gen_smf(x,Q2,'h','fLLn','FL',kind)
                        grid[tar][kind]['fL2p_F2'] = self.gen_smf(x,Q2,'h','fL2p','F2',kind)
                        grid[tar][kind]['fL2n_F2'] = self.gen_smf(x,Q2,'h','fL2n','F2',kind)
                        grid[tar][kind]['f33p_F3'] = self.gen_smf(x,Q2,'h','f33p','F3',kind)
                        grid[tar][kind]['f33n_F3'] = self.gen_smf(x,Q2,'h','f33n','F3',kind)
                        grid[tar][kind]['f22p']    = self.gen_smf(x,Q2,'h','f22p',None,kind)
                        grid[tar][kind]['fLLp']    = self.gen_smf(x,Q2,'h','fLLp',None,kind)
                        grid[tar][kind]['fL2p']    = self.gen_smf(x,Q2,'h','fL2p',None,kind)
                        grid[tar][kind]['f33p']    = self.gen_smf(x,Q2,'h','f33p',None,kind)
                        grid[tar][kind]['f22n']    = self.gen_smf(x,Q2,'h','f22n',None,kind)
                        grid[tar][kind]['fLLn']    = self.gen_smf(x,Q2,'h','fLLn',None,kind)
                        grid[tar][kind]['fL2n']    = self.gen_smf(x,Q2,'h','fL2n',None,kind)
                        grid[tar][kind]['f33n']    = self.gen_smf(x,Q2,'h','f33n',None,kind)
                    #--for tritium, take helium but switch p <--> n
                    if tar in ['t']:
                        grid[tar][kind]['f22p_F2'] = self.gen_smf(x,Q2,'h','f22n','F2',kind)
                        grid[tar][kind]['f22n_F2'] = self.gen_smf(x,Q2,'h','f22p','F2',kind)
                        grid[tar][kind]['fLLp_FL'] = self.gen_smf(x,Q2,'h','fLLn','FL',kind)
                        grid[tar][kind]['fLLn_FL'] = self.gen_smf(x,Q2,'h','fLLp','FL',kind)
                        grid[tar][kind]['fL2p_F2'] = self.gen_smf(x,Q2,'h','fL2n','F2',kind)
                        grid[tar][kind]['fL2n_F2'] = self.gen_smf(x,Q2,'h','fL2p','F2',kind)
                        grid[tar][kind]['f33p_F3'] = self.gen_smf(x,Q2,'h','f33n','F3',kind)
                        grid[tar][kind]['f33n_F3'] = self.gen_smf(x,Q2,'h','f33p','F3',kind)
                        grid[tar][kind]['f22p']    = self.gen_smf(x,Q2,'h','f22n',None,kind)
                        grid[tar][kind]['fLLp']    = self.gen_smf(x,Q2,'h','fLLn',None,kind)
                        grid[tar][kind]['fL2p']    = self.gen_smf(x,Q2,'h','fL2n',None,kind)
                        grid[tar][kind]['f33p']    = self.gen_smf(x,Q2,'h','f33n',None,kind)
                        grid[tar][kind]['f22n']    = self.gen_smf(x,Q2,'h','f22p',None,kind)
                        grid[tar][kind]['fLLn']    = self.gen_smf(x,Q2,'h','fLLp',None,kind)
                        grid[tar][kind]['fL2n']    = self.gen_smf(x,Q2,'h','fL2p',None,kind)
                        grid[tar][kind]['f33n']    = self.gen_smf(x,Q2,'h','f33p',None,kind)

            np.save(filename,grid)
            print('Saving grid for %s dataset %s with ng=%s to %s'%(process,idx,self.ng,filename))

            self.grids[process][idx] = grid

    def load_grid_custom(self,x,Q2,tar):

        if tar=='p' or tar=='n': return

        #--no process
        process = None
        if process not in self.grids: self.grids[process] = {}

        idx = 'tar=%s,x=%s,Q2=%s'%(tar,x,Q2)
        if idx in self.grids[process]: return

        grid = {}
        print('generating custom grid for tar=%s...'%(tar))

        grid[tar] = {}

        #--generate both onshell and offshell components
        for kind in ['onshell','offshell']:
            grid[tar][kind] = {}
            if tar in ['d']:
                grid[tar][kind]['f22_F2']  = self.gen_smf(x,Q2,'d','f22', 'F2',kind)
                grid[tar][kind]['fLL_FL']  = self.gen_smf(x,Q2,'d','fLL', 'FL',kind)
                grid[tar][kind]['fL2_F2']  = self.gen_smf(x,Q2,'d','fL2', 'F2',kind)
                grid[tar][kind]['f33_F3']  = self.gen_smf(x,Q2,'d','f33', 'F3',kind)
                grid[tar][kind]['f22']     = self.gen_smf(x,Q2,'d','f22', None,kind)
                grid[tar][kind]['fLL']     = self.gen_smf(x,Q2,'d','fLL', None,kind)
                grid[tar][kind]['fL2']     = self.gen_smf(x,Q2,'d','fL2', None,kind)
                grid[tar][kind]['f33']     = self.gen_smf(x,Q2,'d','f33', None,kind)
            if tar in ['h']:
                grid[tar][kind]['f22p_F2'] = self.gen_smf(x,Q2,'h','f22p','F2',kind)
                grid[tar][kind]['f22n_F2'] = self.gen_smf(x,Q2,'h','f22n','F2',kind)
                grid[tar][kind]['fLLp_FL'] = self.gen_smf(x,Q2,'h','fLLp','FL',kind)
                grid[tar][kind]['fLLn_FL'] = self.gen_smf(x,Q2,'h','fLLn','FL',kind)
                grid[tar][kind]['fL2p_F2'] = self.gen_smf(x,Q2,'h','fL2p','F2',kind)
                grid[tar][kind]['fL2n_F2'] = self.gen_smf(x,Q2,'h','fL2n','F2',kind)
                grid[tar][kind]['f33p_F3'] = self.gen_smf(x,Q2,'h','f33p','F3',kind)
                grid[tar][kind]['f33n_F3'] = self.gen_smf(x,Q2,'h','f33n','F3',kind)
                grid[tar][kind]['f22p']     = self.gen_smf(x,Q2,'h','f22p',None,kind)
                grid[tar][kind]['fLLp']     = self.gen_smf(x,Q2,'h','fLLp',None,kind)
                grid[tar][kind]['fL2p']     = self.gen_smf(x,Q2,'h','fL2p',None,kind)
                grid[tar][kind]['f33p']     = self.gen_smf(x,Q2,'h','f33p',None,kind)
                grid[tar][kind]['f22n']     = self.gen_smf(x,Q2,'h','f22n',None,kind)
                grid[tar][kind]['fLLn']     = self.gen_smf(x,Q2,'h','fLLn',None,kind)
                grid[tar][kind]['fL2n']     = self.gen_smf(x,Q2,'h','fL2n',None,kind)
                grid[tar][kind]['f33n']     = self.gen_smf(x,Q2,'h','f33n',None,kind)
            #--for tritium, take helium but switch p <--> n
            if tar in ['t']:
                grid[tar][kind]['f22p_F2'] = self.gen_smf(x,Q2,'h','f22n','F2',kind)
                grid[tar][kind]['f22n_F2'] = self.gen_smf(x,Q2,'h','f22p','F2',kind)
                grid[tar][kind]['fLLp_FL'] = self.gen_smf(x,Q2,'h','fLLn','FL',kind)
                grid[tar][kind]['fLLn_FL'] = self.gen_smf(x,Q2,'h','fLLp','FL',kind)
                grid[tar][kind]['fL2p_F2'] = self.gen_smf(x,Q2,'h','fL2n','F2',kind)
                grid[tar][kind]['fL2n_F2'] = self.gen_smf(x,Q2,'h','fL2p','F2',kind)
                grid[tar][kind]['f33p_F3'] = self.gen_smf(x,Q2,'h','f33n','F3',kind)
                grid[tar][kind]['f33n_F3'] = self.gen_smf(x,Q2,'h','f33p','F3',kind)
                grid[tar][kind]['f22p']     = self.gen_smf(x,Q2,'h','f22n',None,kind)
                grid[tar][kind]['fLLp']     = self.gen_smf(x,Q2,'h','fLLn',None,kind)
                grid[tar][kind]['fL2p']     = self.gen_smf(x,Q2,'h','fL2n',None,kind)
                grid[tar][kind]['f33p']     = self.gen_smf(x,Q2,'h','f33n',None,kind)
                grid[tar][kind]['f22n']     = self.gen_smf(x,Q2,'h','f22p',None,kind)
                grid[tar][kind]['fLLn']     = self.gen_smf(x,Q2,'h','fLLp',None,kind)
                grid[tar][kind]['fL2n']     = self.gen_smf(x,Q2,'h','fL2p',None,kind)
                grid[tar][kind]['f33n']     = self.gen_smf(x,Q2,'h','f33p',None,kind)


        self.grids[process][idx] = grid

    def gen_smf(self,X,Q2,tar,fXX,stf,kind):

        #--with TMCs, must evaluate at xi instead of x
        if self.tmc!=False: XI = self.get_xi(X,Q2)
        else:               XI = X

        W = self.mell.W * self.mell.JAC
        nx = X.shape[0]
        #--Gaussian quadrature setup
        XM,  gXM = np.meshgrid(X,  self.gX)
        XIM, gXM = np.meshgrid(XI, self.gX)
        Q2M, gWM = np.meshgrid(Q2, self.gW)
        b = self.ymax[tar]
      
        a = XM
        YM   = 0.5*(b-a)*gXM+0.5*(a+b)
        XM_YM = XM/YM

        a = XIM
        YIM   = 0.5*(b-a)*gXM+0.5*(a+b)
        JIM   = 0.5*(b-a) 
        XIM_YIM = XIM/YIM

        if   tar in ['d']:     SMF=self.dsmf
        elif tar in ['h','t']: SMF=self.hsmf

        #--get TMC, evaluated at xi/y (note that it takes x as the argument and then converts to xi)
        if '22' in fXX or 'LL' in fXX or '33' in fXX:
            proj = fXX[1:3]
            TMC = self.get_TMC(XM_YM,Q2,proj)
            f = np.einsum('gxi,gx->gxi',TMC,SMF.get_fXX2(fXX  ,kind,XIM,Q2M,YIM))
      
        elif 'L2' in fXX:
            TMC_L2 = self.get_TMC(XM_YM,Q2,'L2')
            TMC_22 = self.get_TMC(XM_YM,Q2,'22')
            if tar in ['d']:
                f  = np.einsum('gxi,gx->gxi',TMC_L2,SMF.get_fXX2('fLL',kind,XIM,Q2M,YIM))
                f += np.einsum('gxi,gx->gxi',TMC_22,SMF.get_fXX2('fL2',kind,XIM,Q2M,YIM))
            elif tar in ['h','t']:
                f  = np.einsum('gxi,gx->gxi',TMC_L2,SMF.get_fXX2('fLL%s'%fXX[-1],kind,XIM,Q2M,YIM))
                f += np.einsum('gxi,gx->gxi',TMC_22,SMF.get_fXX2('fL2%s'%fXX[-1],kind,XIM,Q2M,YIM))
 
        XIN  = np.repeat(XIM_YIM,self.N.shape[0],axis=1).reshape(XIM_YIM.shape[0],XIM_YIM.shape[1],self.N.shape[0])

        if stf==None:          func = XIN**(-self.N)
        elif '2' or'L' in stf: func = XIN**(-self.N+1)
        elif '3' in stf:       func = XIN**(-self.N)


        fN = np.einsum('i,gx,gx,gxi,gxi->xi',W,gWM,JIM,f,func) 
        return fN

    ###########################

    # theory calculations
    def get_FX(self,idx,stf,x,Q2,tar,process='idis',channel='all'):

        if 'W' in stf: k = stf[-2]
        else:          k = stf[-1]

        if self.tmc!=False: xi = self.get_xi(x,Q2)
        else:               xi = x 

        # Nucleon structure functions
        if tar in ['p','n']: 
            #--get FXN in mellin space, without TMCs or HTs
            FXN = self.get_FXN_mell(stf,Q2,tar,channel=channel)

            #--multiply by target mass corrections
            FXN = np.einsum('xi,xi->xi',FXN,self.get_TMC(x,Q2,'%s%s'%(k,k)))

            #--extra term for FL/WL
            if 'L' in stf and self.tmc!=False:
                FXN += np.einsum('xi,xi->xi',self.get_FXN_mell('F2',Q2,tar),self.get_TMC(x,Q2,'%s2'%(k)))

            #--add higher twists (does not apply for W)
            #if self.ht_flag and 'W' not in stf:
            #    FXN += self.ht4.get_ht_mell(Q2,tar,stf)
           
            if '2' or 'L' in stf: a = 1
            if '3'        in stf: a = 0

            FXN = self.invert(xi,FXN,a)

            #--add higher twists (does not apply for W), without any extra factors of x
            #--evaluate higher twists at xi
            if self.ht_flag and 'W' not in stf:
                FXN += self.invert(xi,self.ht4.get_ht_mell(Q2,tar,stf),0)

            return FXN


        # Nuclear structure functions
        if tar=='d': p, n = 1,1
        if tar=='h': p, n = 2,1
        if tar=='t': p, n = 1,2

        #--no nuclear smearing for W
        if self.nuc==False or 'W' in stf:
            FXp = self.get_FX(idx,stf,x,Q2,'p') 
            FXn = self.get_FX(idx,stf,x,Q2,'n')
            return (p * FXp + n * FXn)/(p+n)

        return self.get_FXA(idx,stf,x,Q2,tar,process=process,channel=channel)

    #--returns mellin space structure function
    #--note that this returns F2/x, FL/x, F3
    #--note that factor of xi**(-N) is NOT included
    def get_FXN_mell(self,stf,Q2,tar,kind='onshell',nucleus=None,channel='all'):

        C = self.C
 
        alphaS = np.array([self.alphaS(q2) for q2 in Q2])
  
        # Functions to setup PDFs and FFs for einsum
        PDF, e2gl = self.setup_pdfs(stf,Q2,tar,kind=kind,nucleus=nucleus,channel=channel)

        ind_mell_LO  = 'kxi->xi'
        ind_mell_NLO = 'x,kxi,i->xi'
        ind_mell_NLO = 'x,kxi,i->xi'
        
        if 'W' in stf: k = stf[-2]
        else:          k = stf[-1]

        if ind_mell_LO not in self.storage:
            self.storage[ind_mell_LO]  = np.einsum_path(ind_mell_LO,        PDF,           optimize='optimal')[0]
            self.storage[ind_mell_NLO] = np.einsum_path(ind_mell_NLO,alphaS,PDF, C[k]['Q'],optimize='optimal')[0] 

        #--LO term
        if '2' or '3' in stf:
            FX    = np.einsum(ind_mell_LO, PDF, optimize = self.storage[ind_mell_LO])

        if 'L' in stf:
            FX = np.zeros((Q2.shape[0],self.N.shape[0]),dtype=complex) 

        #--NLO term
        if self.order >= 1:
            FX  += np.einsum(ind_mell_NLO,alphaS,PDF, C[k]['Q'],optimize = self.storage[ind_mell_NLO])/4/np.pi\
                 + np.einsum(ind_mell_NLO,alphaS,e2gl,C[k]['G'],optimize = self.storage[ind_mell_NLO])/4/np.pi
 
        if 'W' in stf: factor = 2.0
        if 'F' in stf: factor = 1.0

        FX = factor*FX 
           
        return FX
 
    def get_FXA(self,idx,stf,x,Q2,tar,process='idis',channel='all'):

        C = self.C

        if tar=='d': p, n = 1,1
        if tar=='h': p, n = 2,1
        if tar=='t': p, n = 1,2

        #--if idx is None, generate grid on the fly
        if idx==None: 
            self.load_grid_custom(x,Q2,tar)
            process=None
            idx = 'tar=%s,x=%s,Q2=%s'%(tar,x,Q2)
            grid = self.grids[process][idx]
            i = np.array([k for k in range(len(x))])

        #--load grid and place cuts
        else:
            grid = self.grids[process][idx]
            i = conf['%s tabs'%process][idx]['idx']

        if self.offpdf: kinds = ['onshell','offshell']
        else:           kinds = ['onshell']
        
        if 'W' in stf: k = stf[-2]
        else:          k = stf[-1]

        FXp = np.zeros(Q2.shape[0])
        FXn = np.zeros(Q2.shape[0])

        for kind in kinds:
            if tar=='d':         f = lambda XX, n, stf: grid[tar][kind]['f%s_%s'%(XX,stf)]    [i] 
            if tar in ['h','t']: f = lambda XX, n, stf: grid[tar][kind]['f%s%s_%s'%(XX,n,stf)][i] 

            alphaS = np.array([self.alphaS(q2) for q2 in Q2])
  
            # Functions to setup PDFs and FFs for einsum
            PDF, e2gl = {},{}
            PDF['p'], e2gl['p'] = self.setup_pdfs(stf,Q2,'p',kind=kind,nucleus=tar,channel=channel)
            PDF['n'], e2gl['n'] = self.setup_pdfs(stf,Q2,'n',kind=kind,nucleus=tar,channel=channel)


            ind_LO  = 'xi,kxi->x'
            ind_NLO = 'xi,x,kxi,i->x'
            
            if ind_LO not in self.storage:
                self.storage[ind_LO]  = np.einsum_path(ind_LO, f('22','p','F2'),       PDF['p'],            optimize='optimal')[0]
                self.storage[ind_NLO] = np.einsum_path(ind_NLO,f('22','p','F2'),alphaS,PDF['p'], C[k]['Q'], optimize='optimal')[0] 

            t1 = time.time()

            k = stf[-1]
            #--LO terms
            #--F22,FLL,F33 term
            _FXp    = np.einsum(ind_LO,f('%s%s'%(k,k),'p','F%s'%(k)),PDF['p'],optimize = self.storage[ind_LO])
            _FXn    = np.einsum(ind_LO,f('%s%s'%(k,k),'n','F%s'%(k)),PDF['n'],optimize = self.storage[ind_LO])

            #--FL2 term
            if 'L' in stf:
                _FXp   = np.zeros(Q2.shape[0],dtype=complex) 
                _FXn   = np.zeros(Q2.shape[0],dtype=complex) 


            #--NLO terms
            if self.order >= 1:
                #--F22,FLL,F33 term
                _FXp   += np.einsum(ind_NLO,f('%s%s'%(k,k),'p','F%s'%(k)),alphaS,PDF ['p'],C[k]['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
                        + np.einsum(ind_NLO,f('%s%s'%(k,k),'p','F%s'%(k)),alphaS,e2gl['p'],C[k]['G'], optimize = self.storage[ind_NLO])/4/np.pi


                _FXn   += np.einsum(ind_NLO,f('%s%s'%(k,k),'n','F%s'%(k)),alphaS,PDF ['n'],C[k]['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
                        + np.einsum(ind_NLO,f('%s%s'%(k,k),'n','F%s'%(k)),alphaS,e2gl['n'],C[k]['G'], optimize = self.storage[ind_NLO])/4/np.pi

                #--FL2 terms
                if 'L' in stf:
                    _FXp   += np.einsum(ind_NLO,f('%s2'%(k),'p','F2'),    alphaS,PDF ['p'],C['2']['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
                            + np.einsum(ind_NLO,f('%s2'%(k),'p','F2'),    alphaS,e2gl['p'],C['2']['G'], optimize = self.storage[ind_NLO])/4/np.pi
  

                    _FXn   += np.einsum(ind_NLO,f('%s2'%(k),'n','F2'),    alphaS,PDF ['n'],C['2']['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
                            + np.einsum(ind_NLO,f('%s2'%(k),'n','F2'),    alphaS,e2gl['n'],C['2']['G'], optimize = self.storage[ind_NLO])/4/np.pi


              
            phase = self.mell.phase

            FXp += np.imag(phase*_FXp)/np.pi
            FXn += np.imag(phase*_FXn)/np.pi

        #--add ht
        if self.ht_flag and 'W' not in stf:
            if tar=='d':         f = lambda XX, n: grid[tar]['onshell']['f%s'%(XX)]    [i]
            if tar in ['h','t']: f = lambda XX, n: grid[tar]['onshell']['f%s%s'%(XX,n)][i]

            ht = lambda stf,tar: self.ht4.get_ht_mell(Q2,tar,stf)
            _FXp  = np.einsum('xi,xi->x',f('%s%s'%(k,k),'p'),ht(stf,'p'))
            _FXn  = np.einsum('xi,xi->x',f('%s%s'%(k,k),'n'),ht(stf,'n'))

            if 'L' in stf:
                _FXp  += np.einsum('xi,xi->x',f('%s2'%(k),'p'),ht('F2','p'))
                _FXn  += np.einsum('xi,xi->x',f('%s2'%(k),'n'),ht('F2','n'))

            #--note that extra factor of xi is removed
            FXp += np.imag(phase*_FXp)/np.pi
            FXn += np.imag(phase*_FXn)/np.pi

        return (p*FXp + n*FXn)/(p+n)

    def setup_pdfs(self,stf,Q2,tar,kind='onshell',nucleus=None,channel='all'):
  
        if stf in ['F2', 'FL']:    aX = self.get_aX('aXp',   Q2, channel) 
        if stf in ['F3']:          aX = self.get_aX('aXm',   Q2, channel) 
        if stf in ['W2+','WL+']:   aX = self.get_aX('aXpWp', Q2, channel) 
        if stf in ['W3+']:         aX = self.get_aX('aXmWp', Q2, channel) 
        if stf in ['W2-','WL-']:   aX = self.get_aX('aXpWm', Q2, channel) 
        if stf in ['W3-']:         aX = self.get_aX('aXmWm', Q2, channel) 
       
        if kind=='onshell':  pdf = self.pdf
        if kind=='offshell': pdf = conf['off pdf']
 
        L = Q2.shape[0] 
        for q2 in Q2: pdf.evolve(q2)
 
        PDF  = np.zeros((11,L,self.N.size),dtype=complex)
        e2gl = np.zeros((11,L,self.N.size),dtype=complex)
      
        if kind=='onshell':
            #--isospin relations for free nuclei
            q = {}
            if tar == 'p': flav = {-5:'bb',-4:'cb',-3:'sb',-2:'db',-1:'ub',0:'g',1:'u',2:'d',3:'s',4:'c',5:'b'}
            if tar == 'n': flav = {-5:'bb',-4:'cb',-3:'sb',-2:'ub',-1:'db',0:'g',1:'d',2:'u',3:'s',4:'c',5:'b'}
            for i in flav:
                q[i] = np.array([pdf.storage[q2][flav[i]] for q2 in Q2])
        elif kind=='offshell':
            q = {}
            flav = {-5:'bb',-4:'cb',-3:'sb',-2:'db',-1:'ub',0:'g',1:'u',2:'d',3:'s',4:'c',5:'b'}
            for i in flav:
                q[i] = np.array([pdf.storage[q2][flav[i]] for q2 in Q2])

            q = self.off_model.get_model(q,tar,nucleus)
         
        #--take Nf into account
        Nf_mask = {}
        Nf_mask[0] = np.ones(L) 
        Nf_mask[1] = np.ones(L) 
        Nf_mask[2] = np.ones(L) 
        Nf_mask[3] = np.ones(L) 
        Nf_mask[4] = np.ones(L)*(Q2 > conf['aux'].mc2)
        Nf_mask[5] = np.ones(L)*(Q2 > conf['aux'].mb2) 

        #--no bottom quark for W
        if 'W' in stf: Nf_mask[5] = np.zeros(L)

        for i in [1,2,3,4,5]: Nf_mask[-i] = Nf_mask[i]

 
        for i in range(-5,6):
            if i==0: continue
            PDF[i]  = np.einsum('x,x,xi->xi',Nf_mask[i],aX[i],q[i])
            e2gl[i] = np.einsum('x,x,xi->xi',Nf_mask[i],aX[i],q[0])
     
        return PDF, e2gl  

    #--TMCs in Mellin space
    def get_xi(self,x,Q2):
        mu2=conf['aux'].M**2/Q2
        rho=(1+4*mu2*x**2)**0.5
        xi=2*x/(1+rho)
        return xi

    def get_TMC(self,x,Q2,proj):

        N = self.N
        shape = tuple(x.shape)+tuple(N.shape)

        #--with no TMCS, return zero mass WW
        if self.tmc==False:
            if 'L2' in proj: return np.zeros(shape,dtype=complex)
            else:            return np.ones (shape,dtype=complex)

        xi = self.get_xi(x,Q2)
        mu2=conf['aux'].M**2/Q2
        rho=(1+4*mu2*x**2)**0.5

        XI = np.repeat(xi,self.N.shape[0],axis=-1).reshape(shape)

        if self.tmc=='GP':
            """
            - we follow Eq.(4) of Brady et al., PRD84, 074008 (2011)
            - extra factor of x included for F2 and FL projections
            """

            h2=1/N*(-1+XI**(-N))                          
            g2=(N+XI*(1-N)-XI**(1-N))/N/(1-N)             
            h3=1/(-N)*(1-XI**(-N))

            if '22' in proj:
                C1 = (1+rho)**2/(4*rho**3)
                C2 = 3*x*(rho**2-1)/(2*rho**4)
                C3 = (rho**2-1)/(2*x*rho)
                result  = np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))
                result += np.einsum('...x,...xi->...xi',C2,h2)   /XI**(-N+1)
                result += np.einsum('...x,...xi->...xi',C2*C3,g2)/XI**(-N+1)
  
            elif 'LL' in proj:
                C1= (1+rho)**2/(4*rho)
                result  = np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))
  
            elif 'L2' in proj:
                C2= x*(rho**2-1)/(rho**2)
                C3= (rho**2-1)/(2*x*rho)
                result  = np.einsum('...x,...xi->...xi',C2,h2)   /XI**(-N+1)
                result += np.einsum('...x,...xi->...xi',C2*C3,g2)/XI**(-N+1)
  
            elif '33' in proj:
                C1= (1+rho)/(2*rho**2)
                C2= (rho**2-1)/(2*rho**3)
                result  = np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))
                result += np.einsum('...x,...xi->...xi',C2,h3)   /XI**(-N)

            return result

        """
        - we follow Eq.(17) of "What does kinematical target mass sensitivity in DIS reveal about hadron structure?"
        """
        if self.tmc=='AOT':
            if   '22' in proj: C1 = (1+rho)/(2*rho**2)
            elif 'LL' in proj: C1 = np.ones(tuple(x.shape),dtype=complex) 
            elif 'L2' in proj: C1 = (rho-1)/2
            elif '33' in proj: C1 = 1/rho
            return np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))

    #--inverse Mellin transform
    def invert(self,x,F,a=0):

        N  = self.N
        XN = np.repeat(x,self.N.shape[0],axis=0).reshape(len(x),self.N.shape[0])**(-N+a)
        
        phase = self.mell.phase
        W = self.mell.W * self.mell.JAC

        F = np.einsum('i,xi,xi->x',W,XN,F)
        return np.imag(phase*F)/np.pi


class THEORY_new_but_not_working(AUX):

    def __init__(self):

        self.mell=conf['mellin']
        self.mellnpts = self.mell.density
        self.N=self.mell.N
        self.M=self.mell.N
        self.pdf=conf['pdf']
        if 'ht4' in conf: self.ht4 = conf['ht4']


        if 'dsmf' in conf: self.dsmf=conf['dsmf']
        if 'hsmf' in conf: self.hsmf=conf['hsmf']

        if 'dsmf_type' in conf: self.dsmf_type = conf['dsmf_type']
        if 'hsmf_type' in conf: self.hsmf_type = conf['hsmf_type']
 
        self.setup()
  
        if conf['order']=='LO':  self.order=0
        if conf['order']=='NLO': self.order=1
  
        # nuclear smearing
        if   'idis nuc' in conf and conf['idis nuc']: self.nuc = True
        elif 'nuc'      in conf and conf['nuc']:      self.nuc = True #--for backward compatibility
        else:                                         self.nuc = False
 
        self.Nf     = lambda Q2: conf['alphaS'].get_Nf(Q2)
        self.alphaS = lambda Q2: conf['alphaS'].get_alphaS(Q2)

        #--einsum storage
        self.storage = {}

        #--target mass corrections
        if   'idis tmc' in conf: self.tmc = conf['idis tmc']
        elif 'idis_tmc' in conf: self.tmc = conf['idis_tmc']
        elif 'tmc'      in conf: self.tmc = conf['tmc'] #--for backward compatibility
        else:                    self.tmc = False

        #--higher twist
        if   'idis ht' in conf: self.ht_flag = conf['idis ht']
        elif 'ht'      in conf: self.ht_flag = conf['ht']
        else:                   self.ht_flag = False

        if self.nuc:
            self.setup_nucleus()
            self.grids = {}
            self.load_grids()

        #--offshell corrections
        if   'idis offpdf' in conf: self.offpdf = conf['idis offpdf']
        elif 'idis_offpdf' in conf: self.offpdf = conf['idis_offpdf']
        elif 'offpdf'      in conf: self.offpdf = conf['offpdf'] #--for backward compatibility
        else: self.offpdf = False

        #--cannot have offshell corrections without nuclear corrections
        if self.offpdf: self.nuc = True

        if self.offpdf: self.off_model = OFFSHELL_MODEL()

    ###################################
    #--nuclear smearing grid generation
    ###################################
    #--set up nuclear smearing grids
    def setup_nucleus(self):
        eps_d       = -0.00222
        eps_h       = -0.00772
        eps_t       = -0.00848
        mp          = 0.93827231
        mn          = 0.93956563
        mN          = (mp + mn)/2.0
        mD          =   mp +   mn + eps_d
        mH          = 2*mp +   mn + eps_h
        mT          =   mp + 2*mn + eps_t
        self.ymax = {}
        self.ymax['d']  = mD/mN
        self.ymax['h']  = mH/mN
        self.ymax['t']  = mT/mN
        self.ng = 100
        self.gX,self.gW=np.polynomial.legendre.leggauss(self.ng)

    def load_grids(self,process='idis'):

        if '%s tabs'%process in conf: tabs = conf['%s tabs'%process]
        else: return

        data_dir = '%s/database/%s/expdata' %(os.environ['FITPACK'],process)
        grid_dir = '%s/grids/idis'%(os.environ['FITPACK'])

        self.grids[process] = {}

        checkdir(grid_dir)
        done_grids = os.listdir(grid_dir)

        for idx in tabs:
            #--load full dataset prior to cuts
            tab = pd.read_excel('%s/%s.xlsx'%(data_dir,idx)).to_dict(orient='list')
            x   = np.array(tab['X'])
            Q2  = np.array(tab['Q2'])
            obs = tab['obs'][0]
            tar = tab['target'][0].strip()

            tars = []
            if '/' in tar:
                num,den = tar.split('/')
                tars = [num,den]
            else:
                tars = [tar]

            if tar=='p' or tar=='n': continue

            dsmf_type = self.dsmf_type
            hsmf_type = self.hsmf_type

            name = '%s_%s_order=%s_ng=%s_mellnpts=%s_dsmftype=%s_hsmftype=%s_tmc=%s.npy'%(process,idx,self.order,self.ng,self.mellnpts,dsmf_type,hsmf_type,self.tmc)
            filename = '%s/%s'%(grid_dir,name)
            if name in done_grids:
                self.grids[process][idx] = np.load(filename,allow_pickle=True).item()
                print('Loading nuclear smearing grid %s'%filename)
                continue
            grid = {}
            print('generating nuclear smearing grid for %s dataset %s at %s with mell npts = %s,dsmf type=%s, hsmf type=%s, tmc=%s'\
                %(process,idx,self.order,self.mellnpts,dsmf_type,hsmf_type,self.tmc))

            for tar in tars:
 
                if tar in ['p', 'n']: continue
                grid[tar] = {}

                #--generate both onshell and offshell components
                for kind in ['onshell','offshell']:
                    grid[tar][kind] = {}
                    for channel in ['Q','G']:
                        if tar in ['d']:
                            grid[tar][kind]['F2_%s'%channel]      = self.gen_smf(x,Q2,'d','F2',kind,channel)
                            grid[tar][kind]['FL_%s'%channel]      = self.gen_smf(x,Q2,'d','FL',kind,channel)
                            grid[tar][kind]['F3_%s'%channel]      = self.gen_smf(x,Q2,'d','F3',kind,channel)
                            #grid[tar][kind]['f22']     = self.gen_smf(x,Q2,'d','f22', None,kind)
                            #grid[tar][kind]['fLL']     = self.gen_smf(x,Q2,'d','fLL', None,kind)
                            #grid[tar][kind]['fL2']     = self.gen_smf(x,Q2,'d','fL2', None,kind)
                            #grid[tar][kind]['f33']     = self.gen_smf(x,Q2,'d','f33', None,kind)
                        if tar in ['h']:
                            grid[tar][kind]['F2p_%s'%channel]     = self.gen_smf(x,Q2,'h','F2',kind,'p',channel)
                            grid[tar][kind]['FLp_%s'%channel]     = self.gen_smf(x,Q2,'h','FL',kind,'p',channel)
                            grid[tar][kind]['F3p_%s'%channel]     = self.gen_smf(x,Q2,'h','F3',kind,'p',channel)
                            grid[tar][kind]['F2n_%s'%channel]     = self.gen_smf(x,Q2,'h','F2',kind,'n',channel)
                            grid[tar][kind]['FLn_%s'%channel]     = self.gen_smf(x,Q2,'h','FL',kind,'n',channel)
                            grid[tar][kind]['F3n_%s'%channel]     = self.gen_smf(x,Q2,'h','F3',kind,'n',channel)
                        #--for tritium, take helium but switch p <--> n
                        if tar in ['t']:
                            grid[tar][kind]['F2p_%s'%channel]     = self.gen_smf(x,Q2,'t','F2',kind,'n',channel)
                            grid[tar][kind]['FLp_%s'%channel]     = self.gen_smf(x,Q2,'t','FL',kind,'n',channel)
                            grid[tar][kind]['F3p_%s'%channel]     = self.gen_smf(x,Q2,'t','F3',kind,'n',channel)
                            grid[tar][kind]['F2n_%s'%channel]     = self.gen_smf(x,Q2,'t','F2',kind,'p',channel)
                            grid[tar][kind]['FLn_%s'%channel]     = self.gen_smf(x,Q2,'t','FL',kind,'p',channel)
                            grid[tar][kind]['F3n_%s'%channel]     = self.gen_smf(x,Q2,'t','F3',kind,'p',channel)

            np.save(filename,grid)
            print('Saving grid for %s dataset %s with ng=%s to %s'%(process,idx,self.ng,filename))

            self.grids[process][idx] = grid

    def load_grid_custom(self,x,Q2,tar):

        if tar=='p' or tar=='n': return

        #--no process
        process = None
        if process not in self.grids: self.grids[process] = {}

        idx = 'tar=%s,x=%s,Q2=%s'%(tar,x,Q2)
        if idx in self.grids[process]: return

        grid = {}
        print('generating custom grid for tar=%s...'%(tar))

        grid[tar] = {}

        #--generate both onshell and offshell components
        for kind in ['onshell','offshell']:
            grid[tar][kind] = {}
            if tar in ['d']:
                grid[tar][kind]['f22_F2']  = self.gen_smf(x,Q2,'d','f22', 'F2',kind)
                grid[tar][kind]['fLL_FL']  = self.gen_smf(x,Q2,'d','fLL', 'FL',kind)
                grid[tar][kind]['fL2_F2']  = self.gen_smf(x,Q2,'d','fL2', 'F2',kind)
                grid[tar][kind]['f33_F3']  = self.gen_smf(x,Q2,'d','f33', 'F3',kind)
                grid[tar][kind]['f22']     = self.gen_smf(x,Q2,'d','f22', None,kind)
                grid[tar][kind]['fLL']     = self.gen_smf(x,Q2,'d','fLL', None,kind)
                grid[tar][kind]['fL2']     = self.gen_smf(x,Q2,'d','fL2', None,kind)
                grid[tar][kind]['f33']     = self.gen_smf(x,Q2,'d','f33', None,kind)
            if tar in ['h']:
                grid[tar][kind]['f22p_F2'] = self.gen_smf(x,Q2,'h','f22p','F2',kind)
                grid[tar][kind]['f22n_F2'] = self.gen_smf(x,Q2,'h','f22n','F2',kind)
                grid[tar][kind]['fLLp_FL'] = self.gen_smf(x,Q2,'h','fLLp','FL',kind)
                grid[tar][kind]['fLLn_FL'] = self.gen_smf(x,Q2,'h','fLLn','FL',kind)
                grid[tar][kind]['fL2p_F2'] = self.gen_smf(x,Q2,'h','fL2p','F2',kind)
                grid[tar][kind]['fL2n_F2'] = self.gen_smf(x,Q2,'h','fL2n','F2',kind)
                grid[tar][kind]['f33p_F3'] = self.gen_smf(x,Q2,'h','f33p','F3',kind)
                grid[tar][kind]['f33n_F3'] = self.gen_smf(x,Q2,'h','f33n','F3',kind)
                grid[tar][kind]['f22p']     = self.gen_smf(x,Q2,'h','f22p',None,kind)
                grid[tar][kind]['fLLp']     = self.gen_smf(x,Q2,'h','fLLp',None,kind)
                grid[tar][kind]['fL2p']     = self.gen_smf(x,Q2,'h','fL2p',None,kind)
                grid[tar][kind]['f33p']     = self.gen_smf(x,Q2,'h','f33p',None,kind)
                grid[tar][kind]['f22n']     = self.gen_smf(x,Q2,'h','f22n',None,kind)
                grid[tar][kind]['fLLn']     = self.gen_smf(x,Q2,'h','fLLn',None,kind)
                grid[tar][kind]['fL2n']     = self.gen_smf(x,Q2,'h','fL2n',None,kind)
                grid[tar][kind]['f33n']     = self.gen_smf(x,Q2,'h','f33n',None,kind)
            #--for tritium, take helium but switch p <--> n
            if tar in ['t']:
                grid[tar][kind]['f22p_F2'] = self.gen_smf(x,Q2,'h','f22n','F2',kind)
                grid[tar][kind]['f22n_F2'] = self.gen_smf(x,Q2,'h','f22p','F2',kind)
                grid[tar][kind]['fLLp_FL'] = self.gen_smf(x,Q2,'h','fLLn','FL',kind)
                grid[tar][kind]['fLLn_FL'] = self.gen_smf(x,Q2,'h','fLLp','FL',kind)
                grid[tar][kind]['fL2p_F2'] = self.gen_smf(x,Q2,'h','fL2n','F2',kind)
                grid[tar][kind]['fL2n_F2'] = self.gen_smf(x,Q2,'h','fL2p','F2',kind)
                grid[tar][kind]['f33p_F3'] = self.gen_smf(x,Q2,'h','f33n','F3',kind)
                grid[tar][kind]['f33n_F3'] = self.gen_smf(x,Q2,'h','f33p','F3',kind)
                grid[tar][kind]['f22p']     = self.gen_smf(x,Q2,'h','f22n',None,kind)
                grid[tar][kind]['fLLp']     = self.gen_smf(x,Q2,'h','fLLn',None,kind)
                grid[tar][kind]['fL2p']     = self.gen_smf(x,Q2,'h','fL2n',None,kind)
                grid[tar][kind]['f33p']     = self.gen_smf(x,Q2,'h','f33n',None,kind)
                grid[tar][kind]['f22n']     = self.gen_smf(x,Q2,'h','f22p',None,kind)
                grid[tar][kind]['fLLn']     = self.gen_smf(x,Q2,'h','fLLp',None,kind)
                grid[tar][kind]['fL2n']     = self.gen_smf(x,Q2,'h','fL2p',None,kind)
                grid[tar][kind]['f33n']     = self.gen_smf(x,Q2,'h','f33p',None,kind)


        self.grids[process][idx] = grid

    def gen_smf(self,X,Q2,tar,stf,kind,nucleon='p',channel='Q'):

        #--with TMCs, must evaluate at xi instead of x
        if self.tmc!=False: XI = self.get_xi(X,Q2)
        else:               XI = X

        W = self.mell.W * self.mell.JAC
        nx = X.shape[0]
        #--Gaussian quadrature setup
        XM,  gXM = np.meshgrid(X,  self.gX)
        XIM, gXM = np.meshgrid(XI, self.gX)
        Q2M, gWM = np.meshgrid(Q2, self.gW)
        b = self.ymax[tar]
      
        a = XM
        YM   = 0.5*(b-a)*gXM+0.5*(a+b)
        XM_YM = XM/YM

        a = XIM
        YIM   = 0.5*(b-a)*gXM+0.5*(a+b)
        JIM   = 0.5*(b-a) 
        XIM_YIM = XIM/YIM

        if   tar in ['d']:     SMF=self.dsmf
        elif tar in ['h','t']: SMF=self.hsmf

        if tar in ['d']:     fXX = lambda proj: SMF.get_fXX2('f%s'  %(proj)        ,kind,XIM,Q2M,YIM)
        if tar in ['h','t']: fXX = lambda proj: SMF.get_fXX2('f%s%s'%(proj,nucleon),kind,XIM,Q2M,YIM)

        #--get TMC, evaluated at xi/y (note that it takes x as the argument and then converts to xi)
        TMC = lambda proj: self.get_TMC(XM_YM,Q2,proj)

        N = self.N
        #--get hard kernels
        #--LO
        if 'L' in stf: C = np.zeros((Q2.shape[0],N.shape[0]),dtype=complex)
        else:
            if channel=='Q': C = np.ones ((Q2.shape[0],N.shape[0]),dtype=complex)
            if channel=='G': C = np.zeros((Q2.shape[0],N.shape[0]),dtype=complex)

        #--NLO
        if self.order >= 1:
            alphaS = np.array([self.alphaS(q2) for q2 in Q2])
            C+= np.einsum('x,i->xi',alphaS,self.C[stf[-1]][channel])/4/np.pi

        XIN  = np.repeat(XIM_YIM,self.N.shape[0],axis=1).reshape(XIM_YIM.shape[0],XIM_YIM.shape[1],self.N.shape[0])
        if stf=='F2':
            f = np.einsum('gxi,gxi,gx->gxi',XIN**(-N+1),TMC('22'),fXX('22'))
        elif stf=='F3':
            f = np.einsum('gxi,gxi,gx->gxi',XIN**(-N)  ,TMC('33'),fXX('33'))
        elif stf=='FL':
            f = np.einsum('gxi,gxi,gx->gxi',XIN**(-N+1),TMC('LL'),fXX('LL'))\
              + np.einsum('gxi,gxi,gx->gxi',XIN**(-N+1),TMC('L2'),fXX('LL'))\
              + np.einsum('gxi,gxi,gx->gxi',XIN**(-N+1),TMC('22'),fXX('L2'))
 

        fN = np.einsum('i,gx,gx,gxi,xi->xi',W,gWM,JIM,f,C) 
        return fN

    ###########################

    # theory calculations
    def get_FX(self,idx,stf,x,Q2,tar,process='idis',channel='all'):

        if 'W' in stf: k = stf[-2]
        else:          k = stf[-1]

        if self.tmc!=False: xi = self.get_xi(x,Q2)
        else:               xi = x 

        # Nucleon structure functions
        if tar in ['p','n']: 
            #--get FXN in mellin space, without TMCs or HTs
            FXN = self.get_FXN_mell(stf,Q2,tar,channel=channel)

            #--multiply by target mass corrections
            FXN = np.einsum('xi,xi->xi',FXN,self.get_TMC(x,Q2,'%s%s'%(k,k)))

            #--extra term for FL/WL
            if 'L' in stf and self.tmc!=False:
                FXN += np.einsum('xi,xi->xi',self.get_FXN_mell('F2',Q2,tar),self.get_TMC(x,Q2,'%s2'%(k)))

            #--add higher twists (does not apply for W)
            #if self.ht_flag and 'W' not in stf:
            #    FXN += self.ht4.get_ht_mell(Q2,tar,stf)
           
            if '2' or 'L' in stf: a = 1
            if '3'        in stf: a = 0

            FXN = self.invert(xi,FXN,a)

            #--add higher twists (does not apply for W), without any extra factors of x
            #--evaluate higher twists at xi
            if self.ht_flag and 'W' not in stf:
                FXN += self.invert(xi,self.ht4.get_ht_mell(Q2,tar,stf),0)

            return FXN


        # Nuclear structure functions
        if tar=='d': p, n = 1,1
        if tar=='h': p, n = 2,1
        if tar=='t': p, n = 1,2

        #--no nuclear smearing for W
        if self.nuc==False or 'W' in stf:
            FXp = self.get_FX(idx,stf,x,Q2,'p') 
            FXn = self.get_FX(idx,stf,x,Q2,'n')
            return (p * FXp + n * FXn)/(p+n)

        return self.get_FXA(idx,stf,x,Q2,tar,process=process,channel=channel)

    #--returns mellin space structure function
    #--note that this returns F2/x, FL/x, F3
    #--note that factor of xi**(-N) is NOT included
    def get_FXN_mell(self,stf,Q2,tar,kind='onshell',nucleus=None,channel='all'):

        C = self.C
 
        alphaS = np.array([self.alphaS(q2) for q2 in Q2])
  
        # Functions to setup PDFs and FFs for einsum
        PDF, e2gl = self.setup_pdfs(stf,Q2,tar,kind=kind,nucleus=nucleus,channel=channel)

        ind_mell_LO  = 'kxi->xi'
        ind_mell_NLO = 'x,kxi,i->xi'
        ind_mell_NLO = 'x,kxi,i->xi'
        
        if 'W' in stf: k = stf[-2]
        else:          k = stf[-1]

        if ind_mell_LO not in self.storage:
            self.storage[ind_mell_LO]  = np.einsum_path(ind_mell_LO,        PDF,           optimize='optimal')[0]
            self.storage[ind_mell_NLO] = np.einsum_path(ind_mell_NLO,alphaS,PDF, C[k]['Q'],optimize='optimal')[0] 

        #--LO term
        if '2' or '3' in stf:
            FX    = np.einsum(ind_mell_LO, PDF, optimize = self.storage[ind_mell_LO])

        if 'L' in stf:
            FX = np.zeros((Q2.shape[0],self.N.shape[0]),dtype=complex)

        #--NLO term
        if self.order >= 1:
            FX  += np.einsum(ind_mell_NLO,alphaS,PDF, C[k]['Q'],optimize = self.storage[ind_mell_NLO])/4/np.pi\
                 + np.einsum(ind_mell_NLO,alphaS,e2gl,C[k]['G'],optimize = self.storage[ind_mell_NLO])/4/np.pi
 
        if 'W' in stf: factor = 2.0
        if 'F' in stf: factor = 1.0

        FX = factor*FX 
           
        return FX
 
    def get_FXA(self,idx,stf,x,Q2,tar,process='idis',channel='all'):

        C = self.C

        if tar=='d': p, n = 1,1
        if tar=='h': p, n = 2,1
        if tar=='t': p, n = 1,2

        #--if idx is None, generate grid on the fly
        if idx==None: 
            self.load_grid_custom(x,Q2,tar)
            process=None
            idx = 'tar=%s,x=%s,Q2=%s'%(tar,x,Q2)
            grid = self.grids[process][idx]
            i = np.array([k for k in range(len(x))])

        #--load grid and place cuts
        else:
            grid = self.grids[process][idx]
            i = conf['%s tabs'%process][idx]['idx']

        if self.offpdf: kinds = ['onshell','offshell']
        else:           kinds = ['onshell']
        
        if 'W' in stf: k = stf[-2]
        else:          k = stf[-1]

        FXp = np.zeros(Q2.shape[0])
        FXn = np.zeros(Q2.shape[0])

        invert = lambda F: np.imag(phase*F)/np.pi

        for kind in kinds:
            if tar=='d':         f = lambda nucleon,channel: grid[tar][kind]['%s_%s'  %(stf,        channel)][i] 
            if tar in ['h','t']: f = lambda nucleon,channel: grid[tar][kind]['%s%s_%s'%(stf,nucleon,channel)][i] 

            alphaS = np.array([self.alphaS(q2) for q2 in Q2])

            #for channel in ['Q','G']:
            #    for nucleon in ['p','n']:
  
            # Functions to setup PDFs and FFs for einsum
            PDF, e2gl = {},{}
            PDF['p'], e2gl['p'] = self.setup_pdfs(stf,Q2,'p',kind=kind,nucleus=tar,channel=channel)
            PDF['n'], e2gl['n'] = self.setup_pdfs(stf,Q2,'n',kind=kind,nucleus=tar,channel=channel)


            ind  = 'xi,kxi->x'
            
            if ind not in self.storage:
                self.storage[ind]  = np.einsum_path(ind,f('p','Q'),PDF['p'],optimize='optimal')[0]


            #--quark term
            _FXp    = np.einsum(ind,f('p','Q'),PDF['p'],optimize = self.storage[ind])
            _FXn    = np.einsum(ind,f('n','Q'),PDF['n'],optimize = self.storage[ind])

            #--gluon term
            _FXp   += np.einsum(ind,f('p','G'),e2gl['p'],optimize = self.storage[ind])
            _FXn   += np.einsum(ind,f('n','G'),e2gl['n'],optimize = self.storage[ind])

            #--LO terms
            #--F22,FLL,F33 term
            #_FXp    = np.einsum(ind_LO,f('p'),PDF['p'],optimize = self.storage[ind_LO])
            #_FXn    = np.einsum(ind_LO,f('n'),PDF['n'],optimize = self.storage[ind_LO])

            ##--FL2 term
            #if 'L' in stf:
            #    _FXp   = 0.0
            #    _FXn   = 0.0


            ##--NLO terms
            #if self.order >= 1:
            #    #--F22,FLL,F33 term
            #    _FXp   += np.einsum(ind_NLO,f('p'),alphaS,PDF ['p'],C[k]['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
            #            + np.einsum(ind_NLO,f('p'),alphaS,e2gl['p'],C[k]['G'], optimize = self.storage[ind_NLO])/4/np.pi


            #    _FXn   += np.einsum(ind_NLO,f('n'),alphaS,PDF ['n'],C[k]['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
            #            + np.einsum(ind_NLO,f('n'),alphaS,e2gl['n'],C[k]['G'], optimize = self.storage[ind_NLO])/4/np.pi

            #    #--FL2 terms
            #    #if 'L' in stf:
            #    #    _FXp   += np.einsum(ind_NLO,f('%s2'%(k),'p','F2'),    alphaS,PDF ['p'],C[k]['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
            #    #            + np.einsum(ind_NLO,f('%s2'%(k),'p','F2'),    alphaS,e2gl['p'],C[k]['G'], optimize = self.storage[ind_NLO])/4/np.pi
  

            #    #    _FXn   += np.einsum(ind_NLO,f('%s2'%(k),'n','F2'),    alphaS,PDF ['n'],C[k]['Q'], optimize = self.storage[ind_NLO])/4/np.pi\
            #    #            + np.einsum(ind_NLO,f('%s2'%(k),'n','F2'),    alphaS,e2gl['n'],C[k]['G'], optimize = self.storage[ind_NLO])/4/np.pi


              
            phase = self.mell.phase

            FXp += np.imag(phase*_FXp)/np.pi
            FXn += np.imag(phase*_FXn)/np.pi

        #--add ht
        if self.ht_flag and 'W' not in stf:
            if tar=='d':         f = lambda nucleon: grid[tar]['onshell']['%s'  %(stf)]        [i] 
            if tar in ['h','t']: f = lambda nucleon: grid[tar]['onshell']['%s%s'%(stf,nucleon)][i] 

            ht = lambda stf,tar: self.ht4.get_ht_mell(Q2,tar,stf)
            _FXp  = np.einsum('xi,xi->x',f('p'),ht(stf,'p'))
            _FXn  = np.einsum('xi,xi->x',f('n'),ht(stf,'n'))

            #if 'L' in stf:
            #    _FXp  += np.einsum('xi,xi->x',f('%s2'%(k),'p'),ht('F2','p'))
            #    _FXn  += np.einsum('xi,xi->x',f('%s2'%(k),'n'),ht('F2','n'))

            #--note that extra factor of xi is removed
            FXp += np.imag(phase*_FXp)/np.pi
            FXn += np.imag(phase*_FXn)/np.pi

        return (p*FXp + n*FXn)/(p+n)

    def setup_pdfs(self,stf,Q2,tar,kind='onshell',nucleus=None,channel='all'):
  
        if stf in ['F2', 'FL']:    aX = self.get_aX('aXp',   Q2, channel) 
        if stf in ['F3']:          aX = self.get_aX('aXm',   Q2, channel) 
        if stf in ['W2+','WL+']:   aX = self.get_aX('aXpWp', Q2, channel) 
        if stf in ['W3+']:         aX = self.get_aX('aXmWp', Q2, channel) 
        if stf in ['W2-','WL-']:   aX = self.get_aX('aXpWm', Q2, channel) 
        if stf in ['W3-']:         aX = self.get_aX('aXmWm', Q2, channel) 
       
        if kind=='onshell':  pdf = self.pdf
        if kind=='offshell': pdf = conf['off pdf']
 
        L = Q2.shape[0] 
        for q2 in Q2: pdf.evolve(q2)
 
        PDF  = np.zeros((11,L,self.N.size),dtype=complex)
        e2gl = np.zeros((11,L,self.N.size),dtype=complex)
      
        if kind=='onshell':
            #--isospin relations for free nuclei
            q = {}
            if tar == 'p': flav = {-5:'bb',-4:'cb',-3:'sb',-2:'db',-1:'ub',0:'g',1:'u',2:'d',3:'s',4:'c',5:'b'}
            if tar == 'n': flav = {-5:'bb',-4:'cb',-3:'sb',-2:'ub',-1:'db',0:'g',1:'d',2:'u',3:'s',4:'c',5:'b'}
            for i in flav:
                q[i] = np.array([pdf.storage[q2][flav[i]] for q2 in Q2])
        elif kind=='offshell':
            q = {}
            flav = {-5:'bb',-4:'cb',-3:'sb',-2:'db',-1:'ub',0:'g',1:'u',2:'d',3:'s',4:'c',5:'b'}
            for i in flav:
                q[i] = np.array([pdf.storage[q2][flav[i]] for q2 in Q2])

            q = self.off_model.get_model(q,tar,nucleus)
         
        #--take Nf into account
        Nf_mask = {}
        Nf_mask[0] = np.ones(L) 
        Nf_mask[1] = np.ones(L) 
        Nf_mask[2] = np.ones(L) 
        Nf_mask[3] = np.ones(L) 
        Nf_mask[4] = np.ones(L)*(Q2 > conf['aux'].mc2)
        Nf_mask[5] = np.ones(L)*(Q2 > conf['aux'].mb2) 

        #--no bottom quark for W
        if 'W' in stf: Nf_mask[5] = np.zeros(L)

        for i in [1,2,3,4,5]: Nf_mask[-i] = Nf_mask[i]

 
        for i in range(-5,6):
            if i==0: continue
            PDF[i]  = np.einsum('x,x,xi->xi',Nf_mask[i],aX[i],q[i])
            e2gl[i] = np.einsum('x,x,xi->xi',Nf_mask[i],aX[i],q[0])
     
        return PDF, e2gl  

    #--TMCs in Mellin space
    def get_xi(self,x,Q2):
        mu2=conf['aux'].M**2/Q2
        rho=(1+4*mu2*x**2)**0.5
        xi=2*x/(1+rho)
        return xi

    def get_TMC(self,x,Q2,proj):

        N = self.N
        shape = tuple(x.shape)+tuple(N.shape)

        #--with no TMCS, return zero mass WW
        if self.tmc==False:
            if 'L2' in proj: return np.zeros(shape,dtype=complex)
            else:            return np.ones (shape,dtype=complex)

        xi = self.get_xi(x,Q2)
        mu2=conf['aux'].M**2/Q2
        rho=(1+4*mu2*x**2)**0.5

        XI = np.repeat(xi,self.N.shape[0],axis=-1).reshape(shape)

        if self.tmc=='GP':
            """
            - we follow Eq.(4) of Brady et al., PRD84, 074008 (2011)
            - extra factor of x included for F2 and FL projections
            """

            h2=1/N*(-1+XI**(-N))                          
            g2=(N+XI*(1-N)-XI**(1-N))/N/(1-N)             
            h3=1/(-N)*(1-XI**(-N))

            if '22' in proj:
                C1 = (1+rho)**2/(4*rho**3)
                C2 = 3*x*(rho**2-1)/(2*rho**4)
                C3 = (rho**2-1)/(2*x*rho)
                result  = np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))
                result += np.einsum('...x,...xi->...xi',C2,h2)   /XI**(-N+1)
                result += np.einsum('...x,...xi->...xi',C2*C3,g2)/XI**(-N+1)
  
            elif 'LL' in proj:
                C1= (1+rho)**2/(4*rho)
                result  = np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))
  
            elif 'L2' in proj:
                C2= x*(rho**2-1)/(rho**2)
                C3= (rho**2-1)/(2*x*rho)
                result  = np.einsum('...x,...xi->...xi',C2,h2)   /XI**(-N+1)
                result += np.einsum('...x,...xi->...xi',C2*C3,g2)/XI**(-N+1)
  
            elif '33' in proj:
                C1= (1+rho)/(2*rho**2)
                C2= (rho**2-1)/(2*rho**3)
                result  = np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))
                result += np.einsum('...x,...xi->...xi',C2,h3)   /XI**(-N)

            return result

        """
        - we follow Eq.(17) of "What does kinematical target mass sensitivity in DIS reveal about hadron structure?"
        """
        if self.tmc=='AOT':
            if   '22' in proj: C1 = (1+rho)/(2*rho**2)
            elif 'LL' in proj: C1 = np.ones(tuple(x.shape),dtype=complex) 
            elif 'L2' in proj: C1 = (rho-1)/2
            elif '33' in proj: C1 = 1/rho
            return np.einsum('...x,...xi->...xi',C1,np.ones(shape,dtype=complex))

    #--inverse Mellin transform
    def invert(self,x,F,a=0):

        N  = self.N
        XN = np.repeat(x,self.N.shape[0],axis=0).reshape(len(x),self.N.shape[0])**(-N+a)
        
        phase = self.mell.phase
        W = self.mell.W * self.mell.JAC

        F = np.einsum('i,xi,xi->x',W,XN,F)
        return np.imag(phase*F)/np.pi


 
