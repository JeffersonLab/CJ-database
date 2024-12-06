#!/usr/bin/env python
import sys,os
import numpy as np
import pandas as pd
import time
from scipy.integrate   import quad,fixed_quad
from scipy.interpolate import griddata
from obslib.idis.aux   import AUX 
from tools.config      import conf
from nuclib            import deuterium,helium

class HT:

    def __init__(self):

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

        #--polarized structure functions
        self.params['g1p']=np.array([0.0,0.5,3,0,0])
        self.params['g2p']=np.array([0.0,0.5,3,0,0])

        self.params['g1n']=np.array([0.0,0.5,3,0,0])
        self.params['g2n']=np.array([0.0,0.5,3,0,0])

        self.FLAV = ['F2p','FLp','F3p','F2n','FLn','F3n']
        self.FLAV.extend(['W2p','WLp','W3p','W2n','WLn','W3n'])
        self.FLAV.extend(['g1p','g2p','g1n','g2n'])
        self.PAR = ['N','a','b','c','d']
        
    def setup(self):
        pass

    def get_state(self):
        return (self.params)

    def set_state(self,state):
        self.params = state

    def get_ht(self,x,Q2,nucleon,stf):
        N,a,b,c,d = self.params['%s%s'%(stf,nucleon)]
        return N*x**a*(1-x)**b*(1+c*x**0.5+d*x)/Q2 

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


class THEORY(AUX):
  
    def __init__(self):

        self.mellin  = conf['mellin']
        self.setup()

        if conf['order']=='LO' : self.order=0
        if conf['order']=='NLO': self.order=1

        if 'ht4'    in conf: self.ht4   = conf['ht4']  
        if 'dsmf'   in conf: self.dsmf  = conf['dsmf']
        if 'hsmf'   in conf: self.hsmf  = conf['hsmf']

        #--heavy quark setups
        if 'hq' in conf:
            self.hq=conf['hq']
        else:
            self.hq=False
        self.mc2  = conf['aux'].mc**2
        self.mb2  = conf['aux'].mb**2
   
        #--target mass setups
        if 'tmc' in conf:
            self.tmc=conf['tmc']
        else:
            self.tmc=False
        self.M2=conf['aux'].M2

        #--power corrections setup
        #--flag for higher twist corrections (multiplicative or additive)
        if 'ht'  in conf:
            self.ht=conf['ht']
            #--multiplicative (mult) or additive (add)
            self.ht_type = 'mult'
            if 'ht type' in conf: self.ht_type = conf['ht type']
        else:
            self.ht=False
            self.ht_type=False

        #--nuclear settings
        if 'idis nuc'  in conf:
            self.nuc=conf['idis nuc']
        elif 'nuc'  in conf:
            self.nuc=conf['nuc']
        else:
            self.nuc=False

        if 'offpdf' in conf:
            self.offpdf = conf['offpdf']
            self.model  = OFFSHELL_MODEL()
        else:
            self.offpdf = False

        self.setup_nucleus()

        #--setup interplation grids
        self.setup_interpolation() 

    #--twist 2 unpolarized structure functions 
  
    def get_T2CFX(self,stf,nucleon,Q2,Nf=None,evolve=True,channel='all',offpdf=False,nucleus='d',iveto=[None]):
        """
        CF(2,L,3) = F2/x,FL/x,F3  
        """
        if offpdf: pdf = conf['off pdf']
        else:      pdf = conf['pdf']
 
        if evolve: pdf.evolve(Q2)

        g = pdf.storage[Q2]['g']    
        if Nf==None: Nf=conf['alphaS'].get_Nf(Q2) 
        a=conf['alphaS'].get_a(Q2) 

        if stf=='F2':
            CQ = self.C2Q[0] + a*self.order*self.C2Q[1] 
            CG = self.C2G[0] + a*self.order*self.C2G[1]
            q=np.copy(pdf.storage[Q2]['qp'])   
            aX='ap'
 
        elif stf=='FL':
            CQ = a*self.order*self.CLQ[1]
            CG = a*self.order*self.CLG[1]
            q=np.copy(pdf.storage[Q2]['qp'])
            aX='ap'
  
        elif stf=='F3':
            CQ = self.C3Q[0] + a*self.order*self.C3Q[1]
            CG = 0
            q=np.copy(pdf.storage[Q2]['qm']) 
            aX='am'

        #--isospin symmetry
        if offpdf==False:
            if nucleon=='n':
                qup=np.copy(q[1])
                qdn=np.copy(q[2])
                q[1]=qdn
                q[2]=qup
            elif nucleon=='d':
                qup=np.copy(q[1])
                qdn=np.copy(q[2])
                q[1]=0.5*(qup+qdn)
                q[2]=0.5*(qup+qdn)
        else:
            q = self.model.get_model(q,nucleon,nucleus)

        
        if iveto[0]!=None:
            for i in range(5):
                q[i+1]*=iveto[i]
            
        FX  = np.zeros(self.mellin.N.size,dtype=complex) 
        for i in range(1,Nf+1):
            aXval = self.get_aX(aX,i,Q2,channel)
            FX+=aXval*(CQ*q[i] + 2*CG*g)
 
        return FX

    def get_T2CWX(self,stf,nucleon,Q2,sign,evolve=True):  
        """
        CF(2,L,3) = W2/x,WL/x,W3  
        """ 
        if evolve: conf['pdf'].evolve(Q2)
        g =conf['pdf'].storage[Q2]['g'] 
        Nf=conf['alphaS'].get_Nf(Q2) 
        a=conf['alphaS'].get_a(Q2)  
  
        if stf=='W2':
            CQ = self.C2Q[0] + a*self.order*self.C2Q[1] 
            CG = self.C2G[0] + a*self.order*self.C2G[1]
  
        elif stf=='WL':
            CQ = a*self.order*self.CLQ[1]
            CG = a*self.order*self.CLG[1]
  
        elif stf=='W3':
            CQ = self.C3Q[0] + a*self.order*self.C3Q[1]
            CG = 0
  
        u=conf['pdf'].storage[Q2]['u']
        d=conf['pdf'].storage[Q2]['d']
        s=conf['pdf'].storage[Q2]['s']
        c=conf['pdf'].storage[Q2]['c']
        b=conf['pdf'].storage[Q2]['b']
  
        ub=conf['pdf'].storage[Q2]['ub']
        db=conf['pdf'].storage[Q2]['db']
        sb=conf['pdf'].storage[Q2]['sb']
        cb=conf['pdf'].storage[Q2]['cb']
        bb=conf['pdf'].storage[Q2]['bb']
  
        U =(CQ*u  + CG*g)
        D =(CQ*d  + CG*g) + (CQ*s  + CG*g) 
        UB=(CQ*ub + CG*g)
        DB=(CQ*db + CG*g) + (CQ*sb + CG*g)
  
        if Nf>3:
            U += CQ*c  + CG*g
            UB+= CQ*cb + CG*g
  
        # cannot produce a top
        #if Nf>4:
        #  D+=  CQ*b  + CG*g
        #  DB+= CQ*bb + CG*g 
 
        #--factor of 2 follows PDG definition
        #--pdg.lbl.gov/2019/reviews/rpp2019-rev-structure-functions.pdf
        #--see equation (18.19)
        if sign==+1:
            if   stf=='W2': return 2*(D+UB)
            elif stf=='WL': return 2*(D+UB)
            elif stf=='W3': return 2*(D-UB)
        elif sign==-1:
            if   stf=='W2': return 2*(U+DB)
            elif stf=='WL': return 2*(U+DB)
            elif stf=='W3': return 2*(U-DB)

    #--tmc factors

    def get_TMC_GP(self,N,x,Q2,proj):
        """
        - we follow Eq.(4) of Brady et al., PRD84, 074008 (2011)
        - extra factor of x included for F2 and FL projections
        """
        mu2=self.M2/Q2
        xx=x*x
        rho=(1+4*mu2*xx)**0.5
        xi=2*x/(1+rho)              #kinematical factors
  
        if proj=='22':
  
            C1= (1+rho)**2/(4*rho**3)
            C2= 3*x*(rho**2-1)/(2*rho**4)
            C3= (rho**2-1)/(2*x*rho)                      
            h2=1/N*(-1+xi**(-N))                          
            g2=(N+xi*(1-N)-xi**(1-N))/N/(1-N)             
            return C1*xi**(-N+1) + C2*(h2 + C3*g2)        
  
        elif proj=='LL':
  
            C1= (1+rho)**2/(4*rho)
            return C1*xi**(-N+1) 
  
        elif proj=='L2':
  
            C2= x*(rho**2-1)/(rho**2)
            C3= (rho**2-1)/(2*x*rho)
            h2=1/N*(-1+xi**(-N))
            g2=(N+xi*(1-N)-xi**(1-N))/N/(1-N)
            return C2*(h2 + C3*g2)
  
        elif proj=='33':
  
            C1= (1+rho)/(2*rho**2)
            C2= (rho**2-1)/(2*rho**3)
            h3=1/(-N)*(1-xi**(-N))
            return C1*xi**(-N) + C2*h3 
  
    def get_TMC_AOT(self,N,x,Q2,proj):
        """
        - we follow Eq.(17) of "What does kinematical target mass sensitivity in DIS reveal about hadron structure?"
        - extra factor of x included for F2 and FL projections
        """
        mu2=self.M2/Q2
        xx=x*x
        rho=(1+4*mu2*xx)**0.5
        xi=2*x/(1+rho)
  
        if proj=='22':
  
            C1= (1+rho)/(2*rho**2)
            return C1*xi**(-N+1)    
  
        elif proj=='LL':
            C1= 1 
            return C1*xi**(-N+1) 
  
        elif proj=='L2':
  
            C1= (rho-1)/2
            return C1*xi**(-N+1) 
  
        elif proj=='33':
            C1= 1/rho
            return C1*xi**(-N)

    #--nucleon  structure functions
  
    def _get_FXN(self,x,Q2,stf='F2',nucleon='p',hq=False,tmc=False,ht=False,channel='all',offpdf=False,nucleus='d',evolve=True,iveto=[None]):

        if hq==False: Nf=None
        else:         Nf=3

        if tmc==False:
            if   stf=='F2': FX= x*self.get_T2CFX('F2',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus,iveto=iveto) 
            elif stf=='FL': FX= x*self.get_T2CFX('FL',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus,iveto=iveto)
            elif stf=='F3': FX=   self.get_T2CFX('F3',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus,iveto=iveto)
            FX*=x**(-self.mellin.N) 

        elif tmc=='GP':
            if   stf=='F2': FX= self.get_T2CFX('F2',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_GP(self.mellin.N,x,Q2,'22')
            elif stf=='FL': FX= self.get_T2CFX('FL',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_GP(self.mellin.N,x,Q2,'LL')\
                               +self.get_T2CFX('F2',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_GP(self.mellin.N,x,Q2,'L2')
            elif stf=='F3': FX= self.get_T2CFX('F3',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_GP(self.mellin.N,x,Q2,'33')
  
        elif tmc=='AOT':  
            if   stf=='F2': FX= self.get_T2CFX('F2',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_AOT(self.mellin.N,x,Q2,'22')
            elif stf=='FL': FX= self.get_T2CFX('FL',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_AOT(self.mellin.N,x,Q2,'LL')\
                               +self.get_T2CFX('F2',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_AOT(self.mellin.N,x,Q2,'L2')
            elif stf=='F3': FX= self.get_T2CFX('F3',nucleon,Q2,Nf,evolve,channel,offpdf,nucleus)*self.get_TMC_AOT(self.mellin.N,x,Q2,'33')

        FX=self.mellin.invert(1,FX)

        if ht: 
            if self.ht_type == 'mult': FX*=(1+self.ht4.get_ht(x,Q2,nucleon,stf))
            if self.ht_type == 'add':  
                if channel=='all':     FX+=     self.ht4.get_ht(x,Q2,nucleon,stf)
                elif channel=='gg':    FX+=     self.ht4.get_ht(x,Q2,nucleon,stf)
                else              :    FX+=     0

        return FX

    #--hq routines (ACOT HQ scheme NLO:  Phys. Rev. D 18 (1978) 242)

    def integrand_HQgF2(self,z,x,Q2,g,hq): 
        """
        integrand of the coefficent HG1 (for F2). 
        last update: Nov 22 2019 
        """
        if   hq=='c': mhq = conf['aux'].mc
        elif hq=='b': mhq = conf['aux'].mb
        eps = mhq**2./Q2
        xhq = 1./(1. + 4.*eps)
        beta2 = 1. - 4.*eps*z/(1.-z)
        beta=np.sqrt(beta2)
        a = z**2. + (1.-z)*(1.-z)
        b = 4.*z*(1.-3.*z)
        c = -8.*z*z
        d = 8.*z*(1.-z) - 1.
        e = -4.*z*(1.-z)
        hg1a = (a+ b*eps + c*eps*eps)*np.log((1.+beta)/(1.-beta))
        hg1b = (d + e*eps)*beta
        return (hg1a+hg1b)*g(x/z)/z

    def integrand_HQgFL(self,z,x,Q2,g,hq): 
        """
        integrand of the coefficient HG1 for FL. No nuclear smearing
        last updte: Nov 26 2019 
        """
        if hq=='c':   mhq = conf['aux'].mc
        elif hq=='b': mhq = conf['aux'].mb

        eps   = mhq**2./Q2
        xhq   = 1./(1. + 4.*eps)
        beta2 = 1. - 4.*eps*z/(1.-z)
        beta  = np.sqrt(beta2)
        a     = 4.*z*(1.-z)
        b     = -8.*z**2.
        fac   = np.log((1. + beta)/(1. - beta))
        hgl1a = beta*a + b*eps*fac
        return hgl1a*g(x/z)/z

    def get_HQgFX(self,x,Q2,hq,fX):
        """
        HGF2 (aka HG1) coefficient of F2^hq (enters at order alpha_s).
        """
        g=lambda x: conf['pdf'].get_xF(x,Q2,'g')/x
        N = self.mellin.N
        if   hq=='c': mhq=conf['aux'].mc
        elif hq=='b': mhq=conf['aux'].mb
        eps = mhq**2./Q2
        xhq = 1./(1. + 4.*eps)
        if x>xhq: 
            return 0
        else:
            if fX=='F2':
                return fixed_quad(np.vectorize(lambda z: self.integrand_HQgF2(z,x,Q2,g,hq)),x,xhq,n=20)[0]
            if fX=='FL':
                return fixed_quad(np.vectorize(lambda z: self.integrand_HQgFL(z,x,Q2,g,hq)),x,xhq,n=20)[0]

    def get_HQqFL(self,x,xhq,Q2,N,hq): 
        """
        HQL1 coefficient of FL^hq (enters at order \alpha_s). it's analytic in Mellin. NO smearing
        last update Nov. 26 2019
        """
        return 5*(xhq-0.2)*4/3.*( xhq**(N+1)/(N + 1) -x**(N+1)/(N + 1))/xhq

    def get_T2CFXHQ(self,x,Q2,stf,hq):
        """
        F2c/x and FLc/x for nucleon
        CF(2,L,3) = F2/x,FL/x,F3  
        """

        N = self.mellin.N
        conf['pdf'].evolve(Q2)
        Nf=conf['alphaS'].get_Nf(Q2) # num of active flav
        q = conf['pdf'].storage[Q2]['qp']
        g = conf['pdf'].storage[Q2]['g'] # gluon[N] were N mellin contour
        a = conf['alphaS'].get_a(Q2)*self.order   # a=alphaS/4pi 

        if hq=='c':
            mhq2   = self.mc2
            eps    = self.mc2/Q2
            eq2    = self.get_aX('ap',4,Q2)
            hq_pdf = q[4] 
  
        if hq=='b': 
            mhq2   = self.mb2
            eps    = self.mb2/Q2
            eq2    = self.get_aX('ap',5,Q2)
            hq_pdf = q[5] 

        xhq = 1/(1+4*eps)

        if x>xhq:  return 0

        FX     = np.zeros(self.mellin.N.size,dtype=complex)
        CQ     = self.C2Q[1]
        AHGLOG = self.AHG1LOG[1]

        if stf=='F2':   #--F2/x
            FX+= eq2*(xhq/x)**N*hq_pdf/xhq 
            FX+= eq2*a*(xhq/x)**N*( CQ*hq_pdf - AHGLOG*np.log(Q2/mhq2)*g)/xhq
        elif stf=='FL': #--FL 
            FX+= eq2*a* x**(-N)*self.get_HQqFL(x,xhq,Q2,N,hq)*hq_pdf

        return FX

    def _get_HQFXN_(self,x,Q2,stf,hq): 
        """
        F2L c,b for nucleon. Twist 2 massless 
        """
        if hq=='c':
            eps = self.mc2/Q2
            eq2 = self.get_aX('ap',4,Q2)
        if hq=='b': 
            eps = self.mb2/Q2
            eq2 = self.get_aX('ap',5,Q2)

        #--a=alphaS/4pi 
        a = conf['alphaS'].get_a(Q2)*self.order  

        xhq = 1./(1.+4*eps)
        N = self.mellin.N
 

        if stf=='F2':
            F=0.
            if eps>=1:   # below threshold
                F+= x*eq2*a*2*self.get_HQgFX(x,Q2,hq,'F2') 

            else:        # above threshold
                FX= x*self.get_T2CFXHQ(x,Q2,'F2',hq) 
                F+=self.mellin.invert(1,FX) 
                F+= x*eq2*a*2*self.get_HQgFX(x,Q2,hq,'F2') 
            return F

        elif stf=='FL':
            F=0
            if eps>=1:   # below threshold
                F+= x*eq2*a*2*self.get_HQgFX(x,Q2,hq,'FL') 
            else:        # above threshold
                FX= x*self.get_T2CFXHQ(x,Q2,'FL',hq) 
                F+=self.mellin.invert(1,FX) 
                F+= x*eq2*a*2*self.get_HQgFX(x,Q2,hq,'FL') 
            return F

    def _get_HQFXN(self,x,Q2,stf,hq,tmc=False):

        if tmc==False:
            return self._get_HQFXN_(x,Q2,stf,hq)

        #--works, but is pretty slow since integration is done point by point
        elif tmc=='GP':

            mu2=self.M2/Q2
            xx=x*x
            rho=(1+4*mu2*xx)**0.5
            xi=2*x/(1+rho)
            
            z, w = np.polynomial.legendre.leggauss(3)
            h2int = lambda u: self._get_HQFXN_(u,Q2,'F2',hq)/u**2
            g2int = lambda u: self._get_HQFXN_(u,Q2,'F2',hq)*(u-xi)/u**2
            _u    = 0.5*((1-xi)*z + xi + 1)
            jac   = 0.5*(1-xi)
            h2, g2 = 0, 0
            for i in range(len(z)):
                h2 += w[i]*h2int(_u[i])*jac
                g2 += w[i]*g2int(_u[i])*jac

            if stf=='F2':
                F2=self._get_HQFXN_(xi,Q2,'F2',hq)  
                C1= (1+rho)**2/(4*rho**3)
                C2= 3*x*(rho**2-1)/(2*rho**4)
                C3= (rho**2-1)/(2*x*rho)
                return C1*F2 + C2*(h2 + C3*g2)

            elif stf=='FL':
                F2=self._get_HQFXN_(xi,Q2,'F2',hq)  
                FL=self._get_HQFXN_(xi,Q2,'FL',hq)  
                C1= (1+rho)**2/(4*rho)
                C2= x*(rho**2-1)/(rho**2)
                C3= (rho**2-1)/(2*x*rho)
                return C1*FL + C2*(h2 + C3*g2) 
 
        elif tmc=='AOT':

            mu2=self.M2/Q2
            xx=x*x
            rho=(1+4*mu2*xx)**0.5
            xi=2*x/(1+rho)
  
            if stf=='F2':
                F2=self._get_HQFXN_(xi,Q2,'F2',hq)  
                return (1+rho)/(2*rho**2)*F2
  
            elif stf=='FL':
                F2=self._get_HQFXN_(xi,Q2,'F2',hq)  
                FL=self._get_HQFXN_(xi,Q2,'FL',hq)  
                return FL + (rho-1)/2*F2
  
    #--charge current 
  
    def _get_WXN(self,x,Q2,stf='W2+',nucleon='p',ht=False,evolve=True):
        if '+' in stf:  sign=1
        if '-' in stf:  sign=-1
        if   'W2' in stf: WX= x*self.get_T2CWX('W2',nucleon,Q2,sign,evolve)
        elif 'WL' in stf: WX= x*self.get_T2CWX('WL',nucleon,Q2,sign,evolve)
        elif 'W3' in stf: WX=   self.get_T2CWX('W3',nucleon,Q2,sign,evolve)
        WX=self.mellin.invert(x,WX)
        if ht:  WX+=0
        return WX

    #--setups for nuclear data

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
        self.ymaxD  = mD/mN
        self.ymaxH  = mH/mN
        self.ymaxT  = mT/mN
        if 'idis grid' in conf and conf['idis grid'] == 'prediction': npts = 200
        else:                                                         npts = 100
        self.gX,self.gW=np.polynomial.legendre.leggauss(npts)

    #--nuclear structure functions
        
    def _get_FXA(self,X,Q2,stf='F2',nucleus='d',nuc=False,offpdf=False):


        if   nucleus=='d': p, n = 1,1
        elif nucleus=='h': p, n = 2,1
        elif nucleus=='t': p, n = 1,2

        #--no nuclear smearing for charged current structure functions
        if stf[0] == 'W': nuc=False

        #--no off-shell corrections except for F2
        if stf in ['F2g','F2gZ','F2Z']: offpdf = False

        if nuc==False:
            FXp   = self.get_stf(X,Q2,stf=stf,tar='p')
            FXn   = self.get_stf(X,Q2,stf=stf,tar='n')
            return (p*FXp+n*FXn)/(p+n)

        XM , gXM = np.meshgrid(X , self.gX)
        Q2M, gWM = np.meshgrid(Q2, self.gW)
        a = XM

        #--retrieve smearing class and upper limit of integration 
        if   nucleus=='d': smf, b = self.dsmf,self.ymaxD #--deuterium
        elif nucleus=='h': smf, b = self.hsmf,self.ymaxH #--helium
        elif nucleus=='t': smf, b = self.hsmf,self.ymaxT #--tritium

        switch = False
        #--tritium takes helium and switches p <--> n
        if nucleus=='t': switch = True

        YM   = 0.5*(b-a)*gXM+0.5*(a+b)
        JM   = 0.5*(b-a) 
        XM_YM = XM/YM
       
        if stf in ['F2','F2g','F2gZ','F2Z']:
            if p==n:
                fon22p = smf.get_fXX2('f22','onshell',XM,Q2M,YM)
                fon22n = fon22p
            else:
                fon22p = smf.get_fXX2('f22p','onshell',XM,Q2M,YM)
                fon22n = smf.get_fXX2('f22n','onshell',XM,Q2M,YM)
            if switch: 
                fon22p, fon22n = fon22n[:], fon22p[:]
            F2p   = self.get_stf(XM_YM,Q2M,stf=stf,tar='p')
            F2n   = self.get_stf(XM_YM,Q2M,stf=stf,tar='n')
            integ = p*fon22p*F2p + n*fon22n*F2n

            if offpdf:
                if p==n:
                    fof22p = smf.get_fXX2('f22','offshell',XM,Q2M,YM)
                    fof22n = fof22p
                else:
                    fof22p = smf.get_fXX2('f22p','offshell',XM,Q2M,YM)
                    fof22n = smf.get_fXX2('f22n','offshell',XM,Q2M,YM)
                if switch: 
                    fof22p, fof22n = fof22n[:], fof22p[:]
                F2poff = self.get_stf(XM_YM,Q2M,stf=stf,tar='p',off=True,nucleus=nucleus)
                F2noff = self.get_stf(XM_YM,Q2M,stf=stf,tar='n',off=True,nucleus=nucleus)
                integ+= p*fof22p*F2poff + n*fof22n*F2noff

        elif stf in ['FL','FLg','FLgZ','FLZ']:
            if p==n:
                fonL2p = smf.get_fXX2('fL2','onshell',XM,Q2M,YM)
                fonLLp = smf.get_fXX2('fLL','onshell',XM,Q2M,YM)
                fonL2n = fonL2p
                fonLLn = fonLLp
            else:
                fonL2p = smf.get_fXX2('fL2p','onshell',XM,Q2M,YM)
                fonL2n = smf.get_fXX2('fL2n','onshell',XM,Q2M,YM)
                fonLLp = smf.get_fXX2('fLLp','onshell',XM,Q2M,YM)
                fonLLn = smf.get_fXX2('fLLn','onshell',XM,Q2M,YM)
            if switch: 
                fonL2p, fonL2n = fonL2n[:], fonL2p[:]
                fonLLp, fonLLn = fonLLn[:], fonLLp[:]
            F2p   = self.get_stf(XM_YM,Q2M,stf='F2'+stf[2:],tar='p')
            F2n   = self.get_stf(XM_YM,Q2M,stf='F2'+stf[2:],tar='n')
            FLp   = self.get_stf(XM_YM,Q2M,stf=stf,tar='p')
            FLn   = self.get_stf(XM_YM,Q2M,stf=stf,tar='n')
            integ = p*(fonLLp*FLp + fonL2p*F2p) + n*(fonLLn*FLn + fonL2n*F2n)

        elif stf in ['F3','F3g','F3gZ','F3Z']:
            if p==n:
                fon33p = smf.get_fXX2('f33','onshell',XM,Q2M,YM)
                fon33n = fon33p
            else:
                fon33p = smf.get_fXX2('f33p','onshell',XM,Q2M,YM)
                fon33n = smf.get_fXX2('f33n','onshell',XM,Q2M,YM)
            if switch: 
                fon33p, fon33n = fon33n[:], fon33p[:]
            F3p   = self.get_stf(XM_YM,Q2M,stf=stf,tar='p')
            F3n   = self.get_stf(XM_YM,Q2M,stf=stf,tar='n')
            integ = p*fon33p*F3p + n*fon33n*F3n

        return np.einsum('ij,ij,ij->j',gWM,JM,integ)/(p+n)

    #--interpolation

    def get_grid(self):
       
        Q2min = 1    #--GeV^2
        Q2max = 1e5

        ##################################################################
        #--default fitting grid (sparse and fast, only covers data points)
        ##################################################################
        grid = []
        Q2 = np.geomspace(Q2min,200,10)
        Q2 = np.append(Q2,np.geomspace(300,Q2max,6))


        xmax = 0.99
        for i in range(len(Q2)):
            if i==0:
                X = np.geomspace(1e-5,0.1,10)
                X = np.append(X,np.linspace(0.03,xmax,20))
            if i==1:
                X = np.geomspace(1e-5,0.1,10)
                X = np.append(X,np.linspace(0.03,xmax,25))
            if i==2:
                X = np.geomspace(1e-5,0.1,10)
                X = np.append(X,np.linspace(0.03,xmax,25))
            if i==3:
                X = np.geomspace(1e-5,0.1,10)
                X = np.append(X,np.linspace(0.03,xmax,20))
            if i==4:
                X = np.geomspace(6e-5,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,30))
            if i==5:
                X = np.geomspace(1e-4,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,30))
            if i==6:
                X = np.geomspace(2e-4,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,30))
            if i==7:
                X = np.geomspace(3e-4,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,30))
            if i==8:
                X = np.geomspace(8e-4,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,10))
            if i==9:
                X = np.geomspace(1e-3,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,10))
            if i==10:
                X = np.geomspace(2e-3,0.1,10)
                X = np.append(X,np.linspace(0.07,xmax,10))
            if i==11:
                X = np.geomspace(8e-3,0.1,8)
                X = np.append(X,np.linspace(0.07,xmax,10))
            if i==12:
                X = np.geomspace(3e-2,0.1,5)
                X = np.append(X,np.linspace(0.07,xmax,10))
            if i==13:
                X = np.geomspace(3e-2,0.1,5)
                X = np.append(X,np.linspace(0.07,xmax,10))
            if i==14:
                X = np.linspace(0.11,xmax,5)
            if i==15:
                X = np.linspace(0.11,xmax,3)
            for x in X:
                grid.append([x,Q2[i]])

        #--extra grids for tough region
        Q2 = 1.61
        X  = np.linspace(1e-2,0.11,10)
        for x in X: grid.append([x,Q2])

        Q2 = 2.2
        X  = np.linspace(1e-2,0.11,10)
        for x in X: grid.append([x,Q2])

        #--extra grids for heavy quark
        if self.hq:
            Q2 = 1.61
            X  = 10**np.linspace(-5,-2,10)
            for x in X: grid.append([x,Q2])

            Q2 = 2
            X  = 10**np.linspace(-5,-2,10)
            for x in X: grid.append([x,Q2])

            Q2 = 5.0
            X  = 10**np.linspace(-5,-2,10)
            for x in X: grid.append([x,Q2])

        #--extra points for low Q2 < 1.69 pidis data
        Q2 = 0.9
        X = np.geomspace(2e-3,0.1,10)
        X = np.append(X,np.linspace(0.1,0.3,5))
        for x in X: grid.append([x,Q2])

        #--other grid options
        if 'idis grid' in conf:
            ##################################################################
            #--prediction grid (dense and slow, covers large region)
            ##################################################################
            if conf['idis grid']=='prediction':
                print('Using prediction grid...')
                Q2=[1.30000E+00,1.50159E+00,1.75516E+00,2.07810E+00\
                        ,2.49495E+00,3.04086E+00,3.76715E+00,4.50000E+00\
                        ,4.75000E+00,6.23113E+00,8.37423E+00,1.15549E+01\
                        ,1.64076E+01,2.40380E+01,3.64361E+01,5.73145E+01\
                        ,9.38707E+01,1.60654E+02,2.88438E+02,5.45587E+02\
                        ,1.09231E+03,2.32646E+03,5.30043E+03,1.29956E+04\
                        ,3.45140E+04,1.00000E+05]

                X=[1.00000E-06,1.28121E-06,1.64152E-06,2.10317E-06\
                       ,2.69463E-06,3.45242E-06,4.42329E-06,5.66715E-06\
                       ,7.26076E-06,9.30241E-06,1.19180E-05,1.52689E-05\
                       ,1.95617E-05,2.50609E-05,3.21053E-05,4.11287E-05\
                       ,5.26863E-05,6.74889E-05,8.64459E-05,1.10720E-04\
                       ,1.41800E-04,1.81585E-04,2.32503E-04,2.97652E-04\
                       ,3.80981E-04,4.87518E-04,6.26039E-04,8.00452E-04\
                       ,1.02297E-03,1.30657E-03,1.66759E-03,2.12729E-03\
                       ,2.71054E-03,3.44865E-03,4.37927E-03,5.54908E-03\
                       ,7.01192E-03,8.83064E-03,1.10763E-02,1.38266E-02\
                       ,1.71641E-02,2.11717E-02,2.59364E-02,3.15062E-02\
                       ,3.79623E-02,4.53425E-02,5.36750E-02,6.29705E-02\
                       ,7.32221E-02,8.44039E-02,9.64793E-02,1.09332E-01\
                       ,1.23067E-01,1.37507E-01,1.52639E-01,1.68416E-01\
                       ,1.84794E-01,2.01731E-01,2.19016E-01,2.36948E-01\
                       ,2.55242E-01,2.73927E-01,2.92954E-01,3.12340E-01\
                       ,3.32036E-01,3.52019E-01,3.72282E-01,3.92772E-01\
                       ,4.13533E-01,4.34326E-01,4.55495E-01,4.76836E-01\
                       ,4.98342E-01,5.20006E-01,5.41818E-01,5.63773E-01\
                       ,5.85861E-01,6.08077E-01,6.30459E-01,6.52800E-01\
                       ,6.75387E-01,6.98063E-01,7.20830E-01,7.43683E-01\
                       ,7.66623E-01,7.89636E-01,8.12791E-01,8.35940E-01\
                       ,8.59175E-01,8.82485E-01,9.05866E-01,9.29311E-01\
                       ,9.52817E-01,9.76387E-01,1.00000E+00]
                for q2 in Q2:
                    for x in X:
                        grid.append([x,q2])
            
            ##################################################################
            #--load list of excel files to create grid
            ##################################################################
            elif 'xlsx' in conf['idis grid']:
                if 'overwrite' in conf['idis grid'] and conf['idis grid']['overwrite']:
                    print('Deleting original idis grid...')
                    grid = []
                i = 0
                if 'xlsx' in conf['idis grid']:
                    for fn in conf['idis grid']['xlsx']:
                        if i==0: print('Using custom idis grid %s...'%fn)
                        else:    print('Appending custom idis grid %s...'%fn)
                        data = pd.read_excel(fn)
                        data = data.to_dict(orient='list')
                        if 'X' in data:
                            X  = data['X'] 
                        elif 'W' in data:
                            Q2, W = np.array(data['Q2']), np.array(data['W'])
                            X  = Q2/(W**2 - conf['aux'].M2 + Q2)
                        elif 'W2' in data:
                            Q2, W2 = np.array(data['Q2']), np.array(data['W2'])
                            X  = Q2/(W2 - conf['aux'].M2 + Q2)
                        Q2 = data['Q2']
                        for j in range(len(X)):
                            grid.append([X[j],Q2[j]]) 
                        i+=1

            ##################################################################
            #--create custom grid by hand
            ##################################################################
            elif 'X' in conf['idis grid'] and 'Q2' in conf['idis grid']:
                print('Using custom idis grid...')
                X,Q2 = conf['idis grid']['X'], conf['idis grid']['Q2']
                for q2 in Q2:
                    for x in X:
                        grid.append([x,q2])

            ##################################################################
            #--exact grid (load tables and use points as grid)
            ##################################################################
            elif conf['idis grid']=='exact':
                print('Using exact grid...')
                if 'idis tabs' in conf:
                    tabs = conf['idis tabs']
                    grid = []
                    for idx in tabs:
                        X  = tabs[idx]['X']
                        Q2 = tabs[idx]['Q2']
                        for j in range(len(X)):
                            grid.append([X[j],Q2[j]]) 
                if 'pidis tabs' in conf:
                    tabs = conf['pidis tabs']
                    for idx in tabs:
                        X  = tabs[idx]['X']
                        Q2 = tabs[idx]['Q2']
                        for j in range(len(X)):
                            grid.append([X[j],Q2[j]]) 

        GX  = np.array(grid).T[0]
        GQ2 = np.array(grid).T[1] 
       
        return np.array(GX),np.array(GQ2) 

    def setup_interpolation(self):
        self.X,self.Q2 = self.get_grid()
        self.data={}
        self.data['X']   = self.X
        self.data['Q2']  = self.Q2
        self.data['LX']  = np.log(self.X)
        self.data['LQ2'] = np.log(self.Q2)

        for tar in ['p','n','d','h','t']:
            self.data[tar]={}

        #--look at observables and targets to calculate only necessary structure functions
        if 'idis tabs' in conf:
            tabs = conf['idis tabs']
            for tab in tabs:
                tar = tabs[tab]['target'][0]
                obs = tabs[tab]['obs'][0].strip()

                if   obs == 'F2d/F2p':
                    for _ in ['p','n','d']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'F2n/F2d':
                    for _ in ['p','n','d']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'F2n/F2p':
                    for _ in ['p','n']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'F2d/F2h':
                    for _ in ['p','n','d','h']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'F2h/F2d':
                    for _ in ['p','n','d','h']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'F2t/F2d':
                    for _ in ['p','n','d','t']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'F2h/F2t':
                    for _ in ['p','n','h','t']: self.data[_]['F2']   = np.zeros(self.X.size)

                elif obs == 'sig_r'  :
                    current = tabs[tab]['current'][0].strip()

                    if current=='NC':
                        for _ in ['F2','FL','F3']: self.data[tar][_] = np.zeros(self.X.size)
                        if tar == 'd':
                            for _ in ['F2','FL','F3']: self.data['p'][_] = np.zeros(self.X.size)
                            for _ in ['F2','FL','F3']: self.data['n'][_] = np.zeros(self.X.size)

                    if current=='CC':
                        if '_plus' in tabs[tab]['lepton beam'][0]:
                            for _ in ['W2+','WL+','W3+']: self.data[tar][_] = np.zeros(self.X.size)
                            if tar == 'd':
                                for _ in ['W2+','WL+','W3+']: self.data['p'][_] = np.zeros(self.X.size)
                                for _ in ['W2+','WL+','W3+']: self.data['n'][_] = np.zeros(self.X.size)
                        if '_minus' in tabs[tab]['lepton beam'][0]:
                            for _ in ['W2-','WL-','W3-']: self.data[tar][_] = np.zeros(self.X.size)
                            if tar == 'd':
                                for _ in ['W2-','WL-','W3-']: self.data['p'][_] = np.zeros(self.X.size)
                                for _ in ['W2-','WL-','W3-']: self.data['n'][_] = np.zeros(self.X.size)

                elif obs == 'sigcc_r':
                    for _ in ['F2c','FLc','F2b','FLb']: self.data[tar][_] = np.zeros(self.X.size)

                elif obs == 'sigd/sigp':
                    for tar in ['p','n','d']:
                        for _ in ['F2','FL','F3']: self.data[tar][_] = np.zeros(self.X.size)

                elif obs == 'sigh/sigd':
                    for tar in ['p','n','d','h']:
                        for _ in ['F2','FL','F3']: self.data[tar][_] = np.zeros(self.X.size)

                elif obs == 'sigh/sigt':
                    for tar in ['p','n','t','h']:
                        for _ in ['F2','FL','F3']: self.data[tar][_] = np.zeros(self.X.size)

                elif obs == 'A_PV_e':
                    for _ in ['F2gZ','FLgZ','F3gZ','F2g','FLg']: self.data[tar][_] = np.zeros(self.X.size)
                    if tar == 'd' or tar == 'h':
                        for _ in ['F2gZ','FLgZ','F3gZ','F2g','FLg']: self.data['p'][_] = np.zeros(self.X.size)
                        for _ in ['F2gZ','FLgZ','F3gZ','F2g','FLg']: self.data['n'][_] = np.zeros(self.X.size)

                elif tar == 'd' or tar == 'h':
                    for _ in ['p','n',tar]: self.data[_][obs] = np.zeros(self.X.size)
 
                else:
                    self.data[tar][obs]    = np.zeros(self.X.size)

        #--if tables not available, calculate basic observables
        else:
            self.data['p']['F2']  = np.zeros(self.X.size)
            self.data['p']['FL']  = np.zeros(self.X.size)
            self.data['p']['F3']  = np.zeros(self.X.size)
            self.data['n']['F2']  = np.zeros(self.X.size)
            self.data['d']['F2']  = np.zeros(self.X.size)
            self.data['p']['F2c'] = np.zeros(self.X.size)
            self.data['p']['FLc'] = np.zeros(self.X.size)
            self.data['p']['F2b'] = np.zeros(self.X.size)
            self.data['p']['FLb'] = np.zeros(self.X.size)


        if self.offpdf:
            self.data['p']['F2off d'] = np.zeros(self.X.size)
            self.data['n']['F2off d'] = np.zeros(self.X.size)
            self.data['p']['F2off h'] = np.zeros(self.X.size)
            self.data['n']['F2off h'] = np.zeros(self.X.size)
            self.data['p']['F2off t'] = np.zeros(self.X.size)
            self.data['n']['F2off t'] = np.zeros(self.X.size)

    def _update(self):

        #--in case mc, mb is fitted we need to update
        self.mc2  = conf['aux'].mc**2
        self.mb2  = conf['aux'].mb**2

        n=len(self.X)
        #--nucleons
        for tar in ['p','n']:

            for stf in self.data[tar]:
  
                #--neutral current
                if stf in ['F2','FL','F3','F2g','FLg','F3g','F2gZ','FLgZ','F3gZ','F2Z','FLZ','F3Z']:

                    #--separate channels
                    if stf in   ['F2','FL','F3']:       channel='all'
                    elif stf in ['F2g','FLg','F3g']:    channel='gg'
                    elif stf in ['F2gZ','FLgZ','F3gZ']: channel='gZ'
                    elif stf in ['F2Z','FLZ','F3Z']:    channel='ZZ'

                    self.data[tar][stf]=np.array([self._get_FXN(self.X[i]
                                                               ,self.Q2[i]
                                                               ,stf=stf[:2]
                                                               ,nucleon=tar
                                                               ,hq=self.hq
                                                               ,tmc=self.tmc
                                                               ,ht=self.ht
                                                               ,channel=channel) for i in range(n)])

                #--charged current
                elif stf in ['W2+','WL+','W3+','W2-','WL-','W3-']:
                    self.data[tar][stf]=np.array([self._get_WXN(self.X[i]
                                                               ,self.Q2[i]
                                                               ,stf=stf
                                                               ,nucleon=tar
                                                               ,ht=self.ht)  for i in range(n)])

                #--heavy quarks
                elif stf in ['F2c','FLc','F2b','FLb']:

                    self.data[tar][stf]=np.array([self._get_HQFXN(self.X[i]
                                                                 ,self.Q2[i]
                                                                 ,stf =stf[:2]
                                                                 ,hq  =stf[-1]
                                                                 ,tmc =self.tmc
                                                                 ) for i in range(n)])

                #--offshell structure functions
                elif stf in ['F2off d', 'F2off h', 'F2off t']:
             
                    self.data[tar][stf]=np.array([self._get_FXN(self.X[i]
                                                               ,self.Q2[i]
                                                               ,stf=stf[:2]
                                                               ,nucleon=tar
                                                               ,tmc=self.tmc
                                                               ,offpdf = True
                                                               ,nucleus = stf[-1]
                                                               ,ht=False) for i in range(n)])

        #--composite targets
        for tar in ['d','h','t']:

            for stf in self.data[tar]:
                self.data[tar][stf]=self._get_FXA(self.X
                                                ,self.Q2
                                                ,stf=stf
                                                ,nucleus=tar
                                                ,offpdf=self.offpdf
                                                ,nuc=self.nuc)

    #--master funtion to be called from residuals

    def get_stf(self,X,Q2,stf='F2',tar='p',off=False,nucleus='d'):

        X0   = self.data['X']
        LX0  = self.data['LX']
        LQ20 = self.data['LQ2']

        if self.hq==False:

            if off==False:
                STF  = self.data[tar][stf]
                return griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')
            else:
                STF  = self.data[tar][stf+'off'+' '+nucleus]
                return griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')

        else:

            #--HQ are the same for p and n, so we only compute for p and give the same for n
            if stf=='F2' or stf=='FL':
                STF  = self.data[tar][stf]
                FX=griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')
                STF  = self.data['p'][stf+'c']
                FX+=griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')
                STF  = self.data['p'][stf+'b']
                FX+=griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')
                return FX

            elif stf=='F3':
                STF  = self.data[tar][stf]
                return griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')
  
            elif stf.startswith('W'):
                STF  = self.data[tar][stf]
                return griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')

            elif stf.endswith('c') or stf.endswith('b'):

                STF  = self.data['p'][stf]
                return griddata((LX0,LQ20),STF,(np.log(X) ,np.log(Q2)), fill_value=0, method='cubic')




