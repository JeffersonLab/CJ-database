from mpmath import fp
import numpy as np
from tools.config import conf

class AUX:

    def setup(self):

        ############################################
        # abbreviations
        ############################################
        zeta2=fp.zeta(2)
        zeta3=fp.zeta(3)
        N   = self.mellin.N 
        NP1 = N + 1
        NP2 = N + 2
        NM1 = N - 1
        NS=N*N
        psi=lambda i,N: fp.psi(i,complex(N.real,N.imag))
        S1 = np.array([fp.euler  + psi(0,n+1) for n in N])
        S2 = np.array([zeta2- psi(1,n+1) for n in N])
        S1S=S1**2

        ############################################
        # hard coeffs
        ############################################
        orders=2

        self.C2Q = np.zeros((orders,N.size),dtype=complex)
        self.C2G = np.zeros((orders,N.size),dtype=complex)
        self.CLQ = np.zeros((orders,N.size),dtype=complex)
        self.CLG = np.zeros((orders,N.size),dtype=complex)
        self.C3Q = np.zeros((orders,N.size),dtype=complex)
        self.PC1Q = np.zeros((orders,N.size),dtype=complex)
        self.PC1G = np.zeros((orders,N.size),dtype=complex)
        self.AHG1LOG = np.zeros((orders,N.size),dtype=complex)

        self.C2Q[0]  = np.ones(N.size)  
        self.C3Q[0]  = np.ones(N.size)  
        self.PC1Q[0] = np.ones(N.size)  

        CF=4.0/3.0
        TR=0.5
        # Nucl. Phys. B192 (1981) 417
        # Ellis book
        self.C2Q[1]  = CF*(2.0*S1S - 2.0*S2 + 3.0*S1 - 2.0*S1/N/NP1 + 3.0/N + 4.0/NP1 + 2.0/NS - 9.0)
        self.C2G[1]  = -2*TR*(S1*(NS + N + 2.0)/N/NP1/NP2 + 1.0/N - 1.0/NS - 6.0/NP1 + 6.0/NP2)  

        self.CLQ[1]  = CF*4/NP1  
        self.CLG[1]  = 8.0*TR/NP1/NP2

        self.C3Q[1]  = CF*(2.0*S1S - 2.0*S2 + 3.0*S1 - 2.0*S1/N/NP1 + 1.0/N + 2.0/NP1 + 2.0/NS - 9.0)
        self.AHG1LOG[1] = 2.*(1./N - 2./NP1 + 2./NP2) #NL0 glue coeffcient for F2^{hq}

        # Nucl. Phys. B192 (1981) 417??? need check
        self.PC1Q[1] = CF*(2*S1S-2*S2+2*S1/NP1-2*S1/N+3*S1-2/N/NP1+3/N+2/NS-9)
        self.PC1G[1] = 2*NM1*(1-N-N*S1)/NS/NP1
      
    def get_aX(self,aX,i,Q2,channel='all'):

        ############################################
        # couplings  (see page 18 of Pedro's thesis)
        # {q} = {u,d,s,c,b} 
        ############################################

        alpha = conf['eweak'].get_alpha(Q2)
        sin2w = conf['eweak'].get_sin2w(Q2)
        KQ2 = conf['aux'].GF*conf['aux'].mZ2/(2*2**0.5 * np.pi * alpha) \
            * Q2/(Q2+conf['aux'].mZ2)

        Ve=-0.5+2*sin2w
        Ae=-0.5
        Ve2pAe2=Ve*Ve+Ae*Ae

        eu=2.0/3.0
        Vu=0.5-2*eu*sin2w
        Au=0.5
        Vu2pAu2=Vu*Vu+Au*Au

        ed=-1.0/3.0
        Vd=-0.5-2*ed*sin2w
        Ad=-0.5
        Vd2pAd2=Vd*Vd+Ad*Ad

        #--for individual channels, conventions are taken from PDG
        #--pdg.lbl.gov/2019/reviews/rpp2019-rev-structure-functions.pdf
        #--see Eq. (18.18)
        if aX=='ap':
            if i==1 or i==4:
                if channel=='all':  return eu*eu - 2*eu*Ve*Vu*KQ2 + Ve2pAe2*Vu2pAu2*KQ2**2
                elif channel=='gg': return eu*eu
                elif channel=='ZZ': return Vu2pAu2
                elif channel=='gZ': return 2*eu*Vu
            elif i==2 or i==3 or i==5:
                if channel=='all':  return ed*ed - 2*ed*Ve*Vd*KQ2 + Ve2pAe2*Vd2pAd2*KQ2**2
                elif channel=='gg': return ed*ed
                elif channel=='ZZ': return Vd2pAd2
                elif channel=='gZ': return 2*ed*Vd
        if aX=='am':
            if i==1 or i==4:
                if channel=='all':  return -2*eu*Ae*Au*KQ2 + 4*Ve*Ae*Vu*Au*KQ2**2
                elif channel=='gg': return 0.0 
                elif channel=='ZZ': return 2*Vu*Au
                elif channel=='gZ': return 2*eu*Au
            elif i==2 or i==3 or i==5:
                if channel=='all':  return -2*ed*Ae*Ad*KQ2 + 4*Ve*Ae*Vd*Ad*KQ2**2
                elif channel=='gg': return 0.0 
                elif channel=='ZZ': return 2*Vd*Ad
                elif channel=='gZ': return 2*ed*Ad





