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
        N   = self.mell.N 
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
        self.C = {}
        self.C['2'] = {}
        self.C['L'] = {}
        self.C['3'] = {}

        orders = 2
        self.PC1Q = np.zeros((orders,N.size),dtype=complex)
        self.PC1G = np.zeros((orders,N.size),dtype=complex)
        #self.AHG1LOG = np.zeros((orders,N.size),dtype=complex)

        CF=4.0/3.0
        TR=0.5
        # Nucl. Phys. B192 (1981) 417
        # Ellis book
        self.C['2']['Q']  = CF*(2.0*S1S - 2.0*S2 + 3.0*S1 - 2.0*S1/N/NP1 + 3.0/N + 4.0/NP1 + 2.0/NS - 9.0)
        self.C['2']['G']  = -2*TR*(S1*(NS + N + 2.0)/N/NP1/NP2 + 1.0/N - 1.0/NS - 6.0/NP1 + 6.0/NP2)  

        self.C['L']['Q']  = CF*4/NP1  
        self.C['L']['G']  = 8.0*TR/NP1/NP2

        self.C['3']['Q']  = CF*(2.0*S1S - 2.0*S2 + 3.0*S1 - 2.0*S1/N/NP1 + 1.0/N + 2.0/NP1 + 2.0/NS - 9.0)
        self.C['3']['G']  = np.zeros(N.size,dtype=complex)       
 
        #self.AHG1LOG[1] = 2.*(1./N - 2./NP1 + 2./NP2) #NL0 glue coeffcient for F2^{hq}

        self.PC = {}
        self.PC['1'] = {}
        # Nucl. Phys. B192 (1981) 417??? need check
        self.PC['1']['Q'] = CF*(2*S1S-2*S2+2*S1/NP1-2*S1/N+3*S1-2/N/NP1+3/N+2/NS-9)
        self.PC['1']['G'] = 2*NM1*(1-N-N*S1)/NS/NP1
     
    def get_aX(self,_aX,Q2,channel='all'):
 
        ############################################
        # couplings  (see page 18 of Pedro's thesis)
        # {q} = {u,d,s,c,b} 
        ############################################

        zeros = np.zeros(Q2.shape[0])
        ones  = np.ones (Q2.shape[0])

        alpha = np.array([conf['eweak'].get_alpha(q2) for q2 in Q2])
        sin2w = np.array([conf['eweak'].get_sin2w(q2) for q2 in Q2])
        KQ2   = conf['aux'].GF*conf['aux'].mZ2/(2*2**0.5 * np.pi * alpha) * Q2/(Q2+conf['aux'].mZ2)

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
        aX = {}

        if _aX == 'aXp':
            if channel=='all':
                aX[1] = eu*eu - 2*eu*Ve*Vu*KQ2 + Ve2pAe2*Vu2pAu2*KQ2**2 
                aX[2] = ed*ed - 2*ed*Ve*Vd*KQ2 + Ve2pAe2*Vd2pAd2*KQ2**2 
                aX[3] = ed*ed - 2*ed*Ve*Vd*KQ2 + Ve2pAe2*Vd2pAd2*KQ2**2 
                aX[4] = eu*eu - 2*eu*Ve*Vu*KQ2 + Ve2pAe2*Vu2pAu2*KQ2**2 
                aX[5] = ed*ed - 2*ed*Ve*Vd*KQ2 + Ve2pAe2*Vd2pAd2*KQ2**2 
            if channel=='gg':
                aX[1] = eu*eu*ones
                aX[2] = ed*ed*ones
                aX[3] = ed*ed*ones
                aX[4] = eu*eu*ones
                aX[5] = ed*ed*ones
            if channel=='gZ':
                aX[1] = 2*eu*Vu
                aX[2] = 2*ed*Vd
                aX[3] = 2*ed*Vd
                aX[4] = 2*eu*Vu
                aX[5] = 2*ed*Vd
            if channel=='ZZ':
                aX[1] = Vu2pAu2 
                aX[2] = Vd2pAd2 
                aX[3] = Vd2pAd2 
                aX[4] = Vu2pAu2 
                aX[5] = Vd2pAd2 
            for i in [1,2,3,4,5]: aX[-i] = aX[i]

        if _aX == 'aXm':
            if channel=='all':
                aX[1] = -2*eu*Ae*Au*KQ2 + 4*Ve*Ae*Vu*Au*KQ2**2 
                aX[2] = -2*ed*Ae*Ad*KQ2 + 4*Ve*Ae*Vd*Ad*KQ2**2 
                aX[3] = -2*ed*Ae*Ad*KQ2 + 4*Ve*Ae*Vd*Ad*KQ2**2 
                aX[4] = -2*eu*Ae*Au*KQ2 + 4*Ve*Ae*Vu*Au*KQ2**2 
                aX[5] = -2*ed*Ae*Ad*KQ2 + 4*Ve*Ae*Vd*Ad*KQ2**2 
            if channel=='gg':
                aX[1] = zeros 
                aX[2] = zeros 
                aX[3] = zeros 
                aX[4] = zeros 
                aX[5] = zeros 
            if channel=='gZ':
                aX[1] = 2*eu*Au*ones
                aX[2] = 2*ed*Ad*ones
                aX[3] = 2*ed*Ad*ones
                aX[4] = 2*eu*Au*ones
                aX[5] = 2*ed*Ad*ones
            if channel=='ZZ':
                aX[1] = 2*Vu*Au 
                aX[2] = 2*Vd*Ad 
                aX[3] = 2*Vd*Ad 
                aX[4] = 2*Vu*Au 
                aX[5] = 2*Vd*Ad 
            for i in [1,2,3,4,5]: aX[-i] = -aX[i]

        if _aX == 'aXpWp':
            aX[1]  = np.zeros(Q2.size)
            aX[2]  = np.ones(Q2.size)
            aX[3]  = np.ones(Q2.size)
            aX[4]  = np.zeros(Q2.size)
            aX[5]  = np.ones(Q2.size)
            aX[-1] = np.ones(Q2.size)
            aX[-2] = np.zeros(Q2.size)
            aX[-3] = np.zeros(Q2.size)
            aX[-4] = np.ones(Q2.size)
            aX[-5] = np.zeros(Q2.size)

        if _aX == 'aXmWp':
            aX[1]  = np.zeros(Q2.size)
            aX[2]  = np.ones(Q2.size)
            aX[3]  = np.ones(Q2.size)
            aX[4]  = np.zeros(Q2.size)
            aX[5]  = np.ones(Q2.size)
            aX[-1] = -np.ones(Q2.size)
            aX[-2] = np.zeros(Q2.size)
            aX[-3] = np.zeros(Q2.size)
            aX[-4] = -np.ones(Q2.size)
            aX[-5] = np.zeros(Q2.size)

        if _aX == 'aXpWm':
            aX[1]  = np.ones(Q2.size) 
            aX[2]  = np.zeros(Q2.size)
            aX[3]  = np.zeros(Q2.size)
            aX[4]  = np.ones(Q2.size)
            aX[5]  = np.zeros(Q2.size)
            aX[-1] = np.zeros(Q2.size)
            aX[-2] = np.ones(Q2.size)
            aX[-3] = np.ones(Q2.size)
            aX[-4] = np.zeros(Q2.size)
            aX[-5] = np.ones(Q2.size)

        if _aX == 'aXmWm':
            aX[1]  = np.ones(Q2.size) 
            aX[2]  = np.zeros(Q2.size)
            aX[3]  = np.zeros(Q2.size)
            aX[4]  = np.ones(Q2.size) 
            aX[5]  = np.zeros(Q2.size)
            aX[-1] = np.zeros(Q2.size)
            aX[-2] = -np.ones(Q2.size) 
            aX[-3] = -np.ones(Q2.size) 
            aX[-4] = np.zeros(Q2.size)
            aX[-5] = -np.ones(Q2.size) 


        return aX
 
 

