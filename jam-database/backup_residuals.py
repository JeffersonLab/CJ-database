import sys
import numpy as np
from obslib.idis.reader import READER
from obslib.idis.theory import THEORY
from obslib.idis.aux    import AUX
from tools.residuals    import _RESIDUALS
from tools.config       import conf
import pandas as pd

class RESIDUALS(_RESIDUALS):

    def __init__(self): 
        self.reaction='idis'
        self.tabs=conf['idis tabs']
        self.thy=conf['idis']
        self.setup()

        self.tmc = False
        self.hq  = False
        self.ht  = False
        self.Nf  = None
        self.nuc = False
        self.offshell = False
        if 'hq'   in conf: self.hq  = conf['hq']
        if 'tmc'  in conf: self.tmc = conf['tmc']
        if 'ht'   in conf: self.ht  = conf['ht']
        if 'nuc'  in conf: self.nuc = conf['nuc']
        if 'offshell' in conf: self.offshell = conf['offshell']
        if self.hq : self.Nf  = 3 
 
    def get_F2_ratio(self,x,Q2,num,den):
        F2num=self.thy.get_stf(x,Q2,stf='F2',tar=num)
        F2den=self.thy.get_stf(x,Q2,stf='F2',tar=den)
        rat  = F2num/F2den
        return rat

    def get_sig_ratio(self,X,Q2,Y,num,den,current,lbeam):
        signum=self.get_sigma_r(X,Q2,Y,num,current,lbeam)
        sigden=self.get_sigma_r(X,Q2,Y,den,current,lbeam)
        rat = signum/sigden
        return rat
   
    def get_sigma_r(self,X,Q2,Y,tar,current,lbeam):

        if current=='NC': thy = self._get_sigma_r  (X,Q2,Y,tar,lbeam)
        if current=='CC': thy = self._get_sigma_rcc(X,Q2,Y,tar,lbeam)
  
        return thy
 
    def _get_sigma_r(self,X,Q2,Y,tar,lbeam):
  
        F2  = self.thy.get_stf(X,Q2,stf='F2',tar=tar)
        FL  = self.thy.get_stf(X,Q2,stf='FL',tar=tar)
        F3  = self.thy.get_stf(X,Q2,stf='F3',tar=tar)

        if   '_minus' in lbeam: sign = 1
        elif '_plus'  in lbeam: sign = -1
        elif 'none'   in lbeam: sign = 0
        else: 
            print('ERR in data set %d: lepton id is not identified'%(idx))
            sys.exit()
  
        YP=1+(1-Y)**2
        YM=1-(1-Y)**2

        thy=F2-Y**2/YP*FL+sign*YM/YP*X*F3
        return thy
  
    def _get_sigma_rcc(self,X,Q2,Y,tar,lbeam):
  
        if   '_plus'  in lbeam:
            sign = -1
            W2   = self.thy.get_stf(X,Q2,stf='W2+',tar=tar)
            WL   = self.thy.get_stf(X,Q2,stf='WL+',tar=tar)
            W3   = self.thy.get_stf(X,Q2,stf='W3+',tar=tar)

        elif '_minus' in lbeam:
            sign = 1
            W2   = self.thy.get_stf(X,Q2,stf='W2-',tar=tar)
            WL   = self.thy.get_stf(X,Q2,stf='WL-',tar=tar)
            W3   = self.thy.get_stf(X,Q2,stf='W3-',tar=tar)
  
        YP=1+(1-Y)**2
        YM=1-(1-Y)**2

        thy=YP/4*W2-Y**2/4*WL+sign*YM/4*X*W3
        return thy

    def get_sigma_dxdy(self,X,Q2,Y,tar,lbeam):
  
        F2 =self.thy.get_stf(X,Q2,stf='F2',tar=tar)
        FL =self.thy.get_stf(X,Q2,stf='FL',tar=tar)
        F3 =self.thy.get_stf(X,Q2,stf='F3',tar=tar)
        if    '_plus'  in lbeam: sign=-1
        elif  '_minus' in lbeam: sign=1
        else: 
            print('ERR in data set %d: lepton id is not identified'%(idx))
            sys.exit()
  
        M2 = conf['aux'].M2
        YP=1+(1-Y)**2
        YM=1-(1-Y)**2
  
        thy=2*np.pi*conf['aux'].alfa**2/X/Y/Q2
        thy*=(YP+2*X**2*Y**2*M2/Q2)*F2-Y**2*FL+sign*YM*X*F3
  
        return thy

    def get_sigma_dxdQ2(self,X,Q2,Y,tar,lbeam):
        thy=self._get_sigma_dxdy(X,Q2,Y,tar,lbeam)
        M2=conf['aux'].M2
        s = Q2/X/Y + M2
        return thy/X/(s-M2)

    def get_sigcc_r(self,X,Q2,Y,tar,current,lbeam):

        if current=='NC':

            F2  = self.thy.get_stf(X,Q2,stf='F2c',tar=tar)
            FL  = self.thy.get_stf(X,Q2,stf='FLc',tar=tar)
            F3  = 0#self.thy.get_stf(X,Q2,stf='F3c',tar=tar)
            thy = self._get_sigma_r(X,Q2,Y,tar,lbeam)

        return thy

    def get_A_PV_e(self,X,Q2,Y,tar):
        #--parity violating asymmetry with longitudinally polarized lepton
        #--derived from Eq. (11) in arxiv.org/pdf/hep-ph/9401264.pdf

        #--get structure functions
        F2gZ = self.thy.get_stf(X,Q2,stf='F2gZ',tar=tar)
        FLgZ = self.thy.get_stf(X,Q2,stf='FLgZ',tar=tar)
        F3gZ = self.thy.get_stf(X,Q2,stf='F3gZ',tar=tar)
        F2g  = self.thy.get_stf(X,Q2,stf='F2g', tar=tar)
        FLg  = self.thy.get_stf(X,Q2,stf='FLg', tar=tar)

        M2 = conf['aux'].M2

        rho2 = 1 + 4*x**2*M2/Q2

        YP = Y**2*(rho2+1)/2 - 2*Y +2
        YM = 1-(1-Y)**2

        GF    = conf['aux'].GF

        sin2w = np.array([conf['eweak'].get_sin2w(q2) for q2 in Q2])
        alpha = np.array([conf['eweak'].get_alpha(q2) for q2 in Q2])

        gA = -0.5
        gV = -0.5 + 2*sin2w

        C  = GF*Q2/(2*np.sqrt(2)*np.pi*alpha)

        #--NLO
        #num1 = gA*(YP*F2gZ - y**2*FLgZ)
        #num2 = gV*x*YM*F3gZ
        #num = num1 + num2
        #denom = YP*F2g - y**2*FLg
        #--corrected?
        F1gZ = (rho2*F2gZ - FLgZ)/(2*X)
        F1g  = (rho2*F2g  - FLg )/(2*X)
        
        A = C*(gA*F1gZ/F1g + gV*YM*F3gZ/(2*YP*F1g))

        return A 

    def get_theory(self):

        #--note: for nuclear targets (d, p, h, ...), observables need to be per nucleon (divided by # of protons + neutrons)

        STFs      = ['F2','FL','F3','F2g','FLg','F3g','F2gZ','FLgZ','FLgZ','F3gZ','F2Z','FLZ','F3Z']
        STFs.extend(['W2+','WL+','W3+','W2-','WL-','W3-'])
        F2ratios  = ['F2d/F2p','F2n/F2d','F2d/F2h','F2h/F2d','F2t/F2d','F2h/F2t','F2n/F2p']
        sigratios = ['sigd/sigp','sigh/sigd','sigh/sigt']
        
        for idx in self.tabs:

            X=self.tabs[idx]['X']
            Q2=self.tabs[idx]['Q2']
            target=self.tabs[idx]['target'][0]
            obs=self.tabs[idx]['obs'][0].strip()
            if 'current'     in self.tabs[idx]: current = self.tabs[idx]['current'][0].strip()
            if 'lepton beam' in self.tabs[idx]: lbeam   = self.tabs[idx]['lepton beam'][0]
            if 'Elab' in self.tabs[idx]:
                M2 = conf['aux'].M2
                ELab=self.tabs[idx]['Elab']
                s=M2+2*ELab*M2**0.5
                Y=(Q2/2/X)/((s-M2)/2)
            elif 'Y' in self.tabs[idx]:
                Y=self.tabs[idx]['Y']
            else:
                Y=None

            if   obs in STFs                   :  xsec=self.thy.get_stf    (X,Q2,obs,target)
            elif obs in F2ratios               :  xsec=self.get_F2_ratio   (X,Q2,obs[2],obs[-1])
            elif obs in sigratios              :  xsec=self.get_sig_ratio  (X,Q2,Y,obs[3],obs[-1],current,lbeam)
            elif obs=='sig_r'                  :  xsec=self.get_sigma_r    (X,Q2,Y,target        ,current,lbeam)
            elif obs=='dsig/dxdy'              :  xsec=self.get_sigma_dxdy (X,Q2,Y,target,lbeam)
            elif obs=='dsig/dxdQ2'             :  xsec=self.get_sigma_dxdQ2(X,Q2,Y,target,lbeam)
            #--heavy quark observable
            elif obs=='sigcc_r'                :  xsec=self.get_sigcc_r    (X,Q2,Y,target,current,lbeam)
            #--EIC observable
            elif obs=='A_PV_e'                 :  xsec=self.get_A_PV_e     (X,Q2,Y,target)
            else:
                msg='exp=%d obs=%s and target=%s not implemented'
                raise ValueError(msg%(idx,obs,target))  

            self.tabs[idx]['thy']=xsec
  
    def _gen_report(self,verb=1,level=1):
        """
        verb = 0: Do not print on screen. Only return list of strings
        verv = 1: print on screen the report
        level= 0: only the total chi2s
        level= 1: include point by point 
        """

        L=[]

        if len(self.tabs.keys())!=0:
            nuc      = conf['idis'].nuc
            if nuc:
                if 'dsmf type' in conf: dsmf = conf['dsmf type']
                else:                   dsmf = 'Paris'
                if 'hsmf type' in conf: hsmf = conf['hsmf type']
                else:                   hsmf = 'KPSV'
                nuc = 'True (deuterium:%s, helium:%s)'%(dsmf,hsmf)
            tmc      = conf['idis'].tmc
            offshell = conf['idis'].offpdf
            ht       = conf['idis'].ht
            if ht: ht_type = conf['idis'].ht_type
            else:  ht_type = False
            L.append('reaction: unpol DIS [nuc=%s, tmc=%s, ht=%s, offshell=%s]'%(nuc,tmc,ht_type,offshell))
            for f in conf['datasets']['idis']['filters']:
                L.append('filters: %s'%f)

            L.append('%7s %3s %10s %5s %10s %10s %10s %10s'%('idx','tar','col','npts','chi2','chi2/npts','rchi2','nchi2'))
            for k in self.tabs:
                if len(self.tabs[k])==0: continue 
                res=self.tabs[k]['residuals']

                rres=[]
                for c in conf['rparams']['idis'][k]:
                    rres.append(conf['rparams']['idis'][k][c]['value'])
                rres=np.array(rres)

                if k in conf['datasets']['idis']['norm'] and conf['datasets']['idis']['norm'][k]['fixed']==False:
                    norm=conf['datasets']['idis']['norm'][k]
                    nres=(norm['value']-1)/norm['dN']
                else:
                    nres=0

                
                chi2=np.sum(res**2)
                rchi2=np.sum(rres**2)
                nchi2=nres**2
                tar=self.tabs[k]['target'][0]
                col=self.tabs[k]['col'][0].split()[0]
                npts=res.size
                L.append('%7d %3s %10s %5d %10.2f %10.2f %10.2f %10.2f'%(k,tar,col,npts,chi2,chi2/npts,rchi2,nchi2))

            if level==1:
                L.append('-'*100)  
                for k in self.tabs:
                    if len(self.tabs[k]['value'])==0: continue 
                    if k in conf['datasets']['idis']['norm'] and conf['datasets']['idis']['norm'][k]['fixed']==False:
                        norm=conf['datasets']['idis']['norm'][k]
                        nres=(norm['value']-1)/norm['dN']
                        norm=norm['value']
                    else:
                        norm=1.0
                        nres=0
                    for i in range(len(self.tabs[k]['value'])):
                        x     = self.tabs[k]['X'][i]
                        Q2    = self.tabs[k]['Q2'][i]
                        res   = self.tabs[k]['residuals'][i]
                        thy   = self.tabs[k]['thy'][i]
                        exp   = self.tabs[k]['value'][i]
                        alpha = self.tabs[k]['alpha'][i]
                        rres  = self.tabs[k]['r-residuals'][i]
                        col   = self.tabs[k]['col'][i]
                        shift = self.tabs[k]['shift'][i]
                        tar   = self.tabs[k]['target'][i]
                        msg='col=%7s, tar=%5s, x=%10.3e, Q2=%10.3e, exp=%10.3e, alpha=%10.3e, thy=%10.3e, shift=%10.3e, chi2=%10.3e, res=%10.3e, norm=%10.3e, '
                        L.append(msg%(col,tar,x,Q2,exp,alpha,thy,shift,res**2,res,norm))

        if verb==0:
            return L
        elif verb==1:
            for l in L: print(l)
        return L

    def gen_report(self,verb=1,level=1):

        msg ='reaction: DIS <br />'
        if 'filters' in conf['datasets']['idis']:
            for f in conf['datasets']['idis']['filters']:
                msg+='filters: %s <br />'%f
        msg+='reaction: DIS <br />'


        data={_:[] for _ in ['idx','col','obs','tar','npts','chi2','chi2/npts','rchi2','nchi2']}
        for idx in self.tabs:
            if len(self.tabs[idx])==0: continue
            res=self.tabs[idx]['residuals']

            rres=[]
            for c in conf['rparams']['idis'][idx]:
                rres.append(conf['rparams']['idis'][idx][c]['value'])
            rres=np.array(rres)

            if idx in conf['datasets']['idis']['norm'] and conf['datasets']['idis']['norm'][idx]['fixed']==False:
                norm=conf['datasets']['idis']['norm'][idx]
                nres=(norm['value']-1)/norm['dN']
            else:
                nres=0

            chi2=np.sum(res**2)
            rchi2=np.sum(rres**2)
            nchi2=nres**2
            tar=self.tabs[idx]['target'][0]
            col=self.tabs[idx]['col'][0].split()[0]
            obs=self.tabs[idx]['obs'][0].strip()
            npts=res.size
            data['idx'].append('%7d'%idx)
            data['col'].append('%10s'%col)
            data['obs'].append('%3s'%obs)
            data['tar'].append('%10s'%tar)
            data['npts'].append('%10.2f'%npts)
            data['chi2'].append('%10.2f'%chi2)
            data['chi2/npts'].append('%10.2f'%(chi2/npts))
            data['rchi2'].append('%10.2f'%rchi2)
            data['nchi2'].append('%10.2f'%nchi2)

        data=pd.DataFrame(data)
        msg+=data.to_html(col_space=80,index=False)#,justify='center')
        return msg


class _RESIDUALS(_RESIDUALS):

    def __init__(self): 
        self.reaction='idis'
        self.tabs=conf['idis tabs']
        self.thy=conf['idis']
        self.setup()

        self.tmc = False
        self.hq  = False
        self.ht  = False
        self.Nf  = None
        self.nuc = False
        self.offshell = False
        if 'hq'   in conf: self.hq  = conf['hq']
        if 'tmc'  in conf: self.tmc = conf['tmc']
        if 'ht'   in conf: self.ht  = conf['ht']
        if 'nuc'  in conf: self.nuc = conf['nuc']
        if 'offshell' in conf: self.offshell = conf['offshell']
        if self.hq : self.Nf  = 3 
 
    def get_F2_ratio(self,x,Q2,idx,num,den):
        F2num=self.thy.get_stf(x,Q2,stf='F2',tar=num)
        F2den=self.thy.get_stf(x,Q2,stf='F2',tar=den)
        rat  = F2num/F2den
        return rat

    def get_sig_ratio(self,X,Q2,idx,num,den):
        signum=self.get_sigma_r(X,Q2,idx,num)
        sigden=self.get_sigma_r(X,Q2,idx,den)
        rat = signum/sigden
        return rat
   
    def get_sigma_r(self,X,Q2,idx,tar):

        current=self.tabs[idx]['current'][0].strip()
  
        if current=='NC':

            F2  = self.thy.get_stf(X,Q2,stf='F2',tar=tar)
            FL  = self.thy.get_stf(X,Q2,stf='FL',tar=tar)
            F3  = self.thy.get_stf(X,Q2,stf='F3',tar=tar)
            thy = self._get_sigma_r(idx,F2,FL,F3)
  
        if current=='CC':
  
            if   '_plus'  in self.tabs[idx]['lepton beam'][0]:
                sign = 1
                W2   = self.thy.get_stf(X,Q2,stf='W2+',tar=tar)
                WL   = self.thy.get_stf(X,Q2,stf='WL+',tar=tar)
                W3   = self.thy.get_stf(X,Q2,stf='W3+',tar=tar)

            elif '_minus' in self.tabs[idx]['lepton beam'][0]:
                sign = -1
                W2   = self.thy.get_stf(X,Q2,stf='W2-',tar=tar)
                WL   = self.thy.get_stf(X,Q2,stf='WL-',tar=tar)
                W3   = self.thy.get_stf(X,Q2,stf='W3-',tar=tar)
    
            thy = self._get_sigma_rcc(idx,W2,WL,W3)

        return thy
 
    def _get_sigma_r(self,idx,F2,FL,F3):
  
        if   '_minus' in self.tabs[idx]['lepton beam'][0]: sign = 1
        elif '_plus'  in self.tabs[idx]['lepton beam'][0]: sign = -1
        elif 'none'   in self.tabs[idx]['lepton beam'][0]: sign = 0
  
        else: 
            print('ERR in data set %d: lepton id is not identified'%(idx))
            sys.exit()
  
        x  = self.tabs[idx]['X']
        Q2 = self.tabs[idx]['Q2']
        M2 = conf['aux'].M2
  
        if 'Elab' in self.tabs[idx]:
            ELab=self.tabs[idx]['Elab']
            s=M2+2*ELab*M2**0.5
            y=(Q2/2/x)/((s-M2)/2)
        elif 'Y' in self.tabs[idx]:
            y=self.tabs[idx]['Y']

        YP=1+(1-y)**2
        YM=1-(1-y)**2

        thy=F2-y**2/YP*FL+sign*YM/YP*x*F3
        return thy
  
    def _get_sigma_rcc(self,idx,W2,WL,W3):
  
        if    '_plus'  in self.tabs[idx]['lepton beam'][0]: sign=-1
        elif  '_minus' in self.tabs[idx]['lepton beam'][0]: sign=1
        else: 
            print('ERR in data set %d: lepton id is not identified'%(idx))
            sys.exit()
  
        x=self.tabs[idx]['X']
        Q2=self.tabs[idx]['Q2']
        M2=conf['aux'].M2
  
        if 'Elab' in self.tabs[idx]:
            ELab=self.tabs[idx]['Elab']
            s=M2+2*ELab*M2**0.5
            y=(Q2/2/x)/((s-M2)/2)
        elif 'Y' in self.tabs[idx]:
            y=self.tabs[idx]['Y']

        YP=1+(1-y)**2
        YM=1-(1-y)**2

        thy=YP/4*W2-y**2/4*WL+sign*YM/4*x*W3
        return thy

    def get_sigma_dxdy(self,x,Q2,k,idx,tar):
        F2 =self.thy.get_stf(X,Q2,stf='F2',tar=tar)
        FL =self.thy.get_stf(X,Q2,stf='FL',tar=tar)
        F3 =self.thy.get_stf(X,Q2,stf='F3',tar=tar)
        return self._get_sigma_dxdy(k,idx,F2,FL,F3)
  
    def _get_sigma_dxdy(self,idx,F2,FL,F3):
  
        if    '_plus'  in self.tabs[idx]['lepton beam'][0]: sign=-1
        elif  '_minus' in self.tabs[idx]['lepton beam'][0]: sign=1
        else: 
            print('ERR in data set %d: lepton id is not identified'%(idx))
            sys.exit()
  
        x=self.tabs[idx]['X']
        y=self.tabs[idx]['Y']
        Q2=self.tabs[idx]['Q2'][i]
        M2=conf['aux'].M2
        YP=1+(1-y)**2
        YM=1-(1-y)**2
  
        thy=2*np.pi*conf['aux'].alfa**2/x/y/Q2
        thy*=(YP+2*x**2*y**2*M2/Q2)*F2-y**2*FL+sign*YM*x*F3
  
        return thy

    def _get_sigma_dxdQ2(self,idx,F2,FL,F3):
        thy=self._get_sigma_dxdy(idx,F2,FL,F3)
        x=self.tabs[idx]['X']
        rs=self.tabs[idx]['RS']
        M2=conf['aux'].M2
        return thy/x/(rs**2-M2)

    def get_sigcc_r(self,X,Q2,idx,tar):

        current=self.tabs[idx]['current'][0].strip()
  
        if current=='NC':

            F2  = self.thy.get_stf(X,Q2,stf='F2c',tar=tar)
            FL  = self.thy.get_stf(X,Q2,stf='FLc',tar=tar)
            F3  = 0#self.thy.get_stf(X,Q2,stf='F3c',tar=tar)
            thy = self._get_sigma_r(idx,F2,FL,F3)

        return thy

    def get_A_PV_e(self,x,Q2,idx,tar):
        #--parity violating asymmetry with longitudinally polarized lepton
        #--derived from Eq. (11) in arxiv.org/pdf/hep-ph/9401264.pdf

        #--get structure functions
        F2gZ = self.thy.get_stf(x,Q2,stf='F2gZ',tar=tar)
        FLgZ = self.thy.get_stf(x,Q2,stf='FLgZ',tar=tar)
        F3gZ = self.thy.get_stf(x,Q2,stf='F3gZ',tar=tar)
        F2g  = self.thy.get_stf(x,Q2,stf='F2g', tar=tar)
        FLg  = self.thy.get_stf(x,Q2,stf='FLg', tar=tar)

        M2 = conf['aux'].M2

        rho2 = 1 + 4*x**2*M2/Q2

        #--get y
        if 'Elab' in self.tabs[idx]:
            ELab=self.tabs[idx]['Elab']
            s=M2+2*ELab*M2**0.5
            y=(Q2/2/x)/((s-M2)/2)
        elif 'Y' in self.tabs[idx]:
            y=self.tabs[idx]['Y']

        YP = y**2*(rho2+1)/2 - 2*y +2
        YM = 1-(1-y)**2

        GF    = conf['aux'].GF

        sin2w = np.array([conf['eweak'].get_sin2w(q2) for q2 in Q2])
        alpha = np.array([conf['eweak'].get_alpha(q2) for q2 in Q2])

        gA = -0.5
        gV = -0.5 + 2*sin2w

        C  = GF*Q2/(2*np.sqrt(2)*np.pi*alpha)

        #--NLO
        #num1 = gA*(YP*F2gZ - y**2*FLgZ)
        #num2 = gV*x*YM*F3gZ
        #num = num1 + num2
        #denom = YP*F2g - y**2*FLg
        #--corrected?
        F1gZ = (rho2*F2gZ - FLgZ)/(2*x)
        F1g  = (rho2*F2g  - FLg )/(2*x)
        
        A = C*(gA*F1gZ/F1g + gV*YM*F3gZ/(2*YP*F1g))

        return A 

    def get_theory(self):

        #--note: for nuclear targets (d, p, h, ...), observables need to be per nucleon (divided by # of protons + neutrons)

        STFs      = ['F2','FL','F3','F2g','FLg','F3g','F2gZ','FLgZ','FLgZ','F3gZ','F2Z','FLZ','F3Z']
        STFs.extend(['W2+','WL+','W3+','W2-','WL-','W3-'])
        F2ratios  = ['F2d/F2p','F2n/F2d','F2d/F2h','F2h/F2d','F2t/F2d','F2h/F2t','F2n/F2p']
        sigratios = ['sigd/sigp','sigh/sigd','sigh/sigt']
        for idx in self.tabs:

            X=self.tabs[idx]['X']
            Q2=self.tabs[idx]['Q2']
            target=self.tabs[idx]['target'][0]
            obs=self.tabs[idx]['obs'][0].strip()

            if   obs in STFs                   :  xsec=self.thy.get_stf    (X,Q2,obs,target)
            elif obs in F2ratios               :  xsec=self.get_F2_ratio   (X,Q2,idx,obs[2],obs[-1])
            elif obs in sigratios              :  xsec=self.get_sig_ratio  (X,Q2,idx,obs[3],obs[-1])
            elif obs=='sig_r'                  :  xsec=self.get_sigma_r    (X,Q2,idx,target)
            elif obs=='dsig/dxdy'              :  xsec=self.get_sigma_dxdy (X,Q2,idx,target)
            elif obs=='dsig/dxdQ2'             :  xsec=self.get_sigma_dxdQ2(X,Q2,idx,target)
            #--heavy quark observable
            elif obs=='sigcc_r'                :  xsec=self.get_sigcc_r    (X,Q2,idx,target)
            #--EIC observable
            elif obs=='A_PV_e'                 :  xsec=self.get_A_PV_e     (X,Q2,idx,target)
            else:
                msg='exp=%d obs=%s and target=%s not implemented'
                raise ValueError(msg%(idx,obs,target))  

            self.tabs[idx]['thy']=xsec
  
    def _gen_report(self,verb=1,level=1):
        """
        verb = 0: Do not print on screen. Only return list of strings
        verv = 1: print on screen the report
        level= 0: only the total chi2s
        level= 1: include point by point 
        """

        L=[]

        if len(self.tabs.keys())!=0:
            nuc      = conf['idis'].nuc
            if nuc:
                if 'dsmf type' in conf: dsmf = conf['dsmf type']
                else:                   dsmf = 'Paris'
                if 'hsmf type' in conf: hsmf = conf['hsmf type']
                else:                   hsmf = 'KPSV'
                nuc = 'True (deuterium:%s, helium:%s)'%(dsmf,hsmf)
            tmc      = conf['idis'].tmc
            offshell = conf['idis'].offpdf
            ht       = conf['idis'].ht
            if ht: ht_type = conf['idis'].ht_type
            else:  ht_type = False
            L.append('reaction: unpol DIS [nuc=%s, tmc=%s, ht=%s, offshell=%s]'%(nuc,tmc,ht_type,offshell))
            for f in conf['datasets']['idis']['filters']:
                L.append('filters: %s'%f)

            L.append('%7s %3s %10s %5s %10s %10s %10s %10s'%('idx','tar','col','npts','chi2','chi2/npts','rchi2','nchi2'))
            for k in self.tabs:
                if len(self.tabs[k])==0: continue 
                res=self.tabs[k]['residuals']

                rres=[]
                for c in conf['rparams']['idis'][k]:
                    rres.append(conf['rparams']['idis'][k][c]['value'])
                rres=np.array(rres)

                if k in conf['datasets']['idis']['norm'] and conf['datasets']['idis']['norm'][k]['fixed']==False:
                    norm=conf['datasets']['idis']['norm'][k]
                    nres=(norm['value']-1)/norm['dN']
                else:
                    nres=0

                
                chi2=np.sum(res**2)
                rchi2=np.sum(rres**2)
                nchi2=nres**2
                tar=self.tabs[k]['target'][0]
                col=self.tabs[k]['col'][0].split()[0]
                npts=res.size
                L.append('%7d %3s %10s %5d %10.2f %10.2f %10.2f %10.2f'%(k,tar,col,npts,chi2,chi2/npts,rchi2,nchi2))

            if level==1:
                L.append('-'*100)  
                for k in self.tabs:
                    if len(self.tabs[k]['value'])==0: continue 
                    if k in conf['datasets']['idis']['norm'] and conf['datasets']['idis']['norm'][k]['fixed']==False:
                        norm=conf['datasets']['idis']['norm'][k]
                        nres=(norm['value']-1)/norm['dN']
                        norm=norm['value']
                    else:
                        norm=1.0
                        nres=0
                    for i in range(len(self.tabs[k]['value'])):
                        x     = self.tabs[k]['X'][i]
                        Q2    = self.tabs[k]['Q2'][i]
                        res   = self.tabs[k]['residuals'][i]
                        thy   = self.tabs[k]['thy'][i]
                        exp   = self.tabs[k]['value'][i]
                        alpha = self.tabs[k]['alpha'][i]
                        rres  = self.tabs[k]['r-residuals'][i]
                        col   = self.tabs[k]['col'][i]
                        shift = self.tabs[k]['shift'][i]
                        tar   = self.tabs[k]['target'][i]
                        msg='col=%7s, tar=%5s, x=%10.3e, Q2=%10.3e, exp=%10.3e, alpha=%10.3e, thy=%10.3e, shift=%10.3e, chi2=%10.3e, res=%10.3e, norm=%10.3e, '
                        L.append(msg%(col,tar,x,Q2,exp,alpha,thy,shift,res**2,res,norm))

        if verb==0:
            return L
        elif verb==1:
            for l in L: print(l)
        return L




    def gen_report(self,verb=1,level=1):

        msg ='reaction: DIS <br />'
        if 'filters' in conf['datasets']['idis']:
            for f in conf['datasets']['idis']['filters']:
                msg+='filters: %s <br />'%f
        msg+='reaction: DIS <br />'


        data={_:[] for _ in ['idx','col','obs','tar','npts','chi2','chi2/npts','rchi2','nchi2']}
        for idx in self.tabs:
            if len(self.tabs[idx])==0: continue
            res=self.tabs[idx]['residuals']

            rres=[]
            for c in conf['rparams']['idis'][idx]:
                rres.append(conf['rparams']['idis'][idx][c]['value'])
            rres=np.array(rres)

            if idx in conf['datasets']['idis']['norm'] and conf['datasets']['idis']['norm'][idx]['fixed']==False:
                norm=conf['datasets']['idis']['norm'][idx]
                nres=(norm['value']-1)/norm['dN']
            else:
                nres=0

            chi2=np.sum(res**2)
            rchi2=np.sum(rres**2)
            nchi2=nres**2
            tar=self.tabs[idx]['target'][0]
            col=self.tabs[idx]['col'][0].split()[0]
            obs=self.tabs[idx]['obs'][0].strip()
            npts=res.size
            data['idx'].append('%7d'%idx)
            data['col'].append('%10s'%col)
            data['obs'].append('%3s'%obs)
            data['tar'].append('%10s'%tar)
            data['npts'].append('%10.2f'%npts)
            data['chi2'].append('%10.2f'%chi2)
            data['chi2/npts'].append('%10.2f'%(chi2/npts))
            data['rchi2'].append('%10.2f'%rchi2)
            data['nchi2'].append('%10.2f'%nchi2)

        data=pd.DataFrame(data)
        msg+=data.to_html(col_space=80,index=False)#,justify='center')
        return msg


