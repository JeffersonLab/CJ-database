import sys
import pandas as pd
from qcdlib.aux import AUX
from tools.reader import _READER
from tools.config import conf

class READER(_READER):
  
    def __init__(self):
        self.aux=conf['aux']
  
    def get_X(self,tab):
        cols=tab.columns.values
        if any([c=='X' for c in cols])==False:
            if any([c=='W2' for c in cols]):
                tab['X']=pd.Series(tab['Q2']/(tab['W2']-self.aux.M2+tab['Q2']),index=tab.index)
            elif any([c=='W' for c in cols]):
                tab['X']=pd.Series(tab['Q2']/(tab['W']**2-self.aux.M2+tab['Q2']),index=tab.index)
            else:
                print('cannot retrive X values')
                sys.exit()
        return tab
  
    def get_Y(self,tab):
        cols=tab.columns.values
        if any([c=='Y' for c in cols])==False:
            M=self.aux.M
            if any([c=='Elab' for c in cols]):
                tab['Y']=pd.Series(tab['Q2']/(2*self.aux.M*tab['X']*tab['Elab']),index=tab.index)
            elif any([c=='E' for c in cols]):
                tab['Y']=pd.Series(tab['Q2']/(2*self.aux.M*tab['X']*tab['E']),index=tab.index)
            elif any([c=='RS' for c in cols]):
                tab['Y']=pd.Series(tab['Q2']/(tab['RS']**2-self.aux.M**2)/tab['X'],index=tab.index)
        return tab
  
    def get_W2(self,tab):
        cols=tab.columns.values
        if any([c=='W2' for c in cols])==False: 
            tab['W2'] = pd.Series(self.aux.M2 + tab.Q2/tab.X - tab.Q2,index=tab.index)
        return tab
  
    def get_idx(self,tab):
        tab['idx']=pd.Series(tab.index,index=tab.index)
        return tab
  
    def modify_table(self,tab):
  
        tab=self.get_X(tab)   
        tab=self.get_Y(tab)   
        tab=self.get_W2(tab)   
        tab=self.apply_cuts(tab)
        tab=self.get_idx(tab)
  
        return tab
 

 
