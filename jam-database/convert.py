#!/usr/bin/env python
import sys, os
import matplotlib
matplotlib.use('Agg')
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

import pandas as pd

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

if __name__ == "__main__":

  F = open('f2_data.txt','r')
  L = F.readlines()
  F.close()

  L = [l.strip() for l in L]
  L = [[x for x in l.split()] for l in L]
  #L = np.transpose(L)[0]

  D = {}

  print(np.shape(L))

  print(L[0])

  X, UB, DB, UBUP, UBDO, DBUP, DBDO = [],[],[],[],[],[],[]
  for i in range(len(L)):
      for j in range(len(L[i])):
          if i==0: D[L[i][j]] = []
          else:
              D[L[0][j]].append(L[i][j])

  D = pd.DataFrame(D)
  filename = 'test.xlsx'
  D.to_excel(filename,index=False) 





 
  
        
        
        
        
        
        
        
        
        
        
        
        
