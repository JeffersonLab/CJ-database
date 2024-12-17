#!/usr/bin/env python
import sys, os
import numpy as np
import copy
import pandas as pd


if __name__ == "__main__":

    idxs = [10037, 10038, 10041, 10042, 10050, 10051, 10052, 10053, 10054, 10055, 10056]

    for idx in idxs:
        excel = pd.DataFrame(pd.read_excel('dataframe/%s.xlsx'   %idx))
    
        excel.to_csv  ('dataframe/csv/%s.csv'%idx)

        print('Copying %s.xlsx to %s.csv...'%(idx,idx))
        
        
        
        
        
        
        
        
        