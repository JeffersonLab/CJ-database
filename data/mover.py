#!/usr/bin/env python
import sys, os
import numpy as np
import copy
import pandas as pd


if __name__ == "__main__":

    MOVE = {}
    #--MOVE[old idx] = (new idx, exp, obs, tar, text name)
    MOVE[10037] = (10075, 'SLAC-140x',    'sigma', 'p',   None)
    MOVE[10038] = (10076, 'SLAC-140x',    'sigma', 'd',   None)
    MOVE[10041] = (10077, 'JLab E06-009', 'F2',    'd',   'e06009_d_c')
    MOVE[10042] = (10078, 'JLab E06-009', 'sigma', 'd',   'e06009_sd_c')
    MOVE[10050] = (10079, 'JLab E99-118', 'F2',    'd',   'e99118_d')
    MOVE[10051] = (10080, 'JLab E99-118', 'F2',    'd/p', 'e99118_dp')
    MOVE[10052] = (10081, 'JLab E99-118', 'sigma', 'p',   'e99118_sp_c')
    MOVE[10053] = (10082, 'Jlab E99-118', 'sigma', 'd',   'e99118_sd_c')
    MOVE[10054] = (10083, 'JLab E99-118', 'sigma', 'd/p', 'e99118_sdp')
    MOVE[10055] = (10084, 'JLab JLCEE96', 'sigma', 'p',   'jlcee96_sp_c')
    MOVE[10056] = (10085, 'JLab JLCEE96', 'sigma', 'd',   'jlcee96_sd')

    sys.exit()
    #--DON'T RUN THIS
    #for old_idx in MOVE:
    #    move = MOVE[old_idx]
    #    new_idx, exp, obs, tar, text_name = move[0],move[1],move[2],move[3],move[4]
    #    #--read excel
    #    excel = pd.read_excel('dataframe/%s.xlsx'   %old_idx)
    #
    #    #--read csv
    #    csv   = pd.read_csv  ('dataframe/csv/%s.csv'%old_idx)

    #    #--copy excel file to new location
    #    excel.to_excel('dataframe/%s.xlsx'%new_idx)
    #    excel.to_csv  ('dataframe/csv/%s.csv'%new_idx)

    #    print('Moving %s to %s...'%(old_idx,new_idx))
        
        
        
        
        
        
        
        
        
