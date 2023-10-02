import pandas as pd
import sys
import os
if len(sys.argv)>1:
	excel_files = [sys.argv[1]]
else: 
	# excel_files = ["100%02d.xlsx" %(ii+1) for ii in range(70)]
    excel_files =os.listdir('data/JAM/')

n = 0
for excel_file in excel_files:
    if not ('100' in excel_file):
        continue
# for excel_file in excel_files:
    # if excel_file=="10012.xlsx" or excel_file=="10013.xlsx" or excel_file=="10023.xlsx" or excel_file=="10024.xlsx":
        # continue
    print("Converting '{}'".format(excel_file))
    # try:
    #     df = pd.read_excel("format/"+excel_file,"format")
    #     output = "csv/"+excel_file.split('.')[0]+'.csv'
    #     df.to_csv(output)    
    # except:
    #     print("  Failed to convert")
    df = pd.read_excel("data/JAM/"+excel_file,"format")
    output = "data/JAM/csv/"+excel_file.split('.')[0]+'.csv'
    df.to_csv(output)    
