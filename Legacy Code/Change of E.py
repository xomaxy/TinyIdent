#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 09:46:39 2018

@author: Andrea
"""

#Change of E 

import xlrd
import numpy as np
import glob
import os
from tqdm import tqdm
import pandas as pd 
from pandas import DataFrame
import os
import math
import openpyxl
import glob
import matplotlib.pyplot as plt
from openpyxl import load_workbook

#os.chdir('/Users/Andrea/Desktop/RESULTADOS/Chitosan')
#FileList = glob.glob('**/*.xls', recursive=True)
#print(FileList)
#
DE = []

pos = []
#
#for i in FileList:
# 
path1 = '/Users/Andrea/Desktop/RESULTADOS/Chitosan/CHOR_P%s_%s_0000.ibw_df2.xlsx'
    
#    name = (FileList[1].split('_')[1])
#    save_path = '/Users/Andrea/Desktop/TxtResults '
#    savepath = os.path.join(save_path, name +".txt") 

    
df = pd.read_excel(path1, sheetname='Sheet1')
df.head() 
        
E_nomag= df['E_nomag']
E_mag= df['E_mag']
pos= df['pos']                   
    
de = E_mag -  E_nomag
                 
    
DE.append(de)
De = np.sort(DE)
    
#    np.savetxt('Eccentricity', E)
    
    
#np.histogram (DE,bins='auto')

plt.bar(DE, pos, bins='auto',color='r')

plt.title("Eccentricity of microparticles",fontsize=24)
#plt.xticks

plt.ylabel('Frequency',fontsize = 23)
plt.xlabel("Eccentricity", fontsize=23)

plt.legend(('Ellipsoidal', 'Spherical (commercial)'), loc = 5)


Mean = np.mean (DE)

std = np.std (DE)

print (Mean, std)