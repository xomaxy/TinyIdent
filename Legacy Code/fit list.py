#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 14:09:24 2017

@author: Andrea
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 13:57:28 2017

@author: Andrea
"""

import numpy as np

import os
from tqdm import tqdm
from ForceMetric import (ForceCurve, Multicurve, IndentationFit, ContactPoint,
                         NeoHookeanBead, nearestPoint)
from matplotlib import pyplot as plt
import pandas as pd
#from Tools import reduced_qty
#import brewer2mpl
import seaborn as sns
from pprint import pprint
import itertools
import statistics
import matplotlib.patches as mpatches



path1 = '/Users/andrea/Documents/RESULTADOS/Alginate FD vs normal/Day8/FD_D%s_B%s_0000.ibw'




#        #looping for day of experiment (D%s) and batch number (B%s)
#Htype = 'HDF'
#batch = '_1_'


position = [8]
magnet = [29]
label = ['X%i' %i for i in magnet]
T = [(str(b),str(n)) for b in position for n in magnet]
path1 = [path1 % t for t in T]


MC = [Multicurve(p1) for p1 in path1]


for mc in MC:
    mc.LoadCurves(5)


E = []
c = " "
Mean = []

j = 0
for mc in MC:
    j += 1
    print(j)
    e = []
    for fc in mc.curves[:]:
        try:
            fc.correct()
#            fc.trace = False
            
            ind = 1*fc.indentation.Retrace()
            f= 1*fc.force.Retrace()
            cp = ContactPoint()
            cix = cp.getCP(ind=ind, f=f)
            # move coordinate system for retrace to zero
            ind -= ind[cix]
            f -= f[cix]
            imin = 10e-9
            imax = 400e-9
            e.append(fc.Young(model='h',
                              imin=imin, imax=imax, constant='indentation'))
                #Trace
#            e.append(fc.Young(ind=ind, f=f,
#                               imax=imin, imin=imax, constant='indentation'))
#                #Retrace
        except:
#            c
            np.nan
    E.append(e)     
    mean = np.nanmean(e)
    print ("The mean E* is" , mean)
    print ("The standard deviation is ", np.nanstd(e))
    Mean.append(mean)
    
    rand = np.random.randint(0, len(mc.curves) - 1)
    print('this is curve %i' %rand)
    
    fc = mc.curves[rand]
    
    fig = plt.figure(facecolor='w')
    
#    Line substraction
#    ix = fc.contactidx
#    
#    p = np.polyfit(fc.indentation.Trace()[ix:ix+20], 
#                   fc.force.Trace()[ix:ix+20], 1)
#
#    line = np.polyval(p, fc.indentation.data[ix:])
#
#    fc.force.data[ix:] -= line

#TRACE
    plt.figure(1)
    sns.set_context('paper', font_scale=2)
    sns.set_style('whitegrid')
#    plt.style.use('ggplot')
    ix = fc.contactidx 
    fix = nearestPoint(fc.indentation.Trace(), imin)
    mix = nearestPoint(fc.indentation.Trace(), imax)

    ax = fig.add_subplot(111)

    p3 = ax.plot(fc.indentation.Retrace()*1e6, fc.force.Retrace()*1e9, 'b-', label='Unloading curve')
    p2 = ax.plot(fc.indentation.Trace()*1e6, fc.force.Trace()*1e9, 'r-', label='Loading curve')

    p1 = ax.plot(fc.indentation.Trace()[fix:mix]*1e6, IndentationFit(fc.fit[0], 
                 fc.indentation.Trace()[fix:mix], gamma=1.5)*1e9,'k*', label='Hertz model fit(10-400nm)') 
    #gamma =1, flat
    #punch, 2.5: sneeden, 1.5: Hertz
#    p1 = ax.plot(fc.indentation.Trace()[fix:mix]*1e6, NeoHookeanBead( 
#                 fc.indentation.Trace()[fix:mix], fc.fit[0][0], 10e-6)*1e9)

##RETRACE
#
#    ix = fc.contactidx 
#    fix = nearestPoint(fc.indentation.Retrace(), imin)
#    mix = nearestPoint(fc.indentation.Retrace(), imax)
#
#    ax = fig.add_subplot(111)
#
#    p2 = ax.plot(fc.indentation.Retrace()*1e6, fc.force.Retrace()*1e9, 'ro')
#
#    p1 = ax.plot(fc.indentation.Retrace()[fix:mix]*1e6, 
#                 IndentationFit(fc.fit[0], 
#                 fc.indentation.Retrace()[fix:mix], gamma=1.5)*1e9,'bo')
#
#    fc.plot()
    plt.plot(facecolor='w')
    plt.xlabel(r'$\delta$ in $\mu$m', fontsize=16)
    plt.ylabel(r'F in nN', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc=2, fontsize=14)
    

#plt.clf()
    
plt.figure(2)    
 
for i in range(len(path1)):
    savepath = os.path.abspath(path1[i]).split('0000')[0] + '_Young_fit.xlsx'
    savepath1 = os.path.abspath(path1[i]).split('0000')[0] + '_Young_fit.jpg'
    print('save at: ', savepath)
#    np.savetxt(savepath, E[i])
    df = pd.DataFrame({'E':E[i]})
    df.to_excel(savepath)
    fig = plt.plot(df, 's')
    plt.ylabel("E*(Pa)" )
    plt.xlabel("Bead indentation count")
#    red_patch = mpatches.Patch(color='red', label='South')
    blue_patch = mpatches.Patch(color='blue', label='No magnet')
    green_patch = mpatches.Patch(color='green', label='Magnet')
    plt.legend(handles=[blue_patch, green_patch],
               bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(savepath1)  
#    plt.clf()
    

print (Mean)

  
    
#    fig = sns.swarmplot(y=E[i])     
     
     
     
    ######Average and std.

#    def average(x):
#    
#        d=len(x)
#        av=0
#        b=0
#        
#        #in order to exclude points with null value
#        
#        for e in x:
#            if e:
#                av = av + e
#            else:
#                b += 1
#                
#        significant_d = d - b
#        av = float(av) / significant_d        
#        return av
#    
#    def sig2(x,av=None):
#        ''' sig2(x,av=None) x is a vector '''
#        d=len(x)
#        j=0
#        if not av :
#            av = average(x)
#        s2=0
#        
#        for e in x:
#            if e:
#                dev = e - av
#                s2 = s2 + dev*dev
#            else:
#                j += 1
#        
#        significant_d = d - j
#        
#        if significant_d > 1 :
#            s2 = float(s2) / (significant_d - 1)
#            s = sqrt(s2)
            
#        return s





