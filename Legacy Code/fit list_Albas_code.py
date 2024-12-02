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
import glob
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



if __name__ == "__main__":
    plt.close('all')
    sns.set_style('whitegrid')
    sns.set_context(context = 'talk', font_scale = 1.5)
#path1 = '/Volumes/VERBATIM/chitosan only control/CO%s_%s_0000.ibw'
   
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\17-06-28\\chitosan only control\\CO%s_%s_0000.ibw'
#    delete_points = np.array([3])
    
    
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\17-06-27\\chitosan_wires\\C%s_%s_0000.ibw'

    
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\17-06-13\\C%s_%s_0000.ibw'
   
 
    
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\17-06-27\\Glass control\\G%s_%s_0000.ibw'
#    path_ave_dev = 'C:\\Data\\Chitosan-Fe3O4\\17-06-27\\Glass control\\glass_control_Young_ave.xlsx'
    

    
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\17-09-25\\Glass%s_%s_0000.ibw'
#    path_ave_dev = 'C:\\Data\\Chitosan-Fe3O4\\17-09-25\\glass_Young_ave.xlsx'
    
    
    
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\17-09-25\\C%s_%s_0000.ibw'

###############################################################

# last position = real last position

##########################################################
####################GENEPIN###############################
##########################################################    


#    #genepin - 45
##    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-05_wires_45\\C45_P%s_%s_0000.ibw'
#    path1 = '/Volumes/KINGSTON/Chitosan/18-02-05_wires_45/C45_P%s_%s_0000.ibw'
#    delete_points = np.array([5, 6, 7, 14, 18])
#    last_position = 20
    
    #genepin - HOR    
###    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-05_wires_hor\\CHOR_P%s_%s_0000.ibw'
#    path1 = '/Volumes/KINGSTON/Chitosan/18-02-05_wires_hor/CHOR_P%s_%s_0000.ibw'
#    delete_points = np.array([3, 6, 7, 11, 14])
##  20 is horrible, put last position = 19
#    last_position = 19


#    #genepin - ONLY
##    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-21_chitosan_only_genepin\\chit_only_genepin\\chitosan_only_not_oriented\\CO_G_NO_P%s_%s_0_0000.ibw'
##    path1 = '/Volumes/KINGSTON/Chitosan/18-02-21_chitosan_only_genepin/chit_only_genepin/chitosan_only_not_oriented/CO_G_NO_P%s_%s_0_0000.ibw'
#    path1 = '/Volumes/KINGSTON/Chitosan/18-02-21_chitosan_only_genipin/chit_only_genepin/chitosan_only_not_oriented/CO_G_NO_P%s_%s_0_0000.ibw'
#    delete_points = np.array([2, 3, 12, 14])  # also 16 is horrible, put last position = 15
#    last_position = 15


#    #genepin - random
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-03-15_random\\chitos_wires_RANDOM_genepin\\CG_RAN_P%s_%s_0_0000.ibw'
    path1 = '/Volumes/KINGSTON/Chitosan/18-03-15_random/chitos_wires_RANDOM_genepin/CG_RAN_P%s_%s_0_0000.ibw'
    delete_points = np.array([4, 5, 6, 7, 8, 16, 17, 18, 20, 24])
    last_position =  25

#############################################################
#############################################################
#############################################################
#
   #PDMS
##    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-27_PDMS\\PDMS\\P%s_%s_0_0000.ibw'
#    path1 = '/Volumes/KINGSTON/Chitosan/18-02-27_PDMS/PDMS/P%s_%s_0_0000.ibw'
#    delete_points = np.array([7])
#    # 12 is horrible, put last position = 11
#    last_position = 11

#######################################################
##############GENEPIN AFTER############################
#######################################################

##   #genepin after - ONLY
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-08_chit_only_genepin_after\\CO_NG_NO_P%s_%s_0000.ibw'
#    delete_points = np.array([24])
#    last position = 30
    

#   #genepin after - 45
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-27_PDMS\\chitos_wires_45_genepin_after\\CW_GA_45_P%s_%s_0_0000.ibw'
##    delete_points = np.array([])
#    last_position = 6
    
##   #genepin after - RANDOM
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-03-07_RANDOM_genepin_after\\chitos_wires_RANDOM_genepin_after\\CW_GA_RAN_P%s_%s_0000.ibw'
#    delete_points = np.array([10])
#    last_position =  22


##  #genepin after - HOR
#    path1 = 'C:\\Data\\Chitosan-Fe3O4\\18-02-09_HOR_gp_after\\chit_wires_hori_genepin_after\\CW_NG_HO_P1_0_0_0000.ibw'
##    delete points = np.array([])
#    last_position =  10

###################################################################
###################################################################
##################################################################

    initial_position = 1
    positions = np.arange(initial_position, last_position+1)

######## if no points to be deleted or if you want to try all the positions first ################
######## you need to uncomment the variable: delete_points ################
######## you need to uncomment the variable: new_positions = np.delete(positions, delete_points-1) ################

#    new_positions = positions
    
######## if there are to be deleted ################
######## you need to uncomment the variable:  new_positions = positions ################
    
    new_positions = np.delete(positions, delete_points-1)




    E_no_mag = []
    dev_no_mag = []         
    E_north_mag = []
    dev_north_mag = []
    E_south_mag = []
    dev_south_mag = []

    
    m=1
    l=2
    
    lab_orient = None    
    
    #        #looping
    #Htype = 'HDF'
    #batch = '_1_'

    method_to_use = 'fiv'
    fmin = 30e-9
    fmax = 45e-9  
#    R = 10e-6
    R= 24.605e-6
    
    stop=2


    
    
    path_ave_dev = os.path.abspath(path1).split('%s_0_0000.ibw')[0] + 'Young_ave.xlsx'
    
   
    
    
    for k in new_positions:
        bead = [k]
        b=bead[0]
        
        E_ave = 'E_ave_%i_' %(b) 
        E_dev = 'E_dev_%i_' %(b)    
        
        no = list(range(stop))
        #label = ['B%i' %i for i in no]
        T = [(str(k),str(n)) for k in bead for n in no]
        path2 = [path1 % t for t in T]


        
        
        MC = [Multicurve(p2) for p2 in path2]
        
        
        for mc in MC:
            mc.LoadCurves(100)
        
        
        E = []
      # c = " "
        Mean = [] 
        Dev = []
       
        j = 0
        

        
        for mc in MC:
#            print("for orientation = ", j)    
            e = []
        
        
            for fc in mc.curves[:]:
                try:
                    
                    #print(' i am in the loop')
                    
                    ind = 1*fc.indentation.Trace()
                    f= 1*fc.force.Trace()
                    
                     #detrend the slope in the non-contact region
                    p = np.polyfit(fc.indentation.Trace()[0:1700], 
                           fc.force.Trace()[0:1700], 1)
                    
                    f = f - np.polyval(p, ind)
                    
                    cp = ContactPoint()
                    cix = cp.getCP(ind=ind, f=f, method=method_to_use)
                    
                    
#                    fc.correct(method = method_to_use)
#                    cix = fc.contactidx
                    
                    fc.contactidx = 1*cix
                    
                    fc.force.data -= fc.force.data[cix]
                    fc.indentation.data -= fc.indentation.data[cix]
                    
                    
#                    f_min = nearestPoint(fc.force.Trace(), fmin)
#
#                    f_max = nearestPoint(fc.force.Trace(), fmax)
#                    print('number of points')# %i ' %(f_max - f_min))
    
                    
                    e.append(fc.Young(model='h',
                                      fmin=fmin, fmax=fmax, constant='force', R = R))
                    
        #            fc.trace = False
#                    fc.force.data[:,:fc.contactidx] = 0
                    
#                    ind = 1*fc.indentation.Trace()
#                    f= 1*fc.force.Trace()
#                    
#                    cp = ContactPoint()
#                    cix = cp.getCP(ind=ind, f=f, method=method_to_use)
                    
                                      
                                      
                    #print('cix = %d' %cix)
                    # move coordinate system for retrace to zero
    #                    ind -= ind[cix]
    #                    f -= f[cix]
    
                    
    #                    if cix:
    #                        imin = (ind[-1] - ind[cix])*20/100
    #                        imax = (ind[-1] - ind[cix])*50/100
    #                        print('For imin = %f and imax = %f' %( imin, imax))
    #                    else:                        
    #                       imin = 0.21e-6
    #                       imax = 0.34e-6
    #                       print('caca')
                        
                        #Trace
        #            e.append(fc.Young(ind=ind, f=f,
        #                               imax=imin, imin=imax, constant='indentation'))
        #                #Retrace
                except:
        #            c
                    np.nan
            
            
            E.append(e)         
            mean = np.nanmean(e)
            dev = np.nanstd(e)
            print ("The mean E* is" , mean)
            print ("The standard deviation is ", dev)
                       
            Mean.append(mean)
            Dev.append(dev)
            
           
    #        cols = [E_ave, E_dev]
    #        ave_dev = np.array([Mean, Dev]).T
    #        df_ave_E = pd.DataFrame(data=ave_dev, columns=cols)
    #        df_ave_E.to_excel(savepath_ave_dev)
            
       
            rand= np.random.randint(0, len(mc.curves) - 1)
            print('this is curve %i' %rand)
                    
            fc = mc.curves[rand]
            
        #    Line substraction
#            ix = fc.contactidx
#            
#            p = np.polyfit(fc.indentation.Trace()[ix:ix+40], 
#                           fc.force.Trace()[ix:ix+40], 1)
#        
#            line = np.polyval(p, fc.indentation.data[ix:])
#        
#            fc.force.data[ix:] -= line
        
        #TRACE
            

            ix = fc.contactidx 
            fc.force.data[:fc.contactidx] = 0
            fix = nearestPoint(fc.force.Trace(), fmin)
            mix = nearestPoint(fc.force.Trace(), fmax)
            
            
            save = 0   
            lab_fit = None
            
            
            if j==0:    
                colour = 'ko'
                labell = 'no magnet'
            elif j==1:
                colour = 'ro'
                labell = 'N'
            else:
                colour = 'bo'
                labell = 'S'
                
            if j == stop-1:
                save = 1
                save2 = 1
                lab_fit = 'Hertz Fit'
                
                
            fig1 = plt.figure(m)
            ax = fig1.add_subplot(111)
                 
                
            p2 = ax.plot(fc.indentation.Trace()*1e6, fc.force.Trace()*1e9, colour, linewidth=2, label=labell)
            print('fit =  ', fc.fit)
            p1 = ax.plot(fc.indentation.Trace()[fix:mix]*1e6, IndentationFit(fc.fit[0], 
                             fc.indentation.Trace()[fix:mix], gamma=1.5)*1e9,'y', linewidth=7, label = lab_fit) 
            
            
            plt.title('Force Curves: Trace C%i' %k, size=30)
            plt.xlabel('$\delta$ ($\mu$m)', size=20);
            plt.ylabel('Force (nN)', size=20)
            plt.legend(loc="best", prop={'size':20})
            plt.xticks( color='k', size=20) 
            plt.yticks( color='k', size=20)    
            plt.xlim(left=-1)
            
            
            
            if save == 1:
                savepath0 = os.path.abspath(path2[j]).split('_1_0000')[0] + '_Young_fit.jpg'
                plt.savefig(savepath0, bbox_inches="tight")
            else:
                pass
            
            
    
                
            if k==(last_position):
                lab_orient = labell
            else:
                pass
            
            
            
            mean = np.nanmean(e)
            dev = np.nanstd(e)
            print ("The mean E* is" , mean)
            print ("The standard deviation is ", dev)
            
            if j==0:    
                E_no_mag.append(mean) 
                dev_no_mag.append(dev)
            elif j==1:
                E_north_mag.append(mean) 
                dev_north_mag.append(dev)
            else:
                E_south_mag.append(mean) 
                dev_south_mag.append(dev)
            
            
            
            
            fig0 = plt.figure(0)
            plt.errorbar(k, mean, yerr=3*dev, fmt=colour, linewidth=2, label=lab_orient)
            plt.title('Young\'s Modulus Distribution', size=30)
            plt.xlabel('Position', size=40);
            plt.ylabel("E*(Pa)", size=40)
            plt.legend(loc="best", prop={'size':20})
            plt.xticks( color='k', size=20) 
            plt.yticks( color='k', size=20)   
            plt.semilogy()
            #plt.ylim([0,21000])
            plt.xlim([0,last_position+1])
            
            
            if save == 1 and k==(last_position):
                savepath_averages = os.path.abspath(path2[j]).split('%s_1_0000')[0] + 'Young_averages.jpg'
                plt.savefig(savepath_averages, bbox_inches="tight")
            else:
                pass
            
            j += 1  
            
        
            
            
       
        
        for i in range(len(path2)):
            savepath = os.path.abspath(path2[i]).split('_0000')[0] + '_Young_fit.xlsx'
            
            
            
            print('save at: ', savepath)
                #    np.savetxt(savepath, E[i])
                    
            df = pd.DataFrame({'E':E[i]})
            df.to_excel(savepath)
            
            
            save2 = 0            
            
            if i==0:    
                colour = 'ks'
                labell = 'no magnet'
            elif i==1:
                colour = 'rs'
                labell = 'N'
            else:
                colour = 'bs'
                labell = 'S'
                
            if i == stop-1:
                save2 = 1
                
            plt.figure(l)
             
            fig2 = plt.plot(df.E, colour, label=labell)
            plt.title('Force Curves: Trace C%i' %k, size=30)
            plt.ylabel("E*(Pa)", size=20)
            plt.xlabel("Counts", size=20)
            plt.legend(loc="best", prop={'size':20})
            plt.xticks( color='k', size=20) 
            plt.yticks( color='k', size=20)  
            #plt.ylim([0,30000])
            
            if save2 == 1:
                savepath1 = os.path.abspath(path2[i]).split('_0000')[0] + '_Young_dispersion.jpg'
                plt.savefig(savepath1, bbox_inches="tight")
            else:
                pass
             
            
            
    #        #    plt.clf()
    #            
    #    plt.title('Force Curves', size=30)
    
        m+=2    
        l+=2 
        print (Mean)
        df_ave_E = pd.DataFrame()
        df_ave_E[E_ave] = Mean
        df_ave_E[E_dev] = Dev
        df_ave_E.to_excel(path_ave_dev)     
     

    mag_label = ['magnet' for i in range(len(E_north_mag))]
    no_mag_label = ['no magnet' for i in range(len(E_no_mag))]
    labels = [mag_label, no_mag_label]
    labels = list(itertools.chain.from_iterable(labels))    
    young = np.concatenate([E_north_mag, E_no_mag])
    errs = [dev_north_mag, dev_no_mag]
    errs = list(itertools.chain.from_iterable(errs))
    pos = np.concatenate([np.arange(1,len(E_north_mag)+1), np.arange(1,len(E_no_mag)+1)])
    red_black = ['#ff0000', '000000']
    sns.set_palette(red_black)
    df2 = pd.DataFrame({'pos': pos, 'E': young, 'label': labels, 'error': errs})
    fig, ax = plt.subplots()
    sns.barplot(data=df2, x='pos', y='E', hue='label', ax=ax)
    ax.set_ylabel(r"$\tilde{E}$ [Pa]", size = 20)
    ax.set_xlabel("Position", size = 20)
    ax.legend_.set_title('')

    for i in range(0,len(E_no_mag)):
        ax.plot([i-.2, i-.2], [E_north_mag[i], E_north_mag[i] + dev_north_mag[i]], 'gray')
        ax.plot([i-.2, i-.2], [E_north_mag[i], E_north_mag[i] - dev_north_mag[i]], 'gray')
        ax.plot([i+.2, i+.2], [E_no_mag[i], E_no_mag[i] + dev_no_mag[i]], 'gray')
        ax.plot([i+.2, i+.2], [E_no_mag[i], E_no_mag[i] - dev_no_mag[i]], 'gray')




    if save2 == 1 and k==(last_position):
        savepath_df2 = os.path.abspath(path1).split('%s_0_0000.ibw')[0] + '_df2.xlsx'
        df2.to_excel(savepath_df2)
        savepath_bar_ch = os.path.abspath(path1).split('%s_0_0000.ibw')[0] + 'Young_bar_chart.jpg'
        plt.savefig(savepath_bar_ch, bbox_inches="tight")
    else:
        pass











        
     
#        indexes_bar = np.arange(len(E_no_mag))
#        width = 0.20
#        fig, ax = plt.subplots()
#        ax.set_xlim(-width,len(indexes_bar)+width)
#        ax.set_ylabel('Ẽ (Pa)', size=20)
#        ax.set_xlabel('Position', size=20)
#        #ax.set_title('Chitosan/Fe3O4 stiffness', size=30)
#        ax.set_xticks(indexes_bar + width)
##        ax.legend(loc="best", prop={'size':20})
##        ax.tick_params(labelsize=20)
#        xTickMarks = ['P'+str(i) for i in range(1,last_position+1)]
#        ax.set_xticklabels(xTickMarks)
#        rects1 = ax.bar(indexes_bar, E_no_mag, width, color='k', yerr=dev_no_mag)
#        rects2 = ax.bar(indexes_bar + width, E_north_mag, width, color='r', yerr=dev_north_mag)
#        
#        if stop == 3:
#            rects3 = ax.bar(indexes_bar + 2*width, E_south_mag, width, color='b', yerr=dev_south_mag)
#            ax.legend((rects1[0], rects2[0], rects3[0]), ('No Magnet', 'North', 'South'))
#        else:
#            ax.legend((rects1[0], rects2[0]), ('No Magnet', 'North'))
#        

#            
        
                
     
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
            #   plt.plot()
            #            ,plt.xlabel(r'$\delta$ in um'),
            #            plt.ylabel(r'$F$ in nN'), plt.legend(loc='upper left'))
            
            
                
        
            
                        
                
                
                
                
                
                
                
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
        
        
        
    
    
