#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Daniel 

This script contains all the functions used to analyse the structural data. 

It contains a function to calculate spearmans correlation, cohens d and kruskal test. 
"""

import pandas as pd
from scipy import stats
import math
import pingouin as pin
from statsmodels.stats import multitest
import numpy as np

def post_hoc_mwu(x,y,z, fwe_method=''):
    x_y=pin.mwu(x,y)
    x_z=pin.mwu(x,z)
    y_z=pin.mwu(y,z)
    corrp= multitest.multipletests(np.concatenate((x_y['p-val'].values,x_z['p-val'].values, y_z['p-val'].values )), method=fwe_method)
    return corrp

def correlation(column, volumearray, volume):
    
    '''
    Runs Spearmans correlation.
    
    Input: df with measures that area all used in the correlation. 
    
    Second array/df with thrid variable (string) that states which measure from this second array/df to be used in the correlation.
    
    Returns: pvalues and correlation values in list
    
    '''
    
    pvalues = []
    correlation = []
    for i in column.columns:
        array = pd.concat([volumearray[volume], column[i]],axis=1).dropna()
        c, p = stats.spearmanr(array[volume], array[i])
        pvalues.append(p)
        correlation.append(c)
    return pvalues, correlation



def cohen_d(group1,group2):
    
    '''
    Calculate cohens d.
    
    Input: two series/list
    --------------------------
    
    Output: d
    
    '''
    
    
    diff = group1.mean() - group2.mean()
    pooledstdev = math.sqrt((group1.std()**2 + group2.std())/2 )
    cohend = diff / pooledstdev
    return cohend

def kruskal(x,y,z):
    
    '''
    Runs kruskal-wallis test.
    
    Input: Series from three groups.
    --------------------------------------
    
    Output: df with pvals, degrees of freedom and eta/epsilon effect sizes
    
    
    '''
    
    number_of_observations = len(pd.concat([x,y,z]))
    number_of_groups = len([x, y, z])
    h,p = stats.kruskal(x,y,z)
    degrees_of_freedom = number_of_observations - number_of_groups
    eta = (h-number_of_groups +1)/(number_of_observations-number_of_groups)
    epsilon  = h/((number_of_observations**2 -1) / (number_of_observations +1))
    df= pd.DataFrame(data={'pval': [p], 'kruskal_test_statistic': [h],'df': [degrees_of_freedom],'eta': [eta],'epsilon': [epsilon]})
    return df


def multi_comparisons(dictionary,aan,wr,hc):

    '''
    Function to test if parametric or non-parametric multi comparisons test should be used
    
    Parameters
    ----------- 
    dictionary : dict object of anova/kruskal-wallis results
    aan : DataFrame of aan results
    wr : DataFrame of wr results
    hc : DataFrame of hc results


    Output
    -------
    multi_comp_dict : dict of sorted values into parametric or non-parametric 

    '''
    
    para = []
    non_para = []
    
    for key, val in dictionary.items():
        try:
              if val['PR(>F)'][0] <0.05:
                
                    if stats.normaltest(aan[key].dropna()) and stats.normaltest(wr[key].dropna()) and stats.normaltest(hc[key].dropna()) and stats.levene(aan[key].dropna(), wr[key].dropna(),hc[key].dropna())[1] > 0.05:
                        para.append(key)
                        
                      
                    else:
                        non_para.append(key)

        except KeyError:

            if float(val['pval']) <0.05:
                  
                if stats.normaltest(aan[key].dropna()) and stats.normaltest(wr[key].dropna()) and stats.normaltest(hc[key].dropna()) and stats.levene(aan[key].dropna(), wr[key].dropna(),hc[key].dropna())[1] > 0.05:
                    para.append(key)
                      
                else:
                    non_para.append(key)


    multi_comp_dict = {
      'parametric': para,
      'nonparametric': non_para
      }

    return multi_comp_dict

def mean_values(measures,aan,wr,hc):
    
    '''
    Function to calculate the mean value and std deviation of

    Parameters
    --------------

    measure (str) : list of measures
    aan : DataFrame of aan results
    wr : DataFrame of wr results
    hc : DataFrame of hc results

    Returns
    --------------
    value_dict (dict): dictionary of mean and std of measures



    '''

    aan_mean_values= dict(zip(['aan_'+ measure for measure in measures], [aan[measure].mean() for measure in measures]))
    aan_std = dict(zip(['aan_std_'+ measure for measure in measures], [aan[measure].std() for measure in measures]))
      
    wr_mean_values= dict(zip(['wr_'+ measure for measure in measures], [wr[measure].mean() for measure in measures]))
    wr_std = dict(zip(['wr_std_'+ measure for measure in measures], [wr[measure].std() for measure in measures]))
      
    hc_mean_values= dict(zip(['hc_'+ measure for measure in measures], [hc[measure].mean() for measure in measures]))
    hc_std = dict(zip(['hc_std_'+ measure for measure in measures], [hc[measure].std() for measure in measures]))
      
    value_dict = {**aan_mean_values,**aan_std,**wr_mean_values,**wr_std,**hc_mean_values,**hc_std}

    return value_dict

def ttest_mwu(group1,group2):
    if stats.normaltest(group1.dropna()) and stats.normaltest(group2.dropna()) and stats.levene(group1.dropna(), group2.dropna())[1] > 0.05:
        test= pin.ttest(group1.dropna(),group2.dropna())
    else:
        test = pin.mwu(group1.dropna(),group2.dropna())
    
    return test