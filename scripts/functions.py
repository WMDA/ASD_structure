import pandas as pd
from scipy import stats
import math
import pingouin as pin
from statsmodels.stats import multitest
import numpy as np

"""
This script contains all the functions used to analyse the structural data. 
Function to calculate spearmans correlation, cohens d and kruskal test. 
"""


def post_hoc_mwu(group1:pd.core.series.Series, group2:pd.core.series.Series, 
                 group3:pd.core.series.Series, fwe_method:str='holm-sidak') -> tuple:
    
    '''
    Function that runs man whiteny u tests corrected for multiple comparisons 

    Parameters
    -----------

    group1:pd.core.series.Series As array, Pandas series of data for 1st group
    group2:pd.core.series.Series As array, Pandas series of data for 2nd group
    group3:pd.core.series.Series As array, Pandas series of data for 3rd group
    fwe_method:str method to control fwe rate (optional, default is holm-sidak)
    
    
    Returns
    --------
    corrp : dict like object of corrected pvals
    '''


    group1_group2 = pin.mwu(group1, group2)
    group1_group3 = pin.mwu(group1, group3)
    group2_group3 = pin.mwu(group2 , group3)

    corrp= multitest.multipletests(np.concatenate((group1_group2['p-val'].values, group1_group3['p-val'].values, group2_group3['p-val'].values )), method=fwe_method)
    return corrp

def correlation(behaviour:pd.core.frame.DataFrame, volume:pd.core.frame.DataFrame, volume_name:str) -> list:
    
    '''
    Runs Spearmans correlation.
    
    Parameters
    -----------
    behaviour:pd.core.frame.DataFrame Pandas df of behaviours that wish to be correlated against a volume.
    volume:pd.core.frame.DataFrame Pandas df of volumes that wish to be correlated against behaviours
    volume:pd.core.frame.DataFrame str object, name of volume that behaviours will vbe correlated against.

    
    df with measures that area all used in the correlation. 
    
    Second array/df with thrid variable (string) that states which measure from this second array/df to be used in the correlation.
    
    Returns: 
    ---------
    Pvals:list Lists of pvals
    correlation:List List of rho values 

    '''
    
    pvalues = []
    correlation = []
    for i in behaviour.columns:
        array = pd.concat([volume[volume_name], behaviour[i]], axis=1).dropna()
        corr, pvals = stats.spearmanr(array[volume], array[i])
        pvalues.append(pvals)
        correlation.append(corr)
    return pvalues, correlation



def cohen_d(group1:pd.core.series.Series, group2:pd.core.series.Series) -> float:
    
    '''
    Calculate cohens d.
    
    Parameters: 
    ------------
    group1:pd.core.series.Series array or pandas series to test for effect size.
    group2:pd.core.series.Series array or pandas series to test for effect size.


    Returns
    --------
    Output:float cohen's d value.
    
    '''
    
    
    diff = group1.mean() - group2.mean()
    pooledstdev = math.sqrt((group1.std()**2 + group2.std())/2)
    cohend = diff / pooledstdev
    return cohend

def kruskal(group1:pd.core.series.Series, group2:pd.core.series.Series, group3:pd.core.series.Series) -> pd.core.frame.DataFrame:
    
    '''
    Runs kruskal-wallis test.
    
    
    Parameters
    -----------
    group1:pd.core.series.Series array, Pandas series of data for 1st group
    group2:pd.core.series.Series array, Pandas series of data for 2nd group
    group3:pd.core.series.Series array, Pandas series of data for 3rd group


    Returns
    --------
    df : dataframe of pvals, degrees of freedom and eta/epsilon effect sizes.
    
    
    '''
    
    number_of_observations = len(pd.concat([group1, group2, group3]))
    number_of_groups = len([group1, group2, group3])
    h,p = stats.kruskal(group1, group2, group3)
    degrees_of_freedom = number_of_observations - number_of_groups
    eta = (h-number_of_groups +1)/(number_of_observations-number_of_groups)
    epsilon  = h/((number_of_observations**2 -1) / (number_of_observations +1))
    df= pd.DataFrame(data={'pval': [p], 'kruskal_test_statistic': [h],'df': [degrees_of_freedom],'eta': [eta],'epsilon': [epsilon]})
    return df


def multi_comparisons(dictionary:dict, aan:pd.core.frame.DataFrame, wr:pd.core.frame.DataFrame, hc:pd.core.frame.DataFrame) -> dict:

    '''
    Function to test if parametric or non-parametric multi comparisons test should be used
    
    Parameters
    ----------- 
    dictionary:dict object of anova/kruskal-wallis results
    aan:pd.core.frame.DataFrame DataFrame of aan results
    wr:pd.core.frame.DataFrame DataFrame of wr results
    hc:pd.core.frame.DataFrame DataFrame of hc results


    Returns
    -------
    multi_comp_dict:dict of sorted values into parametric or non-parametric 

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

def mean_values(measures:str, aan:pd.core.frame.DataFrame, wr:pd.core.frame.DataFrame, hc:pd.core.frame.DataFrame) -> dict:
    
    '''
    Function to calculate the mean value and std deviation of

    Parameters
    -----------

    measure:str  list of measures
    aan:pd.core.frame.DataFrame DataFrame of aan results
    wr:pd.core.frame.DataFrame DataFrame of wr results
    hc:pd.core.frame.DataFrame DataFrame of hc results

    Returns
    --------
    value_dict (dict): dictionary of mean and std of measures

    '''

    aan_mean_values= dict(zip(['aan_'+ measure for measure in measures], [aan[measure].mean() for measure in measures]))
    aan_std = dict(zip(['aan_std_'+ measure for measure in measures], [aan[measure].std() for measure in measures]))
      
    wr_mean_values= dict(zip(['wr_'+ measure for measure in measures], [wr[measure].mean() for measure in measures]))
    wr_std = dict(zip(['wr_std_'+ measure for measure in measures], [wr[measure].std() for measure in measures]))
      
    hc_mean_values= dict(zip(['hc_'+ measure for measure in measures], [hc[measure].mean() for measure in measures]))
    hc_std = dict(zip(['hc_std_'+ measure for measure in measures], [hc[measure].std() for measure in measures]))
      
    value_dict = {**aan_mean_values, **aan_std, **wr_mean_values, **wr_std, **hc_mean_values, **hc_std}

    return value_dict

def ttest_mwu(group1:pd.core.series.Series, group2:pd.core.series.Series) -> pd.core.frame.DataFrame:

    '''
    Function to test for group differences. Will run T-test if assumptions met or run Man Whitney U.

    Parameters
    ----------
    group1:pd.core.series.Series array, Pandas series of data for 1st group
    group2:pd.core.series.Series array, Pandas series of data for 2nd group
    
    Returns
    --------
    test:pandas.core.frame.DataFrame df of either mann whiteney U results or T-Test

    '''

    if stats.normaltest(group1.dropna()) and stats.normaltest(group2.dropna()) and stats.levene(group1.dropna(), group2.dropna())[1] > 0.05:
        test= pin.ttest(group1.dropna(),group2.dropna())
    
    else:
        test = pin.mwu(group1.dropna(),group2.dropna())
    
    return test