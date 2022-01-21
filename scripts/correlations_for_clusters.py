import pandas as pd
from statsmodels.stats import multitest
import numpy as np
from colorama import Fore
import os
import functions as fun

'''
This script runs correlations between cluster parameters and behavioural measures using spearmans r

Not very well written (too many variables!!) but does the job.

Prints a dictionary with pvals corrected for multiple comparisons and rho values
'''


os.chdir() #change this to your local file path

global_measures = pd.read_csv('cortical_measures.csv').drop(columns=['mean_curv', 'mean_lgi', 'mean_area','mean_thickness','TotalGrayVol','Total_white_matter'])
behavioural = pd.read_csv('/home/wmda/Documents/Documents_from_old_comp/BEACON/Write_up/Cortical/results/behavioural_results.csv')

corr_df = pd.concat([global_measures,behavioural[['BMI_at_scan','Initial_EDE_Q_Total','IQ', 'ADOS_com_soc', 'ADOS_Creativity',
                                                  'ADOS_sterotyped_and_repetititve','Illness_duration']]], axis=1)

group = corr_df.groupby('age_adjusted_group')
hc = group.get_group('HC')
aan = group.get_group('AAN')
wr = group.get_group('WR')

behavioural_correlations_aan = aan[['Illness_duration', 'Initial_EDE_Q_Total', 'ADOS_com_soc','ADOS_Creativity',
                                   'ADOS_sterotyped_and_repetititve','BMI_at_scan','Age']]

global_measures_correlation_aan = aan[['postcentral','supramarginal']]

behavioural_correlations_wr = wr[['Illness_duration','Initial_EDE_Q_Total', 'ADOS_com_soc','ADOS_Creativity',
                                  'ADOS_sterotyped_and_repetititve','BMI_at_scan','Age']]

global_measures_correlation_wr = wr[['postcentral','supramarginal']]

behavioural_correlations_hc = hc[['Initial_EDE_Q_Total', 'ADOS_com_soc','ADOS_Creativity', 
                                  'ADOS_sterotyped_and_repetititve','BMI_at_scan','Age']]

global_measures_correlation_hc = hc[['postcentral','supramarginal']]

aan_c1p, aan_c1r = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'postcentral')
aan_c2p, aan_c2r = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'supramarginal')

wr_c1p, wr_c1r = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'postcentral')
wr_c2p, wr_c2r = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'supramarginal')

hc_c1p, hc_c1r = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'postcentral')
hc_c2p, hc_c2r = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'supramarginal')

corrp_global_measures = multitest.multipletests(np.concatenate((
                                                       aan_c1p, aan_c2p,
                                                       wr_c1p, wr_c2p,
                                                       hc_c1p, hc_c2p
                                                       )))

aan_keys_measures = ['AAN_'+ measure for measure in behavioural_correlations_aan.columns] 
aan_keys_global = [measure for measure in global_measures_correlation_aan.columns]

wr_keys_measures = ['WR_'+ measure for measure in behavioural_correlations_wr.columns] 
wr_keys_global = [measure for measure in global_measures_correlation_wr.columns]

hc_keys_measures = ['HC_'+ measure for measure in behavioural_correlations_hc.columns] 
hc_keys_global = [measure for measure in global_measures_correlation_hc.columns]

dictionary_keys=[]
for glob in aan_keys_global:
    for measure in aan_keys_measures:
        dictionary_keys.append(measure + '_' + glob)

for glob in wr_keys_global:
    for measure in wr_keys_measures:
        dictionary_keys.append(measure + '_' + glob)

for glob in hc_keys_global:
    for measure in hc_keys_measures:
        dictionary_keys.append(measure + '_' + glob)

p_vals= dict(zip(dictionary_keys,corrp_global_measures[1]))

r2_vals = aan_c1r + aan_c2r + wr_c1r + wr_c2r +  hc_c1r + hc_c2r 
 
r2_val_dict = dict(zip(dictionary_keys, r2_vals))

print(Fore.MAGENTA + '\n Pvals for correlation \n' + Fore.RESET)
for key, val in p_vals.items():
    print(Fore.YELLOW + key + Fore.RESET,':',val)

print(Fore.MAGENTA + '\n Rsquared for correlation \n' + Fore.RESET)
for key, val in r2_val_dict.items():
    print(Fore.YELLOW + key + Fore.RESET, ':', val)