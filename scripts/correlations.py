import pandas as pd
from statsmodels.stats import multitest
import numpy as np
from colorama import Fore
import os
import functions as fun

'''
This script runs correlations between global measures and behavioural measures using spearmans r

Not very well written (too many variables!!) but does the job.

Prints a dictionary with pvals corrected for multiple comparisons and rho values
'''


os.chdir() #change this to your local file path

global_measures = pd.read_csv('cortical_measures.csv')
behavioural = pd.read_csv('behavioural_results.csv')

corr_df = pd.concat([global_measures,behavioural[['BMI_at_scan','Initial_EDE_Q_Total','IQ', 'ADOS_com_soc', 
                                                  'ADOS_Creativity','ADOS_sterotyped_and_repetititve','Illness_duration']]], 
                                                  axis=1)

group = corr_df.groupby('age_adjusted_group')
hc = group.get_group('HC')
aan = group.get_group('AAN')
wr = group.get_group('WR')

behavioural_correlations_aan = aan[['Illness_duration','Initial_EDE_Q_Total', 'ADOS_com_soc','ADOS_Creativity',
                                     'ADOS_sterotyped_and_repetititve','BMI_at_scan','Age']]

global_measures_correlation_aan = aan[['TotalGrayVol','Total_white_matter', 'mean_thickness','mean_curv','mean_lgi','mean_area' ]]

behavioural_correlations_wr = wr[['Illness_duration','Initial_EDE_Q_Total', 'ADOS_com_soc','ADOS_Creativity',
                                  'ADOS_sterotyped_and_repetititve','BMI_at_scan','Age']]

global_measures_correlation_wr = wr[['TotalGrayVol','Total_white_matter', 'mean_thickness','mean_curv','mean_lgi','mean_area' ]]

behavioural_correlations_hc = hc[['Initial_EDE_Q_Total', 'ADOS_com_soc','ADOS_Creativity',
                                  'ADOS_sterotyped_and_repetititve','BMI_at_scan','Age']]

global_measures_correlation_hc = hc[['TotalGrayVol','Total_white_matter', 'mean_thickness','mean_curv','mean_lgi','mean_area']]

aan_gmvp, aan_gmvr = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'TotalGrayVol')
aan_wmvp, aan_wmvr = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'Total_white_matter')

wr_gmvp, wr_gmvr = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'TotalGrayVol')
wr_wmvp, wr_wmvr = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'Total_white_matter')

hc_gmvp, hc_gmvr = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'TotalGrayVol')
hc_wmvp, hc_wmvr = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'Total_white_matter')

aan_ctp, aan_ctr = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'mean_thickness')
aan_ccp, aan_ccr = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'mean_curv')
aan_lgip, aan_lgir = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'mean_lgi')
aan_csap, aan_csar = fun.correlation(behavioural_correlations_aan, global_measures_correlation_aan, 'mean_area')

wr_ctp, wr_ctr = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'mean_thickness')
wr_ccp, wr_ccr = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'mean_curv')
wr_lgip, wr_lgir = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'mean_lgi')
wr_csap, wr_csar = fun.correlation(behavioural_correlations_wr, global_measures_correlation_wr, 'mean_area')

hc_ctp, hc_ctr = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'mean_thickness')
hc_ccp, hc_ccr = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'mean_curv')
hc_lgip, hc_lgir = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'mean_lgi')
hc_csap, hc_csar = fun.correlation(behavioural_correlations_hc, global_measures_correlation_hc, 'mean_area')

corrp_global_measures = multitest.multipletests(np.concatenate((
                                                               aan_gmvp, aan_wmvp, aan_ctp,aan_ccp, aan_lgip, aan_csap, 
                                                               wr_gmvp,wr_wmvp, wr_ctp, wr_ccp, wr_lgip, wr_csap, 
                                                               hc_gmvp, hc_wmvp, hc_ctp, hc_ccp, hc_lgip,hc_csap
                                                               )))

aan_keys_measures = ['AAN_'+ measure for measure in behavioural_correlations_aan.columns] 
aan_keys_global = [measure for measure in global_measures_correlation_aan.columns]

wr_keys_measures = ['WR_'+ measure for measure in behavioural_correlations_wr.columns] 
wr_keys_global = [measure for measure in global_measures_correlation_wr.columns]

hc_keys_measures = ['HC_'+ measure for measure in behavioural_correlations_hc.columns] 
hc_keys_global = [measure for measure in global_measures_correlation_hc.columns]



dictionary_keys = []

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

r2_vals = aan_gmvr + aan_wmvr + aan_ctr + aan_ccr + aan_lgir + aan_csar + wr_gmvr + wr_wmvr + wr_ctr + wr_ccr + wr_lgir + wr_csar +  hc_gmvr + hc_wmvr + hc_ctr + hc_ccr + hc_lgir + hc_csar
 
r2_val_dict = dict(zip(dictionary_keys, r2_vals))

print(Fore.MAGENTA + '\nPvals for correlation\n' + Fore.RESET)
for key, val in p_vals.items():
    print(Fore.YELLOW + key + Fore.RESET,':',val)

print(Fore.MAGENTA + '\nRsquared for correlation\n' + Fore.RESET)
for key, val in r2_val_dict.items():
    print(Fore.YELLOW + key + Fore.RESET,':',val)
