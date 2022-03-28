'''
Script to analyse cluster direction of signficance for clusters  
'''

import pandas as pd
from scipy import stats
import statsmodels.api as sm
import os 
from colorama import Fore
import functions as fun
from decouple import config

data = config('data')

os.chdir(data)

cluster = pd.read_csv('cortical_measures.csv').drop(columns=['mean_curv', 'mean_lgi', 'mean_area','mean_thickness','TotalGrayVol','Total_white_matter'])

#Splits the dataframe into groups
group = cluster.groupby('age_adjusted_group')
hc = group.get_group('HC')
aan= group.get_group('AAN')
wr = group.get_group('WR')


not_cluster =['Unnamed: 0', 'G-Number', 'age_adjusted_group', 'Age','Group']

measures = [measure for measure in cluster.columns if measure not in not_cluster]

mean = fun.mean_values(measures,aan,wr,hc)

for betas in cluster.columns:
    
    if betas not in not_cluster:    
        if stats.normaltest(aan[betas].dropna()) and stats.normaltest(wr[betas].dropna()) and stats.normaltest(hc[betas].dropna()) and stats.levene(aan[betas].dropna(), wr[betas].dropna(),hc[betas].dropna())[1] > 0.05:
            multi_test= sm.stats.multicomp.pairwise_tukeyhsd(cluster[betas], cluster['age_adjusted_group'])
            print(Fore.GREEN + f'\n{betas}' + Fore.RESET, multi_test.summary())

            if multi_test.reject[0] == True:
                if mean['aan_'+ betas] > mean['hc_' + betas]:
                    aan_hc= fun.cohen_d(aan[betas],hc[betas])
                              
                else:
                    aan_hc= fun.cohen_d(hc[betas],aan[betas])
                    print(Fore.GREEN + f'\nCohens d between AAN and HC for {betas}:' + Fore.RESET, aan_hc)

            if multi_test.reject[1] == True:            
                if mean['aan_'+ betas] > mean['wr_' + betas]:
                    aan_wr= fun.cohen_d(aan[betas],wr[betas])
                              
                else:
                    aan_wr= fun.cohen_d(wr[betas],aan[betas])
                print(Fore.GREEN + f'\nCohens d between AAN and wr for {betas}:' + Fore.RESET, aan_wr)

            if multi_test.reject[2] == True:
                if mean['hc_'+ betas] > mean['wr_' + betas]:
                    hc_wr= fun.cohen_d(hc[betas],wr[betas])
                              
                else:
                    hc_wr= fun.cohen_d(wr[betas],aan[betas])

                print(Fore.GREEN + f'\nCohens d between AAN and wr for {betas}:' + Fore.RESET, hc_wr)
    
        else:            
            mwu_corrp= fun.post_hoc_mwu(aan[betas].dropna(), hc[betas].dropna(), wr[betas].dropna())
            print(Fore.GREEN + f'\nBetas for {betas}'+ Fore.RESET,
                '\n AAN-HC:',mwu_corrp[0][0],'p= ',mwu_corrp[1][0], 
                'AAN-WR:', mwu_corrp[0][1], 'p= ',mwu_corrp[1][1], 
                'HC-WR:', mwu_corrp[0][2], 'p= ',mwu_corrp[1][2]
                 )

            if mwu_corrp[0][0] == True:
                if mean['aan_'+ betas] > mean['hc_' + betas]:
                    aan_hc = fun.cohen_d(aan[betas].dropna(), hc[betas].dropna())
                              
                else:
                    aan_hc= fun.cohen_d(hc[betas].dropna(), aan[betas].dropna())
                    print(Fore.GREEN+ f'\nCohens d between AAN and HC for {betas}:' + Fore.RESET, aan_hc)

                if mwu_corrp[0][1] == True:
                    if mean['aan_'+ betas] > mean['wr_' + betas]:
                        aan_wr= fun.cohen_d(aan[betas].dropna(), wr[betas].dropna())
                              
                    else:
                        aan_wr= fun.cohen_d(wr[betas].dropna(), aan[betas].dropna())
                        
                    print(Fore.GREEN + f'\nCohens d between AAN and WR for {betas}:' + Fore.RESET, aan_wr)

                if mwu_corrp[0][2] == True:                              
                    if mean['hc_'+ betas] > mean['wr_' + betas]:
                        hc_wr= fun.cohen_d(hc[betas].dropna(), wr[betas].dropna())
                              
                    else:
                        hc_wr= fun.cohen_d(wr[betas].dropna(), hc[betas].dropna())
                    print(Fore.GREEN+ f'\nCohens d between HC and Wr for {betas}:' + Fore.RESET, hc_wr)

print(Fore.MAGENTA + '\nMean and STD for Betas:\n'+ Fore.RESET)

for key, val in mean.items():
      print(Fore.CYAN + key + Fore.RESET,':',val)