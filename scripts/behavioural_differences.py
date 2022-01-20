import pandas as pd
import functions as fun
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from colorama import Fore
import os

os.chdir() #change this to your local file path

df = pd.read_csv('behavioural_results.csv')

group = df.groupby('age_adjusted_group')
hc = group.get_group('HC')
aan = group.get_group('AAN')
wr = group.get_group('WR')

measures = ['BMI_at_scan','Initial_EDE_Q_Total','Age','IQ', 'ADOS_com_soc', 
            'ADOS_Creativity','ADOS_sterotyped_and_repetititve'] 

mean = fun.mean_values(measures,aan,wr,hc)

anova_list = [measure for measure in measures if stats.normaltest((ols(f'{measure} ~ age_adjusted_group',data = df).fit().resid))[1] > 0.05
              and stats.levene(aan[measure].dropna(), wr[measure].dropna(),hc[measure].dropna())[1] > 0.05] 
              

kruskal_list = [measure for measure in measures if measure not in anova_list]

anova_results_dict = dict(zip(anova_list,[sm.stats.anova_lm(ols(f'{measure} ~ age_adjusted_group', data= df).fit(), typ=1) for measure in anova_list]))

kruskal_results_dict = dict(zip(kruskal_list,[fun.kruskal(aan[measure].dropna(), hc[measure].dropna(), wr[measure].dropna())for measure in kruskal_list]))

for result in anova_results_dict:
    print('\n\n',Fore.BLUE + result + Fore.RESET,'\n',anova_results_dict[result])

for result in kruskal_results_dict:
    print('\n\n',Fore.BLUE + result + Fore.RESET,'\n',kruskal_results_dict[result])

aov_res = fun.multi_comparisons(anova_results_dict,aan,wr,hc)

krus_res = fun.multi_comparisons(kruskal_results_dict,aan,wr,hc)

assumptions_check = [aov_res, krus_res]

for assumption in assumptions_check:

      if len(assumption['parametric'])  > 0:

            for measure in assumption['parametric']:

                  testdf= df[[measure,'age_adjusted_group']].dropna()

                  multi_test = sm.stats.multicomp.pairwise_tukeyhsd(testdf[measure], testdf['age_adjusted_group'])
                  
                  print(Fore.BLUE + f'\n{measure}\n\n'+ Fore.RESET,multi_test.summary())
                  
                  if multi_test.reject[0] == True:

                        if mean['aan_'+ measure] > mean['hc_' + measure]:
                                    
                              aan_hc= fun.cohen_d(aan[measure].dropna(),hc[measure].dropna())
                              
                        else:
                              aan_hc= fun.cohen_d(hc[measure].dropna(),aan[measure].dropna())
                        
                       
                        print(Fore.GREEN + f'\nCohens d between AAN and HC for {measure}:' + Fore.RESET, aan_hc)

                  if multi_test.reject[1] == True:
                              
                        if mean['aan_'+ measure] > mean['wr_' + measure]:
                                    
                              aan_wr= fun.cohen_d(aan[measure].dropna(),wr[measure].dropna())
                              
                        else:
                              aan_wr= fun.cohen_d(wr[measure].dropna(),aan[measure].dropna())
                        
                       
                        print(Fore.GREEN + f'\nCohens d between AAN and WR for {measure}:' + Fore.RESET, aan_wr)

                  
                  if multi_test.reject[2] == True:
                              
                        if mean['hc_'+ measure] > mean['wr_' + measure]:
                                    
                              hc_wr= fun.cohen_d(hc[measure].dropna(),wr[measure].dropna())
                              
                        else:
                              hc_wr= fun.cohen_d(wr[measure].dropna(),aan[measure].dropna())
                        
                       
                        print(Fore.GREEN+ f'\nCohens d between HC and Wr for {measure}:' + Fore.RESET, hc_wr)
                  
      if len(assumption['nonparametric']) > 0:

            for measure in assumption['nonparametric']:

                  mwu_corrp= fun.post_hoc_mwu(aan[measure].dropna(),hc[measure].dropna(),wr[measure].dropna())

                  print(Fore.BLUE + f'\n{measure}\n' + Fore.RESET, 'AAN-HC:',mwu_corrp[0][0],'p= ',mwu_corrp[1][0], 'AAN-WR:', mwu_corrp[0][1], 'p= ',mwu_corrp[1][1], 

                        'HC-WR:', mwu_corrp[0][2], 'p= ',mwu_corrp[1][2])

                  if mwu_corrp[0][0] == True:

                        if mean['aan_'+ measure] > mean['hc_' + measure]:
                                    
                              aan_hc= fun.cohen_d(aan[measure].dropna(),hc[measure].dropna())
                              
                        else:
                              aan_hc= fun.cohen_d(hc[measure].dropna(),aan[measure].dropna())
                        
                       
                        print(Fore.GREEN + f'\nCohens d between AAN and HC for {measure}:' + Fore.RESET, aan_hc)

                  if mwu_corrp[0][1] == True:
                              
                        if mean['aan_'+ measure] > mean['wr_' + measure]:
                                    
                              aan_wr= fun.cohen_d(aan[measure].dropna(),wr[measure].dropna())
                              
                        else:
                              aan_wr= fun.cohen_d(wr[measure].dropna(),aan[measure].dropna())
                        
                       
                        print(Fore.GREEN + f'\nCohens d between AAN and WR for {measure}:' + Fore.RESET, aan_wr)
                  
                  if mwu_corrp[0][2] == True:
                              
                        if mean['hc_'+ measure] > mean['wr_' + measure]:

                              hc_wr= fun.cohen_d(hc[measure].dropna(),wr[measure].dropna())
                              
                        else:
                              hc_wr= fun.cohen_d(wr[measure].dropna(),hc[measure].dropna())
                        
                       
                        print(Fore.GREEN + f'\nCohens d between HC and Wr for {measure}:' + Fore.RESET, hc_wr)

print(Fore.BLUE + '\nIllness Duraion Difference','\n' + Fore.RESET,fun.ttest_mwu(aan['Illness_duration'],wr['Illness_duration']))

print(Fore.MAGENTA + '\nMean and STD for Behavioural Values:\n'+ Fore.RESET)
for key, val in mean.items():
      print(Fore.CYAN + key + Fore.RESET,':',val)
      
print(Fore.CYAN + 'AAN_Illness_duration'+ Fore.RESET, aan['Illness_duration'].mean(),'\n',
      Fore.CYAN + 'WR_Illness_duration'+ Fore.RESET, wr['Illness_duration'].mean(),'\n',
      Fore.CYAN + 'AAN_std_Illness_duration'+ Fore.RESET, aan['Illness_duration'].std(),'\n',
      Fore.CYAN + 'WR_std_Illness_duration'+ Fore.RESET, wr['Illness_duration'].std(),'\n',
)
